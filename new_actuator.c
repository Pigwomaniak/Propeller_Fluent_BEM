#include "udf.h"								//standard library
#include "aerofoilData.h"							//aerofoil data
#include "bladeRotorGeometry.h"							//blade geometry
#include "caseFlowSpecifics.h"							//general specifications for setting up the case

/*******************************************************\
|							|
	ADM model for wind turbines ver. 3.0

	Updates from previous verions:

	v3.0
	- changed architecture to accomodate 2nd rotor
	- NO REVERSE COMPATIBILITY

	v2.1
	- instead of reading omega, now reads tsr
	- handling of tsr added

	v2.0
	- reservation and naming of UDMs at loading
|							|
\*******************************************************/



/*******************************************************\
|							|
|	Declarations					|
|							|
\*******************************************************/

/*maths*/
#define pi		3.1415926535								//pi
#define intwind		5		
								//interpolation window size (number of nodes for interpolation)

#define ACTIVE_PROPS_NUMBER	1

/*global variables*/
double dforcen,dforcet,source;								//actuator forces, source term
double xyz[ND_ND];											//position vector (table of coordinates)
double x,y,z,r;												//coordinates
double vx,vy,vz,vr;											//flow velocity components
double u,vt,w;												//velocity triangle components (axial, tangential, inflow)
double V,A[ND_ND];											//volume, cell area table
double MFR;													//mass flow rate
double pdyn;												//dynamic pressure
double alpha,phi,theta;										//flow angles
double beta,chord;											//blade characteristics
double Cl,Cd,Clbis,Cdbis;									//local force coefficients
double delta;												//precision epsilon
double f1;													//Prandtl correction
double turbepsilon,turbomega,TKE;

double power[Nrotors], powertot[Nrotors];					//variables for summing up the decremental powers
double powerwtun[Nrotors], powerwreal[Nrotors];				//wind power
double Cptun[Nrotors], Cpreal[Nrotors];						//power coefficient
double thrust[Nrotors], thrustot[Nrotors];					//variables for summing up the decremental thrusts
double thrustun[Nrotors], thrustreal[Nrotors];				//wind normal force
double Cttun[Nrotors], Ctreal[Nrotors];						//thrust coefficient
double Rmin[Nrotors], Rmax[Nrotors];
double zmin[Nrotors], zmax[Nrotors];
double diskid[Nrotors], Nbl[Nrotors];
double N_one_own;

double rotDir[Nrotors];										//rotor rotation direction
double vzave=0.;											//average inlet velocity in z direction, initialized to 0
double velref,omega[Nrotors],tsr[Nrotors],rotDir[Nrotors];
int profileNo;												//type of profile at selected radius

/*functions*/
void BEM(cell_t,Thread*,int);
void aitken_0(double*,double*,int,double,double*);
void aitken(double*,double*,int,double,double*,double*);
double prandtl(double,double,double,double);



/*******************************************************\
|							|
|	Code						|
|							|
\*******************************************************/


  
/*various commands to executed at library loading*/
DEFINE_EXECUTE_ON_LOADING(report_version, firestarter)
{
	int k=0;

	Reserve_User_Memory_Vars(6);

	/*set the aerofoil data table*/
	polarDataSetter(k,profilea,Rea,noCoeffa,alfaa,Cla,Cda); k++;	
	polarDataSetter(k,profilea,Reb,noCoeffa,alfaa,Clb,Cdb); k++;
	polarDataSetter(k,profilea,Rec,noCoeffa,alfaa,Clc,Cdc); k++;
	polarDataSetter(k,profilea,Red,noCoeffa,alfaa,Cld,Cdd); k++;
	polarDataSetter(k,profilea,Ree,noCoeffa,alfaa,Cle,Cde); k++;
	polarDataSetter(k,profilea,Ref,noCoeffa,alfaa,Clf,Cdf); k++;
	polarDataSetter(k,profilea,Reg,noCoeffa,alfaa,Clg,Cdg); k++;	
	polarDataSetter(k,profilea,Reh,noCoeffa,alfaa,Clh,Cdh); k++;	
	polarDataSetter(k,profilea,Rei,noCoeffa,alfaa,Cli,Cdi); k++;
	polarDataSetter(k,profilea,Rej,noCoeffa,alfaa,Clj,Cdj); k++;

	/*set case variables*/
	varset();

	/*vectorise variables*/
	diskid[0]=disk_id_one; diskid[1]=disk_id_two;
	Rmin[0]=Rmin_one; Rmin[1]=Rmin_two;
	Rmax[0]=Rmax_one; Rmax[1]=Rmax_two;
	zmin[0]=zmin_one; zmin[1]=zmin_two;
	zmax[0]=zmax_one; zmax[1]=zmax_two;
	Nbl[0]=N_one; Nbl[1]=N_two;
}

/*rename the UDMs at simulation start*/
DEFINE_ON_DEMAND(UDMnamesetter)
{
	Set_User_Memory_Name(0,"LocalReynolds number");
	Set_User_Memory_Name(1,"LocalAoA");
	Set_User_Memory_Name(2,"LocalCl");
	Set_User_Memory_Name(3,"LocalCd");
	Set_User_Memory_Name(4,"Cp");
	Set_User_Memory_Name(5,"Ct");
}

/*reread omega and velref values*/
DEFINE_ON_DEMAND(VARrereader)
{
	varset();
}

/*define axial velocity profile at inlet (see UDF manual p91)*/
DEFINE_PROFILE(z_velocity,t,i)
{
	face_t f;
	double wx,wy,win;
	double multiplier;

	begin_f_loop(f,t)
	{
		F_CENTROID(xyz,f,t);							//stores the cell centroid data in xyz
		x=xyz[0]; y=xyz[1];							//defines x, y
		if (x==0 && y==0)	r=0.000001;
		else			r=sqrt(x*x+y*y);
		phi=asin(fabs(y/r));

		if (x>=0) 	linearInterpolation(xr,xw,17,r,&wx);
		else		linearInterpolation(xr,xw,17,-r,&wx);
		if (y>=0)	linearInterpolation(yr,yw,17,r,&wy);
		else		linearInterpolation(yr,yw,17,-r,&wy);

		win=(wy*phi+wx*(pi/2-phi))/(pi/2);
		multiplier=velref/17.2535725;
		F_PROFILE(f,t,i)=win*multiplier;			
	}
	end_f_loop(f,t)
}

/*define source terms in x direction (UDF manual p125)*/
DEFINE_SOURCE(xmom_source,c,t,dS,eqn)
{
	BEM(c,t,0);											//call BEM, determine source force
	return source;										//returns xmom_source								
}

/*define source terms in y direction (UDF manual p125)*/
DEFINE_SOURCE(ymom_source,c,t,dS,eqn)
{
	BEM(c,t,1);											//call BEM, determine source force
	return source;										//returns xmom_source
}

/*define source terms in z (axial) direction (UDF manual p125)*/
DEFINE_SOURCE(zmom_source,c,t,dS,eqn)
{
	BEM(c,t,2);											//call BEM, determine source force
	return source;										//returns zmom_source
}

/*define source terms for TKE (UDF manual p125)*/
DEFINE_SOURCE(k_source,c,t,dS,eqn)
{
	BEM(c,t,3);											//call BEM, determine TKE
	return source;										//returns TKE
}

/*Vz calculation at the end of an iteration, parallel mode compliant*/
DEFINE_EXECUTE_AT_END(vel_in_MFRave)
{
	#if !RP_HOST										//Host will do nothing in this udf. Serial will
		double VzMFRtotn[1],MFRtotn[1];					//tables for data communication
		double VzMFRtot,MFRtot;							//ut supra
		int i,nelem=1;									//looping variables
		face_t f;										//define face
		cell_t c0;										//define cell
		Domain *domain;									//define domain [pointer to]
		Thread *t, *t0, *tbis;							//define threads [pointers to]
		domain=Get_Domain(1);							//assign domain value
		t=Lookup_Thread(domain,inlet_id);				//assign thread value
		VzMFRtot=0.0;									//zero the variables...
		MFRtot=0.0;										//...before the macro is executed

		begin_f_loop(f,t)								//principal face not applied, as the studied face is boundary
		{
			c0=F_C0(f,t);								//!!! get the 0 cell adjacent to inlet (inside domain)
			t0=THREAD_T0(t);							//!!! get the 0 thread adjacent to inlet (inside domain)
			vz=C_W(c0,t0);								//!!! without the above this will cause error
			MFR=F_FLUX(f,t);							//!!! execute on thread t and NOT t0
		
			VzMFRtotn[0]=VzMFRtotn[0]+vz*MFR;			//store the local data in the tables
			MFRtotn[0]=MFRtotn[0]+MFR;					//ut supra
		}
		end_f_loop(f,t)									//end face looping

		#if !PARALLEL
			vzave=VzMFRtotn[0]/MFRtotn[0];				//if calculation mode is not parallel, calculate the value from single node
		#endif

		#if PARALLEL
			#if RP_NODE
				if(! I_AM_NODE_ZERO_P)							
				{
					PRF_CSEND_REAL(node_zero,VzMFRtotn,nelem,myid);		//send, to node zero, ARRAY of nelem elements with the node id myid
					PRF_CSEND_REAL(node_zero,MFRtotn,nelem,myid);		//ut supra
				}
			#endif
     
			#if RP_NODE
				if(I_AM_NODE_ZERO_P)
				{
					VzMFRtot=VzMFRtotn[0];					//store local node data in global variable
					MFRtot=MFRtotn[0];					//ut supra
				
					compute_node_loop_not_zero(i)
					{
						PRF_CRECV_REAL(i,VzMFRtotn,nelem,i);		//recive from (non-zero) node i an ARRAY of nelem elements
						PRF_CRECV_REAL(i,MFRtotn,nelem,i);		//ut supra
						VzMFRtot=VzMFRtot+VzMFRtotn[0];			//store node data in global variable
						MFRtot=MFRtot+MFRtotn[0];			//ut supra
					}

					vzave=VzMFRtot/MFRtot;					//compute the velocity magnitude from global data
				}
			#endif
		#endif
	#endif
}

/*Cp calculation at the end of an iteration, from BEM data calculated previously, parallel mode compliant*/
DEFINE_EXECUTE_AT_END(power_Cp_ext)
{
	#if !RP_HOST														//Host will do nothing in this udf. Serial will
	int Niter=N_ITER;
	FILE * fp;
	int i,j;

	for (j=0;j<Nrotors;j++)
	{
		powertot[j]=power[j];
		powerwtun[j]=rho*pi*Rmax[j]*Rmax[j]*velref*velref*velref/2;
		powerwreal[j]=rho*pi*Rmax[j]*Rmax[j]*vzave*vzave*vzave/2;

		thrustot[j]=thrust[j];
		thrustun[j]=rho*pi*Rmax[j]*Rmax[j]*velref*velref/2;
		thrustreal[j]=rho*pi*Rmax[j]*Rmax[j]*vzave*vzave/2;
	}

	#if PARALLEL
		#if RP_NODE
			if(! I_AM_NODE_ZERO_P)							
			{
				PRF_CSEND_REAL(node_zero,power,Nrotors,myid);			//send, to node zero, ARRAY of nelem elements with the node id myid
				PRF_CSEND_REAL(node_zero,thrust,Nrotors,myid);			//send, to node zero, ARRAY of nelem elements with the node id myid
			}
		#endif
     
		#if RP_NODE
			if(I_AM_NODE_ZERO_P)
			{
				compute_node_loop_not_zero(i)
				{
					PRF_CRECV_REAL(i,power,Nrotors,i);			//recive from (non-zero) node i an ARRAY of nelem elements
					PRF_CRECV_REAL(i,thrust,Nrotors,i);			//recive from (non-zero) node i an ARRAY of nelem elements
					powertot[0]=powertot[0]+power[0];			//store node data in global variable
					powertot[1]=powertot[1]+power[1];			//store node data in global variable
					thrustot[0]=thrustot[0]+thrust[0];			//store node data in global variable
					thrustot[1]=thrustot[1]+thrust[1];			//store node data in global variable
				}

				for (j=0;j<Nrotors;j++)
				{
					Cptun[j]=powertot[j]/powerwtun[j];		//calculate Cp
					Cpreal[j]=powertot[j]/powerwreal[j];		//calculate Cp
					Cttun[j]=thrustot[j]/thrustun[j];		//calculate Ct
					Ctreal[j]=thrustot[j]/thrustreal[j];		//calculate Ct
				}

				Message("Velref = %f m/s, Velave = %f m/s\nTSR1 = %f [-], Omega1 = %f rad/s, Power1 = %f W, Cptunnel1 = %f, Cpreal1= %f, Thrust1 = %f N, Cttunnel1 = %f, Ctreal1= %f\nTSR2 = %f [-], Omega2 = %f rad/s, Power2 = %f W, Cptunnel2 = %f, Cpreal2= %f, Thrust2 = %f N, Cttunnel2 = %f, Ctreal2= %f\n\n",velref,vzave,tsr[0],omega[0],powertot[0],Cptun[0],Cpreal[0],thrustot[0],Cttun[0],Ctreal[0],tsr[1],omega[1],powertot[1],Cptun[1],Cpreal[1],thrustot[1],Cttun[1],Ctreal[1]);

				fp = fopen ("CpCase.txt", "a");
				fprintf(fp,"iteration %d, Wind speed at -0.2m = %f m/s, Wind speed averaged = %f m/s, TSR1 = %f [-], Omega1 = %f rad/s, Power1 = %f W, Cptunnel1 = %f, Cpreal1= %f, Thrust1 = %f N, Cttunnel1 = %f, Ctreal1= %f, TSR2 = %f [-], Omega2 = %f rad/s, Power2 = %f W, Cptunnel2 = %f, Cpreal2= %f, Thrust2 = %f N, Cttunnel2 = %f, Ctreal2= %f\n",Niter,velref,vzave,tsr[0],omega[0],powertot[0],Cptun[0],Cpreal[0],thrustot[0],Cttun[0],Ctreal[0],tsr[1],omega[1],powertot[1],Cptun[1],Cpreal[1],thrustot[1],Cttun[1],Ctreal[1]);
				fclose (fp);
				fp = fopen ("CpCaseToAve.txt", "a");
				fprintf(fp,"iteration %d Wref %f Wave %f TSR1 %f Power1 %f Pwind1 %f Thrust1 %f Twind1 %f TSR2 %f Power2 %f Pwind2 %f Thrust2 %f Twind2 %f\n",Niter,velref,vzave,tsr[0],powertot[0],powerwtun[0],thrustot[0],thrustun[0],tsr[1],powertot[1],powerwtun[1],thrustot[1],thrustun[1]);
				fclose (fp);
			}
		#endif
	#endif

	/*non-parallel mode not updated*/	
	#if !PARALLEL													//NOT UPDATED
		Cptun[0]=power[0]/powerwtun[0];											//calculate Cp
		Cptun[1]=power[1]/powerwtun[1];
		Cpreal[0]=power[0]/powerwreal[0];	
		Cpreal[1]=power[1]/powerwreal[1];										//calculate Cp

		/*Message("Omega1 = %f rad/s, Power = %f W, Wind power tunnel = %f W, Wind power averaged = %f W, Cptunnel = %f, Cpreal= %f\n",omega,power,powerwtun,powerwreal,Cptun,Cpreal);

		fp = fopen ("CpCase.txt", "a");
		fprintf(fp,"iteration %d, Omega = %f rad/s, Power = %f W, Wind speed at -0.2m = %f m/s, Wind speed averaged = %f m/s, Wind power tunnel = %f W, Wind power averaged = %f W, Cptunnel = %f, Cpreal= %f\n",Niter,omega,power,velref,vzave,powerwtun,powerwreal,Cptun,Cpreal);
		fclose (fp);*/
	#endif

	power[0]=0;						//zero data
	power[1]=0;						//zero data
	thrust[0]=0;						//zero data
	thrust[1]=0;						//zero data
#endif
}

/*BEM function for calculation of source terms*/
void BEM(cell_t c, Thread* t,int a)
{
	face_t f;								//define face
	Thread* tf;								//face thread
	int i, j, n;								//integer for face looping
	double Re;
	double dpower, TKE;
	FILE * fp;

	C_CENTROID(xyz,c,t);							//stores the cell centroid data in xyz
	x=xyz[0]; y=xyz[1]; z=xyz[2];						//defines x, y, z
	r=sqrt(x*x+y*y);							//defines radius r
	source=0;								//sets source to 0

	/*BET*/
	for (j=0;j<Nrotors;j++)
	{
		
		if (j >= ACTIVE_PROPS_NUMBER){
			source = 0;
			return;
		}
		
		if(r<=Rmax[j] && r>=Rmin[j] && (z>=zmin[j]-dz*Rmax[j]) && (z<=zmin[j]+dz*Rmax[j]))
		{
			c_face_loop(c,t,i)					//loop over all faces
			{
				tf=C_FACE_THREAD(c,t,i);			//define face thread
				if(THREAD_ID(tf) == diskid[j])			//check if cell belongs to disk
				{
					f = C_FACE(c,t,i);
				
					/*calculate the axial velocity u*/
					F_AREA(A,f,tf);								//calculate A - the cell area
					V=C_VOLUME(c,t);							//calculate V - the cell volume
					MFR=F_FLUX(f,tf);							//calculate mass flow rate
					phi=atan2(y,x);								//calculate the cell angular coordinate in XY polar
					u=C_W(c,t);
					if(u<0)
					{
						u=0;
					}

					/*calculate the tangential velocity vt*/
					vx=C_U(c,t);								//calculate x flow velocity component
					vy=C_V(c,t);								//calculate y flow velocity component
					vr=vy*cos(phi)-vx*sin(phi);						//calculate tangential flow velocity
					vt=vr-omega[j]*r;							//calculate total tangential velocity

					/*define beta and chord*/
					if(j<0.5)
					{
						linearInterpolation(inRone,inBone,noblel_one,r,&beta);		//interpolate beta (twist angle)
						linearInterpolation(inRone,inCone,noblel_one,r,&chord);		//interpolate chord
					}

					else
					{
						linearInterpolation(inRtwo,inBtwo,noblel_two,r,&beta);		//interpolate beta (twist angle)
						linearInterpolation(inRtwo,inCtwo,noblel_two,r,&chord);		//interpolate chord
					}

					/*calculate velocity triangles etc.*/
					theta=atan2(u,vt)-(1+rotDir[j])*pi/2;					//calculate inflow angle
					w=sqrt(u*u+(vt*vt));							//calculate inflow velocity
					pdyn=rho*w*w/2.0;							//calculate dynamic pressure
					alpha=-rotDir[j]*(180.0*theta/pi-beta);					//calculate AoA in deg
					//N_one_own=5.156*exp(-(tsr[j]-5.03)*(tsr[j]-5.03)/2.119936)+8.844;
					f1=prandtl(r,fabs(theta),N_one,Rmax[j]);				//calculate Prandtl correction
//					f1=prandtl(r,fabs(theta),Nbl[j],Rmax[j]);				//calculate Prandtl correction
					Re=w*chord*rho/(1.831 * 0.00001);					//compute Reynolds number

					/*determine Cl, Cd on the basis of profile type and Re*/
					for (n=0;n<noblelb;n++)
					{
						if((j<0.5)&&(inRbone[n]<=r)&&(r<=inRbone[n+1])) profileNo=5*(inPbone[n]+inPbone[n+1]);		//compute profile number; multiply by 5 to have integers; if r equal to entry of inRb, the value may fall into interpolation region, but will be estimated properly due to interpolation formula (one of the values will be multiplied by 0)
						if((j>0.5)&&(inRbtwo[n]<=r)&&(r<=inRbtwo[n+1])) profileNo=5*(inPbtwo[n]+inPbtwo[n+1]);
					}
					ClCdinterpolator(j,Re);

					/*determine forces, powers etc.*/
					dforcen=f1*Nbl[j]*pdyn*(Cl*cos(theta)-Cd*fabs(sin(theta)))*chord*A[2]/(2.*pi*r);		//calculate normal (axial) aerodynamic force
					dforcet=f1*Nbl[j]*pdyn*(Cl*fabs(sin(theta))+Cd*cos(theta))*chord*A[2]/(2.*pi*r); //*rotDir[j];	//calculate tangential (rotational) aerodynamic force
					dpower=dforcet*r*omega[j];
					TKE=(3*(w*w*Cl*Cl)/8)*(Nbl[j]*chord/(2.*pi*r));							//calculate local TKE

					/*determine source terms*/
					if (z>=zmin[j])
					{
						if (a==0) source= dforcet*sin(phi)/V;								//calculate x component of actuator force from (tangential) aerodynamic force
						if (a==1) source=-dforcet*cos(phi)/V;								//calculate y component of actuator force from (tangential) aerodynamic  force
						if (a==2) source=dforcen/V;									//calculate z component of actuator force from (normal) aerodynamic force
						if (a==3) source= MFR*TKE/V;									//https://www.sharcnet.ca/Software/Fluent6/html/ug/node217.htm
						if (a==0 && PRINCIPAL_FACE_P(f,C_FACE_THREAD(c,t,i))) power[j]=power[j]+dpower;			//power suming, note that power is set to 0 at end of each iteration in power_Cp
						if (a==2 && PRINCIPAL_FACE_P(f,C_FACE_THREAD(c,t,i))) thrust[j]=thrust[j]+dforcen;		//thrust suming, note that thrust is set to 0 at end of each iteration in power_Cp

						/*BOTH ROTORS USE THE SAMR FILE NOW!!! Print something, if desired (unquote variable *fp at declaration section!)*/
						/*if(j>0.5 && a==0)
						{
							fp = fopen ("file.txt", "a");
							fprintf(fp,"r = %f, phi = %f, b = %f, c = %f, u = %f, vx = %f, vy = %f, vr = %f, vt = %f, theta = %f, w = %f, alpha = %f, Re = %f, Cl = %f, Cd = %f, dP = %f, dT = %f\n",r,phi,beta,chord,u,vx,vy,vr,vt,theta,w,alpha,Re,Cl,Cd,dpower,dforcen);
							fclose (fp);
						}*/
					}

					C_UDMI(c,t,0)=Re;
					C_UDMI(c,t,1)=alpha;
					C_UDMI(c,t,2)=Cl;
					C_UDMI(c,t,3)=Cd;
					C_UDMI(c,t,4)=dpower/(rho*pi*Rmax[j]*Rmax[j]*velref*velref*velref/2);
					C_UDMI(c,t,5)=dforcen/(rho*pi*Rmax[j]*Rmax[j]*velref*velref/2);
				}
			}
		}
	}
}

void varset()
{
	velref=RP_Get_Real("velmeasured");
//	tsr[0]=RP_Get_Real("tsrcaseone");
//	tsr[1]=RP_Get_Real("tsrcasetwo");
	rotDir[0]=RP_Get_Real("omegacaseone")/fabs(RP_Get_Real("omegacaseone"));
	rotDir[1]=RP_Get_Real("omegacasetwo")/fabs(RP_Get_Real("omegacasetwo"));
//	omega[0]=rotDir[0]*tsr[0]*velref/Rmax_one;
//	omega[1]=rotDir[1]*tsr[1]*velref/Rmax_two;
	omega[0]=RP_Get_Real("omegacaseone");
//	Message("read value OMEGA1 = %f \n", omega[0]);
//	Message("read raw value OMEGA1 = %f \n", RP_Get_Real("omegacaseone"));
	omega[1]=RP_Get_Real("omegacasetwo");
//	Message("czytanie OMEGA2 = %f \n", omega[1]);
//	Message("read raw value OMEGA2 = %f \n", RP_Get_Real("omegacasetwo"));
//	omega[0]=140.0;
//	omega[1]=145.0;
	Message("Omega1 = %f;  Omega2 = %f;  Velref = %f; \n", omega[0], omega[1], velref);
	tsr[0]=omega[0]*Rmax_one*rotDir[0]/velref;
	tsr[1]=omega[1]*Rmax_two*rotDir[1]/velref;
	host_to_node_double_1(velref);
	host_to_node_double(omega,Nrotors);
	host_to_node_double(tsr,Nrotors);
	host_to_node_double(rotDir,Nrotors);
}

void ClCdinterpolator(int jot, double NRe)
{
	int m;
	//int noDataTwoTemp = 1;
	double coeffReOne[noDataOne],coeffReTwo[1];
	double coeffClOne[noDataOne],coeffClTwo[1];
	double coeffCdOne[noDataOne],coeffCdTwo[1];

	if(profileNo==10*polarData[noDataOne-1][0])
	{
		for (m=0;m<noDataOne;m++)
		{
			aitken_0(polar[m][0],polar[m][1],polarData[m][2],alpha,&Cl);
			aitken_0(polar[m][0],polar[m][2],polarData[m][2],alpha,&Cd);
			coeffReOne[m]=polarData[m][1];
			coeffClOne[m]=Cl;
			coeffCdOne[m]=Cd;
		}

		linearInterpolation(coeffReOne,coeffClOne,m,NRe,&Cl);
		linearInterpolation(coeffReOne,coeffCdOne,m,NRe,&Cd);
	}

	else if(profileNo==10*polarData[noDataOne+noDataTwo-1][0])
	{
		for (m=0;m<noDataTwo;m++)
		{
			aitken_0(polar[m+noDataOne][0],polar[m+noDataOne][1],polarData[m+noDataOne][2],alpha,&Cl);
			aitken_0(polar[m+noDataOne][0],polar[m+noDataOne][2],polarData[m+noDataOne][2],alpha,&Cd);
			coeffReTwo[m]=polarData[m+noDataOne][1];
			coeffClTwo[m]=Cl;
			coeffCdTwo[m]=Cd;
		}

		linearInterpolation(coeffReTwo,coeffClTwo,m,NRe,&Cl);
		linearInterpolation(coeffReTwo,coeffCdTwo,m,NRe,&Cd);
	}

	else
	{
		for (m=0;m<noDataOne;m++)
		{
			aitken_0(polar[m][0],polar[m][1],polarData[m][2],alpha,&Cl);
			aitken_0(polar[m][0],polar[m][2],polarData[m][2],alpha,&Cd);
			coeffReOne[m]=polarData[m][1];
			coeffClOne[m]=Cl;
			coeffCdOne[m]=Cd;
		}

		for (m=0;m<noDataTwo;m++)
		{
			aitken_0(polar[m+noDataOne][0],polar[m+noDataOne][1],polarData[m+noDataOne][2],alpha,&Clbis);
			aitken_0(polar[m+noDataOne][0],polar[m+noDataOne][2],polarData[m+noDataOne][2],alpha,&Cdbis);
			coeffReTwo[m]=polarData[m+noDataOne][1];
			coeffClTwo[m]=Clbis;
			coeffCdTwo[m]=Cdbis;
		}

		linearInterpolation(coeffReOne,coeffClOne,m,NRe,&Cl);
		linearInterpolation(coeffReOne,coeffCdOne,m,NRe,&Cd);
		linearInterpolation(coeffReTwo,coeffClTwo,m,NRe,&Clbis);
		linearInterpolation(coeffReTwo,coeffCdTwo,m,NRe,&Cdbis);

		if(jot<0.5)
		{
			Cl=(Cl*(Rtransup_one-r)+Clbis*(r-Rtransdown_one))/(Rtransup_one-Rtransdown_one);
			Cd=(Cd*(Rtransup_one-r)+Cdbis*(r-Rtransdown_one))/(Rtransup_one-Rtransdown_one);
		}

		else
		{
			Cl=(Cl*(Rtransup_two-r)+Clbis*(r-Rtransdown_two))/(Rtransup_two-Rtransdown_two);
			Cd=(Cd*(Rtransup_two-r)+Cdbis*(r-Rtransdown_two))/(Rtransup_two-Rtransdown_two);
		}
	}
}

/*prepare data for Aitken interpolation*/
void aitken_0(double xin[], double yin[], int n, double x, double *y)

/* xin: pointer, x array
   yin: pointer, y array
   n:   number of elements in array
   x:   point of interpolation
   y:   pointer, interpolated value */
{    
    int istart;											//begining of interpolation interval
	int ipos=-1;										//position of station for interpolation
	int i=0,j=0;										//loops
	double xint[intwind],yint[intwind];					//table with radii for interpolation, table with values to be interpolated
        
	/*ensure that x falls into interpolation interval*/
	if(x<xin[0]) x=xin[0];
	if(x>xin[n-1]) x=xin[n-1];
	
	/*look for position of x in the xin table, store in ipos*/
	for(;i<n-1;i++)
	{
		if((xin[i]<=x)&&(x<=xin[i+1])) ipos=i;
	}

	/*set the first variable for interpolation by ensuring it doesn't extend the data interval*/
	if      (ipos<=(intwind-1)/2)     istart=0;
	else if (ipos>=(n-(intwind-1)/2)) istart=n-intwind;
	else                              istart=ipos-(intwind-1)/2;

	/*create tables for interpolation*/
	for(;j<intwind;j++)
	{
		xint[j]=xin[istart+j];
		yint[j]=yin[istart+j];
	}

	/*call Aitken interpolation*/
	aitken(xint, yint, intwind, x, y, &delta);
	/*	xint:    array with radii
		fint:    array with chords/twist angles
		intwind: interpolation window size
		x:	     point of interpolation
		y:       interpolated value
		delta:   precision epsilon  */
}

/*Aitken interpolation, Copyright (c) Tao Pang 1997*/
void aitken(double xi[], double fi[], int n, double x, double fa[], double dfa[])
/* Subroutine for the interpolation with the Aitken method.
Copyright (c) Tao Pang 1997. */

#define NMAX 5  /* The number of rows in the interpolation/extrapola- */
                /* tion table.                                        */

{
int i, j;
double x1, x2, f1, f2;
double ft[NMAX];

for (i = 0; i <= (n-1); ++i)
{
ft[i] = fi[i];
}
for (i = 0; i <= (n-2); ++i)
{
   for (j = 0; j <= (n-i-2); ++j)
      {
        x1 = xi[j];
        x2 = xi[j+i+1];
        f1 = ft[j];
        f2 = ft[j+1];
        ft[j] = (x-x1)/(x2-x1)*f2+(x-x2)/(x1-x2)*f1;
      }
}
 fa[0] = ft[0];
 dfa[0] = (fabs(fa[0]-f1)+fabs(fa[0]-f2))/2.0;
}

/* linear interpolation*/
void linearInterpolation(double xs[],double ys[],int count, double x, double y[])
{
	int i;
	double dx, dy;

	/*check extrema*/
	if (x < xs[0])		y[0]=ys[0]; /* return minimum element */
	if (x > xs[count-1])	y[0]=ys[count-1]; /* return maximum */

	/* find i, such that xs[i] <= x < xs[i+1] */
	for (i = 0; i < count-1; i++)
	{
		if (xs[i+1] > x) break;
        }

	/* interpolate */
	dx = xs[i+1]-xs[i];	
	dy = ys[i+1]-ys[i];
	y[0]=ys[i]+(x-xs[i])*dy/dx;
}

/*Prandtl correction calculation*/
double prandtl(double r, double t, double N, double Rmax)
{  
   double func,fun;										//Prandtl correction factors
   fun=N*(Rmax-r)/(2.0*r*sin(t));						//!!! 1/sint(theta)
   func=2.0*acos(exp(-fun))/pi;								//Prandtl parameter
   if(func<0.) func=0.0;									//ensure Prandtl parameter is not negative
   if(func>1.) func=1.0;									//ensure Prandtl parameter is not greater than 1
   return func;
}

/*print the contents of table*/
void tablePrinter(double x[],double y[],double z[])
{
	Message("Function output %f,%f,%f\n",x[0],x[1],x[2]);
	Message("Function output %f,%f,%f\n",y[0],y[1],y[2]);
	Message("Function output %f,%f,%f\n",z[0],z[1],z[2]);
}

/*store polar data in tables*/
void polarDataSetter(int a, int d, int e, int f, double *g, double *h, double *i)
{
	polarData[a][0] = d;
	polarData[a][1] = e;
	polarData[a][2] = f;
	polar[a][0] = g;
	polar[a][1] = h;
	polar[a][2] = i;
}





