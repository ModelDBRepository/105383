/*     ----------------------------------------------------

         NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE

         Copyright 2004, The Johns Hopkins University
            School of Medicine. All rights reserved.
			For research use only; commercial use prohibited.
			Distribution without permission of Raimond L. Winslow
			not permitted. rwinslow@bme.jhu.edu

         Name of Program: Guinea Pig C++ (Coupled) (GPC_Coupled)
         Version: Documented Version, version 1.0.1
         Date: February 2004, August 2004

       -----------------------------------------------------  
*/

//Experiment class file

#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>  
#include <iomanip>
#include <float.h>
#include <gpc.h>
#include <states_index.h>
 
//Constructor for Experiment
Experiment::Experiment()
{
	//physical constants
	RT_over_F=(8.314*310.0)/96.5;
	success=true;
	notdone=true;

	/* algebraic membrane potential method**/
	extra_charge_myo=0;
	extra_charge_SS=0;
	extra_charge_NSR=0;
	extra_charge_JSR=0;
	extra_charge_MITO=0;

	extra_q_myo=0;
	extra_q_NSR=0;
	extra_q_JSR=0;
	extra_q_SS=0;
	extra_q_MITO=0;

}

// Setup variables for Experiment
void Experiment::Setup(double t,double step)
{
//	int k;

	//general constants
	FFlag=true;
	h=step;
	start_time=t;
	tstep=t+step;

	success=true;
	notdone=true;

	//weights for integrators
	errweight[index_V] = 1E-2;			//0
	errweight[index_mNa] = 1.0;			//1
	errweight[index_hNa] = 1.0;			//2
	errweight[index_jNa] = 1.0;			//3
	errweight[index_xKs] = 1.0;			//4 
	errweight[index_Nai] = 0.2;			//5
	errweight[index_Ki] = 1.0/132.0;	//6
	errweight[index_Cai] = 1000.0;		//7
	errweight[index_CaNSR] = 0.5;		//8
	errweight[index_CaSS] = 1000.0;		//9
	errweight[index_CaJSR] = 0.05;		//10
	errweight[index_C1_RyR] = 1.0;		//11
	//errweight[index_O1_RyR]	= 1.0;		//12
	errweight[index_O2_RyR] = 1.0;		//12
	errweight[index_C2_RyR] = 1.0;		//13
	errweight[index_C0] = 1.0;			//14
	errweight[index_C1] = 1.0;			//15
	errweight[index_C2] = 1.0;			//16
	errweight[index_C3]	= 1.0;			//17
	errweight[index_C4] = 1.0;			//18
	errweight[index_Open] = 1.0;		//19
	errweight[index_CCa0] = 1.0;		//20
	errweight[index_CCa1] = 1.0;		//21
	errweight[index_CCa2] = 1.0;		//22
	errweight[index_CCa3] = 1.0;		//23
	errweight[index_CCa4] = 1.0;		//24
	errweight[index_OCa] = 1.0;			//25
	errweight[index_yCa] = 1.0;			//26
	errweight[index_LTRPNCa] = 1.0;		//27
	errweight[index_HTRPNCa] = 1.0;		//28
	errweight[index_N0] = 1.0;			//29
	errweight[index_N1] = 1.0;			//30
	errweight[index_P0] = 1.0;			//31
	errweight[index_P1] = 1.0;			//32
	errweight[index_P2] = 1.0;			//33
	errweight[index_P3] = 1.0;			//34
/**after ATP mod**/
	errweight[index_ATPi] = 1.0/8.0;	//35
/**after Mitochondria mod**/
	errweight[index_Cam]= 1.0;			//not sure
	errweight[index_ADPm]= 1.0/10.0;
	errweight[index_Dpsi]= 1.0/200.0;
	errweight[index_NADH]= 1.0/15.0; 
	errweight[index_ISOC]= 1.0;
	errweight[index_AKG]= 1.0;
	errweight[index_SCoA]= 1.0;
	errweight[index_Succ]= 1.0;
	errweight[index_FUM]=1.0;
	errweight[index_MAL]= 1.0;
	errweight[index_Oaa]= 1.0;

}

//sets the current clamp
void Experiment::setCurrent(double start_time)
{
	if (IF_mode){
		//The following code produce the I-F experiment from Rice paper (JT)
		Istim=0;
		
		within_first_pulse = start_time>=time_on_Is1 && start_time<=time_off_Is1;
		within_second_pulse = start_time>=time_on_Is2 && start_time<=time_off_Is2;

		
		if (stimulusFlag) {
				if (within_first_pulse || within_second_pulse)
					Istim=Istim+pulse_amplitude;
		}
		else{
			cerr << "Current Stimulus must be on during I-F experiment!" << endl;
			exit(1);
		}
	}
	else if (BB_mode){
		Istim=0;
		if (start_time-shift>=t1 && start_time-shift<=t2){
			period=(1.0/high_freq)*1000;	//(ms)
		}
		else{ 
			period=(1.0/norm_freq)*1000;	//(ms)
		}
		time_on_Is1=floor((start_time-shift)/period)*period;
		time_off_Is1=time_on_Is1+pulse_duration;

		if (stimulusFlag) {
			if (((start_time-shift)>=time_on_Is1)&&((start_time-shift)<=time_off_Is1))
				Istim=Istim+pulse_amplitude;
		}
		else{
			cerr << "Current Stimulus must be on during BB experiment!" << endl;
			exit(1);
		}
	}
	else{
		//	time_on_Is1=static_cast<int>((start_time-shift)/period)*period;
		time_on_Is1=floor((start_time-shift)/period)*period;
		time_off_Is1=time_on_Is1+pulse_duration;
		Istim=0;

		if (stimulusFlag) {
			if(((start_time-shift)>=time_on_Is1)&&((start_time-shift)<=time_off_Is1))
				Istim=Istim+pulse_amplitude;
		}
	}
}


/*        *************MAIN COMPUTE***************

	This routine calls integrators that do the actual calculations

*/
void Experiment::compute(ofstream& outDataFile,ofstream& outCurrentFile,
				         ofstream& outDerivFile,double states[],double timenow)
{
	// save initial conditions every ssiniperiod ms
	if ((ssiniFlag)&&(fmod(timenow,ssiniperiod)<stepsize))
		create_ini_file(states,statesize,timenow);

#if USE_CVODE==1
	CVode_run(h, states);  // CVode
	update(states);  //assigns local variables with values from the state array	
	printstates(outDataFile,outCurrentFile,outDerivFile,states,timenow);  //Prints the data to an outfile
#else	
	int i;

	//assigns the inital values into the array y0
	for(i=0;i<statesize;i++)
		y0[i]=states[i];

	RK4(h, states);  //Runge Kutta
	update(states);  //assigns local variables with values from the state array	
	printstates(outDataFile,outCurrentFile,outDerivFile,states,timenow);  //Prints the data to an outfile

	for(i=0;i<statesize;i++)    //updating the values in states to their final value y1
		states[i]=y1[i]; 
#endif
}

//Function F, computes values for the differential equations
//and the result is an array F or F1 of the values.  
void Experiment::getF(double start_time,double states[],bool FFlag)
{
	update(states);
	setCurrent(start_time);
	VclampMode(start_time);
	getReversalPotentials();
	getINa();
	getIKs();
	getIK1();
	getINab();
	getIKp();
	getICa();
	getICaK();
	getINaK();
	getINaCa();
	getICab();
	getIpCa();
	getInsCa();
/*after ATP mod**/
	
	getV_AM();
/**after Mitochondria mod**/
	getATPm();
	getDmuH();
	getNAD();
	getVCS();
	getVACO();
	getVIDH();
	getVKGDH();
	getVSL();
	getVSDH();
	getVFH();
	getVMDH();
	getVAAT();
	getVNO_VHNe();
	getVFO_VHFe();
	getVATPase_Vhu();
	getVANT_Vhleak();
//	getVuni();				//moved up before getFNai
//	getVnaCa();

//	getVATP_XB();


	getFNa();
	getFxKs();
	computeInCalFlux();
	get_Force();
	getF_trpmyo();
	getFLTRPNCa();
	getFHTRPNCa();
	computeJtrpn();

	getVuni();			// moved up
	getVnaCa();			// from below (in the Mitene folder)

	CHF();
	getFNai();
	getFKi();
	getFCai();
	getFCaSS();
	getFCaJSR();
	getFCaNSR();
	getFV();
	getFRyR();
	getFCaL();
	getFyCa();
	getFOCa();
/**after ATP mod**/
	getFATPi();
	getF_mitene();


	/**derivatives**/
	if (membranepot_flag==1)
		F[index_V]=0;	//algebraic expression
	else
		F[index_V]=dV;
	
			
	F[index_mNa]=dmNa;		
	F[index_hNa]=dhNa;		
	F[index_jNa]=djNa;		
	F[index_Nai]=dNai;		
	F[index_Ki]=dKi;		
	F[index_Cai]=dCai;
	F[index_CaNSR]=dCaNSR;
	F[index_CaSS]=dCaSS;
	F[index_CaJSR]=dCaJSR;
	F[index_C1_RyR]=dC1_RyR;
	//F[index_O1_RyR]=dO1_RyR;
	F[index_O2_RyR]=dO2_RyR;
	F[index_C2_RyR]=dC2_RyR;
	F[index_xKs]=dxKs;
	F[index_C0]=dC0;
	F[index_C1]=dC1;
	F[index_C2]=dC2;
	F[index_C3]=dC3;
	F[index_C4]=dC4;
	F[index_Open]=dOpen;
	F[index_CCa0]=dCCa0;
	F[index_CCa1]=dCCa1;
	F[index_CCa2]=dCCa2;
	F[index_CCa3]=dCCa3;
	F[index_CCa4]=dCCa4;
	F[index_yCa]=dyCa;
	F[index_OCa]=dOCa;
	F[index_LTRPNCa]=dLTRPNCa;
	F[index_HTRPNCa]=dHTRPNCa;
	F[index_N0] = dN0;
	F[index_P0] = dP0;
	F[index_P1] = dP1;
	F[index_P2] = dP2;
	F[index_P3] = dP3;
	F[index_N1] = dN1;

/**after ATP mod**/
	F[index_ATPi] = dATPi;

/**after Mitochondria mod**/
	F[index_Cam]=dCam;
	F[index_ADPm]=dADPm;
	F[index_Dpsi]=dDpsi;
	F[index_NADH]=dNADH;
	F[index_ISOC]=dISOC;
	F[index_AKG]=dAKG;
	F[index_SCoA]=dSCoA;
	F[index_Succ]=dSucc;
	F[index_FUM]=dFUM;
	F[index_MAL]=dMAL;
	F[index_Oaa]=dOaa;


	// Not really needed with CVode
#if USE_CVODE==0
	if(FFlag==true)
	{
		for(int i=0;i<statesize;i++)
			F1[i]=F[i];
	}
#endif

}

//updates the local state variables to the current values located in the states array 
void Experiment::update(const double states[])
{
	mNa=states[index_mNa];
	hNa=states[index_hNa];
	jNa=states[index_jNa];
	Nai=states[index_Nai];
	Ki=states[index_Ki];
	Cai=states[index_Cai];
	CaNSR=states[index_CaNSR];
	CaSS=states[index_CaSS];
	CaJSR=states[index_CaJSR];

	C1_RyR=states[index_C1_RyR];
	//O1_RyR=states[index_O1_RyR];
	O2_RyR=states[index_O2_RyR];
	C2_RyR=states[index_C2_RyR];
	O1_RyR=1.0-C1_RyR-C2_RyR-O2_RyR;

	xKs=states[index_xKs];

	C0=states[index_C0];
	C1=states[index_C1];
	C2=states[index_C2];
	C3=states[index_C3];
	C4=states[index_C4];
	Open=states[index_Open];
	CCa0=states[index_CCa0];
	CCa1=states[index_CCa1];
	CCa2=states[index_CCa2];
	CCa3=states[index_CCa3];
	CCa4=states[index_CCa4];

	OCa=states[index_OCa];
	yCa=states[index_yCa];
	LTRPNCa=states[index_LTRPNCa];
	HTRPNCa=states[index_HTRPNCa];
	N0=states[index_N0];
	P0=states[index_P0];
	P1=states[index_P1];
	P2=states[index_P2];
	P3=states[index_P3];
	N1=states[index_N1];

/**after ATP mod**/
	ATPi=states[index_ATPi];

/**after Mitochondria mod**/
	Cam=states[index_Cam];
	ADPm=states[index_ADPm];
	Dpsi=states[index_Dpsi];
	NADH=states[index_NADH];
	ISOC=states[index_ISOC];
	AKG=states[index_AKG];
	SCoA=states[index_SCoA];
	Succ=states[index_Succ];
	FUM=states[index_FUM];
	MAL=states[index_MAL];
	Oaa=states[index_Oaa];

	if (membranepot_flag==1)
		CalcMembranePotential();
	else
		V=states[index_V];

}

//******************************************** Start Adding code here

//*************************************************************** CHF
void Experiment::CHF()
{
	if (chf_flag==1)
	{
	    IK1    = chfsc_IK1*IK1;
	    Jup    = chfsc_Jup*Jup;
	    INaCa  = chfsc_INaCa*INaCa;
	}
}

//************************************************************* CVODE
/* 
 * Integrator.C
 *
 * This file implements/calls integrators for solving ordinary differential equations
 *
 */

#if USE_CVODE==1
/*     ----------------------------------------------------

						code to call CVode starts here

       ----------------------------------------------------            */

#include <sundialstypes.h>   /* definitions of types realtype and             */
                             /* integertype, and the constant FALSE           */
#include <cvode.h>           /* prototypes for CVodeMalloc, CVode, and        */
                             /* CVodeFree, constants OPT_SIZE, BDF, NEWTON,   */
                             /* SV, SUCCESS, NST,NFE,NSETUPS, NNI, NCFN, NETF */
#include <cvdense.h>         /* prototype for CVDense, constant DENSE_NJE     */
#include <nvector_serial.h>  /* definitions of type N_Vector and macro        */
                             /* NV_Ith_S, prototypes for N_VNew, N_VFree      */
#include <dense.h>           /* definitions of type DenseMat, macro DENSE_ELEM*/



/*
	CVode specific functions

	CVode is a stiff/non-stiff integrator written in C  (not in C++)
	To work this code needs header files and a cvode library

*/
extern "C" {

static N_Vector cvode_y;
static void *cvode_mem;
static realtype reltol,abstol;
static N_Vector ew;
static double cvode_t;
M_Env MachEnv;

static realtype ropt[OPT_SIZE];
static long int iopt[OPT_SIZE];
}

// Call CVode 
void Experiment::CVode_run(double h,double states[])
{
	int flag=0,i;

	// call CVode to solve ydot(t)=f(y,t) with y(0) given
	flag = CVode(cvode_mem, start_time+h, cvode_y, &cvode_t, NORMAL);

	if (flag != SUCCESS) { 
		cerr<<"CVode failed, flag="<<flag<<"."<<endl; 
	}

	// copy vector y to states which is returned to compute
	for(i=0;i<statesize;i++)
		states[i]=NV_Ith_S(cvode_y,i);

	if (vclamp_flag==1) // set V to correct value
		states[0]=y0[0]; 
}	

// Extra function that acts as a front to CVode
extern "C" void func_f(long N, realtype time, N_Vector y, N_Vector ydot, void *f_data)
{
	Test2->CVode_f(N, time, y, ydot, f_data);
}


// function f(t,y) for CVode, returns derivatives of variables to CVode
void Experiment::CVode_f(int N, double time, N_Vector y, N_Vector ydot, void *f_data)
{
	int z;
	static double states2[MAXSTATES];
	bool flag=true;

	// copy vector y to states2
	for(z=0;z<N;z++)
		states2[z]=NV_Ith_S(y,z);

	// send states2 to getF which uses it to produce derivatives
	getF(time,states2,flag);

	// copy derivatives to vector ydot which is returned to CVode
	for(z=0;z<N;z++) {
		NV_Ith_S(ydot,z)=F[z];
	}
}

// Initializes CVode
void Experiment::CVodeInit(double states[])
{
	int i;
	
	MachEnv = M_EnvInit_Serial((int)statesize); 
	if (MachEnv == NULL) {
		cerr<<"Trouble with MachEnv in CVODE"<<endl;
		exit(-3);
	}

	// Allocate memory for solution vector y(t)
	cvode_y = N_VNew((int)statesize, MachEnv);
	// Allocate memory for solution vector ew(t)
	ew = N_VNew((int)statesize, MachEnv);

	reltol = tolerance_relative; 
	abstol = tolerance_absolute;

	// initialize vector cvode_y with states
	for(i=0;i<statesize;i++)
		NV_Ith_S(cvode_y,i)=states[i];

	// use default values for options
	for(i=0;i<OPT_SIZE;i++) {
		ropt[i]=0.0;
		iopt[i]=0;
	}

	iopt[MXSTEP]=1000000;	//added by JT, taken from integrator.cpp of Canine model

	// except for these
	ropt[HMAX]=step_max;   // Largest step size
	ropt[HMIN]=step_min;  // Smallest step size

	// scale tolerance based on maximum value of state
	for(i=0;i<statesize;i++) {
		NV_Ith_S(ew,i)=abstol/errweight[i];
	}

	/*	CVodeMalloc sets up initial settings for CVode. See 
		Integration method is BDF(Backward Differential Formula)
		Other choice would be ADAMS but it is not as stable */
#if 0
	//	This method of calling CVODE does not pass error weights
	cvode_mem = CVodeMalloc(statesize, func_f, start_time, cvode_y, BDF, NEWTON, SS, &reltol, &abstol,
                          NULL, NULL, TRUE, iopt, ropt, MachEnv);

#else
	// We wish to pass errorweight to CVODE
	cvode_mem = CVodeMalloc(statesize, func_f, start_time, cvode_y, BDF, NEWTON, 
							SV, &reltol, ew, NULL, NULL, TRUE, iopt, ropt, MachEnv);
#endif

	if (cvode_mem == NULL) { 
		cerr<<"CVodeMalloc failed."<<endl;
		exit(1);
	}

	/* CVDense is needed by Newton algorithm for solving linear system
	   The second NULL tells the solver to compute an approximation of the Jacobian. */
	CVDense(cvode_mem, NULL, NULL);
}

void Experiment::CVodeExit(void)
{
	N_VFree(cvode_y);        
	CVodeFree(cvode_mem);
	M_EnvFree_Serial(MachEnv);
}


#else
/*     ----------------------------------------------------

						rk4 code starts here

       ---------------------------------------------------- 
 
 
	Merson Modified Runge-Kutta 4th Order ADAPTIVE Step Algorithm
 	Kubicek, M., Marek, M. (1983). Computational methods in
 	Bifurcation theory and Dissipative Structures. pg. 84.

 	Passed Variables:
 		h			timestep
 		states		Contains initial state at time t on entry.
 					Contains final state at time t + tstep on exit

 	This routine controls the solution of the system of differential
 	equations from time=start_time to time=start_time+h by monitoring 
	the truncation error on each incremental step, and adjusting the 
	step_size based on the error after each attempted step.
*/
void Experiment::RK4(double h,double states[])
{
	int i,z;

	while(notdone)
	{
		if(success)
		{   FFlag=true;//FFlag is used to indicate whether or not to send the derivatives to F1
			getF(start_time,y0,FFlag);
		}
		for(z=0;z<statesize;z++)
		{
			k1[z]=h*F1[z];
			y1[z]=y0[z]+k1[z]/3.0;
			states[z]=y1[z];
		}
		FFlag=false;
		getF(start_time+h/3.0,states,FFlag);
		for(z=0;z<statesize;z++)
		{
			k2[z]=h*F[z];
			y1[z]=y0[z]+(k1[z]+k2[z])/6.0;
			states[z]=y1[z];
		}
			
		getF(start_time+h/3.0,states,FFlag);
		for(z=0;z<statesize;z++)
		{
			k3[z]=h*F[z];
			y1[z]=y0[z]+(k1[z]+3.0*k3[z])/8.0;
			states[z]=y1[z];
		}
		
		getF(start_time+h/2.0,states,FFlag);
		for(z=0;z<statesize;z++)
		{
			k4[z]=h*F[z];
			y4[z]=y0[z]+(.5*k1[z])-(1.5*k3[z])+(2.0*k4[z]);
			states[z]=y4[z];
			y1[z]=y4[z];	
		}
		
		getF(start_time+h,states,FFlag);
		for(z=0;z<statesize;z++)
		{
			k5[z]=h*F[z];
			y1[z]=y0[z]+(k1[z]+4.0*k4[z]+k5[z])/6.0;
			states[z]=y1[z];
		}
		
		// Calculating the error
		/////////////////////////////////////////////////////////////////
		//update(states);
		tr_error = 0;
		for(i=0;i<statesize;i++)
		{
			errtmp = fabs((y4[i]-y1[i])*0.2*errweight[i]);

			if (errtmp<big) 
	 			tr_error = max(errtmp, tr_error);
			else
				tr_error = big;
		}

		if (tr_error< tolerance_absolute) 
		{
			for(i=0;i<statesize;i++) 
		   {
				// If the absolute size of solution is less than 1 thousanth of 
				// error margin, set value to zero
				if (fabs(y1[i])< tolerance_absolute/1000) 
					y1[i] = 0.0;

				states[i]=y1[i];
				y0[i] = y1[i];
			
				
			}

			start_time = start_time+h;
			success = true;
			if (start_time>=tstep)
			{			
				notdone =false;
			}
			else
			{
				h = .85*h*pow((tolerance_absolute/tr_error),.2);
				notdone = true;
			}
		}
		else {
			if (h <= step_min) 
			{
				for(i=0;i<statesize;i++)
				{
					if (fabs(y1[i])<tolerance_absolute/1000) 
					{
						y1[i] = 0.0;
					}
					y0[i] = y1[i];
					states[i]=y1[i];
				}
				start_time = start_time+h;
				success = true;
				if (start_time >= tstep) 
				{
					notdone = false;
				
				}
				else
				{
					h = step_min;
					notdone = true;
				}
			}
			else
			{
				h = .85*h*pow((tolerance_absolute/tr_error),.2);
				success = false;
				notdone = true;
			}
		}
		h = min(step_max,max(step_min,h));
		h = min(h,tstep-start_time);
		if (h < .0000000001) 
		notdone = false;
	
	}
}	
#endif

//******************************************************** MembranPot
//gets the time rate of change for V
void Experiment::getFV()
{	
	if (vclamp_flag==1 || membranepot_flag==1)
	{	
		CalcMembranePotential();
		dV=0;
		y0[0]=V; 
	}
	else
	{
		dV=(-(INa+ICa+ICaK+IKs+IK1+IKp+INaCa+INaK+InsCa+IpCa+ICab+INab+Istim)/C_m);
	} 
}


//sets the voltage clamp
void Experiment::VclampMode(double start_time)
{
	double ramp;
	if (vclamp_flag==1)
	{
		double time_vclamp_on1 = floor((start_time-shift)/period)*period;
		double vclamp_duration = time_vclamp_off - time_vclamp_on;

        if (((start_time-shift) >= time_vclamp_on1+time_vclamp_on) && 
     	      ((start_time-shift) < time_vclamp_on1+time_vclamp_on+vclamp_duration)) 
		{           
			ramp = (((start_time-shift)-time_vclamp_on1-time_vclamp_on)/2.0)
     		      *(vclamp_set-vclamp_hold) + vclamp_hold;
            if (vclamp_hold<=vclamp_set) 
				V = std::min(vclamp_set,ramp); // depol.  steps
            else
                    V = std::max(vclamp_set,ramp); // hyperpol. steps
		}     
        else if ((start_time-shift)<(time_vclamp_on1+time_vclamp_on)) 
		{
            V = vclamp_hold;
		}
		else
		{
			ramp = vclamp_set +((time_vclamp_on1+time_vclamp_on
     			+ vclamp_duration-(start_time-shift))/2.0)*(vclamp_set-vclamp_hold);
            if (vclamp_hold<=vclamp_set) 
                 V = std::max(vclamp_hold,ramp); // depol. step
            else
                 V = std::min(vclamp_hold,ramp); // hyper. step
		}  
	 }
}

/* algebraic membrane potential method */
void Experiment::CalcMembranePotential()
{
	if (membranepot_flag==1) {
			double F_JSR, F_i, Ca_all;
			double Na_all, K_all, extra, Co;		
			double a1;
			
			//Take out the Ca++ buffering in Subspace
			/*a1=CMDNSStot/(CaSS+KmCMDN);
			//there is no EGTA and we need to add CMDNSStot to parameter file
			double F_SS;
			F_SS=1.0+a1;*/
			
			a1=CSQNtot/(CaJSR+KmCSQN);
			F_JSR=1.0+a1;

			a1=CMDNtot/(Cai+KmCMDN);
			F_i=1.0+a1;

			//global
			//Take out the Ca buffereing in subspace
			/*Ca_all=Vmyo*(Cai*F_i+LTRPNCa+HTRPNCa)+VNSR*CaNSR+VSS*(F_SS*CaSS)+VJSR*F_JSR*CaJSR
				+Vmito*Cam;*/
			Ca_all=Vmyo*(Cai*F_i+LTRPNCa+HTRPNCa)+VNSR*CaNSR+VSS*CaSS+VJSR*F_JSR*CaJSR+Vmito*Cam;

			//there is no EGTA
			//note that LTRPNCa is in (mM), unlike the canine version LTRPNCa is %
			
			Na_all=Vmyo*Nai;
			K_all=Vmyo*Ki;
			extra=(extra_charge_myo+extra_q_myo)*Vmyo+(extra_charge_NSR+extra_q_NSR)*VNSR
				+(extra_charge_JSR+extra_q_JSR)*VJSR
				+(extra_charge_SS+extra_q_SS)*VSS
				+(extra_charge_MITO+extra_q_MITO)*Vmito;
			Co=Vtotal*(2*Cao+Ko+Nao);		//volume scaled to match the volume of a cell
/*
			cout << "Nai " << Nai << endl;
			cout << "Ki " << Ki << endl;
			cout << "Na_all " << Na_all << endl;
			cout << "K_all " << K_all << endl;
			cout << "Ca_all " << Ca_all << endl;
			cout << "Co " << Co << endl;
			cout << "extra " << extra << endl;
*/
			V=(Faraday*1000)/(Acap*C_m)*(Na_all+K_all+2*Ca_all-Co+extra);
			y0[0]=V;
//			cout << "calcMembrane: " << V << endl;
	}
}

void Experiment::AlgDiff(double states[])
{
		update(states);	

		double a1=Acap*C_m/(Vmyo*Faraday*1000);
		
		if (vclamp_flag==1){
			states[index_V]=vclamp_hold;
		}
/*
		cout << "stop" << endl;
		//cin >> temp;
		//CalcMembranePotential();
		cout << states[index_V]<< endl;
		cout << V << endl;
		cout << a1 << endl;
		cout << Vtotal << endl;
*/
		extra_charge_myo=(states[index_V]-V)*a1*Vmyo/Vtotal;
		extra_charge_SS+=extra_charge_myo;
		extra_charge_NSR+=extra_charge_myo;
		extra_charge_JSR+=extra_charge_myo;
		extra_charge_MITO+=extra_charge_myo;

		CalcMembranePotential();	//not sure

		cout<<"Difference in concentrations to alg formulation: "<< extra_charge_myo<< " mM" <<endl;
		cout<<"Extra charge used in cytosol: "<<extra_charge_myo<< " (mM)" << endl;
		cout<<"Extra charge used in SS: "<< extra_charge_SS<< " (mM)" << endl;
		cout<<"Extra charge used in NSR: "<< extra_charge_NSR<< " (mM)" << endl;
		cout<<"Extra charge used in JSR: "<< extra_charge_JSR<< " (mM)" << endl;
		cout<<"Extra charge used in MITO: " << extra_charge_MITO<< " (mM) "<< endl;

}

//******************************************************* ReversalPot
// Calculates reversal potentials
void Experiment::getReversalPotentials()
{
	double a1=Ko+0.01833*Nao;
	double a2=Ki+0.01833*Nai;
	
	E_Na=RT_over_F*log(Nao/Nai);
	E_K=RT_over_F*log(Ko/Ki);
	E_Ks=RT_over_F*log(a1/a2);  //IKs is dependent on both Na and K

	if (Cai<1.0e-10)
		Cai=1.0E-10;
	
	E_Ca=0.5*RT_over_F*log(Cao/Cai);
}

//**************************************************************** INa
//getFNa will get the time rate of change of mNa, hNa, and jNa
void Experiment::getFNa()
{
/**paper
	double a1=0.32*(V+47.13);
	double MAlpha=a1/(1.0-exp(-0.1*(V+47.13)));
**/

/**fortran**/
	double a1=0.32*(V+47.13);
	double MAlpha;
	if (V==-47.13)
		MAlpha=3.2;
	else
		MAlpha = a1/(1.0-exp(-0.1*(V+47.13)));			
	// if statements taken out in calculation of MAlpha and dmNa
//

	double MBeta = 0.08*exp(-V/11.0);

/**fortran**/
	if (1.0/(MBeta+MAlpha)<0.03)
		mNa=MAlpha/(MBeta+MAlpha);
//

	dmNa = MAlpha*(1.0-mNa)-MBeta*mNa;
	double HAlpha, HBeta, JAlpha, JBeta;

	if (V<-40)
	{
	    HAlpha = 0.135*exp((80.0+V)/-6.8);
	    HBeta = 3.56*exp(0.079*V)+310000.0*exp(0.35*V);
	    a1 = -127140.0*exp(0.2444*V);
	    double a2 = 3.474E-5*exp(-0.04391*V);
	    double a3 = 1.0+exp(0.311*(V+79.23));
	    JAlpha = (a1-a2)*(V+37.78)/a3;
	    a2 = 1.0+exp(-0.1378*(V+40.14));
	    JBeta = 0.1212*exp(-0.01052*V)/a2;
	}
	else
	{
		HAlpha = 0.0;
	    HBeta = 1.0/(0.13*(1+exp((V+10.66)/-11.1)));
	    JAlpha = 0.0;
	    a1 = 1.0+exp(-0.1*(V+32.0));
	    JBeta = 0.3*exp(-2.535E-7*V)/a1;
	}
	dhNa = HAlpha*(1.0-hNa)-HBeta*hNa;	
	djNa = JAlpha*(1.0-jNa)-JBeta*jNa;	
}


//getINa gets the current value of INa
void Experiment::getINa()
{
	INa=G_Na*(pow(mNa,3.0)*hNa*jNa*(V-E_Na));		
}

//***************************************************** Concentration
//Intercellular Concentration calculations
//These methods get the time rate of change for Nai, Ki, and Cai
void Experiment::getFNai()
{
	double a1=Acap/(Vmyo*Faraday*1000.0);
	/*before Mitochondria mod
	dNai=-( INa+INab+3.0*(INaCa+INaK)+InsNa )*a1;	
	*/
	/**after Mitochondria mod**/
	dNai=-(INa+INab+3.0*(INaCa+INaK)+InsNa)*a1-VnaCa*0.615;
}

void Experiment::getFKi()
{
	double a1=Acap/(Vmyo*Faraday*1000.0);
	dKi = 0;
//	dKi=-(IKs+IK1+IKp+ICaK-2.0*INaK+InsK+Istim)*a1; 
}

void Experiment::getFCai()
{
	double a1=Acap/(2.0*Vmyo*Faraday*1000.0);
/** Commented out by JT, replace with Faraday*1000
	double a1=Acap/(2*Vmyo*Faraday);
**/
	double a3=ICab-2.0*INaCa+IpCa; 
														// scaling factor for volume correction
	/*before Mitochondria mod
	dCai=beta_i*(Jxfer-Jup-Jtrpn-a3*.5*a1);				
	*/
	/**after Mitochondria mod**/
	dCai=beta_i*(Jxfer-Jup-Jtrpn-a3*.5*a1+(-Vuni+VnaCa)*0.615);
}

//************************************************************ TRPNCa
// TROPONIN FRACTION DERIVATIVES
// These methods get the time rate of change for LTRPNCa and HTRPNCa

void Experiment::getFLTRPNCa()
{
/** Fortran **/
	double a1=(kltrpn_minus*LTRPNCa)*(1.0/3.0+2.0/3.0*(1.0-FN_Ca));
	dLTRPNCa=kltrpn_plus*Cai*(LTRPNtot-LTRPNCa)-a1;
	/**/

/**new*
	double a1=kltrpn_minus*(1/3 + (2/3)*(1-Fnorm))*LTRPNCa;
	dLTRPNCa=kltrpn_plus*Cai*(LTRPNtot-LTRPNCa)-a1;
/**JT****/
}

void Experiment::getFHTRPNCa()
{
	double a1=khtrpn_minus*HTRPNCa;
	dHTRPNCa=khtrpn_plus*Cai*(HTRPNtot-HTRPNCa)-a1;
}

//*************************************************************** RyR
/* RyR Receptors */

/* -------COMPUTE DERIVATIVES OF RyR RECEPTOR STATES---------------------------- */


//getFRyR gets the time rate of change for C1_RyR,O2_RyR, C2_RyR, and O1_RyR.
void Experiment::getFRyR()
{
	double a1=pow((CaSS),mcoop);
	double a2=pow((CaSS),ncoop);

	dC1_RyR=-kaplus*a2*C1_RyR+kaminus*O1_RyR;
	dO2_RyR=kbplus*a1*O1_RyR-kbminus*O2_RyR;
	dC2_RyR=kcplus*O1_RyR-kcminus*C2_RyR;
	//dO1_RyR=-(dC1_RyR+dO2_RyR+dC2_RyR);		//reinsert O1_RyR as state variable
}

//************************************************************** ICaL
//L-Type Ca++ current

//L-Type Calcium Channel state transition rates
double C0_to_C1, C1_to_C2, C2_to_C3, C3_to_C4;
double C1_to_C0, C2_to_C1, C3_to_C2, C4_to_C3;
double CCa0_to_CCa1, CCa1_to_CCa2;
double CCa2_to_CCa3, CCa3_to_CCa4;
double CCa1_to_CCa0, CCa2_to_CCa1;
double CCa3_to_CCa2, CCa4_to_CCa3;
double C0_to_CCa0, C1_to_CCa1, C2_to_CCa2;
double C3_to_CCa3, C4_to_CCa4;
double CCa0_to_C0, CCa1_to_C1, CCa2_to_C2;
double CCa3_to_C3, CCa4_to_C4;


//getFCaL gets the time rate of change for C0,C1,C2,C3,C4,Open,CCa0,CCa1,CCa2,CCa3, and CCa4
void Experiment::getFCaL()
{
	double alpha = 0.4*exp((V+12)/10.0);
	double beta = 0.05*exp(-(V+12)/13.0);

	double alpha_prime = aL*alpha;
	double beta_prime = beta/bL;		// change bL to 2.0 in parameters

	C0_to_C1 = 4.0*alpha;
	C1_to_C2 = 3.0*alpha;
	C2_to_C3 = 2.0*alpha;
	C3_to_C4 = alpha;

	CCa0_to_CCa1 = 4.0*alpha_prime;
	CCa1_to_CCa2 = 3.0*alpha_prime;
	CCa2_to_CCa3 = 2.0*alpha_prime;
	CCa3_to_CCa4 = alpha_prime;

	C1_to_C0 =      beta;
	C2_to_C1 = 2.0*beta;
	C3_to_C2 = 3.0*beta;
	C4_to_C3 = 4.0*beta;

	CCa1_to_CCa0 =      beta_prime;
	CCa2_to_CCa1 = 2.0*beta_prime;
	CCa3_to_CCa2 = 3.0*beta_prime;
	CCa4_to_CCa3 = 4.0*beta_prime;
	
	/**paper
		gamma = 0.5625*CaSS;
	**/
	/**fortran code**/
	gamma = 0.140625*CaSS; 

	C0_to_CCa0 = gamma;		
	C1_to_CCa1 = aL*C0_to_CCa0;	// = gamma*aL
	C2_to_CCa2 = aL*C1_to_CCa1;	// = gamma*aL^2
	C3_to_CCa3 = aL*C2_to_CCa2;	// = gamma*aL^3
	C4_to_CCa4 = aL*C3_to_CCa3;	// = gamma*aL^4
		
	CCa0_to_C0 = omega;		// = omega
	CCa1_to_C1 = CCa0_to_C0/bL;	// = omega/bL
	CCa2_to_C2 = CCa1_to_C1/bL;	// = omega/bL^2
	CCa3_to_C3 = CCa2_to_C2/bL;	// = omega/bL^3
	CCa4_to_C4 = CCa3_to_C3/bL;	// = omega/bL^4

	double a1 = (C0_to_C1+C0_to_CCa0)*C0;
	double a2 = C1_to_C0*C1 + CCa0_to_C0*CCa0;
	dC0 = a2 - a1;

	a1 = (C1_to_C0+C1_to_C2+C1_to_CCa1)*C1;
	a2 = C0_to_C1*C0 + C2_to_C1*C2 + CCa1_to_C1*CCa1;
	dC1 = a2 - a1;

	a1 = (C2_to_C1+C2_to_C3+C2_to_CCa2)*C2;
	a2 = C1_to_C2*C1 + C3_to_C2*C3 + CCa2_to_C2*CCa2;
	dC2 = a2 - a1;

	a1 = (C3_to_C2+C3_to_C4+C3_to_CCa3)*C3;
	a2 = C2_to_C3*C2 + C4_to_C3*C4 + CCa3_to_C3*CCa3;
	dC3 = a2 - a1;


	a1 = (C4_to_C3+fL+C4_to_CCa4)*C4;
	a2 = C3_to_C4*C3 + gL*Open + CCa4_to_C4*CCa4;
	dC4 = a2 - a1;

	dOpen =  fL*C4 - gL*Open;		// change fL to 0.3 in parameters

	a1 = (CCa0_to_CCa1+CCa0_to_C0)*CCa0;
	a2 = CCa1_to_CCa0*CCa1 + C0_to_CCa0*C0;
	dCCa0 = a2 - a1;

	a1 = (CCa1_to_CCa0+CCa1_to_CCa2+CCa1_to_C1)*CCa1;
	a2 = CCa0_to_CCa1*CCa0 + CCa2_to_CCa1*CCa2 + C1_to_CCa1*C1;
	dCCa1 = a2 - a1;

	a1 = (CCa2_to_CCa1+CCa2_to_CCa3+CCa2_to_C2)*CCa2;
	a2 = CCa1_to_CCa2*CCa1 + CCa3_to_CCa2*CCa3 + C2_to_CCa2*C2;
	dCCa2 = a2 - a1;

	a1 = (CCa3_to_CCa2+CCa3_to_CCa4+CCa3_to_C3)*CCa3;
	a2 = CCa2_to_CCa3*CCa2 + CCa4_to_CCa3*CCa4 + C3_to_CCa3*C3;
	dCCa3 = a2 - a1;

	// reinsert dCCa4, CCa4 state variable
	a1 = (CCa4_to_CCa3+fprime+CCa4_to_C4)*CCa4;
	a2 = CCa3_to_CCa4*CCa3 + (gprime*OCa) + C4_to_CCa4*C4;		// add fprime and gprime to parameters, add OCa to state vars
	dCCa4 = a2 - a1;

}

void Experiment::getICa()
{
	double VF_over_RT=V/RT_over_F;						//RT_over_F is in (mV)
	double VFsq_over_RT=(1000.0*Faraday)*VF_over_RT;	//VF_over_RT is unitless mV/mV
	double a1=1E-3 *exp(2.0*VF_over_RT)-Cao*.341;
	double a2=exp(2.0*VF_over_RT)-1.0;

// inserted by JT
	if (fabs(V)<LHospitalThreshold) // Use limit V->0, AT
		ICamax=2.0*PCa*1000.0*Faraday*(1.0-340.0*Cao);		//scale by 1000, so current is in uA/uF
	else
		ICamax=PCa*4.0*VFsq_over_RT*(a1/a2);
/** commented out by JT, replace with Limit of V->0 condition statements
	ICamax = PCa*4*VFsq_over_RT*(a1/a2);
JT**/

	/**paper
		ICa=ICamax*yCa*(Open+OCa);
	**/
	/**fortran**/
	ICa=1.50*ICamax*yCa*Open*5.0;	//The 5 factor added to account for low open probability of CAY
									//L-type channel.
	ICa=ICa*(0.5+0.5*0.6);		//modified for calmodulin function
}

//getFyCa gets the time rate of change for yCa
void Experiment::getFyCa()
{
	double yCa_inf=1.0/(1.0+exp((V+55.0)/7.5))+0.5/(1.0+exp((21.0-V)/6.0));
	double tau_yCa=20.0 + 600.0/(1.0+exp((V+30.0)/9.5));

	dyCa = (yCa_inf-yCa)/tau_yCa;
}

//*** New ***
void Experiment::getFOCa()				// add OCa to state variables
{
	dOCa = fprime*CCa4 - gprime*OCa;
}

//******************************************************** Junctional
//These methods get the time rate of change for CaSS, CaJSR, and CaNSR
void Experiment::getFCaSS()
{
	double a2=Acap/(2.0*VSS*Faraday*1000.0);
	double a3=Jrel*VJSR/VSS-Jxfer*Vmyo/VSS;
	dCaSS=beta_SS*(a3-ICa*a2);
}

void Experiment::getFCaJSR()
{
	dCaJSR=beta_JSR*(Jtr-Jrel);
}

void Experiment::getFCaNSR()
{
	dCaNSR=Jup*Vmyo/VNSR-Jtr*VJSR/VNSR;
}

//computes current values for beta_SS, beta_JSR, and beta_i
void Experiment::computeJtrpn()
{
	Jtrpn=dLTRPNCa+dHTRPNCa;
	
	double a1=CMDNtot*KmCMDN/pow((CaSS+KmCMDN),2.0);
	beta_SS=1.0/(1.0+a1);
	
	a1=CSQNtot*KmCSQN/pow((CaJSR+KmCSQN),2.0);
	beta_JSR=(1.0/(1.0+a1));
	
	a1=CMDNtot*KmCMDN/pow((Cai+KmCMDN),2.0);
	beta_i=1.0/(1.0+a1);
}

//computes current values for Jup, Jrel, Jtr, and Jxfer
void Experiment::computeInCalFlux()
{
	double fb=pow((Cai/Kfb),Nfb);
	double rb=pow((CaNSR/Krb),Nrb);
	
/** before ATP mod
	Jup=(vmaxf*fb-vmaxr*rb)/(1+fb+rb);
**/

/**after ATP mod**/
	Jup=KSR*(vmaxf*fb-vmaxr*rb)/(1.0+fb+rb)/((KmATP_SR/ATPi)*(1.0+(8.0-ATPi)/Ki_SR)+(1.0+((8.0-ATPi)/Ki_prime_SR)));
	Jrel = v1*(O1_RyR+O2_RyR)*(CaJSR-CaSS);
	Jtr = (CaNSR - CaJSR)/tautr;
	Jxfer = (CaSS-Cai)/tauxfer;
}

//*************************************************************** IKs

void Experiment::getFxKs()
{
	double xKsAlpha=7.19E-5*(V+30.0)/(1.0-exp(-0.148*(V+30.0)));
	double xKsBeta=1.31E-4*(V+30.0)/(exp(0.0687*(V+30.0))-1.0);

	dxKs = xKsAlpha*(1.0-xKs)-xKsBeta*xKs;
}

//getIKs gets the current value of IKs
void Experiment::getIKs()
{
	/** paper 
		double G_Ks=0.1128*sqrt(Ko/5.4);	//delete G_Ks from parameters
	**/
	/**fortran**/
	double G_Ks=0.282*sqrt(Ko/5.4);

	double Xi=1.0/(1+exp((V-40.0)/40.0));
	IKs=G_Ks*Xi*(pow(xKs,2.0))*(V-E_Ks);
}

//************************************************************** INab
//Background Na+ current
void Experiment::getINab()
{
	INab=G_Nab*(V-E_Na);
}

//************************************************************** ICab
void Experiment::getICab()
{
	ICab = G_Cab*(V-E_Ca);		
}

//*************************************************************** IK1
//IK1 current
void Experiment::getIK1()
{
	double G_K1=0.75*sqrt(Ko/5.4);		
	double K1Alpha=1.02/(1+exp(0.2385*(V-E_K-59.215)));
	double K1Beta=(0.4912*exp(0.08032*(V-E_K+5.476))+exp(0.06175*(V-E_K-594.31)))/(1+exp(-0.5143*(V-E_K+4.753)));
	double K1_inf=K1Alpha/(K1Alpha+K1Beta);

//	IK1=Fcik1*G_K1*K1_inf*(V-E_K);
	IK1=G_K1*K1_inf*(V-E_K);

}

//*************************************************************** IKp
//IKp Current--plateau current
void Experiment::getIKp()
{
	double KpV=1.0/(1.0+exp((7.488-V)/5.98));
	IKp=G_Kp*KpV*(V-E_K);
}

//************************************************************* INaCa
//The sodium-calcium exchanger
void Experiment::getINaCa()
{
	double VF_over_RT=V/RT_over_F;
	double a1 = exp(eta*VF_over_RT)*pow(Nai,3.0)*Cao;
	double a2 = exp((eta-1.0)*VF_over_RT)*pow(Nao,3.0)*Cai;
	double a3 = 1.0+ksat*exp((eta-1.0)*VF_over_RT);
	double a4 = KmCa+Cao;
	double a5 = 1.0/(pow(KmNa,3.0)+pow(Nao,3.0));

	INaCa = kNaCa*a5*(a1-a2)/(a4*a3);
}

//************************************************************** IpCa
//Sarcolemmal pump current
void Experiment::getIpCa()
{
	/**before ATP mod
		IpCa = IpCamax*Cai/(KmpCa+Cai);		
		// value of IpCamax is different in parameters
	**/
	/**after ATP mod**/
	IpCa = IpCamax*Cai/(KmpCa+Cai)*(1.0/(1.0+(Km1ATP_CaP/ATPi)*(1.0+(8.0-ATPi)/KiADP_CaP))+1.0/(1.0+(Km2ATP_CaP/ATPi)));

}

//************************************************************** INaK
//Sodium ATPase pump
void Experiment::getINaK()
{
	double VF_over_RT=V/RT_over_F;
	double sigma = (exp(Nao/67.3)-1.0)/7.0;
	double a1 = 1.0+0.1245*exp(-0.1*VF_over_RT);
	double a2 = 0.0365*sigma*exp(-1*VF_over_RT);
	double fNaK = 1.0/(a1+a2);
	a1 = Ko/(Ko+KmKo);
    a2 = 1.0+pow((KmNai/Nai),1.5);
	/**before ATP mod
		INaK = INaKmax*fNaK*(a1/a2);		
		// value of INaKmax is different in parameters
	**/
	/**after ATP mod**/
	INaK = INaKmax*fNaK*(a1/a2)/(1.0+(Km1AT_NaK/ATPi)*(1.0+(8.0-ATPi)/Ki1AD_NaK));
}

//************************************************************* InsCa
// InsCa added to Guinea pig project (JT)
void Experiment::getInsCa()
{	
	double VF_over_RT=V/RT_over_F;					
	double VFsq_over_RT=(1000.0*Faraday)*VF_over_RT;

	//use Faraday*1000 instead of Faraday (JT), so current is in uA/uF
	double InsNa_bar = PnsNa*VFsq_over_RT*(0.75*Nai*exp(V/RT_over_F)-0.75*Nao)/(exp(V/RT_over_F)-1.0);
	InsNa = InsNa_bar/(1.0+pow(KmnsCa/Cai,3.0));

	double InsK_bar = PnsK*VFsq_over_RT*(0.75*Ki*exp(V/RT_over_F)-0.75*Ko)/(exp(V/RT_over_F)-1.0);
	InsK = InsK_bar/(1+pow(KmnsCa/Cai,3.0));

	InsCa = InsNa + InsK;
}

//**************************************************** Tropomyosin_XB
//tropomyo.cpp added to project (JT)
void Experiment::getF_trpmyo()
{
	
	/** Model on Paper
	kTrop_np = kTrop_pn * pow(((LTRPNCa/LTRPNtot)/Ktrop_half), Ntrop);	//corrected version

	dN0 = kTrop_pn*P0 - kTrop_np*N0 + g_01*N1;
	dP0 = -(kTrop_pn + f_01)*P0 + kTrop_np*N0 + g_01*P1;
	dP1 = -(kTrop_pn + f_12 + g_01)*P1 + kTrop_np*N1 + f_01*P0 + g_12*P2;
	dP2 = -(f_23 + g_12)*P2 + f_12*P1 + g_23*P3;
	dP3 = -(g_23*P3) + f_23*P2;
	dN1 = kTrop_pn*P1 - (kTrop_np + g_01)*N1;
	**/
	/** Real model in Fortran**/
	kTrop_np = kTrop_pn * pow(((LTRPNCa/LTRPNtot)/Ktrop_half), Ntrop);

	//dN0 = kTrop_pn*P0 - kTrop_np*N0 + g_01_mod*N1;		
	//Rice did not include this in fortran model
	dP0 = -(kTrop_pn + f_01)*P0 + kTrop_np*N0 + g_01_mod*P1;
	dP1 = -(kTrop_pn + f_12 + g_01_mod)*P1 + kTrop_np*N1 + f_01*P0 + g_12_mod*P2;
	dP2 = -(f_23 + g_12_mod)*P2 + f_12*P1 + g_23_mod*P3;
	dP3 = -(g_23_mod*P3) + f_23*P2;
	dN1 = kTrop_pn*P1 - (kTrop_np + g_01_off_mod)*N1;
	dN0 = -dP0-dP1-dP2-dP3-dN1;							//added by JT
}

//************************************************************* Force
void Experiment::get_Force()
{
	FN_Ca = alpha_SL * (P1 + N1 + P2 + P3)/fnormmax2;
	Fnorm = alpha_SL * ((P1 + N1 + 2*P2+3*P3)/3.0)/(fnormmax);
	Force = zeta * Fnorm;
}

//************************************************************** ICaK
//L-Type Ca channel permeable to K+
void Experiment::getICaK()
{
	double VF_over_RT=V/RT_over_F;
	/**paper
		double VFsq_over_RT=(1000.0*Faraday)*VF_over_RT;
		double PKprime = PK/(1.0+ICamax/ICahalf);
		double a1 = Ki*exp(2.0*VF_over_RT)-Ko;
   		double a2 = exp(2.0*VF_over_RT)-1.0; // singular
		ICaK = PKprime*(Open+OCa)*yCa*VFsq_over_RT*(a1/a2); //modified
		//OCa added to keep consistency, this eq is now the same as 
		//what is in the Rice 2000 paper (JT)
	**/
	/**fortran**/
	double PKprime = PK/(1.0+std::min(ICamax,0.0)/ICahalf)*4.0;
   	double a1 = Ki*exp(VF_over_RT)-Ko;
   	double a2 = exp(VF_over_RT)-1.0; // singular
	ICaK = PKprime*yCa*Open*(a1/a2);
}

//************************************************************** V_AM
void Experiment::getV_AM()
{
/**before ATP mod
double V_XB = P0*g_01 + P1*g_12 + P2*g_23;
VATP_XB = No*V_XB;
**/
/**after ATP mod**/
V_AM=V_AM_scaler*V_AM_max*((f_01*P0+f_12*P1+f_23*P2)/(f_01+f_12+f_23))*(1.0/(1.0+(KmATP_AM/ATPi)*(1.0+(8.0-ATPi)/Ki_AM)));
}


void Experiment::getFATPi()
{
//dATPi = VANT*0.615 - Jup - VATP_XB - INaK - IpCa;
// double Vant_on=1.0;
// if (ATPi>7.99 | ATPi<0.01)
//         Vant_on=0.0;
dATPi = 0.615*VANT-V_AM-Jup-(6.371e-5*(INaK+IpCa));

/* if (ATPi+dATPi>8.0 | ATPi+dATPi<0.0){
         cout << "VANT: " << VANT << endl;
         cout << "V_AM: " << V_AM << endl;
         cout << "Jup: " << Jup << endl;
         cout << "INaK: " << 6.371e-5*INaK << endl;
         cout << "IpCa: " << 6.371e-5*IpCa << endl;
         cout << "ATPi: " <<  ATPi << endl;
         cout << "ADPm: " <<  ADPm << endl;
         cout << "ATPm: " << ATPm << endl;
         cout << "Dpsi: " << Dpsi << endl;
         cout << "Cm: " << Cm << endl;
         while (true){
         }
}
*/
}

//*******************************************************************
//***************************** Mitene ******************************
//*******************************************************************

//******************************************************* getF_mitene
void Experiment::getF_mitene()
{
	//equation (12)
	dCam = fm * (Vuni - VnaCa);
	//equation (1)
	dADPm = VANT - VATPase - VSL;
	//equation (2)
	dDpsi = -(-VHNe - VHFe + Vhu + VANT + Vhleak + 2.0*b*VnaCa + 2.0*Vuni) / Cmito;
	//equation (3)
	dNADH = -VNO + VIDH + VKGDH + VMDH;
	//equation (4)
	dISOC = VACO - VIDH;
	//equation (5)
	dAKG = VIDH + VAAT - VKGDH;
	//equation (6)
	dSCoA = VKGDH - VSL;
	//equation (7)
	dSucc = VSL - VSDH;
	//equation (8)
	dFUM = VSDH - VFH;
	//equation (9)
	dMAL = VFH - VMDH;
	//equation (10)
	dOaa = VMDH - VCS - VAAT;
}

//************************************************************** ATPm
void Experiment::getATPm()
{
	//atp = pme(13) - v(2);
	ATPm=Cm-ADPm;
	//mM = mM - mM
}

//************************************************************** Dmuh
void Experiment::getDmuH()
{
	//change in proton motive force (mV)
	//dmu = -2.303*pme(4)*pme(5)*pme(3)/pme(6) + v(3);

	DmuH = -2.303*RToverF*DpH+Dpsi;
	//mV = (Unitls)*mV*(Unitls)+mV
}

//*************************************************************** NAD
void Experiment::getNAD()
{
	//CPN --> total sum of mito pyridine nucleotides
	NAD = CPN-NADH;
	//mM = mM-mM	
}

//*************************************************************** VCS
void Experiment::getVCS()
{
	//vc = pme(31)*pme(32)/(1 +pme(29)/pme(24) +pme(33)/v(11)+ pme(29)/pme(24) *pme(33)/v(11));

	VCS=KCS*EtCS/(1.0+KmAcCoA/AcCoA+KmOaa/Oaa+KmAcCoA/AcCoA*KmOaa/Oaa);
	//mM/ms=(1/ms)*(mM)/(1+ mM/mM + mM/mM + mM/mM*mM/mM)

}

//************************************************************** VACO
void Experiment::getVACO()
{
	//Aconitase
	//ci = pme(28)-v(5)-v(6)-v(7)-v(8)-v(9)-v(10)-v(11);
	//vac = pme(34)*(CIT(v,pme) - v(5)/pme(35));

	CIT=CIK-ISOC-AKG-SCoA-Succ-FUM-MAL-Oaa;
	//mM

	VACO=kfACO*(CIT-ISOC/KACOeq);
	//mM/ms = 1/(ms) * (mM-mM/const)
}

//************************************************************** VIDH
void Experiment::getVIDH()
{
	//a = (1+v(2)/pme(37))*(1+v(1)/pme(43));
	//i = (1+v(4)/pme(36));
	//vid = pme(38)*pme(39)/(((1 + pme(40)/8.1*10^(-5) + 5.98*10^(-5)/pme(40)) + (pme(41)/v(5))^pme(44)/Fa(v,pme))+    ...
    //    (pme(42)/NAD(v,pme))*Fi(v,pme) +((pme(41)/v(5))^pme(44)*pme(42)/NAD(v,pme))*Fi(v,pme)/Fa(v,pme));


//	NAD();

	Fa=(1+ADPm/KADP)*(1+Cam/KaCa);
	// unitless

	Fi=(1+NADH/KidhNADH);
	// unitless

	//equation (16)
	//VIDH=kIDH*EtID/(((1.0+H/kh_1+kh_2/H)+pow((Kmiso/ISOC),nID)/Fa)+(KmIDNAD/NAD)*Fi+(pow((Kmiso/ISOC),nID)*(KmIDNAD/NAD)*Fi/Fa));
	VIDH=kIDH*EtID/(((1.0+H/kh_1+kh_2/H)+pow(Kmiso/ISOC,nID)/Fa)+(KmIDNAD/NAD)*Fi+(pow(Kmiso/ISOC,nID)*KmIDNAD/NAD)*Fi/Fa);
	//mM/ms = mM/ms

	/*
	kh_1=8.1*10^-5
	kh_2=5.98*10^-5
	*/
}

//************************************************************* VKGDH
void Experiment::getVKGDH()
{
	//vkgd = pme(45)*pme(46)/(1+(pme(47)/v(6))/((1+pme(52)/pme(49))*(1+v(1)/pme(50)))+(pme(48)/NAD(v,pme))^pme(51)/((1+pme(52)   ...
    //    /pme(49))*(1+v(1)/pme(50))));


	//NAD();
	//equation (17) 
	//There is a typo on the paper, the matlab version is correct
	//VKGDH=kKGDH*EtKG/(1.0+pow((KmKG/AKG),nKG)/((1.0+Mg/Kmg)*(1.0+Cam/Kca))+(KmKGNAD/NAD)/((1.0+Mg/Kmg)*(1.0+Cam/Kca)));
	VKGDH = kKGDH*EtKG/(1.0+(KmKG/AKG)/((1.0+Mg/Kmg)*(1.0+Cam/Kca))+pow((KmKGNAD/NAD),nKG)/((1.0+Mg/Kmg)*(1.0+Cam/Kca)));
	//mM/ms = mM/ms

}

//*************************************************************** VSL
void Experiment::getVSL()
{
	//vs = pme(53)*(v(7)*v(2)-v(8)*ATPm(v,pme)*pme(55)/pme(54));


	//ATPm

	//equation (18)
	VSL=kfSL*(SCoA*ADPm-Succ*ATPm*CoA/KSLeq);
	//mM/ms=mM/ms

}

//************************************************************** VSDH
void Experiment::getVSDH()
{
	//vsd = pme(56)*pme(57)/(1+pme(58)/v(8)*(1+v(9)/pme(59))*(1+v(11)/pme(60)));

	//equation (19)
	VSDH=kSDH*EtSDH/(1.0+KmSucc/Succ*(1.0+FUM/KiFUM)*(1.0+Oaa/KiOxaa));
	//mV/ms=(1/ms)*mM

}

//*************************************************************** VFH
void Experiment::getVFH()
{
	//vf = pme(61)*(v(9)-v(10)/pme(62));

	VFH=kfFH*(FUM-MAL/KFHeq);
	//mM/ms=(1/ms)*(mM-mM/const)
}

//************************************************************** VMDH
void Experiment::getVMDH()
{
	//fi = (1/(1+pme(63)/pme(40)+pme(63)*pme(64)/(pme(40)^2)))^2;
	//fa = 1/(1+pme(40)/pme(65)+(pme(40)^2)/(pme(65)*pme(66)))+pme(67);
	//f = Fhi(pme)*Fha(pme);
	//vmd = pme(68)*Fh(pme)*pme(69)/(1+pme(70)/v(10)*(1+v(11)/pme(71))+pme(72)/NAD(v,pme)+pme(70)/v(10)*(1+v(11)/pme(71))*pme(72)/NAD(v,pme));


	//NAD();
	
	//equation (23)
	Fhi=pow((1.0/(1.0+Kh3/H+Kh3*Kh4/pow(H,2))),2);
	//Fhi=pow((1.0/(1.0+Kh3/H+Kh3*Kh4/pow(H,2))),2);

	//equation (22)
	Fha=1.0/(1.0+H/Kh1+pow(H,2)/(Kh1*Kh2))+Koff;
	//Fha=1.0/(1.0+H/Kh1+pow(H,2)/(Kh1*Kh2))+Koff;
	
	Fh=Fhi*Fha;
	//unitless

	//VMDH=kMDH*Fh*EtMD/(1.0+Kmal/MAL*(1.0+Oaa/Kioaa)+KmmNAD/NAD+Kmal/MAL*(1.0+Oaa/Kioaa)*KmmNAD/NAD);
	VMDH=kMDH*Fh*EtMD/(1.0+Kmal/MAL*(1.0+Oaa/Kioaa)+KmmNAD/NAD+Kmal/MAL*(1.0+Oaa/Kioaa)*KmmNAD/NAD);
	//mM/ms=mM/ms
}

//************************************************************** VAAT
void Experiment::getVAAT()
{
	//vasp = pme(75)*v(12);
	//vaa = pme(73)*(pme(26)*v(11) - v(12)*v(6)/pme(74));
	//matlab version:
	//vaa = pme(73)*pme(26)*v(11)*pme(75)*pme(74)/(pme(75)*pme(74) + v(6)*pme(73));

	/**JT
	VAAT=kfAAT*(GLU*Oaa-ASP*AKG/KAATeq);
	**/

	VAAT=kfAAT*GLU*Oaa*kcnsASP*KAATeq/(kcnsASP*KAATeq+AKG*kfAAT);
	//mM = (1/(ms*mM)*(mM)*(mM)*(1/ms)*const)/((1/ms)*const+mM*(1/(ms*mM)))

}

//********************************************************** VNO_VHNe
void Experiment::getVNO_VHNe()
{
	//arn = (pme(4)*pme(5)/pme(6))*log((1.35*10^18)*sqrt(v(4)/(NAD(v,pme))));

	//vn = 30*pme(1)*((6.394*10^(-10)+2.656*10^(-19)*exp(6*pme(6)*pme(7)/pme(4)/pme(5))) *exp(pme(6)*AREN(v,pme)/pme(4)/pme(5))-6.394*   ...
    //    10^(-10)*exp(6*0.85*pme(6)*DmuH(v,pme)/pme(4)/pme(5))+8.632*10^(-27)*exp(pme(6)*AREN(v,pme)/pme(4)/pme(5)) *exp(6*0.85*       ...
    //    pme(6)*DmuH(v,pme)/pme(4)/pme(5))) / (((1+2.077*10^(-18)*exp(pme(6)*AREN(v,pme)/pme(4)/pme(5))) *exp(6*pme(6)*pme(7)/         ...
    //    pme(4)/pme(5))+(1.728*10^(-9)+1.059*10^(-26)*exp(pme(6)*AREN(v,pme)/pme(4)/pme(5))) *exp(6*0.85*pme(6)*DmuH(v,pme)/pme(4)/pme(5))));

	//vhn = 360*pme(1)*(6.394*10^(-10)*exp(pme(6)*AREN(v,pme)/pme(4)/pme(5))-(6.394*10^(-10)+1.762*10^(-13))*exp(6*0.85*pme(6)*           ...
    //   DmuH(v,pme)/pme(5)/pme(4)))/(((1+2.077*10^(-18)*exp(pme(6)*AREN(v,pme)/pme(4)/pme(5))) *exp(6*pme(6)*pme(7)/pme(5)/pme(4))+    ...
    //   (1.728*10^(-9)+1.059*10^(-26)*exp(pme(6)*AREN(v,pme)/pme(4)/pme(5))) *exp(6*0.85*pme(6)*DmuH(v,pme)/pme(4)/pme(5))));




	//NAD();

	//equation (28)
	AREN=RToverF*log(kres*sqrt(NADH/NAD));
	
	//DmuH();

	VNO=0.5*rhoREN*((ra+rc1*exp(6.0*Dpsio/RToverF))*exp(AREN/RToverF)-ra*exp(6.0*g*DmuH/RToverF)+rc2*exp(AREN/RToverF)*exp(6.0*g*DmuH/RToverF))/((1.0+r1*exp(AREN/RToverF))*exp(6.0*Dpsio/RToverF)+(r2+r3*exp(AREN/RToverF))*exp(6.0*g*DmuH/RToverF));
	//VNO= 30*rhoREN*((ra+rc1*exp(6*Dpsio/RToverF))*exp(AREN/RToverF)-ra*exp(6*g*DmuH/RToverF)+rc2*exp(AREN/RToverF)*exp(6*g*DmuH/RToverF))/((1+r1*exp(AREN/RToverF))*exp(6*Dpsio/RToverF)+(r2+r3*exp(AREN/RToverF))*exp(6*g*DmuH/RToverF));
	

	VHNe=6.0*rhoREN*(ra*exp(AREN/RToverF)-(ra+rb)*exp(6.0*g*DmuH/RToverF))/((1.0+r1*exp(AREN/RToverF))*exp(6.0*Dpsio/RToverF)+(r2+r3*exp(AREN/RToverF))*exp(6.0*g*DmuH/RToverF));
	//VHNe=360*rhoREN*(ra*exp(AREN/RToverF)-(ra+rb)*exp(6*g*DmuH/RToverF))/((1+r1*exp(AREN/RToverF))*exp(6*Dpsio/RToverF)+(r2+r3*exp(AREN/RToverF))*exp(6*g*DmuH/RToverF));

	/*
	ra=6.394*10^(-10)
	rb=1.762*10^-13
	rc1=2.656*10^(-19)
	g=0.85
	rc2=8.632*10^(-27)
	r1=2.077*10^(-18)
	r2=1.728*10^(-9)
	r3=1.059*10^(-26)
	kres=1.35*10^18
	*/
}

//********************************************************** VSDH
void Experiment::getVFO_VHFe()
{
	//arf = (pme(4)*pme(5)/pme(6))*log((5.497*10^13)*sqrt(pme(10)/(FAD(pme))));

	//vf = 30*pme(2)*((6.394*10^(-10)+2.656*10^(-19)*exp(6*pme(6)*pme(7)/pme(4)/pme(5))) *exp(pme(6)*AREF(pme)/pme(4)/pme(5))-6.394*   ...
    //    10^(-10)*exp(6*0.85*pme(6)*DmuH(v,pme)/pme(4)/pme(5))+8.632*10^(-27)*exp(pme(6)*AREF(pme)/pme(4)/pme(5)) *exp(6*0.85*       ...
    //    pme(6)*DmuH(v,pme)/pme(4)/pme(5))) / (((1+2.077*10^(-18)*exp(pme(6)*AREF(pme)/pme(4)/pme(5))) *exp(6*pme(6)*pme(7)/         ...
    //    pme(4)/pme(5))+(1.728*10^(-9)+1.059*10^(-26)*exp(pme(6)*AREF(pme)/pme(4)/pme(5))) *exp(6*0.85*pme(6)*DmuH(v,pme)/pme(4)/    ...
    //    pme(5))));

	//vhf = 240*pme(2)*(6.394*10^(-10)*exp(pme(6)*AREF(pme)/pme(4)/pme(5))-(6.394*10^(-10)+1.762*10^(-13))*exp(6*0.85*pme(6)*          ...
    //   DmuH(v,pme)/pme(5)/pme(4)))/(((1+2.077*10^(-18)*exp(pme(6)*AREF(pme)/pme(4)/pme(5))) *exp(6*pme(6)*pme(7)/pme(5)/pme(4))+    ...
    //   (1.728*10^(-9)+1.059*10^(-26)*exp(pme(6)*AREF(pme)/pme(4)/pme(5))) *exp(6*0.85*pme(6)*DmuH(v,pme)/pme(4)/pme(5))));



	//equation (29)
	AREF=RToverF*log(kresf*sqrt(FADH2/FAD));		
	//AREF=RToverF*log(kresf*sqrt(FDAH2/FAD));

	//DmuH();
	
	//equation (26)
	VFO=0.5*rhoREF*((ra+rc1*exp(6.0*Dpsio/RToverF))*exp(AREF/RToverF)-ra*exp(6.0*g*DmuH/RToverF)+rc2*exp(AREF/RToverF)*exp(6.0*g*DmuH/RToverF))/((1.0+r1*exp(AREF/RToverF))*exp(6.0*Dpsio/RToverF)+(r2+r3*exp(AREF/RToverF))*exp(6.0*g*DmuH/RToverF));
	//VFO=0.5*rhoREF*((ra+rc1*exp(6.0*Dpsio/RToverF))*exp(AREF/RToverF)-ra*exp(6.0*g*DmuH/RToverF)+rc2*exp(AREF/RToverF)*exp(6.0*g*DmuH/RToverF))/((1.0+r1*exp(AREF/RToverF))*exp(6.0*Dpsio/RToverF)+(r2+r3*exp(AREF/RToverF))*exp(6.0*g*DmuH/RToverF));


	//equation (27) (same equations from VNO_VHNE, but 4 instead of 6)
	VHFe=4.0*rhoREF*(ra*exp(AREF/RToverF)-(ra+rb)*exp(6.0*g*DmuH/RToverF))/((1.0+r1*exp(AREF/RToverF))*exp(6.0*Dpsio/RToverF)+(r2+r3*exp(AREF/RToverF))*exp(6.0*g*DmuH/RToverF));
	//VHFe=240*rhoREF*(ra*exp(AREF/RToverF)-(ra+rb)*exp(6.0*g*DmuH/RToverF))/((1.0+r1*exp(AREF/RToverF))*exp(6.0*Dpsio/RToverF)+(r2+r3*exp(AREF/RToverF))*exp(6.0*g*DmuH/RToverF));



	/*
	kresf=5.497*10^13	??? slightly different
	ra=6.394*10^(-10)
	rb=1.762*10^(-13)
	rc1=2.656*10^(-19)
	g=0.85
	rc2=8.632*10^(-27)
	r1=2.077*10^(-18)
	r2=1.728*10^(-9)
	r3=1.059*10^(-26)
	*/



}

//******************************************************* VATPase_Vhu
void Experiment::getVATPase_Vhu()
{
	//af = pme(4)*pme(5)/pme(6)*log(1.71*10^6*ATPm(v,pme)/(v(2)*pme(11)));
	
	//vatp = -60*pme(12)*((1.656*10^(-3)+9.651*10^(-14) *exp(3*pme(6)*pme(7)/pme(4)/pme(5))) *exp(pme(6)*AF1(v,pme)/pme(4)/  ...
    //    pme(5))-1.656*10^(-5)*exp(3*pme(6)*DmuH(v,pme)/pme(4)/pme(5))+4.845*10^(-14) *exp(pme(6)*AF1(v,pme)/pme(4)/pme(5))...
    //    *exp(3*pme(6)*DmuH(v,pme)/pme(4)/pme(5))) / (((1+1.346*10^(-8)*exp(pme(6)*AF1(v,pme)/pme(4)/pme(5))) *exp(3*pme(6)...
    //    *pme(7)/pme(4)/pme(5))+(7.739*10^(-7)+6.65*10^(-15)*exp(pme(6)*AF1(v,pme)/pme(4)/pme(5))) *exp(3*pme(6)*           ...
    //    DmuH(v,pme)/pme(4)/pme(5))));

	//vh = -180*pme(12)*(1.656*10^(-3)*(1+exp(pme(6)*AF1(v,pme)/pme(4)/pme(5)))-(1.656*10^(-5)+3.373*10^(-7))*exp(3*pme(6)*  ...
    //    DmuH(v,pme)/pme(4)/pme(5))) / (((1+1.346*10^(-8)*exp(pme(6)*AF1(v,pme)/pme(4)/pme(5)))*exp(3*pme(6)*pme(7)/pme(4) ...
    //    /pme(5)) +(7.739*10^(-7)+6.65*10^(-15)*exp(pme(6)*AF1(v,pme)/pme(4)/pme(5)))*exp(3*pme(6)*DmuH(v,pme)/pme(4)/pme(5))));




	AF1=RToverF*log(kf1*ATPm/(ADPm*Pi));
	//AF1=RToverF*log(kf1*ATPm/(ADPm*Pi));
	//mV

	//DmuH();

	//VATPase=-rhoF1*((100*pa+pc1*exp(3*Dpsio/RToverF))*exp(AF1/RToverF)+(-pa*exp(3*DmuH/RToverF)+pc2*exp(AF1/RToverF)*exp(3*DmuH/RToverF)))/((1+p1*exp(AF1/RToverF))*exp(3*Dpsio/RToverF)+(p2+p3*exp(AF1/RToverF))*exp(3*DmuH/RToverF));
	//matlab:
	VATPase=-rhoF1*((100.0*pa+pc1*exp(3.0*Dpsio/RToverF))*exp(AF1/RToverF)-pa*exp(3.0*DmuH/RToverF)+pc2*exp(AF1/RToverF)*exp(3.0*DmuH/RToverF))/(((1.0+p1*exp(AF1/RToverF))*exp(3.0*Dpsio/RToverF)+(p2+p3*exp(AF1/RToverF))*exp(3.0*DmuH/RToverF)));
	//VATPase=-60*rhoF1*((100.0*pa+pc1*exp(3.0*Dpsio/RToverF))*exp(AF1/RToverF)-pa*exp(3.0*DmuH/RToverF)+pc2*exp(AF1/RToverF)*exp(3.0*DmuH/RToverF))/(((1.0+p1*exp(AF1/RToverF))*exp(3.0*Dpsio/RToverF)+(p2+p3*exp(AF1/RToverF))*exp(3.0*DmuH/RToverF)));
	//is pa 1.656*10^-5 or -3?

	//mM/ms = mM


	//Vhu=-3*rhoF1*(100*pa*(1+exp(AF1/RToverF))-(pa+pb)*exp(3*DmuH/RToverF))/((1+p1*exp(AF1/RToverF))*exp(3*Dpsio/RToverF)+(p2+p3*exp(AF1/RToverF))*exp(3*DmuH/RToverF));
	//matlab:
	Vhu = -3*rhoF1*(100.0*pa*(1.0+exp(AF1/RToverF))-(pa+pb)*exp(3.0*DmuH/RToverF))/(((1.0+p1*exp(AF1/RToverF))*exp(3.0*Dpsio/RToverF) +(p2+p3*exp(AF1/RToverF))*exp(3.0*DmuH/RToverF)));
	//Vhu = -180*rhoF1*(100.0*pa*(1.0+exp(AF1/RToverF))-(pa+pb)*exp(3.0*DmuH/RToverF))/(((1.0+p1*exp(AF1/RToverF))*exp(3.0*Dpsio/RToverF) +(p2+p3*exp(AF1/RToverF))*exp(3.0*DmuH/RToverF)));

	//same questions as above

	/*
	pa=1.656*10^(-3)		for VATPase ** unit?? also there is some discrepancies for the value
	pb=3.373*10^-7
	pc1=9.651*10^(-14) 
	pc2=4.845*10^(-14)		??? slight difference
	p1=1.346*10^(-8)
	p2=7.739*10^(-7)
	p3=6.65*10^(-15)
	kf1=1.71*10^6

	pa=1.656*10^-5			for Vhu
	*/

}

//******************************************************* VANT_Vhleak
void Experiment::getVANT_Vhleak()
{
//vad_t = pme(14)*(1-(0.05*pme(16)*0.45*0.8*v(2))/(0.45*pme(15)*0.05*ATPm(v,pme))*exp(-pme(6)*v(3)/pme(4)/pme(5)))/
//         ((1+0.05*pme(16)/(0.45*pme(15))*exp(-0.5*pme(6)*v(3)/pme(4)/pme(5)))*(1+0.45*0.8*v(2)/(0.05*ATPm(v,pme))));

//vleak = pme(22)*20.83*DmuH(v,pme);

//ATPm();

//the following expression is different from the paper. (35, extra exp. term)
VANT=VmDT*(0.75-(0.25*ATPi*0.45*1.0*ADPm)/(0.225*(8.0-ATPi)*0.025*ATPm)*exp(-Dpsi/RToverF))/((1.0+0.25*ATPi/(0.225*(8.0-ATPi))*exp(-hm*Dpsi/RToverF))*(1.0+0.45*1.0*ADPm/(0.025*ATPm)));
// original VANT=VmDT*(1.0-(0.05*ATPi*0.45*0.8*ADPm)/(0.45*(8.0-ATPi)*0.05*ATPm)*exp(-Dpsi/RToverF))/((1.0+0.05*ATPi/(0.45*(8.0-ATPi))*exp(-hm*Dpsi/RToverF))*(1.0+0.45*0.8*ADPm/(0.05*ATPm)));
//VANT=VmDT*(1.0-(0.05*ATPi*0.45*0.8*ADPm)/(0.45*(8.0-ATPi)*0.05*ATPm)*exp(-Dpsi/RToverF))/((1.0+0.05*ATPi/(0.45*(8.0-ATPi))*exp(-hm*Dpsi/RToverF))*(1.0+0.45*0.8*ADPm/(0.05*ATPm)));
//vad_t = pme(14)*(1-(0.05*pme(16)*0.45*0.8*v(2))/(0.45*pme(15)*0.05*ATPm(v,pme))*exp(-pme(6)*v(3)/pme(4)/pme(5)))/((1+0.05*pme(16)/(0.45*pme(15))*exp(-0.5*pme(6)*v(3)/pme(4)/pme(5)))*(1+0.45*0.8*v(2)/(0.05*ATPm(v,pme))));
//mM/ms= mM/ms

//DmuH();
Vhleak=gh*DmuH;
//mM/ms= mM/(ms*mV) * mV


//h=0.5

}

//************************************************************** Vuni
void Experiment::getVuni()
{
	//vun = pme(19)*(pme(20)/0.019)*(1 + pme(20)/0.019)^3/((1+pme(20)/0.019)^4 + 110.0/(1 + pme(20)/0.00038)^2.8)*    ...
    //    (2*pme(6)*(v(3)-0.091)/pme(4)/pme(5)/(1-exp(-2*pme(6)*(v(3)-0.091)/pme(4)/pme(5))));

	//equation (38)
	Vuni=Vmuni*(Cai/ktrans)*pow((1+Cai/ktrans),3)/(pow((1+Cai/ktrans),4)+L/pow((1+Cai/kact),na))*(2*(Dpsi-91.0)/RToverF/(1-exp(-2*(Dpsi-91.0)/RToverF)));
	//Vuni=Vmuni*(Cai/ktrans)*pow((1+Cai/ktrans),3)/(pow((1+Cai/ktrans),4)+L/pow((1+Cai/kact),na))*(2*(Dpsi-91.0)/RToverF/(1-exp(-2*(Dpsi-91.0)/RToverF)));


	/*
	ktrans=0.019
	L=110.0
	kact=0.00038
	na=2.8
	*/

}

//************************************************************* VnaCa
void Experiment::getVnaCa()
{
	//vnac = pme(17)*exp(pme(18)*pme(6)/pme(4)/pme(5)*(v(3)-0.091))*v(1)*exp(-log(pme(20)/v(1)))/((1+9.4/10)^(2+2*pme(18))*(v(1)+0.00375));

	//equation (39)
	//VnaCa=VmNC*exp(b*(1/RToverF)*(Dpsi-91.0))*exp(-log(Cai/Cam))/(pow((1+Kna/Nai),n)*(1+Knca/Cam)); //not sure about the negative sign in front of log!
	//matlab version:
	//VnaCa=VmNC*exp(b*(1/RToverF)*(Dpsi-91.0))*Cam*exp(-log(Cai/Cam))/(pow((1+Kna/Nai),n))*(Cam+Knca));
	VnaCa=VmNC*exp(b*(1/RToverF)*(Dpsi-91.0))*exp(-log(Cai/Cam))/(pow((1+Kna/Nai),n)*(1+Knca/Cam));
	//ask sonia
	//mM/ms=mM/ms


	/*
	Kna=9.4
	Kca=3.75*10^-3
	*/

}




