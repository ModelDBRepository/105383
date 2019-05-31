
/*     ----------------------------------------------------

         NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE

         Copyright 2004, The Johns Hopkins University
            School of Medicine. All rights reserved.
			For research use only; commercial use prohibited.
			Distribution without permission of Raimond L. Winslow
			not permitted. rwinslow@bme.jhu.edu

         Name of Program: Guinea Pig C++: Coupled, Algbraic, BB, MCA
         Version: Documented Version, version 1.0.1
         Date: August 2004

       -----------------------------------------------------  
*/

#include "integratorcvode.h"

//Define CVode F function
//Note: direct access of y, ydot is not constant
void func_f(long N, realtype time, N_Vector y, N_Vector ydot, void *f_data) { 
	((Model*)f_data)->F(time, N_VGetData(ydot), N_VGetData(y)); 
}

IntegratorCVode::IntegratorCVode(void) {
	ropt = new double[OPT_SIZE];
	iopt = new long[OPT_SIZE];
	usingErrorWeights = true;
}

IntegratorCVode::~IntegratorCVode(void) {
	N_VFree(cvode_y);   
	N_VFree(cvode_ic);   
	//N_VDispose(cvode_y);
	N_VFree(ew);  
	CVodeFree(cvode_mem);
	delete[] ropt;
	delete[] iopt;
}

//Integrates the current problem, attempting to reach time
//State variables and actual time recorded in cvode_y, &cvode_t.
void IntegratorCVode::iterateToTime(double time) {
	// call CVode to solve ydot(t)=f(y,t) with y(0) given
	//Normal Mode
	int flag = CVode(cvode_mem, time, cvode_y, &cvode_t, NORMAL);
	//*
	//One step mode: ONE_STEP
//	int flag;
//	while(cvode_t < time) {
//		flag = CVode(cvode_mem, time, cvode_y, &cvode_t, ONE_STEP);
		if (flag != SUCCESS) { 
			cout << "CVode failed, flag=" << flag << "." << endl; 
//			char c;
//			cin >> c;
		}
//	}
}

//Initializes 
void IntegratorCVode::setupCVode(Model *model) {
	//Parameters: step_max, step_min, usingErrorWeights, start_time
	N_Vector errweight;

	machEnv = NULL;
	machEnv = M_EnvInit_Serial(model->getProblemSize()); 
	if (machEnv == NULL) {
		cerr << "Trouble with MachEnv in CVODE" << endl;
		exit(-3);
	}

	// Allocate memory for solution vector y(t)
	cvode_y = N_VNew( model->getProblemSize(), machEnv);
//	cvode_y = N_VMake( model->getProblemSize(), model->getStates(), machEnv);

	// Allocate memory for solution vector ew(t)
	ew = N_VNew( model->getProblemSize(), machEnv);
	errweight = N_VMake( model->getProblemSize(), model->getErrorWeights(), machEnv);
	// scale tolerance based on maximum value of state
	N_VInv(errweight, ew);
	N_VScale(abstol, ew, ew);

	//Allocate and set cvode_ic
	cvode_ic = N_VNew( model->getProblemSize(), machEnv);
	N_Vector ic = N_VMake( model->getProblemSize(), model->getStates(), machEnv);
	N_VAddConst(ic, 0.0, cvode_ic);
	N_VDispose(ic);

	// use default values for options
	for(int i = 0; i < OPT_SIZE; i++) {
		ropt[i]=0.0;
		iopt[i]=0L;
	}
	//Set optional parameters
	iopt[MXSTEP] = 1000000;	//added by JT, taken from integrator.cpp of Canine model
	ropt[HMAX] = step_max;   // Largest step size
	ropt[HMIN] = step_min;  // Smallest step size


	/*	CVodeMalloc sets up initial settings for CVode. See 
		Integration method is BDF(Backward Differential Formula)
		Other choice would be ADAMS but it is not as stable */
	if (usingErrorWeights) {
		// We wish to pass errorweight to CVODE
		cvode_mem = CVodeMalloc(model->getProblemSize(), func_f, 0, cvode_ic, BDF,
			NEWTON, SV, &reltol, ew     , model, NULL, TRUE, iopt, ropt, machEnv);
	} else {
		//	This method of calling CVODE does not pass error weights
		cvode_mem = CVodeMalloc(model->getProblemSize(), func_f, 0, cvode_ic, BDF,
			NEWTON, SS, &reltol, &abstol, model, NULL, TRUE, iopt, ropt, machEnv);
	}

	if ( cvode_mem == NULL ) { 
		cerr << "CVodeMalloc failed." << endl;
		exit(1);
	}

	/* CVDense is needed by Newton algorithm for solving linear system
	   The second NULL tells the solver to compute an approximation of the Jacobian. */
	CVDense(cvode_mem, NULL, NULL);

	N_VDispose(errweight);
}
//Refresh the CVode integrator for another run
void IntegratorCVode::refreshCVode(Model *model)  {
	//Parameters: step_max, step_min, usingErrorWeights, start_time
	int error;

	N_Vector ic = N_VMake( model->getProblemSize(), model->getStates(), machEnv);
	N_VAddConst(ic, 0.0, cvode_ic);
	N_VDispose(ic);

	/*	CVodeMalloc sets up initial settings for CVode. See 
		Integration method is BDF(Backward Differential Formula)
		Other choice would be ADAMS but it is not as stable */
	if (usingErrorWeights) {
		// We wish to pass errorweight to CVODE
		error = CVReInit(cvode_mem, func_f, 0, cvode_ic, BDF,
			NEWTON, SV, &reltol, ew     , model, NULL, TRUE, iopt, ropt, machEnv);
	} else {
		//	This method of calling CVODE does not pass error weights
		error = CVReInit(cvode_mem, func_f, 0, cvode_ic, BDF,
			NEWTON, SS, &reltol, &abstol, model, NULL, TRUE, iopt, ropt, machEnv);
	}
	if (SUCCESS != error) {
		cout << "Error in IntegratorCVode::refreshCVode(), code :" << error << "." << endl;
		exit(-1);
	}
	/* CVDense is needed by Newton algorithm for solving linear system
	   The second NULL tells the solver to compute an approximation of the Jacobian. */
	CVDense(cvode_mem, NULL, NULL);
}
//Set all the parameters for the integrator
void IntegratorCVode::setParameters(fileIO& data) {
	reltol		= data["reltol"];
	abstol		= data["abstol"];
	step_max	= data["step_max"];
	step_min	= data["step_min"];
	start_time	= data["start_time"];
	stepsize	= data["stepsize"];
	cvode_t		= start_time;
	usingErrorWeights = data.getBoolean("usingErrorWeights");
}
//Do a single integration of the model
void IntegratorCVode::integrateModel(Model *model, fileIO *out, int runNumber) {
	//Some stimulas mode parameters
//	model->setInitialConditions(*out);
	model->setupIFmode(runNumber);

	//CVode Refresh -> Initial conditions:
	refreshCVode(model);

	//during each iteration of this loop, an experiment object is set up.  A set of statevalues
	//is also computed.
	double stopTime = model->getStopTime();
	
	//Since we don't want to max out intermediate steps getting to time t
	for(double t = stepsize; t < start_time; ) {
		//Periodic output to let you know how it's going
		for(int i = 0; (t < start_time) && (i < 1000); i++, t = cvode_t + stepsize) {
			iterateToTime(t);
		}
		cout<<'.';
	}
	for(double t = start_time  + stepsize; t < stopTime; ) {
		//Periodic output to let you know how it's going
		for(int i = 0; (t < stopTime) && (i < 1000); i++, t = cvode_t + stepsize) {
			iterateToTime(t);
			out->writeLnOutputFiles( N_VGetData(cvode_y), cvode_t, model);
		}
		cout<<'.';
	}
	cout << "Complete." << endl;
}
