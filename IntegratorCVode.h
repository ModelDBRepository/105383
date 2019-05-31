#pragma once

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

#include "integrator.h"
#include "model.h"
#include "fileIO.h"

#include <sundialstypes.h>	/* definitions of types realtype and             */
							/* integertype, and the constant FALSE           */
#include <cvode.h>			/* prototypes for CVodeMalloc, CVode, and CVodeFree, */
							/* constants OPT_SIZE, BDF, NEWTON, SV, SUCCESS,     */
							/* NST, NFE, NSETUPS, NNI, NCFN, NETF                */
#include <cvdense.h>		/* prototype for CVDense, constant DENSE_NJE         */
#include <nvector.h>		/* definitions of type N_Vector and macro N_VIth,    */
#include <nvector_serial.h>	/* definitions of type N_Vector and macro N_VIth,    */
							/* prototypes for N_VNew, N_VFree                    */
#include <dense.h>			/* definitions of type DenseMat, macro DENSE_ELEM    */

#include <stdlib.h>  
#include <iostream>
using namespace std;

class IntegratorCVode :
	public Integrator
{
public:
	IntegratorCVode(void);
	virtual ~IntegratorCVode(void);

	void setParameters(fileIO& data);

	void setupCVode(Model *model);

	void integrateModel(Model *model, fileIO *out, int runNumber);

protected:
	void iterateToTime(double time);
	void refreshCVode(Model *model);

	//Internal variables
	double cvode_t;
	bool usingErrorWeights;
	N_Vector ew;
	N_Vector cvode_y;
	N_Vector cvode_ic;
	void* cvode_mem;
	double* ropt;
	long* iopt;
	M_Env machEnv;

	//Parameters
	double reltol;
	double abstol;
	double step_max;
	double step_min;
	double start_time;
	double stepsize;
};
