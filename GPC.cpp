// GPC.cpp : Defines the entry point for the console application.
//

#include "Integrator.h"
#include "integratorcvode.h"
#include "fileIO.h"
#include "Model.h"

#include <time.h>
#include <iostream>
using namespace std;

void displayCopyright() {
	// Outputting Introduction and Copyright information
	cout << "-------------------------------------------------------\n\n";
	cout << "      NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE    \n";
	cout << "       Copyright 2004, The Johns Hopkins University    \n";
	cout << "        School of Medicine. All rights reserved.       \n\n";
	cout << "  For research use only; commercial use prohibited.    \n";
	cout << "  Distribution without permission of Raimond L. Winslow\n";
	cout << "  not permitted. rwinslow@bme.jhu.edu					\n\n";
	cout << "  Name of Program: Guinea Pig C++ (Coupled)			\n";
    cout << "  Version: Documented Version, version 1.0.8			\n";
	cout << "  Date: November 2004									\n";
    cout << "-------------------------------------------------------\n";  
}

int main(int argc, char* argv[])
{
	//Variable Declarations
	clock_t tic = clock();
	Model gpc;
	fileIO data;
	IntegratorCVode cvode;

	displayCopyright();

	//Load parameters in classes
	data.loadFile("parameters.txt");
	data.loadFile("initial_conditions.txt");
	data.loadFile("control.txt");
	gpc.setParameters(data);
	gpc.setInitialConditions(data);
	cvode.setParameters(data);

	//Init CVode structures
	cvode.setupCVode(&gpc);

	//Initialize the control variables
	data.setupFileOut();

	//Open files
	data.openOutputFiles(&gpc, "States.txt", "Currents.txt", "Derivatives.txt" );

	//Setup MCA and averaging
	data.MCA_Setup(&gpc);

	//Should we perform an MCA
	if ( data.MCA_isEnabled() ) {
		for (int param = 0; param < data.MCA_GetNumberOfParameters(); param++) {
			for (int percent = 0; percent < data.MCA_GetNumberOfPercents(); percent++) {
				data.MCA_Update(&gpc,param,percent);
				cvode.integrateModel(&gpc, &data, 0);
				data.MCA_Analyze(&gpc,param,percent);
			}
		}
	} else {
		//IF mode loop
		for(int i = 0; i < gpc.getNumRun(); i++) {
			data.setupFileOut();
			cvode.integrateModel(&gpc, &data, i);
		}
	}

	data.finalizeFileOut();
	data.closeOutputFiles();
	data.writeInitialConditionsFile("SSConditions.txt", gpc);
	// calculate time spent
	cout << "Time used: "<< (clock() - tic) / (double)CLOCKS_PER_SEC<<"s\a\a\a"<< endl;
	return 0;
}