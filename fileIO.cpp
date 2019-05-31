#include "fileio.h"

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

fileIO::fileIO(void) {
	writeDependent = true;
	writeDerivatives = false;
	writestates = true;

	lastState = NULL;

	separator = char(9);

	sample_size = 0;
}

fileIO::~fileIO(void) {
}
//New data storage stuff
const double& fileIO::operator[](const string &param) {
	realDataType::iterator p = realData.find(param);
	if (p != realData.end())
		return p->second;
	//Might want to throw an exception in the future, but for now error and and pause
	cout << "An error has occured in 'fileIO::operator[]', numeric data: " << param << " parameter was not found." << endl;
	exit(-1);
}
//New data accessors
double fileIO::getReal(const string &param) { 	return operator[](param); }
int fileIO::getInteger(const string &param) { 
	double d = operator[](param);
	if ( (floor(d) - d) == 0 )
        return (int)(operator[](param)); 
	//Might want to throw an exception in the future, but for now error and and pause
	//Could also right round code using modf
	cout << "An error has occured in 'fileIO::operator[]', numeric data: " << param << " was: " << d << " not integer." << endl;
	exit(-1);
}
const string& fileIO::getKeyword(const string &param) {
	keywordDataType::iterator p = keywordData.find(param);
	if (p != keywordData.end())
		return p->second;
	//Might want to throw an exception in the future, but for now error and and pause
	cout << "An error has occured in 'fileIO::operator[]', keyword data: " << param << " parameter was not found." << endl;
	exit(-1);
}
bool fileIO::getBoolean(const string &param) {
	keywordDataType::iterator p = keywordData.find(param);
	if (p == keywordData.end()) {
		//Might want to throw an exception in the future, but for now error and and pause
		cout << "An error has occured in 'fileIO::operator[]', keyword data: " << param << " parameter was not found." << endl;
		exit(-1);
	}
	if( _strnicmp(p->second.c_str(), "true", 4) == 0 )
		return true;
	if( _strnicmp(p->second.c_str(), "false", 4) == 0 )
		return false;
	//Again with exceptions
	cout << "An error has occured in 'fileIO::operator[]', keyword data: " << param << " was: " << p->second << " not boolean." << endl;
	exit(-1);
}
//Test to see if a value is loaded
bool fileIO::parameterExists(const string &param) {
	keywordDataType::iterator p = keywordData.find(param);
	if (p == keywordData.end()) {
		realDataType::iterator r = realData.find(param);
		if (r == realData.end()) {
			return false;
		}
	}
	return true;
}

//readFile reads all the initial conditions in from an input file called "ini.txt"
void fileIO::loadFile(const char *ParameterFileName) {
	//Open the file and set up the iterators for the parsers
	file_iterator<> start(ParameterFileName);
    if (!start) {
		cerr << "No such file as" << ParameterFileName << "." << endl;
		exit(-1);
    }
    file_iterator<> end = start.make_end();

	parameterFileGrammer pfg(this);

	if (!(parse(start,end, pfg)).full) {//True if good parse, so the !
		cout << "Invalid file format for file: " << ParameterFileName << endl;
		exit(-1);
	}
}
// Writing Functions **************************************************
//Write the steadystate conditions file
//Might be useful: Form sci8; sci8.scientific().precision(8), instead of set with, see page 633
void fileIO::writeInitialConditionsFile(const char *filename, Model &model) {
	if (lastState == NULL) return;
	return writeInitialConditionsFile(filename, model, lastState);
}
void fileIO::writeInitialConditionsFile(const char *filename, Model &model, const double *states) {
	//Write file only if the flag is set
	if (getBoolean("write_end_initial_conditions")) {
		ofstream out(filename, ios::out);

		if( !out.is_open() ) {
				cout << "Error opening: " << filename << " for output." << endl;
				exit(-10);
		}

		time_t rawtime;
		tm *timenow;
		time(&rawtime);
		timenow = localtime(&rawtime);
		int year  = 1900 + timenow->tm_year;
		int month = timenow->tm_mon;
		int day   = timenow->tm_mday;
		int hour  = timenow->tm_hour;
		int min   = timenow->tm_min;

		out <<"/*"<<endl
			<< "Initial conditions file." << endl
			<< "Date: " << year << "." << month << "." << day << " at " << hour << ":" << min << " (ISO format: Year.Month.Date)" << endl
			<< endl
			<< "Format of the parameter file:" << endl
			<< char(9) << "Comment lines start with //, !, %, Or block" << endl
			<< char(9) << "Initial conditions are in the following format:" << endl
			<< char(9) << "[variable] = [value]" << endl
			<<"*/"<<endl;
		//setw(14)


		//ios_base:: scientific, ios_base:: floatfield
		//ios_base::boolalpha
		// ios_base::dec, ios_base::adjustfield
		out.setf(ios_base::scientific, ios_base:: floatfield);
		out.setf(ios_base::boolalpha);
		out.setf(ios_base::dec, ios_base::adjustfield);
		out.precision(16);

		const double *s = model.getInitialConditions(states);
		for(int i = 0; i < model.getInitialConditionsSize(); i++) {
			out << model.getInitialConditionLabel(i) << char(9) << "= " << s[i] << endl;
		}

		out.close();
	}
}
//printstates outputs to 3 separate files depending on the input file.
//Each file contains either the states, currents, or derivatives.
void fileIO::openOutputFiles(Model* model, const char * statesFilename, const char * currentsFilename, const char * derivativeFilename) {
	if (writestates) {
		if( outStates.is_open() )
			outStates.close();
		outStates.open(statesFilename, ios::out );
		if( !outStates.is_open() ) {
			cout << "Error opening: " << statesFilename << " for output." << endl;
			exit(-1);
		}
		outStates.setf(ios_base::scientific, ios_base:: floatfield);
		outStates.setf(ios_base::boolalpha);
		outStates.setf(ios_base::dec, ios_base::adjustfield);
		outStates.precision(getInteger("States_Precision"));

		keywordList *stateList = &(keywordListData["Ouput_State_List"]);
		keywordList::iterator b = stateList->begin();
		stateIndexes.resize(stateList->size());
		if(outputMode == MCA) {
			outStates << "Sample_Size";
			outStates << separator << "Parameter";
			outStates << separator << "Percent";
		} else if (outputMode == Average) {
			outStates << "Metric";
			outStates << separator << "Sample_Size";
		} else {
			outStates << "Time";
		}
		for ( unsigned int i = 0; i < stateList->size(); i++) {
			outStates << separator << *b;
			stateIndexes[i] = model->getStateIndex( (*b).c_str() );
			b++;
		}			
		outStates << endl;
	}
	
	if (writeDependent) {
		if( outCurrents.is_open() )
			outCurrents.close();
		outCurrents.open(currentsFilename, ios::out);
		if( !outCurrents.is_open() ) {
			cout << "Error opening: " << currentsFilename << " for output." << endl;
			exit(-1);
		}
		outCurrents.setf(ios_base::scientific, ios_base:: floatfield);
		outCurrents.setf(ios_base::boolalpha);
		outCurrents.setf(ios_base::dec, ios_base::adjustfield);
		outCurrents.precision(getInteger("Dependent_Precision"));

		keywordList *dependentList = &(keywordListData["Ouput_Dependent_List"]);
		keywordList::iterator b = dependentList->begin();
		dependentIndexes.resize(dependentList->size());
		if(outputMode == MCA) {
			outCurrents << "Sample_Size";
			outCurrents << separator << "Parameter";
			outCurrents << separator << "Percent";
		} else if (outputMode == Average) {
			outCurrents << "Metric" << separator;
			outCurrents << "Sample_Size";
		} else {
			outCurrents << "Time";
		}
		for ( unsigned int i = 0; i < dependentList->size(); i++) {
			outCurrents << separator << *b;
			dependentIndexes[i] = model->getDependentVariableIndex( (*b).c_str() );
			b++;
		}			
		outCurrents << endl;
	}

	if (writeDerivatives) {
		if( outDerivatives.is_open() )
			outDerivatives.close();
		outDerivatives.open(derivativeFilename, ios::out);
		if( !outDerivatives.is_open() ) {
			cout << "Error opening: " << derivativeFilename << " for output." << endl;
			exit(-1);
		}
		outDerivatives.setf(ios_base::scientific, ios_base:: floatfield);
		outDerivatives.setf(ios_base::boolalpha);
		outDerivatives.setf(ios_base::dec, ios_base::adjustfield);
		outDerivatives.precision(getInteger("Derivatives_Precision"));

		keywordList *derivativeList = &(keywordListData["Ouput_Derivative_List"]);
		keywordList::iterator b = derivativeList->begin();
		derivativeIndexes.resize(derivativeList->size());
		outDerivatives << "Time";
		for ( unsigned int i = 0; i < derivativeList->size(); i++) {
			outDerivatives << separator << *b;
			derivativeIndexes[i] = model->getStateIndex( (*b).c_str() );
			b++;
		}			
		outDerivatives << endl;
	}

	
}
//Setup the various output file variables before writting
void fileIO::setupFileOut() {

	//Setup file output
	string out_mode = getKeyword("Output_Mode");
	if (out_mode == "numeric")
		outputMode = Numeric;
	if (out_mode == "average")
		outputMode = Average;
	if (out_mode == "MCA")
		outputMode = MCA;
	
	writeDependent = getBoolean("Write_Dependent");
	writestates = getBoolean("Write_States");
	writeDerivatives = getBoolean("Write_Derivatives");
}
//Ouput averages etc
void fileIO::finalizeFileOut() {
	if(outputMode==Average) {
		int s_size = sample_size;	
		MCA_AverageData();
		if (!outStates.is_open() && writestates) {
			cout << "Attempting to output to the states file without opening." << endl;
			exit(-1);
		} else if (writestates) { //Actually output all the states
			outStates << "Mean" << separator;
			outStates << s_size;
			for ( unsigned int i = 0; i < stateAverage.size(); i++) { outStates << separator << stateAverage[ i ]; }			
			outStates << endl << "Max" << separator << s_size;
			for ( unsigned int i = 0; i < stateMax.size(); i++) { outStates << separator << stateMax[ i ]; }			
			outStates << endl << "Min" << separator << s_size;
			for ( unsigned int i = 0; i < stateMin.size(); i++) { outStates << separator << stateMin[ i ]; }
			outStates << endl;
		}
		if (!outCurrents.is_open() && writeDependent) {
			cout << "Attempting to output to the currents file without opening." << endl;
			exit(-1);
		} else if (writeDependent) { //Actually output all the states
			outCurrents << "Mean" << separator;
			outCurrents << s_size;
			for ( unsigned int i = 0; i < dependentAverage.size(); i++) { outCurrents << separator << dependentAverage[ i ]; }			
			outCurrents << endl << "Max" << separator << s_size;
			for ( unsigned int i = 0; i < dependentMax.size(); i++) { outCurrents << separator << dependentMax[ i ]; }			
			outCurrents << endl << "Min" << separator << s_size;
			for ( unsigned int i = 0; i < dependentMin.size(); i++) { outCurrents << separator << dependentMin[ i ]; }
			outCurrents << endl;
		}
		if (!outDerivatives.is_open() && writeDerivatives) {
			cout << "Attempting to output to the Derivatives file without opening." << endl;
			exit(-1);
		} else if (writeDerivatives) { //Actually output all the states
			outDerivatives << time;
			for ( unsigned int i = 0; i < derivativeIndexes.size(); i++) {
				outDerivatives << separator << dependentAverage[ i ];
			}			
			outDerivatives << endl;
		}	
	}
}
//Close any open file handles
void fileIO::closeOutputFiles() {
	if (outStates.is_open())
		outStates.close();
	if (outCurrents.is_open())
		outCurrents.close();
	if (outDerivatives.is_open())
		outDerivatives.close();
}
// write incremental output
//Fundamental flaw: CAN NOT use the model for output variables, CVode interpolates back to get y.
//CVode has a built in function for derivatives
void fileIO::writeLnOutputFiles(double * states, double time, Model* model) {
	//Save the last state for IC
	lastState = states;
	sample_size++;

	model->F(time, model->getStatesDerivatives(), states);
	const double *derivs = model->getStatesDerivatives();
	//const double *vars = model.getDependentVariables(); //Need to write this, old currents
	
	//Output the states
	//[ V, Nai, Ki, Cai, CaNSR, CaSS, CaJSR, ATPi, Cam, ADPm, Dpsi, NADH, ISOC, AKG, SCoA, Succ, FUM, MAL, Oaa ]
	if(writestates) {
		if (!outStates.is_open() ) {
			cout << "Attempting to output to the states file without opening." << endl;
			exit(-1);
		} //Actually output all the states
		switch(outputMode) {
		case Numeric:
			outStates << time;
			//outStates.write( (char*)(&time), sizeof(double) );
			for ( unsigned int i = 0; i < stateIndexes.size(); i++) {
				outStates << separator << states[ stateIndexes[i] ];
				//outStates.write( (char*)(&states[ stateIndexes[i] ]), sizeof(double) );
			}			
			outStates << endl;
			break;
		case MCA:
		case Average:
			for ( unsigned int i = 0; i < stateIndexes.size(); i++) {
				 stateAverage[i] += states[ stateIndexes[i] ];
				 stateMin[i] = __min( stateMin[i], states[ stateIndexes[i] ] );
				 stateMax[i] = __max( stateMax[i], states[ stateIndexes[i] ] );
			}		
			break;
		}
	}

	if (writeDependent) {
		if (!outCurrents.is_open() ) {
			cout << "Attempting to output to the currents file without opening." << endl;
			exit(-1);
		} //Actually output all the states
		switch(outputMode) {
		case Numeric:
			outCurrents << time;
			for ( unsigned int i = 0; i < dependentIndexes.size(); i++) {
				outCurrents << separator << model->getDependentVariable( dependentIndexes[i] );
			}			
			outCurrents << endl;
			break;
		case MCA:
		case Average:
			for ( unsigned int i = 0; i < dependentIndexes.size(); i++) {
				 dependentAverage[i] += model->getDependentVariable( dependentIndexes[i] );
				 dependentMin[i] = __min( dependentMin[i], model->getDependentVariable( dependentIndexes[i] ) );
				 dependentMax[i] = __max( dependentMax[i], model->getDependentVariable( dependentIndexes[i] ) );
			}		
			break;
		}
	}

	//Output the derivatives
	//[ V, Nai, Ki, Cai, CaNSR, CaSS, CaJSR, ATPi, Cam, ADPm, Dpsi, NADH, ISOC, AKG, SCoA, Succ, FUM, MAL, Oaa ]
	if (writeDerivatives) {
		if (!outDerivatives.is_open()) {
			cout << "Attempting to output to the Derivatives file without opening." << endl;
			exit(-1);
		} //Actually output all the states
		switch(outputMode) {
		case Numeric:
			outDerivatives << time;
			for ( unsigned int i = 0; i < derivativeIndexes.size(); i++) {
				outDerivatives << separator << derivs[ derivativeIndexes[i] ];
			}			
			outDerivatives << endl;
			break;
		case MCA:
		case Average:
			for ( unsigned int i = 0; i < derivativeIndexes.size(); i++) {
				 derivativeAverage[i] += derivs[ stateIndexes[i] ];
			}		
			break;
		}
	}
}

// ******************************************************************
// ****************** MCA Specific Functions ************************
// ******************************************************************
void fileIO::MCA_Setup(Model* model) {
	//Make sure the sample size is 0 to start
	sample_size = 0;


	stateAverage.resize(stateIndexes.size());
	stateMin.resize(stateIndexes.size());
	stateMax.resize(stateIndexes.size());
	derivativeAverage.resize(derivativeIndexes.size());
	dependentMin.resize(dependentIndexes.size());
	dependentMax.resize(dependentIndexes.size());
	dependentAverage.resize(dependentIndexes.size());

	for ( unsigned int i = 0; i < stateAverage.size(); i++) { stateAverage[i] = 0; }
 	for ( unsigned int i = 0; i < stateMin.size(); i++) { stateMin[i] = DBL_MAX; }
	for ( unsigned int i = 0; i < stateMax.size(); i++) { stateMax[i] = DBL_MIN; }
	for ( unsigned int i = 0; i < derivativeAverage.size(); i++) { derivativeAverage[i] = 0; }
 	for ( unsigned int i = 0; i < dependentMin.size(); i++) { dependentMin[i] = DBL_MAX; }
	for ( unsigned int i = 0; i < dependentMax.size(); i++) { dependentMax[i] = DBL_MIN; }
	for ( unsigned int i = 0; i < dependentAverage.size(); i++) { dependentAverage[i] = 0; }

	writestates		= getBoolean("Write_States");
	writeDependent	= getBoolean("Write_Dependent");
	writeDerivatives= getBoolean("Write_Derivatives");

	stateOldAverage.resize(stateIndexes.size());
	derivativeOldAverage.resize(derivativeIndexes.size());
	dependentOldAverage.resize(dependentIndexes.size());

	// Load the percent and parameter lists
	//realList MCA_percents;
	//keywordList MCA_Params;
	keywordList *paramList = &(keywordListData["MCA_Params"]);
	MCA_Params.resize( paramList->size() );
	for ( unsigned int i = 0; i < paramList->size(); i++ ) {
		MCA_Params[i] = paramList->at(i);
	}
	
	realList *percentList = &(realListData["MCA_percents"]);
	MCA_percents.resize(percentList->size());
	for ( unsigned int i = 0; i < percentList->size(); i++ ) {
		MCA_percents[i] = percentList->at(i);
	}
}

//Average the data over the sample size
void fileIO::MCA_AverageData() {
	//Calculate Means
	for(unsigned int i = 0; i < stateAverage.size(); i++) {
		stateAverage[i] /= sample_size;
	}
	for(unsigned int i = 0; i < derivativeAverage.size(); i++) {
		derivativeAverage[i] /= sample_size;
	}
	for(unsigned int i = 0; i < dependentAverage.size(); i++) {
		dependentAverage[i] /= sample_size;
	}
	//Reset sample size
	sample_size = 0;
}
void fileIO::MCA_Analyze(Model* model, int param, int percent) {
    //Is this not the first trail of this parameter?
	if( percent != 0 ) {
		//Perform the MCA analysis and output
		double inv_flux = 1.0 / ( log(MCA_percent)-log(MCA_percent_old) );
		//Ouput header info stuff including sample size
		int n = sample_size;
		//Average the current data set
		MCA_AverageData();
		//Ouput data
		if(writestates) {
			outStates << n;
			outStates << separator << MCA_Params[param].c_str();
			outStates << separator << MCA_percents[percent];
			cout << "writestates " << stateAverage.size() << endl;
			for ( unsigned int i = 0; i < stateAverage.size(); i++) {
				outStates << separator << ( log(stateAverage[i]) - log(stateOldAverage[i]) ) * inv_flux ;
			}			
			outStates << endl;
		}
		if(writeDerivatives) {
			outDerivatives << n;
			outDerivatives << separator << MCA_Params[param].c_str();
			outDerivatives << separator << MCA_percents[percent];
			for ( unsigned int i = 0; i < derivativeAverage.size(); i++) {
				outDerivatives << separator << ( log(derivativeAverage[i]) - log(derivativeOldAverage[i]) ) * inv_flux ;
			}			
			outDerivatives << endl;
		}
		if(writeDependent) {
			outCurrents << n;
			outCurrents << separator << MCA_Params[param].c_str();
			outCurrents << separator << MCA_percents[percent];
			for ( unsigned int i = 0; i < dependentAverage.size(); i++) {
				outCurrents << separator << ( log(dependentAverage[i]) - log(dependentOldAverage[i]) ) * inv_flux ;
			}			
			outCurrents << endl;
		}
	} else {
		MCA_AverageData();
		if (param != 0) {
			if(writestates)
				outStates << endl;
			if(writeDependent)
				outCurrents << endl;
			if(writeDerivatives)
				outDerivatives << endl;
		}
	}
	//Update Old average and reset average
	MCA_percent_old = MCA_percent;
	for ( unsigned int i = 0; i < stateAverage.size(); i++) {
		stateOldAverage[i] = stateAverage[i];
		stateAverage[i] = 0;
	}
	for ( unsigned int i = 0; i < derivativeAverage.size(); i++) {
		derivativeOldAverage[i] = derivativeAverage[i];
		derivativeAverage[i] = 0;
	}
	for ( unsigned int i = 0; i < dependentAverage.size(); i++) {
		dependentOldAverage[i] = dependentAverage[i];
		dependentAverage[i] = 0;
	}
	//Reset parameters from saved data
	model->setParamater(MCA_Params[param].c_str(), MCA_param_value );
}
int fileIO::MCA_GetNumberOfParameters() { return MCA_Params.size(); }
int fileIO::MCA_GetNumberOfPercents() { return MCA_percents.size(); }
void fileIO::MCA_Update(Model* model, int param, int percent) {	
	// Determine the percentage change
	MCA_percent = (100.0+MCA_percents[percent])/100.0;
	//Find the default value
	MCA_param_value = getReal(MCA_Params[param].c_str());
	//Update the model with a modified parameter
	model->setParamater(MCA_Params[param].c_str(), MCA_percent * MCA_param_value );
}
bool fileIO::MCA_isEnabled() { return outputMode == MCA; }
