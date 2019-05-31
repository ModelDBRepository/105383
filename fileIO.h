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

#include <string>
#include <map>
#include <ctype.h>
#include <fstream>
#include <stdlib.h>  
#include <errno.h>  
#include <iomanip>
#include <math.h>
#include <iostream>
#include <time.h>
#include <iterator>
#include <vector>
#include <limits>

//Boost Parsing Library
#include <boost/spirit/core.hpp>
#include <boost/spirit/utility/confix.hpp>
#include <boost/spirit/actor/assign_actor.hpp>
#include <boost/spirit/actor/insert_key_actor.hpp>
#include <boost/spirit/actor/clear_actor.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#include <boost/spirit/iterator/file_iterator.hpp> 

#include "Model.h"
#include "states_index.h"

class Model;

using namespace std;
using namespace boost::spirit;

//TODO use boolalpha set on streamsa to autoconvert true/false to boolean and viceversa
class fileIO
{
public:
	fileIO(void);
	virtual ~fileIO(void);

	void loadFile(const char *ParameterFileName);

	void writeInitialConditionsFile(const char *filename, Model &model);
	void writeInitialConditionsFile(const char *filename, Model &model, const double * states);

	void openOutputFiles(Model* model, const char * statesFilename, const char * currentsFilename, const char * derivativeFilename);
	void writeLnOutputFiles(double * states, double time, Model* model);

	//New data accessors
	const double& operator[](const string &param);
	double getReal(const string &param);
	int getInteger(const string &param);
	const string& getKeyword(const string &param);
	bool getBoolean(const string &param);
	bool parameterExists(const string &param);

	void setupFileOut();
	void finalizeFileOut();
	void closeOutputFiles();

	int MCA_GetNumberOfParameters();
	int MCA_GetNumberOfPercents();
	void MCA_Analyze(Model* model, int param, int percent);
	void MCA_Setup(Model* model);
	void MCA_Update(Model* model, int param, int percent);
	bool MCA_isEnabled();
private:
	typedef map< string, double > realDataType;
	typedef map< string, string > keywordDataType;
	typedef pair< string, double > realDataValueType;
	typedef pair< string, string > keywordDataValueType;
	typedef vector< string > keywordList;
	typedef vector< double > realList;
	typedef map< string, keywordList > keywordListDataType;
	typedef map< string, realList > realListDataType;
	typedef pair< string, keywordList > keywordListDataValueType;
	typedef pair< string, realList > realListDataValueType;
	typedef vector< int > indexList;


	//MCA Stuff
	void MCA_AverageData();

	int MCA_Flux_index;
	bool MCA_Flux_isState;
	int sample_size; //The size of the average/MCA
	realList MCA_percents;
	keywordList MCA_Params;
	double MCA_percent_old;
	double MCA_percent;
	double MCA_param_value;

	//Output mode
	enum OutputMode {Numeric, Average, MCA};
	OutputMode outputMode;

	//File handles:
	ofstream outStates;
	ofstream outCurrents;
	ofstream outDerivatives;
	bool writeDependent;
	bool writestates;
	bool writeDerivatives;

	//New data Storage
	realDataType realData;
	keywordDataType keywordData;

	//List storage types
	keywordListDataType keywordListData;
	realListDataType realListData;

	realList stateMin;
	realList stateMax;
	realList stateAverage;
	realList dependentMin;
	realList dependentMax;
	realList derivativeAverage;
	realList dependentAverage;

	//MCA Only
	realList stateOldAverage;
	realList derivativeOldAverage;
	realList dependentOldAverage;


	//Range storage types, special list limited to 2 or 3 datum
	realListDataType realRangeData;

	//Output variables
	char separator; //tab
	indexList stateIndexes;
	indexList dependentIndexes;
	indexList derivativeIndexes;

	const double * lastState;	

	struct parameterFileGrammer : public grammar<parameterFileGrammer>
    {
        template <typename ScannerT> struct definition
        {
			realDataValueType n;
			keywordDataValueType s;
			keywordListDataValueType sl;
			realListDataValueType nl;

			rule<ScannerT> numPair, strPair, comments, allParams, label, equals;
			rule<ScannerT> rightBrace, leftBrace, seperator, numList, strList, numListPair, strListPair;
			rule<ScannerT> numRange, numRangePair;

			definition(parameterFileGrammer const &self)
            {
				//Shared Grammers
				label = (alpha_p | "_") >> *( alpha_p | digit_p | "_");//Allow underscores in names
				equals = *blank_p >> '=' >> *blank_p;
				comments = comment_p("!") | comment_p("//") | comment_p("/*","*/") | comment_p("%") | space_p;
				//List based grammers
				rightBrace = *comments >> (ch_p(')') | ch_p('}') | ch_p(']'));
				leftBrace  = (ch_p('(') | ch_p('{') | ch_p('[')) >> *comments;
				seperator  =  (*comments >> ch_p(',') >> *comments) | +comments;
				numList = confix_p( 
					leftBrace, real_p[push_back_a(nl.second)] >> *( seperator >> real_p[push_back_a(nl.second)] ), rightBrace );
				strList = confix_p( 
					leftBrace, label[push_back_a(sl.second)] >> *( seperator >> label[push_back_a(sl.second)] ), rightBrace );
				numListPair = confix_p(
					label[assign_a(nl.first)] , equals, numList >> !ch_p(';')
					)[insert_key_a(self.data->realListData,nl)] >> epsilon_p[clear_a(nl.second)];
				strListPair = confix_p(
					label[assign_a(sl.first)] , equals, strList >> !ch_p(';')
					)[insert_key_a(self.data->keywordListData,sl)] >> epsilon_p[clear_a(sl.second)];
				//Range Grammer
				numRange = confix_p( 
					leftBrace, 
					*comments >> real_p[push_back_a(nl.second)] >> *comments >>  ':'  >> *comments >>  real_p[push_back_a(nl.second)] >> *comments >> !(':'  >> *comments >>  real_p[push_back_a(nl.second)] >> *comments),
					rightBrace );
				numRangePair = confix_p(
					label[assign_a(nl.first)] , equals, numRange >> !ch_p(';')
					)[insert_key_a(self.data->realRangeData,nl)] >> epsilon_p[clear_a(nl.second)];
				//Single Parameter Grammers
				numPair = confix_p(
					label[assign_a(n.first)] , equals, real_p[assign_a(n.second)] >> !ch_p(';')
					)[insert_key_a(self.data->realData,n)];
				strPair = confix_p(
					label[assign_a(s.first)] , equals, label[assign_a(s.second)] >> !ch_p(';') 
					)[insert_key_a(self.data->keywordData,s)];
				//Complete Grammer
				allParams = *(numPair | strPair | comments | strListPair | numListPair | numRangePair );
            }
            rule<ScannerT> const& start() const { return allParams; }
        };
		parameterFileGrammer(fileIO *dataStore) {
			data = dataStore;
		}
		fileIO * data;
    };
	friend parameterFileGrammer;
};
