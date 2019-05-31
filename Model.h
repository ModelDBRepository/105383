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
// CVODE Types are included in sundialstypes.h. Default is double

#include "states_index.h"
#include "fileIO.h"

#include <math.h>
#include <stdlib.h>  
#include <iostream>
#include <map>
#include <string>
using namespace std;
//#include <float.h>
//#include <algorithm>
class fileIO;

class Model 
{
public:
	Model(void);
	~Model(void);

	//Check main to see if this is the proper usage
	//Something is funny about start_time usage
	void F(double start_time, double *ydot, double *y);

	double* getErrorWeights();
	double* getStates();
	double* getStatesDerivatives();
	int getProblemSize();

	int getStateIndex(const char * label);
	int getDependentVariableIndex(const char * label);
	double getDependentVariable(int index);

	void setParameters(fileIO& data);
	void setParamater(const string &name, double value);
	void setInitialConditions(fileIO& data);
	
	//These are only for final state reporting, e.g. SS_conditions.txt
	const char * getStateLabel(int index);
	const char * getInitialConditionLabel(int index);
	const double* getInitialConditions(const double *y); //Use ONLY for SS/IC output and not integrator
	int getInitialConditionsSize();

	//Should get moved out of here into a "Control" class, todo in future
	//Checks and only runs if in IF mode, run for each trial.
	void setupIFmode(int num_run);
	double getStopTime();
	double getNumRun();
	
private:
	enum StimulasMode {	BB, IF, periodic, none };
	enum ModelType { Coupled, Mitochondria, Force };

	//New data handling stuff
	typedef map< string, double* > parameterLookupType;
	typedef map< string, StimulasMode > stimulasModeLookupType;
	typedef map< string, ModelType > modelTypeLookupType;
	typedef map< int, string > stateLabelLookupType;
	typedef map< string, int > stateIndexLookupType;
	parameterLookupType parameterLookup;
	stimulasModeLookupType stimulasModeLookup;
	modelTypeLookupType modelTypeLookup;
	stateLabelLookupType stateLabelLookup;
	stateIndexLookupType stateIndexLookup;
	stateIndexLookupType dependentVariableIndexLookup;

	void setStateWithLink(fileIO& data, const string &name, int index);
	void setParamaterWithLink(fileIO& data,const string &name, double &param);
	void linkParameter(const string &name, double &param);

	//Basic model Properties
	void initializeModel(bool algebraicMode = false);

	double *s;		//States data
	double *dS;	//Derivative States data
	double *errweight;
	int getIndexOffset();

	//remember with flags: 1 = true, 0 = false
	bool usingAlgebraicMembranePotential;
	bool usingCHF;
	bool clampVoltage;
	bool usingASP;
	bool usingCK;

	StimulasMode stimulasMode;
	ModelType modelType;

	void F_GPC(double start_time);
	void F_Mitochondria();
	void F_Force();

	void linkStatesReferences_GPC( double *ydot, double *y );
	void linkStatesReferences_CK_Only( double *ydot, double *y , int offset);
	void linkStatesReferences_Full( double *ydot, double *y );
	void linkStatesReferences_Alg( double *ydot, double *y );
	void linkStatesReferences_Force( double *ydot, double *y );
	void linkStatesReferences_Mito( double *ydot, double *y );
	void linkStatesReferences_Mito_ASP( double *ydot, double *y );

	//Called by CVode_f as does most of the work
	//Funtion Prototypes
	
	void calculateConstantModelParameters();
	void sharedVariableUpdate_GPC_1(double start_time);
	void sharedVariableUpdate_GPC_2();
	void sharedVariableUpdate_Mitochondria();
	void sharedVariableUpdate_Force();

	void calculateAlgebraicMembranePotential();
	void calculateStimulatedCurrent(double start_time);
	void calculateClampedVoltage();
	void calculateReversalPotentials();
	void calculateINa();
	void calculateIKs();
	void calculateIK1();
	void calculateINab();
	void calculateIKp();
	void calculateICa();
	void calculateICaK();
	void calculateINaK();
	void calculateINaCa();
	void calculateICab();
	void calculateIpCa();
	void calculateInsCa();
	void calculateV_AM();
	void calculateATPm();			//Mitochondria Model +
	void calculateDmuH();			//Mitochondria Model +
	void calculateNAD();			//Mitochondria Model +
	void calculateVCS();			//Mitochondria Model +
	void calculateVACO();			//Mitochondria Model +
	void calculateVIDH();			//Mitochondria Model +
	void calculateVKGDH();			//Mitochondria Model +
	void calculateVSL();			//Mitochondria Model +
	void calculateVSDH();			//Mitochondria Model +
	void calculateVFH();			//Mitochondria Model +
	void calculateVMDH();			//Mitochondria Model +
	void calculateVAAT();			//Mitochondria Model +
	void calculateVAAT_ASP();		//Mitochondria Model +
	void calculateVNO_VHNe();		//Mitochondria Model - This is the funny one, close enough for jazz
	void calculateVFO_VHFe();		//Mitochondria Model +
	void calculateVATPase_Vhu();	//Mitochondria Model +
	void calculateVANT_Vhleak();	//Mitochondria Model +
	void calculateFNa();
	void calculateFxKs();
	void calculateInCalFlux();
	void calculateForce();			//Force Model
	void calculateF_trpmyo();		//Force Model
	void calculateFLTRPNCa();		//Force Model
	void calculateFHTRPNCa();
	void calculateJtrpn();
	void calculateVuni();			//Mitochondria Model +
	void calculateVnaCa();			//Mitochondria Model +
	void calculateFNai();
	void calculateFKi();
	void calculateFCai();
	void calculateFCaSS();
	void calculateFCaJSR();
	void calculateFCaNSR();
	void calculateFV();
	void calculateFRyR();
	void calculateFCaL();
	void calculateFyCa();
	void calculateFOCa();
	void calculateFATPi();
	void calculateCK();				//CK Mod
	void calculateF_mitene();		//Mitochondria Model +
	void calculateFCrPi();			

	//Variable Definitions
	//Direct refences to the states array
	double *V;	
	double *mNa;   
	double *hNa;
	double *jNa;
	double *xKs;
	double *Nai;
	double *Ki;
	double *Cai;
	double *CaNSR;
	double *CaSS;
	double *CaJSR;
	double *C1_RyR; 
	double *O2_RyR;
	double *C2_RyR;
	double *C0;
	double *C1;
	double *C2;
	double *C3;
	double *C4;
	double *Open; 
	double *CCa0;
	double *CCa1;
	double *CCa2;
	double *CCa3;
	double *CCa4;
	double *OCa;
	double *yCa;
	double *LTRPNCa;
	double *HTRPNCa;
	double *N0;
	double *N1;
	double *P0;
	double *P1;
	double *P2;
	double *P3;
	double *ATPi;
	double *Cam;
	double *ADPm;
	double *Dpsi;
	double *NADH;
	double *ISOC;
	double *AKG;
	double *SCoA;
	double *Succ;
	double *FUM;
	double *MAL;
	double *Oaa;
	double *ASP;
	double *ATPi_cyto;
	double *CrPi_mito;
	double *CrPi_cyto;

	//Direct refences to the derivative states array
	double *dV;	
	double *dmNa;   
	double *dhNa;
	double *djNa;
	double *dxKs;
	double *dNai;
	double *dKi;
	double *dCai;
	double *dCaNSR;
	double *dCaSS;
	double *dCaJSR;
	double *dC1_RyR; 
	double *dO2_RyR;
	double *dC2_RyR;
	double *dC0;
	double *dC1;
	double *dC2;
	double *dC3;
	double *dC4;
	double *dOpen; 
	double *dCCa0;
	double *dCCa1;
	double *dCCa2;
	double *dCCa3;
	double *dCCa4;
	double *dOCa;
	double *dyCa;
	double *dLTRPNCa;
	double *dHTRPNCa;
	double *dN0;
	double *dN1;
	double *dP0;
	double *dP1;
	double *dP2;
	double *dP3;
	double *dATPi;
	double *dCam;
	double *dADPm;
	double *dDpsi;
	double *dNADH;
	double *dISOC;
	double *dAKG;
	double *dSCoA;
	double *dSucc;
	double *dFUM;
	double *dMAL;
	double *dOaa;
	double *dASP;
	double *dATPi_cyto;
	double *dCrPi_mito;
	double *dCrPi_cyto;

	//Loaded Model Parameters
	double kt_2;		//CK Mod
	double kf_2;		//CK Mod
	double kf_3;		//CK Mod
	double keq;			//CK Mod
	double CRT_cyto;	//CK Mod
	double CRT_mito;	//CK Mod
	double VATPase_cyto;//CK Mod
	double Acap;
	double AcCoA;		//Mitochondria Model
	double aL;
	double b;			//Mitochondria Model
	double bL;
	double C_m;
	double Cao;
	double chfsc_IK1;
	double chfsc_INaCa;
	double chfsc_Jup;
	double CIK;			//Mitochondria Model
	double Cm;			//Mitochondria Model
	double CMDNtot;
	double Cmito;		//Mitochondria Model
	double CPN;			//Mitochondria Model
	double CoA;			//Mitochondria Model
	double CSQNtot;
	double DpH;			//Mitochondria Model
	double Dpsio;		//Mitochondria Model
	double ESI_increment;	//IFmode Setup
	double eta;
	double EtCS;		//Mitochondria Model
	double EtID;		//Mitochondria Model
	double EtKG;		//Mitochondria Model
	double EtMD;		//Mitochondria Model
	double EtSDH;		//Mitochondria Model
	double FAD;			//Mitochondria Model
	double FADH2;		//Mitochondria Model
	double fL;
	double fm;			//Mitochondria Model
	double fprime;
	double g;			//Mitochondria Model
	double G_Cab;
	double G_Kp;
	double G_Na;
	double G_Nab;
	double gh;			//Mitochondria Model
	double gL;
	double GLU;			//Mitochondria Model
	double gprime;
	double H;			//Mitochondria Model
	double high_freq;	//Brandes and Bers parameters
	double hm;			//Mitochondria Model
	double HTRPNtot;
	double ICahalf;
	double INaKmax;
	double IpCamax;
	double KAATeq;		//Mitochondria Model
	double KaCa;		//Mitochondria Model
	double KACOeq;		//Mitochondria Model
	double kact;		//Mitochondria Model
	double KADP;		//Mitochondria Model
	double kaminus;
	double kaplus;
	double kbminus;
	double kbplus;
	double Kca;			//Mitochondria Model
	double kcminus;
	double kcnsASP;		//Mitochondria Model
	double kcplus;
	double KCS;			//Mitochondria Model
	double kf1;			//Mitochondria Model
	double kfAAT;		//Mitochondria Model
	double kfACO;		//Mitochondria Model
	double Kfb;
	double kfFH;		//Mitochondria Model
	double KFHeq;		//Mitochondria Model
	double kfSL;		//Mitochondria Model
	double kh_1;		//Mitochondria Model
	double kh_2;		//Mitochondria Model
	double Kh1;			//Mitochondria Model
	double Kh2;			//Mitochondria Model
	double Kh3;			//Mitochondria Model
	double Kh4;			//Mitochondria Model
	double khtrpn_minus;
	double khtrpn_plus;
	double Ki_AM;
	double Ki_prime_SR;
	double Ki_SR;
	double Ki1AD_NaK;
	double KiADP_CaP;
	double kIDH;		//Mitochondria Model
	double KidhNADH;	//Mitochondria Model
	double KiFUM;		//Mitochondria Model
	double Kioaa;		//Mitochondria Model
	double KiOxaa;		//Mitochondria Model
	double kKGDH;		//Mitochondria Model
	double kltrpn_plus;	//Force Model
	double kltrpn_minus;//Force Model
	double Km1AT_NaK;
	double Km1ATP_CaP;
	double Km2ATP_CaP;
	double KmAcCoA;		//Mitochondria Model
	double Kmal;		//Mitochondria Model
	double KmATP_AM;
	double KmATP_SR;
	double KmCMDN;
	double KmCa;
	double KmCSQN;
	double kMDH;		//Mitochondria Model
	double Kmg;			//Mitochondria Model
	double KmIDNAD;		//Mitochondria Model
	double Kmiso;		//Mitochondria Model
	double KmKG;		//Mitochondria Model
	double KmKGNAD;		//Mitochondria Model
	double KmKo;
	double KmmNAD;		//Mitochondria Model
	double KmNa;
	double KmNai;
	double KmnsCa;
	double KmOaa;		//Mitochondria Model
	double KmpCa;
	double KmSucc;		//Mitochondria Model
	double Kna;			//Mitochondria Model
	double kNaCa;
	double Knca;		//Mitochondria Model
	double Ko;
	double Koff;		//Mitochondria Model
	double Krb;
	double kres;		//Mitochondria Model
	double kresf;		//Mitochondria Model
	double ksat;
//	double Fcik1;
	double kSDH;		//Mitochondria Model
	double KSLeq;		//Mitochondria Model
	double KSR;
	double ktrans;		//Mitochondria Model
	double kTrop_pn;	//Force Model
	double L;			//Mitochondria Model
	double LTRPNtot;	//Force Model
	double mcoop;
	double Mg;			//Mitochondria Model
	double n;			//Mitochondria Model
	double na;			//Mitochondria Model
	double Nao;
	double ncoop;
	double Nfb;
	double nID;			//Mitochondria Model
	double nKG;			//Mitochondria Model
	double norm_freq;	//Brandes and Bers parameters
	double Nrb;
	double omega;
	double p1;			//Mitochondria Model
	double p2;			//Mitochondria Model
	double p3;			//Mitochondria Model
	double pa;			//Mitochondria Model
	double pb;			//Mitochondria Model
	double pc1;			//Mitochondria Model
	double pc2;			//Mitochondria Model
	double PCa;
	double period;		//Assigned in some run modes
	double PESI;			//IFmode Setup
	double Pi;			//Mitochondria Model
	double PK;
	double PnsK;
	double PnsNa;
	double pulse_amplitude;
	double pulse_duration;
	double r1;			//Mitochondria Model
	double r2;			//Mitochondria Model
	double r3;			//Mitochondria Model
	double ra;			//Mitochondria Model
	double rb;			//Mitochondria Model
	double rc1;			//Mitochondria Model
	double rc2;			//Mitochondria Model
	double refrac_buffer;	//IFmode Setup
	double rhoF1;		//Mitochondria Model
	double rhoREF;		//Mitochondria Model
	double rhoREN;		//Mitochondria Model
	double shift;
	double t1;			//Brandes and Bers parameters
	double t2;			//Brandes and Bers parameters
	double t3;			//Brandes and Bers parameters, stop time
	double tautr;
	double tauxfer;
	double time_off_Is1;//Assigned in some run modes
	double time_off_Is2;
	double time_on_Is1;	//Assigned in some run modes
	double time_on_Is2;
	double time_vclamp_off;
	double time_vclamp_on;
	double V_AM_scaler;
	double V_AM_max;
	double v1;
	double vclamp_hold;
	double vclamp_set;
	double VJSR;
	double vmaxf;
	double vmaxr;
	double VmDT;		//Mitochondria Model
	double VmNC;		//Mitochondria Model
	double Vmuni;		//Mitochondria Model
	double Vmyo;
	double VNSR;
	double VSS;
	double zeta;		//Force Model

	double f_xb;
	double SL;
	double gmin_xb;

	double stopTime;	//From parameter file, original stop_time value;
	double numRun;			//IF model parameter

	//Calculated Model Parameters
	double inv_keq;							//CK Mod
	double Vmito;
	double Vtotal;
	double g_01_off_mod;//Force Model
	double g_01_mod;	//Force Model
	double g_12_mod;	//Force Model
	double g_23_mod;	//Force Model
	double Ktrop_half;	//Force Model
	double Ntrop;		//Force Model
	double fnormmax;	//Force Model
	double fnormmax2;	//Force Model
	double alpha_SL;	//Force Model
	double f_01;		//Force Model
	double f_12;		//Force Model
	double f_23;		//Force Model
	double zeta_alpha_SL_fnormmax;	//Dependent Variables

	double one_inv_KACOeq;					//Mitochondria Model
	double inv_11;
	double inv_11p1;
	double neg_inv_13;
	double inv_5p98;
	double inv_6;
	double inv_6p8;
	double inv_7p5;
	double inv_9p5;
	double inv_ATPi;
	double inv_bL;
	double inv_C_m;
	double inv_Cmito;						//Mitochondria Model
	double inv_ICahalf;
	double inv_KADP;						//Mitochondria Model
	double inv_KaCa;						//Mitochondria Model
	double inv_kact;						//Mitochondria Model
	double inv_Kfb;
	double inv_Ki_prime_SR;
	double inv_Ki1AD_NaK;
	double inv_KiADP_CaP;
	double inv_KidhNADH;					//Mitochondria Model
	double inv_KiOxaa;						//Mitochondria Model
	double inv_KmNai;
	double KmNaiP1p5;
	double inv_Krb;
	double inv_ktrans;						//Mitochondria Model
	double inv_LTRPNtot_Ktrop_half;			//Force Model
	double inv_tautr;
	double inv_tauxfer;
	double tenDiv9;							//Mitochondria Model
	double FaradayE3;
	double twoThirds;								//Force Model
	double two_b;								//Mitochondria Model
	double Pca_4En3;
	double FRT2;							//Mitochondria Model
	double PKFe3;
	double Acap_Vmyo_F;
	double Acap_VSS_F;
	double ADP;
	double alpha_SL_fnormmax;				//Force Model
	double alpha_SL_fnormmax2;				//Force Model
	double AREF;							//Mitochondria Model
	double b_05;							//Mitochondria Model
	double Cao_341;
	double CMDNtot_KmCMDN;
	double Co;
	double CoA_KSLeq;						//Mitochondria Model
	double CSQNtot_KmCSQN;
	double DmuH_Constant;					//Mitochondria Model
	double E_Ca_Cai_Min;
	double eta_1;
	double exp_AREF_FRT;					//Mitochondria Model
	double exp_3_FRT_Dpsio;					//Mitochondria Model
	double exp_6_FRT_Dpsio;					//Mitochondria Model
	double f_23_g_12_mod;					//Force Model
	double F_over_RT;
	double FRT_3;							//Mitochondria Model
	double FRT_6_g;							//Mitochondria Model
	double G_K1;		
	double G_Ks;
	double high_freq_hz;
	double hm_F_over_RT;					//Mitochondria Model
	double hNa_HAlpha_C1;
	double hNa_HBeta_C1;
	double ICamax_LHospital;
	double INaKmax_Ko_Ko_KmKo;
	double kcnsASP_KAATeq_kfAAT;			//Mitochondria Model
	double kf1_Pi;							//Mitochondria Model
	double KfAAT_GLU;						//Mitochondria Model
	double KfAAT_KAATeq;					//Mitochondria Model
	double kfFH_KFHeq;						//Mitochondria Model
	double kIDH_EtID;						//Mitochondria Model
	double kKGDH_EtKG;						//Mitochondria Model
	double Kmal_Kioaa;						//Mitochondria Model
	double KmATP_AM_Ki_AM;
	double KmATP_SR_Ki_SR;
	double KmCa_Cao;
	double KmCa_Cao_ksat;
	double kMDH_Fh_EtMD;					//Mitochondria Model
	double KmKGNAD_KmIDNAD;					//Mitochondria Model
	double KmnsCa_p3;
	double KmSucc_KiFUM;					//Mitochondria Model
	double kres_sq_KmIDNAD;					//Mitochondria Model
	double kSDH_EtSDH;						//Mitochondria Model
	double kTrop_pn_f_01;					//Force Model
	double kTrop_pn_f_12_g_01_mod;			//Force Model
	double Mg_Kmg_1;						//Mitochondria Model
	double Mg_Kmg_1_Kca;					//Mitochondria Model
	double Nao_p3;
	double norm_freq_hz;
	double p1_exp_3_FRT_Dpsio;				//Mitochondria Model
	double pa_pb_3;							//Mitochondria Model
	double pa_300;							//Mitochondria Model
	double r1_exp_6_FRT_Dpsio;				//Mitochondria Model
	double r2_r3_exp_AREF_FRT;				//Mitochondria Model
	double ra_rc1_exp_6_FRT_Dpsio;			//Mitochondria Model
	double ra_exp_AREF_FRT;					//Mitochondria Model
	double ra_rb;							//Mitochondria Model
	double ra_rc2_exp_AREF_FRT;				//Mitochondria Model
	double rhoRen_6_ra;						//Mitochondria Model
	double rhoRen_6_ra_rb;					//Mitochondria Model
	double rhoREN_ra;						//Mitochondria Model
	double rhoREN_ra_rc1_exp_6_FRT_Dpsio;	//Mitochondria Model
	double rhoREN_rc2;						//Mitochondria Model
	double RTF_05;
	double RTF_05_log_Cao;
	double RTF_log_Ko;
	double RTF_log_Nao;
	double RTF_log_Nao_Ko;
	double sigma;
	double V_30;
	double V_AM_scaler_max_1_f_01_12_23;
	double VAAT_Constant;					//Mitochondria Model
	double VATPase_C1;						//Mitochondria Model
	double VCS_C1;							//Mitochondria Model
	double VFO_C1;							//Mitochondria Model
	double VFO_VHFe_C1;						//Mitochondria Model
	double VIDH_Constant;					//Mitochondria Model
	double VJSR_VNSR;
	double VJSR_VSS;
	double VmDT_75;							//Mitochondria Model
	double VmDT_20;							//Mitochondria Model
	double Vmuni_ktrans;					//Mitochondria Model
	double Vmyo_1000_F_Acap_C_m;
	double Vmyo_VNSR;
	double Vmyo_VSS;

	//Shared variables
	double KmIDNAD_NAD;			//Mitochondria Model
	double exp_FRT_6_g_DmuH;	//Mitochondria Model
	double FRT2_Dpsi;			//Mitochondria Model
	double kltrpn_plus_Cai;		//Force Model
	double start_time_shift;
	double start_time_shift_time_on;
	double V_E_K;
	double VFsq_over_RT;
	double exp2VFRT;
	double VF_over_RT;
	double exp_VF_over_RT;
	double O1_RyR;
	double P1_N1_P2_P3;	//Force Model + precalcs

	//Function Assignments
	double Istim;
	double E_Na;
	double E_K;
	double E_Ks;
	double E_Ca;
	double INa;
	double IKs;
	double IK1;
	double INab;
	double IKp;
	double ICamax;
	double ICa;
	double ICaK;
	double INaK;
	double INaCa;
	double ICab;
	double IpCa;
	double InsK;
	double InsNa;
	double InsCa;
	double V_AM;
	double ATPm;		//Mitochondria Model
	double DmuH;		//Mitochondria Model
	double NAD;			//Mitochondria Model
	double VCS;			//Mitochondria Model
	double VACO;		//Mitochondria Model
	double VIDH;		//Mitochondria Model
	double VKGDH;		//Mitochondria Model
	double VSL;			//Mitochondria Model
	double VSDH;		//Mitochondria Model
	double VFH;			//Mitochondria Model
	double VMDH;		//Mitochondria Model
	double VAAT;		//Mitochondria Model
	double VNO;			//Mitochondria Model
	double VHNe;		//Mitochondria Model
	double VHFe;		//Mitochondria Model
	//double VFO;		//Mitochondria Model
	double VATPase;		//Mitochondria Model
	double Vhu;			//Mitochondria Model
	double VANT;		//Mitochondria Model
	double Vhleak;		//Mitochondria Model
	double Jup;
	double Jrel;
	double Jtr;
	double Jxfer;
	double FN_Ca;		//Force Model
	//double Fnorm;		//Force Model
	//double Force;		//Force Model
	double Jtrpn;
	double beta_SS;
	double beta_JSR;
	double beta_i;
	double Vuni;		//Mitochondria Model
	double VnaCa;		//Mitochondria Model
	double stop_time;	//BB / IF mode setup
};