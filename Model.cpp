
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

#include "model.h"

Model::Model(void) { 
	errweight = NULL;
	s = NULL;
	dS = NULL;

	initializeModel();
}

// So the model can be reset midway through
void Model::initializeModel(bool algebraicMode) {
	if ( errweight != NULL ) delete[] errweight;
	if ( s != NULL ) delete[] s;
	if ( dS != NULL ) delete[] dS;

	modelTypeLookup["Coupled"] = Coupled;
	modelTypeLookup["Mitochondria"] = Mitochondria;
	modelTypeLookup["Force"] = Force;
	stimulasModeLookup["BB"] = BB;
	stimulasModeLookup["IF"] = IF;
	stimulasModeLookup["periodic"] = periodic;
	stimulasModeLookup["none"] = none;

    errweight	= new double[MAXSTATES];
    s			= new double[MAXSTATES];
    dS			= new double[MAXSTATES];
	
	linkStatesReferences_Full(dS, s);

	usingAlgebraicMembranePotential = algebraicMode;

	//Initial/none Current stimulation level
	Istim = 0;

	//Error weights for integrators
	errweight[index_V] = 1E-2;			//46
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
	errweight[index_ASP]= 1.0;
	errweight[index_ASP]= 1.0;
	errweight[index_ASP]= 1.0;
	errweight[index_ASP]= 1.0;
	//CK mod
	errweight[index_ATPi_cyto] = 1.0/8.0;
	errweight[index_CrPi_cyto] = 1.0/40.0;
	errweight[index_CrPi_mito] = 1.0/40.0;

	//Setup the Dependent (i.e. current, outputs)
	dependentVariableIndexLookup["INa"]		= 1;
	dependentVariableIndexLookup["IKs"]		= 2;
	dependentVariableIndexLookup["IK1"]		= 3;
	dependentVariableIndexLookup["INab"]	= 4;
	dependentVariableIndexLookup["IKp"]		= 5;
	dependentVariableIndexLookup["ICa"]		= 6;
	dependentVariableIndexLookup["ICaK"]	= 7;
	dependentVariableIndexLookup["INaK"]	= 8;
	dependentVariableIndexLookup["RNaK"]	= 9;
	dependentVariableIndexLookup["InsCa"]	= 10;
	dependentVariableIndexLookup["INaCa"]	= 11;
	dependentVariableIndexLookup["ICab"]	= 12;
	dependentVariableIndexLookup["IpCa"]	= 13;
	dependentVariableIndexLookup["RpCa"]	= 14;
	dependentVariableIndexLookup["Jup"]		= 15;
	dependentVariableIndexLookup["Jrel"]	= 16;
	dependentVariableIndexLookup["Jtr"]		= 17;
	dependentVariableIndexLookup["Jxfer"]	= 18;
	dependentVariableIndexLookup["V_AM"]	= 19;
	dependentVariableIndexLookup["VNO"]		= 20;
	dependentVariableIndexLookup["VANT"]	= 21;
	dependentVariableIndexLookup["VATPase"]	= 22;
	dependentVariableIndexLookup["VnaCa"]	= 23;
	dependentVariableIndexLookup["Vuni"]	= 24;
	dependentVariableIndexLookup["VIDH"]	= 25;
	dependentVariableIndexLookup["VKGDH"]	= 26;
	dependentVariableIndexLookup["FN_Ca"]	= 27;
	dependentVariableIndexLookup["Fnorm"]	= 28;
	dependentVariableIndexLookup["Force"]	= 29;
}

Model::~Model(void) {
	delete[] errweight;
	delete[] s;
	delete[] dS;
}

//Should be a simpler way to do all this, i.e. too many function calls.
//This is a wrapper to perform getF with arrays and not N_Vectors
//In theory (check documentation to see if change is made) y, and ydot correspond to s, dS
void Model::F_GPC(double start_time) {
	
	//Because the code base changes this (a hack by another person)
	//We have to save and restore the value
	sharedVariableUpdate_Force();
	
	if (usingAlgebraicMembranePotential)
		calculateAlgebraicMembranePotential();

	sharedVariableUpdate_GPC_1(start_time);

	calculateStimulatedCurrent(start_time);
	if (clampVoltage) 
		calculateClampedVoltage();
	calculateReversalPotentials();

	sharedVariableUpdate_GPC_2();

	calculateINa();
	calculateIKs();
	calculateIK1();
	calculateINab();
	calculateIKp();
	calculateICa();
	calculateICaK();
	calculateINaK();
	calculateINaCa();
	calculateICab();
	calculateIpCa();
	calculateInsCa();
	/*after ATP mod**/
	calculateV_AM();
	/**after Mitochondria mod**/
	calculateATPm();
	calculateDmuH();
	calculateNAD();

	sharedVariableUpdate_Mitochondria();

	calculateVCS();
	calculateVACO();
	calculateVIDH();
	calculateVKGDH();
	calculateVSL();
	calculateVSDH();
	calculateVFH();
	calculateVMDH();
	calculateVAAT();
	calculateVNO_VHNe();
	calculateVFO_VHFe();
	calculateVATPase_Vhu();
	calculateVANT_Vhleak();

	//Starting Derivatives
	calculateFNa();
	calculateFxKs();
	calculateInCalFlux(); //Not derivative
	calculateForce();
	calculateF_trpmyo();
	calculateFLTRPNCa();
	calculateFHTRPNCa();
	calculateJtrpn();

	calculateVuni();			// moved up
	calculateVnaCa();			// from below (in the Mitene folder)

	//CHF()
	if (usingCHF) {
	    IK1    = chfsc_IK1 * IK1;
	    Jup    = chfsc_Jup * Jup;
	    INaCa  = chfsc_INaCa * INaCa;
	}

	calculateFNai();
	calculateFKi();
	calculateFCai();
	calculateFCaSS();
	calculateFCaJSR();
	calculateFCaNSR();
	calculateFV();		//Don't need to do in algebriac case
	calculateFRyR();
	calculateFCaL();
	calculateFyCa();
	calculateFOCa();
	/**after ATP mod**/
	if(usingCK) {
		calculateCK();		//CK Mod
	} else {
		calculateFATPi();	//No CK mod
	}
	calculateF_mitene();

}
//Model F function for the mitochondria model
void Model::F_Mitochondria() {
	//Shared updated variables

	//Calculate Intemediates
	calculateATPm();
	calculateDmuH();
	calculateNAD();

	sharedVariableUpdate_Mitochondria();

	calculateVCS();
	calculateVACO();
	calculateVIDH();
	calculateVKGDH();
	calculateVSL();
	calculateVSDH();
	calculateVFH();
	calculateVMDH();
	if (usingASP)
		calculateVAAT_ASP();
	else
		calculateVAAT();
	calculateVNO_VHNe();
	calculateVFO_VHFe();
	calculateVATPase_Vhu();
	calculateVANT_Vhleak();
	calculateVuni();
	calculateVnaCa();
	
	//Put them all together
	calculateF_mitene();
}
//Model F for the Force Model
void Model::F_Force() {
	sharedVariableUpdate_Force();
	calculateForce();
	calculateF_trpmyo();
	calculateFLTRPNCa();
}
//Common variable computations for the GPC Model Part 1
//Must be done after calculateAlgebraicMembranePotential()
void Model::sharedVariableUpdate_GPC_1(double start_time) {
	//Before IpCa
	ADP = 8.0 - (*ATPi);
	inv_ATPi = 1.0 / (*ATPi);

	//Done after V -> Algebraic Potential, before ICA
	VF_over_RT = (*V) * F_over_RT;
	exp_VF_over_RT = exp( VF_over_RT );
	VFsq_over_RT = FaradayE3 * VF_over_RT;	//VF_over_RT is unitless mV/mV
	exp2VFRT = exp_VF_over_RT * exp_VF_over_RT;

	//Done before calculateStimulatedCurrent
	start_time_shift = start_time - shift;
	start_time_shift_time_on = floor( (start_time_shift) / period ) * period;

	//Intermediate Variables
	O1_RyR = 1.0 - *C1_RyR - *C2_RyR - *O2_RyR;

	V_30 = (*V) + 30.0;
}
//Common variable computations for the GPC Model Part 2
//Done after calculateReversalPotentials and before calculateIk1
void Model::sharedVariableUpdate_GPC_2() {
	//Done After calculateReversalPotentials, before Ik1
	V_E_K = (*V) - E_K;
}

//Common variable computations for the Mitochondria Model
void Model::sharedVariableUpdate_Mitochondria() {
	//Must be done after calculateNAD is Called
	KmIDNAD_NAD = KmIDNAD / NAD;
	//After calculateDmuH
	exp_FRT_6_g_DmuH = exp(FRT_6_g * DmuH);
	//Anytime
	FRT2_Dpsi = FRT2 * ( (*Dpsi) - 91.0 );
}
//Common variable computations for the Force Model, well if there were any
void Model::sharedVariableUpdate_Force() {}
//Perform pre-Computations
void Model::calculateConstantModelParameters() {
	//Moved from reading.cpp
	Vtotal = (Vmyo + VJSR + VNSR + VSS) / 0.64;
	Vmito  = Vtotal * 0.36;	//assuming Vmito is 36% of the cellular volume
	
	f_01 = 3.0  * f_xb;
	f_12 = 10.0 * f_xb;
	f_23 = 7.0  * f_xb;

	//g_xbSL=gmin_xb*(1-pow((1-SLnorm),1.6));
	double g0_01 = 1.0 * gmin_xb;
	double g0_12 = 2.0 * gmin_xb;
	double g0_23 = 3.0 * gmin_xb;

	double paths = g0_01 * g0_12 * g0_23 + f_01 * g0_12 * g0_23 + f_01 * f_12 * g0_23 + f_01 * f_12 * f_23;
	double P1max = (f_01 * (2.0 * gmin_xb) * (3.0 * gmin_xb)) / paths;
	double P2max = (f_01 * f_12 * (3.0 * gmin_xb)) / paths;
	double P3max = (f_01 * f_12 * f_23) / paths;
	double Fmax  = P1max + 2.0 * P2max + 3 * P3max;
	fnormmax = Fmax / 3.0;
	double SLnorm = (SL - 1.7) / 0.6;

	double Ktrop_Ca = kltrpn_minus / kltrpn_plus;
	
	Ktrop_half = 1.0 / (1.0 + (Ktrop_Ca / (1.7 / 1000.0 + ((0.9 / 1000.0 - 1.7 / 1000.0) / (2.3 - 1.7)) * (SL - 1.7))));
	Ntrop = 3.5 * SL - 2.0;

	 //Real model in Fortran
	fnormmax2 = (P1max + P2max + P3max);
	
	double La = 1;			 // um
	double Lm_prime = 1.5;	 // um
	double Lz = 0.1;			 // um
	double Lb = 0.1;			 // um
	double Lm = Lm_prime - Lb;	 // um

	double mod_factor = 1 + (2.3 - SL) / pow((2.3 - 1.7),1.6);  // dimensionless

	if (SL <2.2)
			alpha_SL = __min(1.0,(SL - 2.0 * La + Lm_prime - Lz) / Lm);	 // dimensionless
	else
			alpha_SL = 1 - (SL - 2.2) / Lm;

	g_01_mod = g0_01 * mod_factor;	 // 1 / ms
	g_12_mod = g0_12 * mod_factor;	 // 1 / ms
	g_23_mod = g0_23 * mod_factor;	 // 1 / ms
	double g_01_off = 30.0 / 1000.0;			 // New variable -  - >unit: (1 / ms)
	g_01_off_mod = g_01_off * mod_factor;	 // 1 / ms

	//Optimization Precalcs
	Vmyo_1000_F_Acap_C_m = Vmyo * ( Faraday * 1000) / (Acap * C_m);
	Co = Vtotal * ( 2 * Cao + Ko + Nao) / Vmyo;
	high_freq_hz = ( 1.0 / high_freq ) * 1000; //(ms)
	norm_freq_hz = ( 1.0 / norm_freq ) * 1000; //(ms)
	RTF_05			= RT_over_F * 0.5;
	RTF_log_Nao		= RT_over_F * log(Nao);
	RTF_log_Ko		= RT_over_F * log(Ko);
	RTF_log_Nao_Ko	= RT_over_F * log(Ko + 0.01833 * Nao);
	RTF_05_log_Cao	= RTF_05 * log(Cao);
	//When Cai is clamped to __min 1e-10
	E_Ca_Cai_Min	= RTF_05_log_Cao + 10 * RTF_05; 
	G_Ks = 0.282 * sqrt( Ko / 5.4 );
	G_K1 = 0.75 * sqrt( Ko / 5.4 );		
	inv_5p98 = 1.0 / 5.98;
	FaradayE3 = (1000.0*Faraday);
	Cao_341 = Cao * 341;
	//scale by 1000, so current is in uA/uf
	ICamax_LHospital = 2.0 * PCa * 1000.0 * Faraday *(1.0 - 341.0 * Cao);
	Pca_4En3 =  4.0E-3 * PCa;
	F_over_RT = 1.0 / RT_over_F;
	inv_ICahalf = 1.0 / ICahalf;
	PKFe3 = FaradayE3 * PK;
	sigma = 0.0365 * ( exp(Nao / 67.3)-1.0)/7.0;
	inv_KmNai = 1.0 / KmNai;
	KmNaiP1p5 = sqrt(KmNai*KmNai*KmNai); //KmNai^1.5
	INaKmax_Ko_Ko_KmKo = INaKmax * Ko / (Ko + KmKo);
	inv_Ki1AD_NaK = 1.0 / Ki1AD_NaK;
	eta_1 = eta-1.0;
	Nao_p3 = pow( Nao, 3.0) / Cao;
	KmCa_Cao = (KmCa + Cao) * (pow(KmNa, 3.0) + pow(Nao, 3.0)) / kNaCa /Cao;
	KmCa_Cao_ksat = KmCa_Cao * ksat;
	inv_KiADP_CaP = 1.0 / KiADP_CaP;
	KmnsCa_p3 = pow(KmnsCa, 3.0);
	V_AM_scaler_max_1_f_01_12_23 = V_AM_scaler * V_AM_max / (f_01 + f_12 + f_23);
	KmATP_AM_Ki_AM = KmATP_AM / Ki_AM;

	//Mitochondira Constants : ATPm
	DmuH_Constant = -2.303 * RT_over_F * DpH;
	VCS_C1 = (KCS * EtCS * AcCoA) / (KmAcCoA + AcCoA);
	one_inv_KACOeq = 1.0 + 1.0 / KACOeq;
	VIDH_Constant = 1.0 + H / kh_1 + kh_2 / H;
	kIDH_EtID = kIDH * EtID;
	inv_KADP = 1.0 / KADP;
	inv_KaCa = 1.0 / KaCa;
	inv_KidhNADH = 1.0 / KidhNADH;
	KmKGNAD_KmIDNAD = KmKGNAD / KmIDNAD;
	Mg_Kmg_1 = Mg / Kmg + 1.0;
	Mg_Kmg_1_Kca = Mg_Kmg_1 / Kca;
	kKGDH_EtKG = kKGDH * EtKG;
	CoA_KSLeq = CoA / KSLeq;
	kSDH_EtSDH = kSDH * EtSDH;
	KmSucc_KiFUM = KmSucc / KiFUM;
	inv_KiOxaa = 1.0 / KiOxaa;
	kfFH_KFHeq = kfFH / KFHeq;
	kMDH_Fh_EtMD = pow( ( 1.0 / ( 1.0 + Kh3 / H + Kh3 * Kh4 / pow(H,2) ) ) ,2) *
		                ( 1.0 / ( 1.0 + H / Kh1 + pow(H,2) / ( Kh1 * Kh2 ) ) + Koff) *
						  kMDH * EtMD;
	Kmal_Kioaa = Kmal / Kioaa;
	VAAT_Constant = kfAAT * GLU * kcnsASP * KAATeq / kfAAT;
	kcnsASP_KAATeq_kfAAT = kcnsASP * KAATeq / kfAAT;
	KfAAT_GLU = kfAAT * GLU;
	KfAAT_KAATeq = kfAAT/KAATeq;
	kres_sq_KmIDNAD = kres * kres / KmIDNAD;
	exp_6_FRT_Dpsio = exp( 6.0 * Dpsio * F_over_RT);
	FRT_6_g = 6.0 * g * F_over_RT;
	ra_rc1_exp_6_FRT_Dpsio = ra + rc1 * exp_6_FRT_Dpsio;
	r1_exp_6_FRT_Dpsio = r1 * exp_6_FRT_Dpsio;
	rhoREN_ra_rc1_exp_6_FRT_Dpsio = 0.5 * rhoREN * ra_rc1_exp_6_FRT_Dpsio;
	rhoREN_rc2 = 0.5 * rhoREN * rc2;
	rhoREN_ra = 0.5 * rhoREN * ra;
	rhoRen_6_ra = 6.0 * rhoREN * ra;
	rhoRen_6_ra_rb = 6.0 * rhoREN * ( ra + rb );
	AREF = RT_over_F * log ( kresf * sqrt ( FADH2 / FAD ) );
	exp_AREF_FRT = exp( AREF * F_over_RT );
	ra_rc2_exp_AREF_FRT = 0.5 * (ra + rc2 * exp_AREF_FRT);
	VFO_C1 = ( ra + rc1 * exp_6_FRT_Dpsio ) * exp_AREF_FRT * 0.5;
	ra_exp_AREF_FRT = 4.0 * ra * exp_AREF_FRT;
	ra_rb = 4.0 * (ra + rb);
	VFO_VHFe_C1 = ( 1.0 + r1 * exp_AREF_FRT) * exp_6_FRT_Dpsio;
	r2_r3_exp_AREF_FRT = r2 + r3 * exp_AREF_FRT;
	exp_3_FRT_Dpsio = exp( 3.0 * Dpsio * F_over_RT);
	FRT_3 = 3.0 * F_over_RT;
	kf1_Pi = kf1 / Pi;
	VATPase_C1 = ( 100.0 * pa + pc1 * exp_3_FRT_Dpsio);
	pa_pb_3 = 3.0 * (pa + pb);
	pa_300 = 300.0 * pa;
	p1_exp_3_FRT_Dpsio = p1 * exp_3_FRT_Dpsio;
	hm_F_over_RT = hm * F_over_RT;
	VmDT_75 = 0.75 * VmDT;
	VmDT_20 = 20.0 * VmDT;
	tenDiv9 = 10.0 / 9.0;
	//Pause doing Mitochondria : VHLeak
	
	inv_11 = 1.0 / 11.0;
	inv_11p1 = -1.0 / 11.1;
	inv_6p8 = -1.0 / 6.8;
	hNa_HAlpha_C1 =  0.135*exp(-80.0/6.8);
	hNa_HBeta_C1 = 0.13*exp( - 10.66 / 11.1);
	inv_Kfb = 1.0 / Kfb;
	inv_Krb = 1.0 / Krb;
	inv_tautr = 1.0 / tautr;
	inv_tauxfer = 1.0 / tauxfer;
	KmATP_SR_Ki_SR = KmATP_SR / Ki_SR;
	inv_Ki_prime_SR = 1.0 / Ki_prime_SR;

	//Start Force Model : Force
	alpha_SL_fnormmax2 = alpha_SL / fnormmax2;
	alpha_SL_fnormmax  = alpha_SL / fnormmax / 3.0;
	inv_LTRPNtot_Ktrop_half = 1.0 / ( LTRPNtot * Ktrop_half );
	kTrop_pn_f_01 = -kTrop_pn - f_01;
	kTrop_pn_f_12_g_01_mod = -(kTrop_pn + f_12 + g_01_mod);
	f_23_g_12_mod = -(f_23 + g_12_mod);
	twoThirds = 2.0 / 3.0;
	//End Force Model : FLTRPNCa

	CMDNtot_KmCMDN = CMDNtot * KmCMDN;
	CSQNtot_KmCSQN = CSQNtot * KmCSQN;

	//Start Mitochondria Again : Vuni
	inv_ktrans = 1.0 / ktrans;
	inv_kact = 1.0 / kact;
	Vmuni_ktrans = Vmuni / ktrans;
	FRT2 = 2.0 * F_over_RT;
	b_05 = b * 0.5;
	//Pause Mitochondria Again : Vanca

	Acap_Vmyo_F = Acap / ( Vmyo * Faraday * 1000.0);
	Acap_VSS_F = Acap / ( 2.0 * VSS * Faraday * 1000.0);
	VJSR_VSS = VJSR / VSS;
	Vmyo_VSS = Vmyo / VSS;
	Vmyo_VNSR = Vmyo / VNSR;
	VJSR_VNSR = VJSR / VNSR;
	inv_C_m = 1.0 / C_m;
	neg_inv_13 = -1.0 / 13;
	inv_bL = 1.0 / bL;
	inv_7p5 = 1.0 / 7.5;
	inv_6 = 1.0 / 6.0;
	inv_9p5 = 1.0 / 9.5;

	//Mitochondria Last Time : F_Mitene
	inv_Cmito = 1.0 / Cmito;
	two_b = 2.0 * b;

	//CK Mod
	inv_keq = 1.0 / keq;

	//Dependent Variables
	zeta_alpha_SL_fnormmax = zeta * alpha_SL_fnormmax;
}
// algebraic membrane potential method; Currently has some issues 
void Model::calculateAlgebraicMembranePotential() {
	//Parameters: Faraday, Acap, C_m, Vtotal, Cao, Ko, Nao, Vmyo, CMDNtot, KmCMDN, VJSR, VNSR, CSQNtot, KmCSQN, VSS, Vmito
	double Ca_all = 2.0 * ( (*Cai) + (*Cai) * CMDNtot / ( (*Cai) + KmCMDN ) + (*LTRPNCa) + (*HTRPNCa) ) + 
					VJSR * (1.0 + CSQNtot / ( (*CaJSR) + KmCSQN )) * (*CaJSR) + 
					VNSR * (*CaNSR) + VSS * (*CaSS) + Vmito * (*Cam);
	(*V) = Vmyo_1000_F_Acap_C_m * ( (*Nai) + (*Ki) + Ca_all  - Co ); // + extra);
}
// Cell Stimulation, BB, IF, periodic and none modes
void Model::calculateStimulatedCurrent(double start_time)
{
	//Parameters: high_freq, norm_freq, start_time, shift, time_on_Is2, time_off_Is2, pulse_amplitude, t1, t2
	//Assigned: Istim, [period, time_on_Is1, time_off_Is1]
	switch (stimulasMode)
	{
		case IF:
			if ( ((start_time >= time_on_Is1) && (start_time <= time_off_Is1)) ||
				( (start_time >= time_on_Is2) && (start_time <= time_off_Is2) ) )
				Istim = pulse_amplitude;
			else
				Istim = 0;
			break;
		case BB:
			start_time_shift = start_time - shift;
			if ( (start_time_shift >= t1) && (start_time_shift <= t2) ) 
				period = high_freq_hz;
			else
				period = norm_freq_hz;

			time_on_Is1  = start_time_shift_time_on;
			time_off_Is1 = time_on_Is1 + pulse_duration;
			if ( ( (start_time_shift) >= time_on_Is1) && 
				 ( (start_time_shift) <= time_off_Is1))
				Istim = pulse_amplitude;
			else
				Istim = 0;
			break;
		case periodic:
			time_on_Is1  = start_time_shift_time_on;
			time_off_Is1 = time_on_Is1 + pulse_duration;
			Istim = 0;
			if ( ( (start_time_shift) >= time_on_Is1) && 
				 ( (start_time_shift) <= time_off_Is1))
				Istim = pulse_amplitude;
			else
				Istim = 0;
			break;
		case none:
			//Change this if mode changes
//			time_on_Is1=floor((start_time-shift)/period)*period;
//			time_off_Is1=time_on_Is1+pulse_duration;
//			Istim = 0;
			break;
	}
}
//sets the voltage clamp
//Optimizations can be done if 2 or more of these are parameters:
//time_vclamp_on, vclamp_set, vclamp_hold
void Model::calculateClampedVoltage() {
	//Parameters: time_vclamp_on, time_vclamp_off, vclamp_set, vclamp_hold
	bool test1 = start_time_shift >= (start_time_shift_time_on + time_vclamp_on);
    if ( test1 && (start_time_shift <  (start_time_shift_time_on + time_vclamp_off)) ) {           
		double ramp = ( (start_time_shift - start_time_shift_time_on - time_vclamp_on)*0.5 )
     		    *(vclamp_set - vclamp_hold) + vclamp_hold;
        if ( vclamp_hold <= vclamp_set ) 
			(*V) = __min(vclamp_set, ramp); // depol.  steps
        else
			(*V) = __max(vclamp_set, ramp); // hyperpol. steps
	} else if ( !test1 ) {
        (*V) = vclamp_hold;
	} else {
		double ramp = (start_time_shift_time_on + time_vclamp_off - start_time_shift) *
						0.5 * (vclamp_set - vclamp_hold) + vclamp_set;
        if (vclamp_hold <= vclamp_set) 
                (*V) = __max(vclamp_hold, ramp); // depol. step
        else
                (*V) = __min(vclamp_hold, ramp); // hyper. step
	}  
}
// Calculates reversal potentials
//Calcs log(*Nai), log(*Ki), log(Cai),log(*Ki + 0.01833 * *Nai): Optomize?
void Model::calculateReversalPotentials() {
	//Parameters: Nao, Ko, Cao
	//Fixed Precomputes

	E_Na = RTF_log_Nao    - RT_over_F * log( (*Nai) );
	E_K  = RTF_log_Ko     - RT_over_F * log( (*Ki) );
	E_Ks = RTF_log_Nao_Ko - RT_over_F * log( (*Ki) + 0.01833 * (*Nai));

	if ((*Cai) < 1.0e-10) {
//		(*Cai) = 1.0E-10;
		E_Ca = E_Ca_Cai_Min;
		cerr << "calculateReversalPotentials: Below minimum calcium levels." <<endl;
	} else {
		E_Ca = RTF_05_log_Cao - RTF_05 * log( (*Cai) );
	}
}
//The current value of INa, an intermediate variable;E_Na
void Model::calculateINa() {
	//Parameters: G_Na
	INa = G_Na *( (*mNa)*(*mNa)*(*mNa)*(*hNa)*(*jNa)*((*V) - E_Na) );	
	//	INa=G_Na*(pow(mNa,3.0)*hNa*jNa*(V-E_Na));		
}	
//The current value of IKs, Fortran Code; exp((v-40)/40); intermediate Variable
void Model::calculateIKs() {
	//Parameters: Ko
	IKs = G_Ks * (*xKs) * (*xKs) * ( (*V) - E_Ks ) / ( 1.0 + exp(( (*V) - 40.0) / 40.0) );
}
//IK1 current; some math optimizations possible; V-E_K(V_E_K )
void Model::calculateIK1() {
	//Parameters: Ko
	double K1Alpha = 1.02 / ( 1.0 + exp( 0.2385* ( V_E_K - 59.215)));
	double K1Beta  = ( 0.4912 * exp( 0.08032 * ( V_E_K + 5.476)) +
					   exp( 0.06175 * ( V_E_K - 594.31))            ) / 
					   ( 1.0 + exp( -0.5143 * ( (*V) - E_K + 4.753 ) ) );
	double K1_inf = K1Alpha / ( K1Alpha + K1Beta );
	IK1 = G_K1 * K1_inf * (V_E_K);
//	IK1 = 0.68* G_K1 * K1_inf * (V_E_K);
}
//Background Na+ current; intermediate Variable; V-E_Na
void Model::calculateINab() {
	//Parameters: G_Nab
	INab = G_Nab * ( (*V) - E_Na);
}
//IKp Current--plateau current, Math optimizations possible; Intermediate Variable; V_E_K
void Model::calculateIKp() {
	//Parameters: G_Kp
	IKp = G_Kp * (V_E_K) / ( 1.0 + exp( ( 7.488 - (*V) ) * inv_5p98) );
}
//Fortran and JT;VF_over_RT;VFsq_over_RT;exp_VF_over_RT
void Model::calculateICa() {
	//Parameters Pca, Cao
	//Taylor series: x/exp(2x) => 0.5 + -0.5 x, x = VFRT, => 0.5 + -0.02V
	if (fabs(*V) < LHospitalThreshold) { // Use limit V->0, AT
		ICamax = Pca_4En3 * (exp2VFRT - Cao_341) * ( 0.5 - 0.02 * (*V) );
	} else {
		ICamax = Pca_4En3 * VFsq_over_RT * (exp2VFRT - Cao_341)/(exp2VFRT - 1.0);
	}
	//The 5 factor added to account for low open probability of CAY L-type channel.
	ICa = 6.0 * ICamax * (*yCa) * (*Open);	
}
//L-Type Ca channel permeable to K+; Fortran Code;F_over_RT;exp_VF_over_RT;VF_over_RT;PKprime
void Model::calculateICaK() {
	//Parameters: PK, ICahalf, Ko
	//Taylor series: x/exp(2x) => 0.5 + -0.5 x, x = VFRT, => 0.5 + -0.02V
	if (fabs(*V) < LHospitalThreshold) { // Use limit V->0, AT
		ICaK = PKFe3 * ( (*Open) + (*OCa) ) * (*yCa) * ( (*Ki) * exp2VFRT - Ko ) * ( 0.5 - 0.02 * (*V) )
			/ ( 1.0 + ICamax * inv_ICahalf );
	} else {
		ICaK = PKFe3 * ( (*Open) + (*OCa) ) * (*yCa) * ( (*Ki) * exp2VFRT - Ko ) * VF_over_RT
			/ ( ( exp2VFRT - 1.0 ) * ( 1.0 + ICamax * inv_ICahalf ) );
	}
}
//Sodium ATPase pump;ATP Mod;VF_over_RT;exp_VF_over_RT;ADP;inv_ATPi
void Model::calculateINaK() {
	//Parameters: Ko, KmKo, KmNai, INaKmax, Ki1AD_NaK, Km1AT_NaK
	double NaiP1p5 = sqrt((*Nai)*(*Nai)*(*Nai));
	INaK = INaKmax_Ko_Ko_KmKo * NaiP1p5 / 
		( ( NaiP1p5 + KmNaiP1p5)*
		  ( 1.0 + 0.1245 * exp ( -0.1 * VF_over_RT) + sigma / exp_VF_over_RT)*
		  ( 1.0 + (Km1AT_NaK * inv_ATPi) * ( 1.0 + ADP * inv_Ki1AD_NaK ) ));
}
//The sodium-calcium exchanger;VF_over_RT;exp_VF_over_RT
void Model::calculateINaCa() {
	//Parameters eta, Nao, Cao, ksat, KmCa
	double exp_eta_VF_over_RT = exp( eta * VF_over_RT);
	double exp_eta1_VF_over_RT = exp_eta_VF_over_RT / exp_VF_over_RT;
	INaCa =( exp_eta_VF_over_RT * (*Nai) * (*Nai) * (*Nai) - 
		     exp_eta1_VF_over_RT * Nao_p3 * (*Cai) )/ 
		   ( KmCa_Cao + KmCa_Cao_ksat * exp_eta1_VF_over_RT);
}
//ICab;E_Ca
void Model::calculateICab() {
	//G_Cab is a parameter
	ICab = G_Cab * ((*V) - E_Ca);		
}
//Sarcolemmal pump current; After ATP mod;ADP;inv_ATPi
void Model::calculateIpCa() {
	//Parameters: Km2ATP_CaP, KiADP_CaP, Km1ATP_CaP, KmpCa, IpCamax
	IpCa = IpCamax * (*Cai) / ( KmpCa + (*Cai) ) * 
		( 1.0 / ( 1.0+( Km1ATP_CaP * inv_ATPi)*( 1.0 + ADP * inv_KiADP_CaP) ) 
		  + 1.0 / ( 1.0 + ( Km2ATP_CaP*inv_ATPi )));
}
// InsCa added to Guinea pig project (JT);VF_over_RT;exp_VF_over_RT;VFsq_over_RT;exp_VF_over_RT_1
//Do a taylor series expansion on this as soon as possible
void Model::calculateInsCa() {	
	//VFsq_over_RT, exp_VF_over_RT_1, Cai_p3, VFsq_over_RT, KmnsCa_p3, PnsK_75, PnsNa_75, PnsK_Nao_PnsNa_Ko_75
	//Parameters:Nao, PnsNa, KmnsCa, PnsK
	//Taylor series: x/exp(x) => 1 + -0.5 x, x = VFRT, => 1 + -0.02V
	double common;
	double CaiP3 = (*Cai)*(*Cai)*(*Cai);
	if (fabs(*V) < LHospitalThreshold) {// Use limit V->0, AT
		common = 0.75 * CaiP3 * ( 1.0 - 0.02 * (*V) ) / ( CaiP3 + KmnsCa_p3 );
	} else {
		common = 0.75 * CaiP3 * VFsq_over_RT / ( (exp_VF_over_RT - 1.0) * (CaiP3 + KmnsCa_p3) );
	}
		InsNa = PnsNa * common * ((*Nai) * exp_VF_over_RT - Nao);
		InsK  = PnsK  * common * ((*Ki)  * exp_VF_over_RT - Ko);
		InsCa = InsNa + InsK;
}

//After ATP mod;Uses P1,2,3; ADP; inv_ATPi
void Model::calculateV_AM() {
	//Parameters: V_AM_scaler, V_AM_max, f_01, f_12, f_23, KmATP_AM, Ki_AM
	V_AM = V_AM_scaler_max_1_f_01_12_23 * (f_01 * (*P0) + f_12 * (*P1) + f_23 * (*P2))
		 / (1.0 + inv_ATPi * ( KmATP_AM + KmATP_AM_Ki_AM * ADP ));
}
//Mitochondia; Intemediate Variable; Optimize Things that use this
void Model::calculateATPm() {
	ATPm = Cm - (*ADPm);
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this
void Model::calculateDmuH() {
	DmuH = DmuH_Constant + (*Dpsi);
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this
void Model::calculateNAD() {
	NAD = CPN - (*NADH);
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this
void Model::calculateVCS() {
	//Parameter: KmOaa, KCS, EtCS, AcCoA, KmAcCoA
	VCS = VCS_C1 * (*Oaa) / ((*Oaa) +  KmOaa);
}
//Mitochondia;Also CIT; Intemediate Variable;  Optimize Things that use this
void Model::calculateVACO() {
	VACO = kfACO * ( CIK - (*AKG) - (*SCoA) - (*Succ) - (*FUM) - (*MAL) - (*Oaa) - (*ISOC)* one_inv_KACOeq );
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this;1/Isoc;Dependent on NAD;KmIDNAD_NAD
void Model::calculateVIDH() {
	double Fa = 1.0 / (( 1.0 + (*ADPm) * inv_KADP) * (1.0 + (*Cam) * inv_KaCa));
	double Fi = 1.0 + (*NADH) * inv_KidhNADH;
	VIDH = kIDH_EtID /
		  (VIDH_Constant + KmIDNAD_NAD * Fi + pow( Kmiso/(*ISOC), nID) * Fa * (1.0 + KmIDNAD_NAD  * Fi));
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this; Dependent on Calculate NAD(KmIDNAD_NAD); 1/akg
void Model::calculateVKGDH() {		//corrected on 07/26)
	double a = ( (Mg_Kmg_1 + Mg_Kmg_1_Kca*(*Cam)) );
	VKGDH = kKGDH_EtKG * a / ( a + pow(KmKG/(*AKG),nKG) + KmKGNAD_KmIDNAD * KmIDNAD_NAD );//corregir
}
//Mitochondia;Also CIT; Intemediate Variable;  Optimize Things that use this; Dependent on calculateATPm
void Model::calculateVSL() {
	VSL = kfSL * ( (*SCoA) * (*ADPm) - CoA_KSLeq * (*Succ) * ATPm);
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this;1/Succ
void Model::calculateVSDH() {
	//Parameters: kSDH, EtSDH, KmSucc, KiFUM, KiOxaa
	VSDH = kSDH_EtSDH * (*Succ) / ((*Succ) + (KmSucc + KmSucc_KiFUM*(*FUM)) * (1.0 + inv_KiOxaa*(*Oaa)) );
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this;Dependent on NAD(KmIDNAD_NAD)
void Model::calculateVFH() {
	//Parameters: kfFH, KFHeq
	VFH = kfFH * (*FUM) - kfFH_KFHeq * (*MAL);
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this;
void Model::calculateVMDH() {
	//Parameters: Kh1, Kh2, Kh3, Kh4, Koff, H, Kioaa, KmmNAD, Kmal, EtMD, kMDH
	VMDH = kMDH_Fh_EtMD * (*MAL) * NAD / 
		( ( (*MAL) + Kmal + (*Oaa) * Kmal_Kioaa ) * ( KmmNAD + NAD ) );
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this;
void Model::calculateVAAT() {
	//Parameters: kfAAT, GLU, kcnsASP, KAATeq, 
	VAAT = VAAT_Constant * (*Oaa) / ( kcnsASP_KAATeq_kfAAT + (*AKG) );
}
//Mitochondia; Matlab model Calculation of VAAT. Uses ASP
void Model::calculateVAAT_ASP() {
	//Parameters: kfAAT, GLU, KAATeq
    VAAT = KfAAT_GLU * (*Oaa) - KfAAT_KAATeq * (*ASP) * (*AKG);
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this;RT_over_F; F_over_RT;KmIDNAD_NAD;DmuH;exp_6_FRT_Dpsio;exp_FRT_6_g_DmuH;
void Model::calculateVNO_VHNe() {
	//Parameters: kres, rhoREN, Dpsio, g, ra, rc1, r1, rc2, rb
	double AREN =  sqrt ( (*NADH) * kres_sq_KmIDNAD * KmIDNAD_NAD) ;
	double denominator = 1.0 / ( (exp_6_FRT_Dpsio + r1_exp_6_FRT_Dpsio *  AREN ) + ( r2 + r3 *  AREN ) * exp_FRT_6_g_DmuH );
	
	VNO  = ( (rhoREN_ra_rc1_exp_6_FRT_Dpsio + rhoREN_rc2 * exp_FRT_6_g_DmuH) *  AREN  - rhoREN_ra * exp_FRT_6_g_DmuH ) * denominator;
	VHNe = (rhoRen_6_ra * AREN  - rhoRen_6_ra_rb * exp_FRT_6_g_DmuH ) *	denominator;
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this;RT_over_F;exp_6_FRT_Dpsio;exp_FRT_6_g_DmuH;
void Model::calculateVFO_VHFe() {
	//Parameters: kresf, FADH2, FAD, r2, r3, rhoREF
	double denominator = rhoREF / (VFO_VHFe_C1 + r2_r3_exp_AREF_FRT * exp_FRT_6_g_DmuH);
	//Unused
	//VFO  = ( VFO_C1 - ra_rc2_exp_AREF_FRT * exp_FRT_6_g_DmuH) * denominator;
	VHFe = ( ra_exp_AREF_FRT - ra_rb*exp_FRT_6_g_DmuH ) * denominator;
}
//Mitochondia; Intemediate Variable;  Optimize Things that use this; RT_over_F; exp_3_FRT_Dpsio;F_over_RT; 1/ADPm; ATPm
void Model::calculateVATPase_Vhu() {
	//This has a numeric instability due to round off (I guess) dependent on the order of evaluation of variables in VATPase
	// DO NOT combine the exp_3FRT_DmuH statements, breaks the model
	//Parameters: kf1, Pi, pa, pc1, pc2, p1, p2, p3 ,rhoF1
	//Bad naming on some of the parameters
	double exp_3FRT_DmuH = exp( FRT_3 * DmuH );
	double AF1 =  kf1_Pi * ATPm / (*ADPm);
	double denominator = -rhoF1 / ( exp_3_FRT_Dpsio + p1_exp_3_FRT_Dpsio * AF1 + ( p2 + p3 * AF1 ) * exp_3FRT_DmuH);
//	VATPase = ( VATPase_C1 * AF1 - (pa + pc2 * AF1) * exp_3FRT_DmuH ) * denominator;
	VATPase = ( (VATPase_C1 + pc2 * exp_3FRT_DmuH) * AF1 - pa * exp_3FRT_DmuH ) * denominator;
	Vhu     = ( pa_300 + pa_300 * AF1  - pa_pb_3 * exp_3FRT_DmuH ) * denominator;
}
//Mitochondia;Different from paper; Intemediate Variable;  Optimize Things that use this; ATPm;ADP
void Model::calculateVANT_Vhleak() {
	//Parameters: gh, VmDT; hm
	double ATPi_ADP = (*ATPi) / ADP;
	double ADPm_ATPm = (*ADPm) / ATPm;
	VANT = ( VmDT_75 - VmDT_20 * ATPi_ADP * ADPm_ATPm * exp(-F_over_RT * (*Dpsi)) ) /
		( ( 1.0 + tenDiv9 * ATPi_ADP * exp( -hm_F_over_RT * (*Dpsi) ) ) * (1.0 + 18 * ADPm_ATPm ) );
	Vhleak = gh * DmuH;
}
//getFNa will get the time rate of change of mNa, hNa, and jNa; fortran;
void Model::calculateFNa() {
	//Parameters:
	//Verbose check done
	double MAlpha;
	double MBeta;
	double inv_MBeta_MAlpha;
	double HAlpha;
	double HBeta;
	double JAlpha;
	double JBeta;
	double tmNa;

	if ( (*V) == -47.13)
		MAlpha = 3.2;
	else
		MAlpha = 0.32 * ( (*V) + 47.13 ) / ( 1.0 - exp(-0.1 * ( (*V) + 47.13) ) );

	MBeta = 0.08*exp (  -(*V) * inv_11);
	inv_MBeta_MAlpha = 1.0 / ( MBeta + MAlpha);
	if ( inv_MBeta_MAlpha < 0.03) {
		tmNa = MAlpha * inv_MBeta_MAlpha;
	} else {
		tmNa = (*mNa);
	}

	if ( (*V) < -40 ) {
	    HAlpha = hNa_HAlpha_C1*exp( inv_6p8 * (*V) );
	    HBeta =  3.56 * exp(0.079 * (*V) ) + 310000.0 * exp( 0.35 * (*V) );
	    JAlpha = (-127140.0 * exp(0.2444 * (*V)) - 3.474E-5 * exp ( - 0.04391 * (*V) ) ) *
			     ( (*V) + 37.78 ) / (1.0 + exp( 0.311 * ( (*V) + 79.23) ));
	    JBeta = 0.1212 * exp( -0.01052 * (*V) ) / ( 1.0 + exp( -0.1378 * ( (*V) + 40.14 ) ) );
	} else {
		HAlpha = 0.0;
	    JAlpha = 0.0;
	    HBeta = 1.0 / (  0.13 + hNa_HBeta_C1 * exp( (*V) * inv_11p1));
	    JBeta = 0.3*exp( - 2.535E-7 * (*V)) / (1.0 + exp( -0.1 * (*V) - 3.2) );
	}
	(*dmNa) = MAlpha * ( 1.0 - (tmNa) ) - MBeta * (tmNa);
	(*dhNa) = HAlpha * (1.0 - (*hNa)) - HBeta * (*hNa);	
	(*djNa) = JAlpha * (1.0 - (*jNa)) - JBeta * (*jNa);	
}
//xKs
void Model::calculateFxKs() {
	(*dxKs) = 7.19E-5 * V_30 / (1.0 - exp( -0.148 * V_30)) * (1.0-(*xKs)) - 
		   1.31E-4 * V_30 /  (exp(0.0687 * V_30) - 1.0) * (*xKs);
}
//computes current values for Jup, Jrel, Jtr, and Jxfer; After ATPmod;O1_RyR;ADP;inv_ATPi
void Model::calculateInCalFlux(){
	//Paramerters: Kfb, Krb, Nfb,Nrb, KSR, vmaxf, vmaxr, KmATP_SR, Ki_SR, Ki_prime_SR, v1, tauxfer, tautr
	double fb = pow( (*Cai) * inv_Kfb, Nfb);
	double rb = pow( (*CaNSR) * inv_Krb, Nrb);
	Jup   = KSR * ( vmaxf * fb - vmaxr * rb) / 
		( (1.0+fb+rb)*( inv_ATPi*( KmATP_SR + ADP * KmATP_SR_Ki_SR) + (1.0 + ADP * inv_Ki_prime_SR) ) );
	Jrel  = v1 * ( O1_RyR + (*O2_RyR) ) * ( (*CaJSR) - (*CaSS) );
	Jtr   = ( (*CaNSR) - (*CaJSR) ) * inv_tautr;
	Jxfer = ( (*CaSS)  - (*Cai) )   * inv_tauxfer;
}
//Force Model; FN_Ca;Fnorm; Force
void Model::calculateForce() {
	//Parameters: alpha_SL, zeta
	P1_N1_P2_P3 = (*P1) + (*N1) + (*P2) + (*P3);
	FN_Ca = alpha_SL_fnormmax2 * P1_N1_P2_P3;
}
//Force Model; Tropomyosin_XB added to project (JT);LTRPNCa;kTrop_np;dP0;dP1;dP2;dP3;dN0;dN1;Fortran;
void Model::calculateF_trpmyo() {
	//Parameters: kTrop_pn; LTRPNtot; Ktrop_half; Ntrop; f_01; g_01_mod; f_12; g_12_mod; g_23_mod; f_23
	double kTrop_np = kTrop_pn * pow( (*LTRPNCa) * inv_LTRPNtot_Ktrop_half, Ntrop);
	(*dP0) = kTrop_pn_f_01 * (*P0) + kTrop_np * (*N0) + g_01_mod * (*P1);
	(*dP1) = kTrop_pn_f_12_g_01_mod * (*P1) + kTrop_np * (*N1) + f_01 * (*P0) + g_12_mod * (*P2);
	(*dP2) = f_23_g_12_mod * (*P2) + f_12 * (*P1) + g_23_mod * (*P3);
	(*dP3) = -g_23_mod * (*P3) + f_23 * (*P2);
	(*dN1) = kTrop_pn * (*P1) - (kTrop_np + g_01_off_mod) * (*N1);
	(*dN0) = -(*dP0)-(*dP1)-(*dP2)-(*dP3)-(*dN1);							//added by JT
}
// TROPONIN FRACTION DERIVATIVES; Fortran; FN_Ca; kltrpn_plus_Cai
// These methods get the time rate of change for LTRPNCa and HTRPNCa
//Force Model;
void Model::calculateFLTRPNCa() {
	//Parameters:kltrpn_plus;LTRPNtot;kltrpn_minus
	(*dLTRPNCa) = kltrpn_plus * (*Cai) * (LTRPNtot - (*LTRPNCa))- ( kltrpn_minus * (*LTRPNCa) ) * ( 1.0 - twoThirds * FN_Ca );
}
//dHTRPNCa; kltrpn_plus_Cai
void Model::calculateFHTRPNCa() {
	//Parameters:kltrpn_plus;HTRPNtot;kltrpn_minus
	(*dHTRPNCa) = khtrpn_plus * (*Cai) * (HTRPNtot - (*HTRPNCa)) - khtrpn_minus * (*HTRPNCa);
}
//computes current values for beta_SS, beta_JSR, and beta_i
void Model::calculateJtrpn() {
	//Parameters: CMDNtot; KmCMDN; CSQNtot;KmCSQN
	Jtrpn = (*dLTRPNCa) + (*dHTRPNCa);

	beta_SS  = 1.0 / ( 1.0 + CMDNtot_KmCMDN / ( ((*CaSS)  + KmCMDN)*((*CaSS)  + KmCMDN) ) );
	beta_JSR = 1.0 / ( 1.0 + CSQNtot_KmCSQN / ( ((*CaJSR) + KmCSQN)*((*CaJSR) + KmCSQN) ) );
	beta_i   = 1.0 / ( 1.0 + CMDNtot_KmCMDN / ( ((*Cai)   + KmCMDN)*((*Cai)   + KmCMDN) ) );
}
//Mitochondia;Vuni; F_over_RT; FRT2_Dpsi; Treat Cai as parameter in Mitochondia only model
void Model::calculateVuni() {
	//Parameters: Vmuni; ktrans; L; kact; na
	double Cai_ktrans_plus1 = 1.0 + (*Cai) * inv_ktrans;
	double Cai_ktrans_plus1_p3 = Cai_ktrans_plus1 * Cai_ktrans_plus1 * Cai_ktrans_plus1;
	Vuni =  Vmuni_ktrans * (*Cai) * FRT2_Dpsi * Cai_ktrans_plus1_p3  / 
		( ( Cai_ktrans_plus1_p3 * Cai_ktrans_plus1 + L / pow( 1.0 + (*Cai) * inv_kact, na) ) * ( 1.0 - exp(-FRT2_Dpsi) ) );
}
///Mitochondia; F_over_RT; FRT2_Dpsi; 1/Nai; 1/Cam; 1/Cai ; Cai Nai as params in Mito only
void Model::calculateVnaCa() {
	//Parameters: b; VmNC; Kna; n; Knca;
	VnaCa = VmNC * exp( b_05 * FRT2_Dpsi ) * (*Cam) / ( (*Cai) * pow( ( 1.0 + Kna / (*Nai) ) , n) * ( 1.0 + Knca / (*Cam) ) );
}
//Intercellular Concentration calculations
//These methods get the time rate of change for Nai, Ki, and Cai; Post-Mitochondrial mod
//INa;INab;INaCa;INaK;InsNa;VnaCa | Acap_Vmyo_F
void Model::calculateFNai() {
	//Parameters: InsNa; Vmyo; Acap
	(*dNai) = - ( INa + INab + InsNa + 3.0 * ( INaCa + INaK ) ) * Acap_Vmyo_F - VnaCa * 0.615;
	//(*dNai) = 0.0;
}
//Ki; Acap_Vmyo_F | INaK, Istim, ICaK, IKp, IK1, IKs
void Model::calculateFKi() {
	//Parameters: InsK
	//temporarily change in lines to calculate elasticities!
	(*dKi) = -(InsK + IKs + IK1 + IKp + ICaK + Istim - 2.0 * INaK) * Acap_Vmyo_F; 
	//(*dKi) = 0.0;
}
//Cai; Acap_Vmyo_F; ICab; INaCa; IpCa; beta_i; Jxfer; Jup; Jtrpn
void Model::calculateFCai() {
	//temporarily change in line 959 to calculate elasticities!
	(*dCai) = beta_i * ( Jxfer - Jup - Jtrpn - 0.25 * Acap_Vmyo_F * (ICab - 2.0 * INaCa + IpCa) + ( VnaCa - Vuni ) * 0.615);
	//(*dCai) = 0.0;
}
//These methods get the time rate of change for CaSS, CaJSR, and CaNSR
//CaSS; ICa, beta_SS
void Model::calculateFCaSS() {
	//Parameters: VSS
	(*dCaSS) = beta_SS * ( Jrel * VJSR_VSS - Jxfer * Vmyo_VSS - ICa * Acap_VSS_F);
	//(*dCaSS) = 0.0;
}
//CaJSR;beta_JSR; Jtr; Jrel
void Model::calculateFCaJSR() {
	(*dCaJSR) = beta_JSR * ( Jtr - Jrel );
	//(*dCaJSR) = 0.0;
}
//CaNSR; Jup; Jtr
void Model::calculateFCaNSR() {
	//parameter: Vmyo,| VNSR; VJSR
	(*dCaNSR) = Jup * Vmyo_VNSR - Jtr * VJSR_VNSR;
	//(*dCaNSR) = 0.0;
}
//gets the time rate of change for V; INa; INaCa; INab
void Model::calculateFV() {	
	//Parameters: C_m
	if ( clampVoltage || usingAlgebraicMembranePotential ) {	
		(*dV) = 0;
	} else {
		(*dV)= -inv_C_m *( INa + ICa + ICaK + IKs + IK1 + IKp + INaCa + INaK + InsCa + IpCa + ICab + INab + Istim );
	}
}
//getFRyR gets the time rate of change for C1_RyR,O2_RyR, C2_RyR.
void Model::calculateFRyR() {
	//Parameters: kaplus; kbplus; kcplus; kaminus; kbminus; kcminus; ncoop; mcoop
	(*dC1_RyR) = -kaplus * pow( (*CaSS), ncoop) * (*C1_RyR) + kaminus * O1_RyR;
	(*dO2_RyR) =  kbplus * pow( (*CaSS), mcoop) * O1_RyR    - kbminus * (*O2_RyR);
	(*dC2_RyR) =  kcplus *                        O1_RyR    - kcminus * (*C2_RyR);
}
//L-Type Ca++ current; fortran; Some Variables could be optimized( Gamma|Omega)
//getFCaL gets the time rate of change for C0,C1,C2,C3,C4,Open,CCa0,CCa1,CCa2,CCa3, and CCa4
//Long and barely worked on. Check this out in detail when their's some more time.
void Model::calculateFCaL() {
	//Parameters: omega;gprime;fprime Other*?*, gL, fL, aL
	//LCa: L-Type Calcium Channel state transition rates
	double gamma;
	double C0_to_C1;
	double C1_to_C2; 
	double C2_to_C3; 
	double C3_to_C4;
	double C1_to_C0;
	double C2_to_C1; 
	double C3_to_C2;
	double C4_to_C3;
	double CCa0_to_CCa1; 
	double CCa1_to_CCa2;
	double CCa2_to_CCa3;
	double CCa3_to_CCa4;
	double CCa1_to_CCa0;
	double CCa2_to_CCa1;
	double CCa3_to_CCa2;
	double CCa4_to_CCa3;
	double C0_to_CCa0;
	double C1_to_CCa1;
	double C2_to_CCa2;
	double C3_to_CCa3;
	double C4_to_CCa4;
	double CCa0_to_C0; 
	double CCa1_to_C1;
	double CCa2_to_C2;
	double CCa3_to_C3;
	double CCa4_to_C4;

	//A whole bunch of Cx_to_Cx Variables
	//On April 04 2006 I modified the expression for the gamma =0.5625 instead of 0.140625
	// and changed the expression of the alpha and beta substituting the "12" by a "2" as in the canine code. S.C.
	double alpha = 0.4 *  exp( ( (*V) + 2 ) * 0.1);
	double beta  = 0.05 * exp(( (*V) + 2 ) * neg_inv_13); 

	double alpha_prime = aL * alpha;
	double beta_prime  = beta * inv_bL;		// change bL to 2.0 in parameters

	C0_to_C1 = 4.0 * alpha;
	C1_to_C2 = 3.0 * alpha;
	C2_to_C3 = 2.0 * alpha;
	C3_to_C4 = alpha;

	CCa0_to_CCa1 = 4.0 * alpha_prime;
	CCa1_to_CCa2 = 3.0 * alpha_prime;
	CCa2_to_CCa3 = 2.0 * alpha_prime;
	CCa3_to_CCa4 = alpha_prime;

	C1_to_C0 = beta;
	C2_to_C1 = 2.0 * beta;
	C3_to_C2 = 3.0 * beta;
	C4_to_C3 = 4.0 * beta;

	CCa1_to_CCa0 = beta_prime;
	CCa2_to_CCa1 = 2.0 * beta_prime;
	CCa3_to_CCa2 = 3.0 * beta_prime;
	CCa4_to_CCa3 = 4.0 * beta_prime;
	
	gamma = 0.1875 * (*CaSS); // = 0.140625

	C0_to_CCa0 = gamma;		
	C1_to_CCa1 = aL * C0_to_CCa0;	// = gamma*aL
	C2_to_CCa2 = aL * C1_to_CCa1;	// = gamma*aL^2
	C3_to_CCa3 = aL * C2_to_CCa2;	// = gamma*aL^3
	C4_to_CCa4 = aL * C3_to_CCa3;	// = gamma*aL^4
		
	CCa0_to_C0 = omega;		// = omega
	CCa1_to_C1 = CCa0_to_C0 * inv_bL;	// = omega/bL
	CCa2_to_C2 = CCa1_to_C1 * inv_bL;	// = omega/bL^2
	CCa3_to_C3 = CCa2_to_C2 * inv_bL;	// = omega/bL^3
	CCa4_to_C4 = CCa3_to_C3 * inv_bL;	// = omega/bL^4

	(*dC0) = C1_to_C0 * (*C1) + CCa0_to_C0 * (*CCa0) - 
			( C0_to_C1 + C0_to_CCa0 ) * (*C0); 
	(*dC1) = C0_to_C1 * (*C0) + C2_to_C1 * (*C2) + CCa1_to_C1 * (*CCa1) - 
			( C1_to_C0 + C1_to_C2 + C1_to_CCa1 ) * (*C1);
	(*dC2) = C1_to_C2 * (*C1) + C3_to_C2 * (*C3) + CCa2_to_C2 * (*CCa2) - 
			( C2_to_C1 + C2_to_C3 + C2_to_CCa2) * (*C2);
	(*dC3) = C2_to_C3 * (*C2) + C4_to_C3 * (*C4) + CCa3_to_C3 * (*CCa3) - 
			( C3_to_C2 + C3_to_C4 + C3_to_CCa3) * (*C3);
	(*dC4) = C3_to_C4 * (*C3) + gL * (*Open) + CCa4_to_C4 * (*CCa4) - 
			( C4_to_C3 + fL + C4_to_CCa4 ) * (*C4);
	(*dOpen) =  fL * (*C4) - gL * (*Open);
	(*dCCa0) = CCa1_to_CCa0 * (*CCa1) + C0_to_CCa0 * (*C0) - 
			( CCa0_to_CCa1 + CCa0_to_C0 ) * (*CCa0);
	(*dCCa1) = CCa0_to_CCa1 * (*CCa0) + CCa2_to_CCa1 * (*CCa2) + C1_to_CCa1 * (*C1) - 
			( CCa1_to_CCa0 + CCa1_to_CCa2 + CCa1_to_C1) * (*CCa1);
	(*dCCa2) = CCa1_to_CCa2 * (*CCa1) + CCa3_to_CCa2 * (*CCa3) + C2_to_CCa2 * (*C2) - 
			( CCa2_to_CCa1 + CCa2_to_CCa3 + CCa2_to_C2) * (*CCa2);
	(*dCCa3) = CCa2_to_CCa3 * (*CCa2) + CCa4_to_CCa3 * (*CCa4) + C3_to_CCa3 * (*C3) - 
			( CCa3_to_CCa2 + CCa3_to_CCa4 + CCa3_to_C3 ) * (*CCa3);
	(*dCCa4) = CCa3_to_CCa4 * (*CCa3) + gprime * (*OCa) + C4_to_CCa4 * (*C4) - 
			( CCa4_to_CCa3 + fprime + CCa4_to_C4 ) * (*CCa4);
}
//getFyCa gets the time rate of change for yCa
void Model::calculateFyCa() {
	(*dyCa) = ( 1.0 / ( 1.0 + exp( ( (*V) + 55.0 ) * inv_7p5 ) ) + 0.5 / ( 1.0 + exp( ( 21.0 - (*V) ) * inv_6 ) ) - (*yCa) ) / 
		      ( 20.0 + 600.0 / ( 1.0 + exp( ( (*V) + 30.0 ) * inv_9p5 ) ) );
}
//OCa|VANT, V_AM, Jup, INaK, IpCa
void Model::calculateFOCa(){
	//Parameters: fprime, gprime
	(*dOCa) = fprime * (*CCa4) - gprime * (*OCa);
} 
//C_K Mod, Comes before calculateFATPi. Combination of getFCrPi(),getCK()
void Model::calculateCK(){
	//Parameters: kt_2, kf_2, kf_3, keq, CRT_cyto, CRT_cyto, CRT_mito, VATPase_cyto
	//Precomputes: inv_keq = 1.0 / keq;
	double Vt_CRP2  = kt_2 * ( (*CrPi_mito) - (*CrPi_cyto) );
	double VCK_cyto = kf_2 * ( ( CRT_cyto - (*CrPi_cyto) ) * (*ATPi_cyto) -  (*CrPi_cyto) * ( 8.0 - (*ATPi_cyto) ) * inv_keq);
	double VCK_mito = kf_3 * ( ( CRT_mito - (*CrPi_mito) ) * (*ATPi) - (*CrPi_mito) * ADP * inv_keq);

	(*dCrPi_mito) = VCK_mito - Vt_CRP2;
	(*dCrPi_cyto) = Vt_CRP2 + VCK_cyto;
	(*dATPi) = 0.615 * VANT - V_AM - 0.5*Jup - ( 6.371e-5 * ( INaK + IpCa ) ) - VCK_mito; // defined normally in calculateFATPi
	(*dATPi_cyto) = - VCK_cyto - VATPase_cyto;
}
//ATP
void Model::calculateFATPi() {
	(*dATPi) = 0.615 * VANT - V_AM - 0.5* Jup - ( 6.371e-5 * ( INaK + IpCa ) );
}
//All the mitochondria stuff
//Missing ASP (due to change in VAAT calc), and some other slight modifications
void Model::calculateF_mitene() {
	//Parameters:Cmito, fm, b, kcnsASP
	(*dCam)  = fm * (Vuni - VnaCa);
	(*dADPm) = VANT - VATPase - VSL;
	(*dDpsi) = -(-VHNe - VHFe + Vhu + VANT + Vhleak + two_b * VnaCa + 2.0 * Vuni ) * inv_Cmito;
	(*dNADH) = -VNO + VIDH + VKGDH + VMDH;
	(*dISOC) = VACO - VIDH;
	(*dAKG)  = VIDH + VAAT - VKGDH;
	(*dSCoA) = VKGDH - VSL;
	(*dSucc) = VSL - VSDH;
	(*dFUM)  = VSDH - VFH;
	(*dMAL)  = VFH - VMDH;
	(*dOaa)  = VMDH - VCS - VAAT;
	(*dASP)  = VAAT - kcnsASP * (*ASP);
}
int Model::getIndexOffset() {
	switch(modelType) {
		case Coupled:
			if(usingAlgebraicMembranePotential)
				return index_Alg_Offset;
			else
				return 0;
			break;
		case Mitochondria:
			return  index_Mito_Offset;
			break;
		case Force:
			return index_Force_Offset;
			break;
	}
	cout << "Something wrong with Model::getIndexOffset()." << endl;
	exit(-1);
}
//Model Access functions, Skeletons put in but will have to be changed to handle multiple models
//The errorweight matrix
double* Model::getErrorWeights() {
	return errweight + getIndexOffset();
}
//Return the s array
double* Model::getStates() {
	return s + getIndexOffset();
}
//Return the dS array
double* Model::getStatesDerivatives() {
	return dS + getIndexOffset();
}
//Returns the number of state variable in this particular problem mode
int Model::getProblemSize() {
	switch(modelType) {
		case Coupled:
			if(usingAlgebraicMembranePotential)
				if (usingCK)
					return PROBLEM_SIZE_GPC_CK_ALG;
				else
					return PROBLEM_SIZE_GPC_ALG;
			else
				if (usingCK)
					return PROBLEM_SIZE_GPC_CK;
				else
					return PROBLEM_SIZE_GPC;
			break;
		case Mitochondria:
			if(usingASP)
				return PROBLEM_SIZE_MITO_ASP;
			else
				return PROBLEM_SIZE_MITO;
			break;
		case Force:
			return PROBLEM_SIZE_FORCE;
			break;
	}
	cout << "Something wrong with Model::getProblemSize()." << endl;
	exit(-1);
	return -1;
}
void Model::F(double start_time, double *ydot, double *y) {
	//Setup pointer references and shortcuts (Coupled, Mitochondria, Force), usingAlgebraicMembranePotential, usingASP
	//And run the correct F function.
	switch(modelType) {
		case Coupled:
			if(usingAlgebraicMembranePotential) {
				linkStatesReferences_Alg(ydot, y);
				linkStatesReferences_CK_Only(ydot, y, index_Alg_Offset);
			} else {
				linkStatesReferences_GPC(ydot, y);
				linkStatesReferences_CK_Only(ydot, y, 0);
			}
			F_GPC(start_time);
			break;
		case Mitochondria:
			if(usingASP)
				linkStatesReferences_Mito_ASP(ydot, y);
			else
				linkStatesReferences_Mito(ydot, y);
			F_Mitochondria();
			break;
		case Force:
			linkStatesReferences_Force(ydot, y);
			F_Force();
			break;
	}

}
//Setup params for the IF mode experiment
void Model::setupIFmode(int num_run){
	//Vars stop_time
	//Parameters: ESI_increment, PESI, refrac_buffer
	if(stimulasMode == IF) {
		time_on_Is1 = num_run * ESI_increment + refrac_buffer;
		time_off_Is1 = time_on_Is1 + pulse_duration;
		time_on_Is2 = time_off_Is1 + PESI;
		time_off_Is2 = time_on_Is2 + pulse_duration;
		stop_time = time_off_Is2 + 800.0;
	}
}
//This set of functions assign shortcut pointers to all data;
void Model::linkStatesReferences_GPC( double *ydot, double *y ) {
	V = &y[index_V];	
	mNa = &y[index_mNa];
	hNa = &y[index_hNa];
	jNa = &y[index_jNa];
	xKs = &y[index_xKs];
	Nai = &y[index_Nai];
	Ki = &y[index_Ki];
	Cai = &y[index_Cai];
	CaNSR = &y[index_CaNSR];
	CaSS = &y[index_CaSS];
	CaJSR = &y[index_CaJSR];
	C1_RyR = &y[index_C1_RyR];
	O2_RyR = &y[index_O2_RyR];
	C2_RyR = &y[index_C2_RyR];
	C0 = &y[index_C0];
	C1 = &y[index_C1];
	C2 = &y[index_C2];
	C3 = &y[index_C3];
	C4 = &y[index_C4];
	Open = &y[index_Open];
	CCa0 = &y[index_CCa0];
	CCa1 = &y[index_CCa1];
	CCa2 = &y[index_CCa2];
	CCa3 = &y[index_CCa3];
	CCa4 = &y[index_CCa4];
	OCa = &y[index_OCa];
	yCa = &y[index_yCa];
	LTRPNCa = &y[index_LTRPNCa];
	HTRPNCa = &y[index_HTRPNCa];
	N0 = &y[index_N0];
	N1 = &y[index_N1];
	P0 = &y[index_P0];
	P1 = &y[index_P1];
	P2 = &y[index_P2];
	P3 = &y[index_P3];
	ATPi = &y[index_ATPi];
	Cam = &y[index_Cam];
	ADPm = &y[index_ADPm];
	Dpsi = &y[index_Dpsi];
	NADH = &y[index_NADH];
	ISOC = &y[index_ISOC];
	AKG = &y[index_AKG];
	SCoA = &y[index_SCoA];
	Succ = &y[index_Succ];
	FUM = &y[index_FUM];
	MAL = &y[index_MAL];
	Oaa = &y[index_Oaa];

	//Derivative refences
	dV = &ydot[index_V];	
	dmNa = &ydot[index_mNa];
	dhNa = &ydot[index_hNa];
	djNa = &ydot[index_jNa];
	dxKs = &ydot[index_xKs];
	dNai = &ydot[index_Nai];
	dKi = &ydot[index_Ki];
	dCai = &ydot[index_Cai];
	dCaNSR = &ydot[index_CaNSR];
	dCaSS = &ydot[index_CaSS];
	dCaJSR = &ydot[index_CaJSR];
	dC1_RyR = &ydot[index_C1_RyR];
	dO2_RyR = &ydot[index_O2_RyR];
	dC2_RyR = &ydot[index_C2_RyR];
	dC0 = &ydot[index_C0];
	dC1 = &ydot[index_C1];
	dC2 = &ydot[index_C2];
	dC3 = &ydot[index_C3];
	dC4 = &ydot[index_C4];
	dOpen = &ydot[index_Open];
	dCCa0 = &ydot[index_CCa0];
	dCCa1 = &ydot[index_CCa1];
	dCCa2 = &ydot[index_CCa2];
	dCCa3 = &ydot[index_CCa3];
	dCCa4 = &ydot[index_CCa4];
	dOCa = &ydot[index_OCa];
	dyCa = &ydot[index_yCa];
	dLTRPNCa = &ydot[index_LTRPNCa];
	dHTRPNCa = &ydot[index_HTRPNCa];
	dN0 = &ydot[index_N0];
	dN1 = &ydot[index_N1];
	dP0 = &ydot[index_P0];
	dP1 = &ydot[index_P1];
	dP2 = &ydot[index_P2];
	dP3 = &ydot[index_P3];
	dATPi = &ydot[index_ATPi];
	dCam = &ydot[index_Cam];
	dADPm = &ydot[index_ADPm];
	dDpsi = &ydot[index_Dpsi];
	dNADH = &ydot[index_NADH];
	dISOC = &ydot[index_ISOC];
	dAKG = &ydot[index_AKG];
	dSCoA = &ydot[index_SCoA];
	dSucc = &ydot[index_Succ];
	dFUM = &ydot[index_FUM];
	dMAL = &ydot[index_MAL];
	dOaa = &ydot[index_Oaa];
}
void Model::linkStatesReferences_CK_Only( double *ydot, double *y , int offset) {
	ATPi_cyto = &y[index_ATPi_cyto - offset ];
	CrPi_cyto = &y[index_CrPi_cyto - offset ];
	CrPi_mito = &y[index_CrPi_mito - offset ];
		
	dATPi_cyto = &ydot[index_ATPi_cyto - offset ];
	dCrPi_cyto = &ydot[index_CrPi_cyto - offset ];
	dCrPi_mito = &ydot[index_CrPi_mito - offset ];
}
void Model::linkStatesReferences_Full( double *ydot, double *y ) {
	linkStatesReferences_GPC( ydot, y );
	linkStatesReferences_CK_Only( ydot, y , 0);
	ASP = &y[index_ASP];
	dASP = &ydot[index_ASP];
}
void Model::linkStatesReferences_Alg( double *ydot, double *y ) {
	mNa = &y[index_mNa - index_Alg_Offset];
	hNa = &y[index_hNa - index_Alg_Offset];
	jNa = &y[index_jNa - index_Alg_Offset];
	xKs = &y[index_xKs - index_Alg_Offset];
	Nai = &y[index_Nai - index_Alg_Offset];
	Ki = &y[index_Ki - index_Alg_Offset];
	Cai = &y[index_Cai - index_Alg_Offset];
	CaNSR = &y[index_CaNSR - index_Alg_Offset];
	CaSS = &y[index_CaSS - index_Alg_Offset];
	CaJSR = &y[index_CaJSR - index_Alg_Offset];
	C1_RyR = &y[index_C1_RyR - index_Alg_Offset];
	O2_RyR = &y[index_O2_RyR - index_Alg_Offset];
	C2_RyR = &y[index_C2_RyR - index_Alg_Offset];
	C0 = &y[index_C0 - index_Alg_Offset];
	C1 = &y[index_C1 - index_Alg_Offset];
	C2 = &y[index_C2 - index_Alg_Offset];
	C3 = &y[index_C3 - index_Alg_Offset];
	C4 = &y[index_C4 - index_Alg_Offset];
	Open = &y[index_Open - index_Alg_Offset];
	CCa0 = &y[index_CCa0 - index_Alg_Offset];
	CCa1 = &y[index_CCa1 - index_Alg_Offset];
	CCa2 = &y[index_CCa2 - index_Alg_Offset];
	CCa3 = &y[index_CCa3 - index_Alg_Offset];
	CCa4 = &y[index_CCa4 - index_Alg_Offset];
	OCa = &y[index_OCa - index_Alg_Offset];
	yCa = &y[index_yCa - index_Alg_Offset];
	LTRPNCa = &y[index_LTRPNCa - index_Alg_Offset];
	HTRPNCa = &y[index_HTRPNCa - index_Alg_Offset];
	N0 = &y[index_N0 - index_Alg_Offset];
	N1 = &y[index_N1 - index_Alg_Offset];
	P0 = &y[index_P0 - index_Alg_Offset];
	P1 = &y[index_P1 - index_Alg_Offset];
	P2 = &y[index_P2 - index_Alg_Offset];
	P3 = &y[index_P3 - index_Alg_Offset];
	ATPi = &y[index_ATPi - index_Alg_Offset];
	Cam = &y[index_Cam - index_Alg_Offset];
	ADPm = &y[index_ADPm - index_Alg_Offset];
	Dpsi = &y[index_Dpsi - index_Alg_Offset];
	NADH = &y[index_NADH - index_Alg_Offset];
	ISOC = &y[index_ISOC - index_Alg_Offset];
	AKG = &y[index_AKG - index_Alg_Offset];
	SCoA = &y[index_SCoA - index_Alg_Offset];
	Succ = &y[index_Succ - index_Alg_Offset];
	FUM = &y[index_FUM - index_Alg_Offset];
	MAL = &y[index_MAL - index_Alg_Offset];
	Oaa = &y[index_Oaa - index_Alg_Offset];
	ASP = &y[index_ASP - index_Alg_Offset];

	//Derivative refences
	dV = &ydot[index_V - index_Alg_Offset];	
	dmNa = &ydot[index_mNa - index_Alg_Offset];
	dhNa = &ydot[index_hNa - index_Alg_Offset];
	djNa = &ydot[index_jNa - index_Alg_Offset];
	dxKs = &ydot[index_xKs - index_Alg_Offset];
	dNai = &ydot[index_Nai - index_Alg_Offset];
	dKi = &ydot[index_Ki - index_Alg_Offset];
	dCai = &ydot[index_Cai - index_Alg_Offset];
	dCaNSR = &ydot[index_CaNSR - index_Alg_Offset];
	dCaSS = &ydot[index_CaSS - index_Alg_Offset];
	dCaJSR = &ydot[index_CaJSR - index_Alg_Offset];
	dC1_RyR = &ydot[index_C1_RyR - index_Alg_Offset];
	dO2_RyR = &ydot[index_O2_RyR - index_Alg_Offset];
	dC2_RyR = &ydot[index_C2_RyR - index_Alg_Offset];
	dC0 = &ydot[index_C0 - index_Alg_Offset];
	dC1 = &ydot[index_C1 - index_Alg_Offset];
	dC2 = &ydot[index_C2 - index_Alg_Offset];
	dC3 = &ydot[index_C3 - index_Alg_Offset];
	dC4 = &ydot[index_C4 - index_Alg_Offset];
	dOpen = &ydot[index_Open - index_Alg_Offset];
	dCCa0 = &ydot[index_CCa0 - index_Alg_Offset];
	dCCa1 = &ydot[index_CCa1 - index_Alg_Offset];
	dCCa2 = &ydot[index_CCa2 - index_Alg_Offset];
	dCCa3 = &ydot[index_CCa3 - index_Alg_Offset];
	dCCa4 = &ydot[index_CCa4 - index_Alg_Offset];
	dOCa = &ydot[index_OCa - index_Alg_Offset];
	dyCa = &ydot[index_yCa - index_Alg_Offset];
	dLTRPNCa = &ydot[index_LTRPNCa - index_Alg_Offset];
	dHTRPNCa = &ydot[index_HTRPNCa - index_Alg_Offset];
	dN0 = &ydot[index_N0 - index_Alg_Offset];
	dN1 = &ydot[index_N1 - index_Alg_Offset];
	dP0 = &ydot[index_P0 - index_Alg_Offset];
	dP1 = &ydot[index_P1 - index_Alg_Offset];
	dP2 = &ydot[index_P2 - index_Alg_Offset];
	dP3 = &ydot[index_P3 - index_Alg_Offset];
	dATPi = &ydot[index_ATPi - index_Alg_Offset];
	dCam = &ydot[index_Cam - index_Alg_Offset];
	dADPm = &ydot[index_ADPm - index_Alg_Offset];
	dDpsi = &ydot[index_Dpsi - index_Alg_Offset];
	dNADH = &ydot[index_NADH - index_Alg_Offset];
	dISOC = &ydot[index_ISOC - index_Alg_Offset];
	dAKG = &ydot[index_AKG - index_Alg_Offset];
	dSCoA = &ydot[index_SCoA - index_Alg_Offset];
	dSucc = &ydot[index_Succ - index_Alg_Offset];
	dFUM = &ydot[index_FUM - index_Alg_Offset];
	dMAL = &ydot[index_MAL - index_Alg_Offset];
	dOaa = &ydot[index_Oaa - index_Alg_Offset];
}
void Model::linkStatesReferences_Force( double *ydot, double *y ) {
	LTRPNCa = &y[index_LTRPNCa - index_Force_Offset];
	N0 = &y[index_N0 - index_Force_Offset];
	N1 = &y[index_N1 - index_Force_Offset];
	P0 = &y[index_P0 - index_Force_Offset];
	P1 = &y[index_P1 - index_Force_Offset];
	P2 = &y[index_P2 - index_Force_Offset];
	P3 = &y[index_P3 - index_Force_Offset];

	//Derivative refences
	dLTRPNCa = &ydot[index_LTRPNCa - index_Force_Offset];
	dN0 = &ydot[index_N0 - index_Force_Offset];
	dN1 = &ydot[index_N1 - index_Force_Offset];
	dP0 = &ydot[index_P0 - index_Force_Offset];
	dP1 = &ydot[index_P1 - index_Force_Offset];
	dP2 = &ydot[index_P2 - index_Force_Offset];
	dP3 = &ydot[index_P3 - index_Force_Offset];
}
void Model::linkStatesReferences_Mito( double *ydot, double *y ) {
	Cam = &y[index_Cam - index_Mito_Offset];
	ADPm = &y[index_ADPm - index_Mito_Offset];
	Dpsi = &y[index_Dpsi - index_Mito_Offset];
	NADH = &y[index_NADH - index_Mito_Offset];
	ISOC = &y[index_ISOC - index_Mito_Offset];
	AKG = &y[index_AKG - index_Mito_Offset];
	SCoA = &y[index_SCoA - index_Mito_Offset];
	Succ = &y[index_Succ - index_Mito_Offset];
	FUM = &y[index_FUM - index_Mito_Offset];
	MAL = &y[index_MAL - index_Mito_Offset];
	Oaa = &y[index_Oaa - index_Mito_Offset];

	//Derivative refences
	dCam = &ydot[index_Cam - index_Mito_Offset];
	dADPm = &ydot[index_ADPm - index_Mito_Offset];
	dDpsi = &ydot[index_Dpsi - index_Mito_Offset];
	dNADH = &ydot[index_NADH - index_Mito_Offset];
	dISOC = &ydot[index_ISOC - index_Mito_Offset];
	dAKG = &ydot[index_AKG - index_Mito_Offset];
	dSCoA = &ydot[index_SCoA - index_Mito_Offset];
	dSucc = &ydot[index_Succ - index_Mito_Offset];
	dFUM = &ydot[index_FUM - index_Mito_Offset];
	dMAL = &ydot[index_MAL - index_Mito_Offset];
	dOaa = &ydot[index_Oaa - index_Mito_Offset];
}
void Model::linkStatesReferences_Mito_ASP( double *ydot, double *y ) {
	linkStatesReferences_Mito( ydot, y );
	ASP = &y[index_ASP - index_Mito_Offset];
	dASP = &ydot[index_ASP - index_Mito_Offset];
}
//Set all the model parameters
void Model::setParameters(fileIO& data) {
	setParamaterWithLink(data,"Acap",Acap);
	setParamaterWithLink(data,"AcCoA",AcCoA);
	setParamaterWithLink(data,"aL",aL);
	setParamaterWithLink(data,"b",b);
	setParamaterWithLink(data,"bL",bL);
	setParamaterWithLink(data,"C_m",C_m);
	setParamaterWithLink(data,"Cao",Cao);
	setParamaterWithLink(data,"chfsc_IK1",chfsc_IK1);
	setParamaterWithLink(data,"chfsc_INaCa",chfsc_INaCa);
	setParamaterWithLink(data,"chfsc_Jup",chfsc_Jup);
	setParamaterWithLink(data,"CIK",CIK);
	setParamaterWithLink(data,"Cm",Cm);
	setParamaterWithLink(data,"CMDNtot",CMDNtot);
	setParamaterWithLink(data,"Cmito",Cmito);
	setParamaterWithLink(data,"CPN",CPN);
	setParamaterWithLink(data,"CoA",CoA);
	setParamaterWithLink(data,"CSQNtot",CSQNtot);
	setParamaterWithLink(data,"DpH",DpH);
	setParamaterWithLink(data,"Dpsio",Dpsio);
	setParamaterWithLink(data,"ESI_increment",ESI_increment);
	setParamaterWithLink(data,"eta",eta);
	setParamaterWithLink(data,"EtCS",EtCS);
	setParamaterWithLink(data,"EtID",EtID);
	setParamaterWithLink(data,"EtKG",EtKG);
	setParamaterWithLink(data,"EtMD",EtMD);
	setParamaterWithLink(data,"EtSDH",EtSDH);
	setParamaterWithLink(data,"FAD",FAD);
	setParamaterWithLink(data,"FADH2",FADH2);
	setParamaterWithLink(data,"fL",fL);
	setParamaterWithLink(data,"fm",fm);
	setParamaterWithLink(data,"fprime",fprime);
	setParamaterWithLink(data,"g",g);
	setParamaterWithLink(data,"G_Cab",G_Cab);
	setParamaterWithLink(data,"G_Kp",G_Kp);
	setParamaterWithLink(data,"G_Na",G_Na);
	setParamaterWithLink(data,"G_Nab",G_Nab);
	setParamaterWithLink(data,"gh",gh);
	setParamaterWithLink(data,"gL",gL);
	setParamaterWithLink(data,"GLU",GLU);
	setParamaterWithLink(data,"gprime",gprime);
	setParamaterWithLink(data,"H",H);
	setParamaterWithLink(data,"high_freq",high_freq);
	setParamaterWithLink(data,"hm",hm);
	setParamaterWithLink(data,"HTRPNtot",HTRPNtot);
	setParamaterWithLink(data,"ICahalf",ICahalf);
	setParamaterWithLink(data,"INaKmax",INaKmax);
	setParamaterWithLink(data,"IpCamax",IpCamax);
	setParamaterWithLink(data,"KAATeq",KAATeq);
	setParamaterWithLink(data,"KaCa",KaCa);
	setParamaterWithLink(data,"KACOeq",KACOeq);
	setParamaterWithLink(data,"kact",kact);
	setParamaterWithLink(data,"KADP",KADP);
	setParamaterWithLink(data,"kaminus",kaminus);
	setParamaterWithLink(data,"kaplus",kaplus);
	setParamaterWithLink(data,"kbminus",kbminus);
	setParamaterWithLink(data,"kbplus",kbplus);
	setParamaterWithLink(data,"Kca",Kca);
	setParamaterWithLink(data,"kcminus",kcminus);
	setParamaterWithLink(data,"kcnsASP",kcnsASP);
	setParamaterWithLink(data,"kcplus",kcplus);
	setParamaterWithLink(data,"KCS",KCS);
	setParamaterWithLink(data,"kf1",kf1);
	setParamaterWithLink(data,"kfAAT",kfAAT);
	setParamaterWithLink(data,"kfACO",kfACO);
	setParamaterWithLink(data,"Kfb",Kfb);
	setParamaterWithLink(data,"kfFH",kfFH);
	setParamaterWithLink(data,"KFHeq",KFHeq);
	setParamaterWithLink(data,"kfSL",kfSL);
	setParamaterWithLink(data,"kh_1",kh_1);
	setParamaterWithLink(data,"kh_2",kh_2);
	setParamaterWithLink(data,"Kh1",Kh1);
	setParamaterWithLink(data,"Kh2",Kh2);
	setParamaterWithLink(data,"Kh3",Kh3);
	setParamaterWithLink(data,"Kh4",Kh4);
	setParamaterWithLink(data,"khtrpn_minus",khtrpn_minus);
	setParamaterWithLink(data,"khtrpn_plus",khtrpn_plus);
	setParamaterWithLink(data,"Ki_AM",Ki_AM);
	setParamaterWithLink(data,"Ki_prime_SR",Ki_prime_SR);
	setParamaterWithLink(data,"Ki_SR",Ki_SR);
	setParamaterWithLink(data,"Ki1AD_NaK",Ki1AD_NaK);
	setParamaterWithLink(data,"KiADP_CaP",KiADP_CaP);
	setParamaterWithLink(data,"kIDH",kIDH);
	setParamaterWithLink(data,"KidhNADH",KidhNADH);
	setParamaterWithLink(data,"KiFUM",KiFUM);
	setParamaterWithLink(data,"Kioaa",Kioaa);
	setParamaterWithLink(data,"KiOxaa",KiOxaa);
	setParamaterWithLink(data,"kKGDH",kKGDH);
	setParamaterWithLink(data,"kltrpn_plus",kltrpn_plus);
	setParamaterWithLink(data,"kltrpn_minus",kltrpn_minus);
	setParamaterWithLink(data,"Km1AT_NaK",Km1AT_NaK);
	setParamaterWithLink(data,"Km1ATP_CaP",Km1ATP_CaP);
	setParamaterWithLink(data,"Km2ATP_CaP",Km2ATP_CaP);
	setParamaterWithLink(data,"KmAcCoA",KmAcCoA);
	setParamaterWithLink(data,"Kmal",Kmal);
	setParamaterWithLink(data,"KmATP_AM",KmATP_AM);
	setParamaterWithLink(data,"KmATP_SR",KmATP_SR);
	setParamaterWithLink(data,"KmCMDN",KmCMDN);
	setParamaterWithLink(data,"KmCa",KmCa);
	setParamaterWithLink(data,"KmCSQN",KmCSQN);
	setParamaterWithLink(data,"kMDH",kMDH);
	setParamaterWithLink(data,"Kmg",Kmg);
	setParamaterWithLink(data,"KmIDNAD",KmIDNAD);
	setParamaterWithLink(data,"Kmiso",Kmiso);
	setParamaterWithLink(data,"KmKG",KmKG);
	setParamaterWithLink(data,"KmKGNAD",KmKGNAD);
	setParamaterWithLink(data,"KmKo",KmKo);
	setParamaterWithLink(data,"KmmNAD",KmmNAD);
	setParamaterWithLink(data,"KmNa",KmNa);
	setParamaterWithLink(data,"KmNai",KmNai);
	setParamaterWithLink(data,"KmnsCa",KmnsCa);
	setParamaterWithLink(data,"KmOaa",KmOaa);
	setParamaterWithLink(data,"KmpCa",KmpCa);
	setParamaterWithLink(data,"KmSucc",KmSucc);
	setParamaterWithLink(data,"Kna",Kna);
	setParamaterWithLink(data,"kNaCa",kNaCa);
	setParamaterWithLink(data,"Knca",Knca);
	setParamaterWithLink(data,"Ko",Ko);
	setParamaterWithLink(data,"Koff",Koff);
	setParamaterWithLink(data,"Krb",Krb);
	setParamaterWithLink(data,"kres",kres);
	setParamaterWithLink(data,"kresf",kresf);
	setParamaterWithLink(data,"ksat",ksat);
//	setParamaterWithLink(data,"Fcik1",Fcik1);
	setParamaterWithLink(data,"kSDH",kSDH);
	setParamaterWithLink(data,"KSLeq",KSLeq);
	setParamaterWithLink(data,"KSR",KSR);
	setParamaterWithLink(data,"ktrans",ktrans);
	setParamaterWithLink(data,"kTrop_pn",kTrop_pn);
	setParamaterWithLink(data,"L",L);
	setParamaterWithLink(data,"LTRPNtot",LTRPNtot);
	setParamaterWithLink(data,"mcoop",mcoop);
	setParamaterWithLink(data,"Mg",Mg);
	setParamaterWithLink(data,"n",n);
	setParamaterWithLink(data,"na",na);
	setParamaterWithLink(data,"Nao",Nao);
	setParamaterWithLink(data,"ncoop",ncoop);
	setParamaterWithLink(data,"Nfb",Nfb);
	setParamaterWithLink(data,"nID",nID);
	setParamaterWithLink(data,"nKG",nKG);
	setParamaterWithLink(data,"norm_freq",norm_freq);
	setParamaterWithLink(data,"Nrb",Nrb);
	setParamaterWithLink(data,"omega",omega);
	setParamaterWithLink(data,"p1",p1);
	setParamaterWithLink(data,"p2",p2);
	setParamaterWithLink(data,"p3",p3);
	setParamaterWithLink(data,"pa",pa);
	setParamaterWithLink(data,"pb",pb);
	setParamaterWithLink(data,"pc1",pc1);
	setParamaterWithLink(data,"pc2",pc2);
	setParamaterWithLink(data,"PCa",PCa);
	setParamaterWithLink(data,"period",period);
	setParamaterWithLink(data,"PESI",PESI);
	setParamaterWithLink(data,"Pi",Pi);
	setParamaterWithLink(data,"PK",PK);
	setParamaterWithLink(data,"PnsK",PnsK);
	setParamaterWithLink(data,"PnsNa",PnsNa);
	setParamaterWithLink(data,"pulse_amplitude",pulse_amplitude);
	setParamaterWithLink(data,"pulse_duration",pulse_duration);
	setParamaterWithLink(data,"r1",r1);
	setParamaterWithLink(data,"r2",r2);
	setParamaterWithLink(data,"r3",r3);
	setParamaterWithLink(data,"ra",ra);
	setParamaterWithLink(data,"rb",rb);
	setParamaterWithLink(data,"rc1",rc1);
	setParamaterWithLink(data,"rc2",rc2);
	setParamaterWithLink(data,"refrac_buffer",refrac_buffer);
	setParamaterWithLink(data,"rhoF1",rhoF1);
	setParamaterWithLink(data,"rhoREF",rhoREF);
	setParamaterWithLink(data,"rhoREN",rhoREN);
	setParamaterWithLink(data,"shift",shift);
	setParamaterWithLink(data,"t1",t1);
	setParamaterWithLink(data,"t2",t2);
	setParamaterWithLink(data,"tautr",tautr);
	setParamaterWithLink(data,"tauxfer",tauxfer);
	setParamaterWithLink(data,"time_vclamp_off",time_vclamp_off);
	setParamaterWithLink(data,"time_vclamp_on",time_vclamp_on);
	setParamaterWithLink(data,"V_AM_scaler",V_AM_scaler);
	setParamaterWithLink(data,"V_AM_max",V_AM_max);
	setParamaterWithLink(data,"v1",v1);
	setParamaterWithLink(data,"vclamp_hold",vclamp_hold);
	setParamaterWithLink(data,"vclamp_set",vclamp_set);
	setParamaterWithLink(data,"VJSR",VJSR);
	setParamaterWithLink(data,"vmaxf",vmaxf);
	setParamaterWithLink(data,"vmaxr",vmaxr);
	setParamaterWithLink(data,"VmDT",VmDT);
	setParamaterWithLink(data,"VmNC",VmNC);
	setParamaterWithLink(data,"Vmuni",Vmuni);
	setParamaterWithLink(data,"Vmyo",Vmyo);
	setParamaterWithLink(data,"VNSR",VNSR);
	setParamaterWithLink(data,"VSS",VSS);
	setParamaterWithLink(data,"zeta",zeta);
	setParamaterWithLink(data,"stopTime",stopTime);
	setParamaterWithLink(data,"numRun",numRun);
	setParamaterWithLink(data,"t3",t3);
	//CK Mod
	setParamaterWithLink(data,"kt_2",kt_2);
	setParamaterWithLink(data,"kf_2",kf_2);
	setParamaterWithLink(data,"kf_3",kf_3);
	setParamaterWithLink(data,"keq",keq);
	setParamaterWithLink(data,"CRT_cyto",CRT_cyto);
	setParamaterWithLink(data,"CRT_mito",CRT_mito);
	setParamaterWithLink(data,"VATPase_cyto",VATPase_cyto);

	setParamaterWithLink(data,"f_xb",f_xb);
	setParamaterWithLink(data,"SL",SL);
	setParamaterWithLink(data,"gmin_xb",gmin_xb);

	//Precalculate varibales
	calculateConstantModelParameters();

	//Store Flags
	usingAlgebraicMembranePotential	= data.getBoolean("usingAlgebraicMembranePotential");
	clampVoltage					= data.getBoolean("clampVoltage");
	usingCHF						= data.getBoolean("usingCHF");
	usingASP						= data.getBoolean("usingASP");
	usingCK							= data.getBoolean("usingCK");
	
	stimulasMode = stimulasModeLookup[ data.getKeyword("StimulasMode") ];
	modelType    =    modelTypeLookup[ data.getKeyword("ModelType")    ];
}
//The veriable stop time values given the mode
double Model::getStopTime() {
	switch(stimulasMode) {
		case BB:
			return t3;
			break;
		case IF:
			return stop_time;
			break;
		default:
			return stopTime;
			break;
	}
	return stopTime;
}
//IF mode parameter, the number of trial runs; 1 in other cases
double Model::getNumRun() {
	if(stimulasMode == IF)
		return numRun;
	return 1;
}
//Loads data and sets the lookup
void Model::setParamaterWithLink(fileIO& data, const string &name, double &param) {
	parameterLookup[name] = &param;
	param = data[name];
}
//Simply set a parameter
void Model::setParamater(const string &name, double value) {
	parameterLookupType::iterator p = parameterLookup.find(name);
	if (p == parameterLookup.end()) {
		cout << "An error has occured in 'Model::setParamater()': " << name << " parameter was not found." << endl;
		exit(-1);
	}
	*(p->second) = value;
	calculateConstantModelParameters();
	return;
}
//Create an entry in parameter lookup
void Model::linkParameter(const string &name, double &param) {
	parameterLookup[name] = &param;
}
const char * Model::getStateLabel(int index) {
	index -= getIndexOffset(); //Raw index
	return getInitialConditionLabel(index);
}
const char * Model::getInitialConditionLabel(int index) {
    if((index == index_ASP)&&usingASP) //The overloaded one
		return "ASP";
	//Otherwise use lookup:
	if( (index < 0) || (index > MAXSTATES ) ) {
		cerr << "Error in Model::getInitialConditionLabel: index was:" << index <<"."<<endl;
		return "";
	}
	return stateLabelLookup[index].c_str();
}
void Model::setStateWithLink(fileIO& data, const string &name, int index) {
	stateLabelLookup[index] = name;
	stateIndexLookup[name] = index;
    s[index] = data[name];
}
void Model::setInitialConditions(fileIO& data) {
	setStateWithLink(data,"V",index_V);
	setStateWithLink(data,"mNa",index_mNa);
	setStateWithLink(data,"hNa",index_hNa);
	setStateWithLink(data,"jNa",index_jNa);
	setStateWithLink(data,"xKs",index_xKs);
	setStateWithLink(data,"Nai",index_Nai);
	setStateWithLink(data,"Ki",index_Ki);
	setStateWithLink(data,"Cai",index_Cai);
	setStateWithLink(data,"CaNSR",index_CaNSR);
	setStateWithLink(data,"CaSS",index_CaSS);
	setStateWithLink(data,"CaJSR",index_CaJSR);
	setStateWithLink(data,"C1_RyR",index_C1_RyR);
	setStateWithLink(data,"O2_RyR",index_O2_RyR);
	setStateWithLink(data,"C2_RyR",index_C2_RyR);
	setStateWithLink(data,"C0",index_C0);
	setStateWithLink(data,"C1",index_C1);
	setStateWithLink(data,"C2",index_C2);
	setStateWithLink(data,"C3",index_C3);
	setStateWithLink(data,"C4",index_C4);
	setStateWithLink(data,"Open",index_Open);
	setStateWithLink(data,"CCa0",index_CCa0);
	setStateWithLink(data,"CCa1",index_CCa1);
	setStateWithLink(data,"CCa2",index_CCa2);
	setStateWithLink(data,"CCa3",index_CCa3);
	setStateWithLink(data,"CCa4",index_CCa4);
	setStateWithLink(data,"OCa",index_OCa);
	setStateWithLink(data,"yCa",index_yCa);
	setStateWithLink(data,"HTRPNCa",index_HTRPNCa);
	//Force model
	setStateWithLink(data,"LTRPNCa",index_LTRPNCa);
	setStateWithLink(data,"N0",index_N0);
	setStateWithLink(data,"N1",index_N1);
	setStateWithLink(data,"P0",index_P0);
	setStateWithLink(data,"P1",index_P1);
	setStateWithLink(data,"P2",index_P2);
	setStateWithLink(data,"P3",index_P3);
	//after ATP mod
	setStateWithLink(data,"ATPi",index_ATPi);
	//after Mitochondria mod
	//Mitochondria model
	setStateWithLink(data,"Cam",index_Cam);
	setStateWithLink(data,"ADPm",index_ADPm);
	setStateWithLink(data,"Dpsi",index_Dpsi);
	setStateWithLink(data,"NADH",index_NADH);
	setStateWithLink(data,"ISOC",index_ISOC);
	setStateWithLink(data,"AKG",index_AKG);
	setStateWithLink(data,"SCoA",index_SCoA);
	setStateWithLink(data,"Succ",index_Succ);
	setStateWithLink(data,"FUM",index_FUM);
	setStateWithLink(data,"MAL",index_MAL);
	setStateWithLink(data,"Oaa",index_Oaa);
	//after ATP & creatine kinase mod
	setStateWithLink(data,"ATPi_cyto",index_ATPi_cyto);
	setStateWithLink(data,"CrPi_cyto",index_CrPi_cyto);
	setStateWithLink(data,"CrPi_mito",index_CrPi_mito);
	if (usingASP) {
		setStateWithLink(data,"ASP",index_ASP);
	}
}
const double* Model::getInitialConditions(const double *y) {
	//Set all "active" parameters to the state array in y
	for (int i = 0; i < getProblemSize(); i++) {
		s[i + getIndexOffset()] = y[i];
	}
	return s;
}
int Model::getInitialConditionsSize() {
	return MAXSTATES;
}
int Model::getStateIndex(const char * label) { return stateIndexLookup[label] - getIndexOffset(); }
int Model::getDependentVariableIndex(const char * label) {
	return dependentVariableIndexLookup[label];
}
double Model::getDependentVariable(int index) {
	switch(index) {
		case 1:		return INa;
		case 2:		return IKs;
		case 3:		return IK1;
		case 4:		return INab;
		case 5:		return IKp;
		case 6:		return ICa;
		case 7:		return ICaK;
		case 8:		return INaK;
		case 9:		return (INaK*6.3710E-5);	//RNak
		case 10:	return InsCa;
		case 11:	return INaCa;
		case 12:	return ICab;
		case 13:	return IpCa;
		case 14:	return (IpCa*6.3710E-5);	//RpCa
		case 15:	return Jup;
		case 16:	return Jrel;
		case 17:	return Jtr;
		case 18:	return Jxfer;
		case 19:	return V_AM;
		case 20:	return VNO;
		case 21:	return VHNe;
		case 22:	return VHFe;
		case 23:	return VANT;
		case 24:    return Vhu;
		case 25:	return VATPase;
		case 26:	return VnaCa;
		case 27:	return Vuni;
		case 28:	return VIDH;
		case 29:	return VKGDH;
		case 30:	return Vhleak;
		case 31:	return alpha_SL_fnormmax2 * P1_N1_P2_P3;								//FN_Ca
		case 32:	return alpha_SL_fnormmax * ( P1_N1_P2_P3 + (*P2) + (*P3) + (*P3));		//Fnorm
		case 33:	return zeta_alpha_SL_fnormmax * ( P1_N1_P2_P3 + (*P2) + (*P3) + (*P3));	//Fnorm * Zeta = Force;
		default:	return 0;
	}
}
