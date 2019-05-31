#pragma once

/*     ----------------------------------------------------

         NOTICE OF COPYRIGHT AND OWNERSHIP OF SOFTWARE

         Copyright 2004, The Johns Hopkins University
            School of Medicine. All rights reserved.
			For research use only; commercial use prohibited.
			Distribution without permission of Raimond L. Winslow
			not permitted. rwinslow@bme.jhu.edu

         Name of Program: Guinea Pig C++ (Coupled) (GPC_Coupled)
         Version: Documented Version, version 1.0
         Date: February 2004

       -----------------------------------------------------  
*/
/*	   This file contains Index definitions 
	   that should correspond with states variables and
	   initial conditions.

	   The reason for using indexing like this, is to make
	   changing order of and adding of states easier. This
	   should be the only place where state numbers are referred
	   to explicitly.

	   It also contains fixed defined physical constants and any hard 
	   coded parameters.

*/
//Switch to parameter reading for these where possible
#define	RT_over_F 26.708186528497409326424870466321
#define Faraday	96.5
#define MAXSTATES 50	//Was 200, redifined as the total number of states
#define LHospitalThreshold 1e-7

//Problem size definitions
#define PROBLEM_SIZE_GPC 47
#define PROBLEM_SIZE_GPC_ALG 46
#define PROBLEM_SIZE_GPC_CK 50
#define PROBLEM_SIZE_GPC_CK_ALG 49
#define PROBLEM_SIZE_FORCE 7
#define PROBLEM_SIZE_MITO 11
#define PROBLEM_SIZE_MITO_ASP 12

//Must be the State maximum/minimum, for algebraic mode
#define index_V     0	
#define index_Alg_Offset	1
#define index_mNa   1    
#define index_hNa   2  
#define index_jNa   3  
#define index_xKs   4 
#define index_Nai   5  
#define index_Ki    6    
#define index_Cai   7  
#define index_CaNSR  8  
#define index_CaSS   9    
#define index_CaJSR  10  
#define index_C1_RyR   11  
#define index_O2_RyR   12
#define index_C2_RyR   13  
#define index_C0   14  
#define index_C1   15  
#define index_C2   16  
#define index_C3   17  
#define index_C4   18 
#define index_Open  19  
#define index_CCa0  20  
#define index_CCa1  21  
#define index_CCa2  22  
#define index_CCa3  23  
#define index_CCa4  24  
#define index_OCa   25
#define index_yCa	26  
#define index_HTRPNCa  27
//Force model
#define index_Force_Offset	28
#define index_LTRPNCa  28  
#define index_N0       29
#define index_N1       30
#define index_P0       31
#define index_P1       32
#define index_P2       33
#define index_P3       34

/**after ATP mod**/
#define index_ATPi		35

/**after Mitochondria mod**/
//Mitochondria model
#define index_Mito_Offset	36
#define index_Cam      36
#define index_ADPm     37
#define index_Dpsi     38
#define index_NADH     39
#define index_ISOC     40
#define index_AKG      41
#define index_SCoA     42
#define index_Succ     43
#define index_FUM      44
#define index_MAL      45
#define index_Oaa      46
#define index_ASP      47

/**after ATP & creatine kinase mod**/
#define index_ATPi_cyto	47	//Uses ASP's index since ASP isn't in the coupled model
#define index_CrPi_cyto	48
#define index_CrPi_mito	49

