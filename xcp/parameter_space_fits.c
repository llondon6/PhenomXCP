
/*
* Copyright (C) 2021 The University of Amsterdam
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifdef __cplusplus
extern "C" {
#endif

/**
* \author Lionel London
*/

//
#include "LALSimIMRPhenomX_PNR_deviations.h"

#ifndef _OPENMP
#define omp ignore
#endif

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#else
#define DEBUG 1
#endif


// Import usefuls
#include <math.h>



// MU1 fit implementation 
// Header formatting for MU1 of (l,m)=(2,2) multipole
double IMRPhenomXCP_MU1_l2m2( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	UNUSED double eta3 = eta2*eta;
	double mu1;

	// Evaluate fit for this parameter
	mu1 = -2.46767869e+01*a1 + 9.30143130e-01*u + -4.98153241e+01*eta + 2.44125029e+00 + 6.25481470e+01*a12 + -7.52831252e+00*(u*a1) + 3.21253969e+02*eta2 + 4.94048965e+02*(eta*a1) + -1.72954088e+01*(u*eta*a1) + -3.14001680e+03*(eta2*a1) + -4.74367603e+01*(a12*a1) + 2.90350868e+01*(u2*a1) + 2.67057418e+01*(u*a12) + -4.62899701e+01*(u*eta2) + 6.08834111e+00*(u2*eta) + -6.37477667e+02*(eta2*eta) + -1.22651507e+03*(eta*a12) + -1.55135335e+00*(u2*u) + -9.64898066e+01*(u2*a12) + -1.06500864e+02*(u2*eta2) + -1.53653214e+02*(u*eta*a12) + 7.69572698e+03*(eta2*a12) + -5.99225072e+02*(u2*eta*a1) + 9.12018348e+02*(eta*a13) + -2.52797705e+01*(u*a13) + 1.24873540e+02*(u*eta2*eta) + 1.43241326e+01*(u3*a1) + -7.57957572e+00*(u3*u) + 6.20269493e+03*(eta3*a1) + 5.23430627e+02*(u*eta2*a1) + -6.68469672e+00*(u3*eta) + -1.50864536e+04*(eta3*a12) + 1.79920464e+03*(u2*eta*a12) + -1.34498926e+03*(u*eta3*a1) + 2.86180007e+02*(u2*eta2*eta) + 4.19593228e+03*(u2*eta2*a1) + 1.50214819e+02*(u3*eta2) + 2.06993585e+01*(u4*a1) + -4.00346718e+01*(u3*a12) + 2.58944677e+02*(u*eta*a13) + 7.91531312e+01*(u2*a13) + -5.63322785e+03*(eta2*a13) + 1.37367328e+02*(u4*eta) + 3.32395606e+01*(u3*a13) + -3.21891497e+02*(u4*eta*a1) + 8.17990485e+02*(u*eta3*a12) + -8.35580724e+02*(u3*eta2*a1) + 1.09227160e+04*(eta3*a13) + -1.15721758e+04*(u2*eta2*a12) + -1.52755998e+01*(u4*a12) + -9.68995987e+02*(u*eta2*a13) + -8.66270395e+03*(u2*eta3*a1) + -3.92763193e+02*(u3*eta2*eta) + -7.39243221e+02*(u4*eta2) + -1.38972015e+03*(u2*eta*a13) + 1.94435826e+02*(u3*eta*a12) + 2.56930940e+02*(u4*eta*a12) + 8.47562593e+03*(u2*eta2*a13) + 5.65531115e+02*(u3*eta2*a12) + 1.31159567e+03*(u4*eta2*eta) + -2.82575295e+02*(u3*eta*a13) + 2.26358931e+04*(u2*eta3*a12) + 2.44134813e+03*(u3*eta3*a1) + 1.18153425e+03*(u4*eta2*a1) + 1.28448387e+03*(u*eta3*a13) + -1.59463281e+04*(u2*eta3*a13) + -4.53004012e+02*(u4*eta2*a12) + -1.43416481e+03*(u4*eta3*a1) + 6.26541854e+02*(u3*eta2*a13) + -2.92340145e+03*(u3*eta3*a12) + -6.06105289e+01*(u4*eta*a13);

	// Return answer
	return mu1;

} // END of MU1 (2,2) fit implementation



// MU2 fit implementation 
// Header formatting for MU2 of (l,m)=(2,2) multipole
double IMRPhenomXCP_MU2_l2m2( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	UNUSED double eta3 = eta2*eta;
	double mu2;

	// Evaluate fit for this parameter
	mu2 = 1.43416530e+01*a1 + 2.28304527e+01*u + 4.57882049e+01*eta + -9.62461763e-01 + -2.74706947e+01*a12 + -4.36879236e+02*(u*eta) + -1.29887511e+02*(u*a1) + -4.50197958e+02*eta2 + -2.64724113e+01*(u*u) + -4.75325985e+02*(eta*a1) + 2.37673156e+03*(u*eta*a1) + 4.07761971e+03*(eta2*a1) + 1.61940577e+01*(a12*a1) + 6.86600976e+01*(u2*a1) + 2.39004469e+02*(u*a12) + 2.73992366e+03*(u*eta2) + 4.41571032e+02*(u2*eta) + 1.13407678e+03*(eta2*eta) + 1.06598199e+03*(eta*a12) + -1.57151064e+01*(u2*u) + -7.37190314e+01*(u2*a12) + -2.29838597e+03*(u2*eta2) + -4.29040928e+03*(u*eta*a12) + -9.28314015e+03*(eta2*a12) + -7.90178586e+02*(u2*eta*a1) + -7.24996452e+02*(eta*a13) + -1.45993899e+02*(u*a13) + -5.57352428e+03*(u*eta2*eta) + 4.54480348e+01*(u3*a1) + 2.00417205e+01*(u3*u) + -9.67310407e+03*(eta3*a1) + -1.45414110e+04*(u*eta2*a1) + 3.32576830e+02*(u3*eta) + 2.59041167e+04*(u*eta2*a12) + 2.19510296e+04*(eta3*a12) + -9.50019417e+02*(u3*eta*a1) + 2.92390378e+04*(u*eta3*a1) + 4.00320120e+03*(u2*eta2*eta) + 1.83726216e+03*(u2*eta2*a1) + -2.33279776e+03*(u3*eta2) + 2.61142846e+03*(u*eta*a13) + 3.78821488e+01*(u2*a13) + 6.42782524e+03*(eta2*a13) + -4.31776717e+02*(u4*eta) + -2.87049188e+01*(u3*a13) + 4.63991957e+02*(u4*eta*a1) + -5.16941385e+04*(u*eta3*a12) + 7.42193306e+03*(u3*eta2*a1) + -1.52191121e+04*(eta3*a13) + 7.09659691e+03*(u2*eta2*a12) + -6.60180515e+01*(u4*a12) + -1.56050805e+04*(u*eta2*a13) + 5.23257920e+03*(u3*eta2*eta) + 2.61274224e+03*(u4*eta2) + 4.73200982e+02*(u2*eta*a13) + 2.94841216e+02*(u4*eta*a12) + -8.24417220e+03*(u2*eta2*a13) + -2.90495497e+03*(u3*eta2*a12) + 2.50526397e+01*(u4*a13) + -5.04700806e+03*(u4*eta2*eta) + 6.44153442e+02*(u3*eta*a13) + -2.21661627e+04*(u2*eta3*a12) + -1.85972927e+04*(u3*eta3*a1) + -3.31677967e+03*(u4*eta2*a1) + 3.08296315e+04*(u*eta3*a13) + 2.21734954e+04*(u2*eta3*a13) + -5.52307258e+02*(u4*eta2*a12) + 6.90574273e+03*(u4*eta3*a1) + -2.44634268e+03*(u3*eta2*a13) + 1.38173234e+04*(u3*eta3*a12) + -4.49106888e+01*(u4*eta*a13);

	// Return answer
	return mu2;

} // END of MU2 (2,2) fit implementation



// MU3 fit implementation 
// Header formatting for MU3 of (l,m)=(2,2) multipole
double IMRPhenomXCP_MU3_l2m2( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	UNUSED double eta3 = eta2*eta;
	double mu3;

	// Evaluate fit for this parameter
	mu3 = -1.58314450e-01*a1 + -3.58498212e-01*u + -2.09565821e-01*eta + 2.35714985e-02 + 5.58441939e-01*a12 + 4.48990294e+00*(u*eta) + 2.92076718e+00*(u*a1) + -3.73241355e-01*(u*u) + 8.75234093e-01*(eta*a1) + -3.50118215e+01*(u*eta*a1) + 5.02871430e+00*(eta2*a1) + -4.03290110e-01*(a12*a1) + 1.79899846e+00*(u2*a1) + -6.84519742e+00*(u*a12) + -1.65435023e+01*(u*eta2) + 5.08799936e+00*(u2*eta) + 1.57541095e+00*(eta2*eta) + -4.53874143e+00*(eta*a12) + 3.95180784e-01*(u2*u) + -3.53443719e+00*(u2*a12) + -2.45438501e+01*(u2*eta2) + 8.04952265e+01*(u*eta*a12) + -2.04817851e+01*(u2*eta*a1) + 3.29476877e+00*(eta*a13) + 4.76512744e+00*(u*a13) + 1.67765212e+01*(u*eta2*eta) + -3.24203757e+00*(u3*a1) + 3.21514884e-01*(u3*u) + -2.28585808e+01*(eta3*a1) + 1.20514568e+02*(u*eta2*a1) + -4.26953734e+00*(u3*eta) + -2.69029088e+02*(u*eta2*a12) + 3.50432971e+01*(eta3*a12) + 3.83050763e+01*(u2*eta*a12) + 3.42998970e+01*(u3*eta*a1) + -1.04653092e+02*(u*eta3*a1) + 4.09622583e+01*(u2*eta2*eta) + 8.46359010e+01*(u2*eta2*a1) + 1.09459412e+01*(u3*eta2) + -1.18448266e+00*(u4*a1) + 7.76016965e+00*(u3*a12) + -5.53075006e+01*(u*eta*a13) + 2.40863061e+00*(u2*a13) + -4.28757831e+00*(u4*eta) + -5.49260181e+00*(u3*a13) + 9.97725717e+00*(u4*eta*a1) + 2.16457898e+02*(u*eta3*a12) + -8.44258024e+01*(u3*eta2*a1) + -2.58895176e+01*(eta3*a13) + -1.49027956e+02*(u2*eta2*a12) + 1.69781306e+00*(u4*a12) + 1.81513728e+02*(u*eta2*a13) + -1.28167221e+02*(u2*eta3*a1) + 2.11598683e+01*(u4*eta2) + -2.64604066e+01*(u2*eta*a13) + -8.21518275e+01*(u3*eta*a12) + -7.19724029e+00*(u4*eta*a12) + 1.03094926e+02*(u2*eta2*a13) + 2.04878729e+02*(u3*eta2*a12) + -1.11571122e+00*(u4*a13) + -3.73469809e+01*(u4*eta2*eta) + 5.82632419e+01*(u3*eta*a13) + 2.13992763e+02*(u2*eta3*a12) + -1.00893701e+01*(u3*eta3*a1) + -3.50245564e+01*(u4*eta2*a1) + -1.39324365e+02*(u*eta3*a13) + -1.45832151e+02*(u2*eta3*a13) + 5.98872767e+01*(u4*eta3*a1) + -1.47805007e+02*(u3*eta2*a13) + 1.35506328e+01*(u3*eta3*a12) + 4.68529851e+00*(u4*eta*a13);

	// Return answer
	return mu3;

} // END of MU3 (2,2) fit implementation



// NU4 fit implementation 
// Header formatting for NU4 of (l,m)=(2,2) multipole
double IMRPhenomXCP_NU4_l2m2( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	UNUSED double eta3 = eta2*eta;
	double nu4;

	// Evaluate fit for this parameter
	nu4 = -1.38848648e-01*a1 + -2.15222655e-02*u + -2.69106927e-01*eta + 1.15967006e-02 + 3.55249490e-01*a12 + 2.12351256e-01*(u*eta) + 2.02065381e-01*(u*a1) + 1.79212381e+00*eta2 + 1.05083328e-02*(u*u) + 3.03588902e+00*(eta*a1) + -2.34938640e+00*(u*eta*a1) + -1.96632904e+01*(eta2*a1) + -2.61775484e-01*(a12*a1) + 2.19689058e-01*(u2*a1) + -5.71111734e-01*(u*a12) + -4.97195206e-01*(u*eta2) + -3.69311682e+00*(eta2*eta) + -7.80838651e+00*(eta*a12) + 1.85062854e-02*(u2*u) + -6.23447554e-01*(u2*a12) + -6.76752295e-01*(u2*eta2) + 7.75122855e+00*(u*eta*a12) + 5.05379625e+01*(eta2*a12) + -5.66023676e+00*(u2*eta*a1) + 5.82645485e+00*(eta*a13) + 4.31820204e-01*(u*a13) + -1.83488568e-01*(u3*a1) + -4.81004336e-02*(u3*u) + 3.93928460e+01*(eta3*a1) + 8.66276258e+00*(u*eta2*a1) + -3.60491772e+01*(u*eta2*a12) + -1.00717143e+02*(eta3*a12) + 1.57373683e+01*(u2*eta*a12) + 7.96260402e-01*(u3*eta*a1) + -1.03064945e+01*(u*eta3*a1) + 2.20973581e+00*(u2*eta2*eta) + 3.77807023e+01*(u2*eta2*a1) + -1.44112349e+00*(u3*eta2) + 7.40700826e-02*(u4*a1) + 4.71349707e-01*(u3*a12) + -6.37180005e+00*(u*eta*a13) + 4.59549617e-01*(u2*a13) + -3.79701270e+01*(eta2*a13) + 7.56381075e-01*(u4*eta) + -3.37436659e-01*(u3*a13) + -6.64137423e-01*(u4*eta*a1) + 5.76029282e+01*(u*eta3*a12) + 5.55511068e+00*(u3*eta2*a1) + 7.58580099e+01*(eta3*a13) + -1.03981743e+02*(u2*eta2*a12) + -8.44071341e-02*(u4*a12) + 3.23212792e+01*(u*eta2*a13) + -7.63331591e+01*(u2*eta3*a1) + 4.55267993e+00*(u3*eta2*eta) + -4.16084724e+00*(u4*eta2) + -1.16921696e+01*(u2*eta*a13) + -3.34980622e+00*(u3*eta*a12) + 7.76004667e+01*(u2*eta2*a13) + 5.73720971e-02*(u4*a13) + 7.30707580e+00*(u4*eta2*eta) + 2.99664221e+00*(u3*eta*a13) + 2.08041779e+02*(u2*eta3*a12) + -2.28271919e+01*(u3*eta3*a1) + 3.68426289e+00*(u4*eta2*a1) + -5.56216015e+01*(u*eta3*a13) + -1.55508103e+02*(u2*eta3*a13) + -6.44341734e+00*(u4*eta3*a1) + -6.36112798e+00*(u3*eta2*a13) + 2.22524904e+01*(u3*eta3*a12);

	// Return answer
	return nu4;

} // END of NU4 (2,2) fit implementation



// NU5 fit implementation 
// Header formatting for NU5 of (l,m)=(2,2) multipole
double IMRPhenomXCP_NU5_l2m2( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	UNUSED double eta3 = eta2*eta;
	double nu5;

	// Evaluate fit for this parameter
	nu5 = 1.45981410e-01*u + 5.74729648e-02*eta + -5.81258283e-03 + -2.33415978e-01*a12 + -1.21584615e+00*(u*eta) + -1.36069478e+00*(u*a1) + -3.65415485e-01*eta2 + 6.17467603e-02*(u*u) + -1.43432438e-02*(eta*a1) + 1.22089540e+01*(u*eta*a1) + 1.93805509e-01*(a12*a1) + 3.43252084e+00*(u*a12) + -1.00499733e+00*(u2*eta) + 8.95502244e-01*(eta2*eta) + 3.07410256e+00*(eta*a12) + -1.99305986e-01*(u2*u) + 5.31839109e+00*(u2*eta2) + -3.19569789e+01*(u*eta*a12) + -1.49896010e+01*(eta2*a12) + -2.58206987e+00*(eta*a13) + -2.51537037e+00*(u*a13) + 1.05975326e+01*(u*eta2*eta) + 2.02014136e+00*(u3*a1) + -1.09701049e-01*(u3*u) + -1.05859921e+01*(u*eta2*a1) + 1.70425675e+00*(u3*eta) + 4.00782555e+01*(u*eta2*a12) + 2.61484129e+01*(eta3*a12) + -2.03539191e+01*(u3*eta*a1) + -6.96745845e+01*(u*eta3*a1) + -9.19965335e+00*(u2*eta2*eta) + 1.97782716e-01*(u4*a1) + -5.24360162e+00*(u3*a12) + 2.37938807e+01*(u*eta*a13) + -1.06968335e-01*(u2*a13) + 1.16149831e+01*(eta2*a13) + 1.81665708e+00*(u4*eta) + 3.91740860e+00*(u3*a13) + -3.09258495e+00*(u4*eta*a1) + 1.39765189e+02*(u*eta3*a12) + 3.49558250e+01*(u3*eta2*a1) + -1.79445989e+01*(eta3*a13) + -3.54945028e+01*(u*eta2*a13) + -1.51848855e+01*(u3*eta2*eta) + -9.78982098e+00*(u4*eta2) + 1.50085783e+00*(u2*eta*a13) + 5.62030646e+01*(u3*eta*a12) + -6.50085555e+00*(u2*eta2*a13) + -1.27369903e+02*(u3*eta2*a12) + 1.71834912e+01*(u4*eta2*eta) + -4.36278264e+01*(u3*eta*a13) + 6.18664375e+01*(u3*eta3*a1) + 1.66651908e+01*(u4*eta2*a1) + -8.34091320e+01*(u*eta3*a13) + 8.84594982e+00*(u2*eta3*a13) + -2.97830232e+01*(u4*eta3*a1) + 1.13887582e+02*(u3*eta2*a13) + -6.63416139e+01*(u3*eta3*a12);

	// Return answer
	return nu5;

} // END of NU5 (2,2) fit implementation



// NU6 fit implementation 
// Header formatting for NU6 of (l,m)=(2,2) multipole
double IMRPhenomXCP_NU6_l2m2( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	UNUSED double eta3 = eta2*eta;
	double nu6;

	// Evaluate fit for this parameter
	nu6 = -3.73708021e-01*a1 + -2.33908226e-01*u + -1.04751149e+00*eta + 4.44310872e-02 + 7.01827004e-01*a12 + 4.46812265e+00*(u*eta) + 1.33151778e+00*(u*a1) + 7.26057425e+00*eta2 + 2.32268400e-01*(u*u) + 8.67059989e+00*(eta*a1) + -2.42825702e+01*(u*eta*a1) + -5.95412509e+01*(eta2*a1) + -3.87831714e-01*(a12*a1) + -1.05824303e+00*(u2*a1) + -2.41614099e+00*(u*a12) + -2.77210808e+01*(u*eta2) + -3.79356363e+00*(u2*eta) + -1.52576222e+01*(eta2*eta) + -1.76114743e+01*(eta*a12) + 1.07299519e-01*(u2*u) + 2.42563592e+00*(u2*a12) + 2.07197233e+01*(u2*eta2) + 4.31317947e+01*(u*eta*a12) + 1.24643868e+02*(eta2*a12) + 1.30566192e+01*(u2*eta*a1) + 1.06551764e+01*(eta*a13) + 1.41030562e+00*(u*a13) + 5.58594485e+01*(u*eta2*eta) + -1.05279296e-01*(u3*u) + 1.24308261e+02*(eta3*a1) + 1.45843407e+02*(u*eta2*a1) + -2.87146711e+00*(u3*eta) + -2.54614075e+02*(u*eta2*a12) + -2.63866348e+02*(eta3*a12) + -2.54870849e+01*(u2*eta*a12) + 5.25128813e+00*(u3*eta*a1) + -2.88154658e+02*(u*eta3*a1) + -3.89601066e+01*(u2*eta2*eta) + -5.11056091e+01*(u2*eta2*a1) + 2.17894497e+01*(u3*eta2) + -2.29553344e-01*(u4*a1) + -1.13534789e+00*(u3*a12) + -2.52360001e+01*(u*eta*a13) + -1.88735183e+00*(u2*a13) + -7.84463911e+01*(eta2*a13) + 2.71212603e+00*(u4*eta) + 1.04799358e+00*(u3*a13) + 4.97329529e+02*(u*eta3*a12) + -5.91969038e+01*(u3*eta2*a1) + 1.69517099e+02*(eta3*a13) + 6.96314126e+01*(u2*eta2*a12) + 7.16852536e-01*(u4*a12) + 1.48002235e+02*(u*eta2*a13) + 7.18661282e+01*(u2*eta3*a1) + -5.03076038e+01*(u3*eta2*eta) + -1.95211585e+01*(u4*eta2) + 1.89465601e+01*(u2*eta*a13) + 9.84749436e+00*(u3*eta*a12) + -6.60985495e+00*(u4*eta*a12) + -4.43484038e+01*(u2*eta2*a13) + -9.84327592e-02*(u4*a13) + 4.38398261e+01*(u4*eta2*eta) + -1.20835825e+01*(u3*eta*a13) + -3.79845120e+01*(u2*eta3*a12) + 1.64851484e+02*(u3*eta3*a1) + 1.69725918e+01*(u4*eta2*a1) + -2.86904520e+02*(u*eta3*a13) + 1.88005054e+01*(u4*eta2*a12) + -6.41232204e+01*(u4*eta3*a1) + 3.47421666e+01*(u3*eta2*a13) + -1.06500824e+02*(u3*eta3*a12);

	// Return answer
	return nu6;

} // END of NU6 (2,2) fit implementation



// ZETA1 fit implementation 
// Header formatting for ZETA1 of (l,m)=(2,2) multipole
double IMRPhenomXCP_ZETA1_l2m2( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	UNUSED double eta3 = eta2*eta;
	double zeta1;

	// Evaluate fit for this parameter
	zeta1 = 7.19903828e-05*a1 + 3.68830495e-05*u + -3.67068197e-04*eta + 3.57245526e-05 + -4.19803606e-04*a12 + -7.50721516e-04*(u*eta) + -1.82997391e-04*(u*a1) + 9.15091305e-04*eta2 + -4.32569180e-04*(u*u) + -3.77190605e-03*(eta*a1) + 4.57008622e-03*(u*eta*a1) + 3.33902142e-02*(eta2*a1) + 4.06275340e-04*(a12*a1) + 1.25643527e-03*(u2*a1) + 2.40692269e-04*(u*a12) + 4.72198789e-03*(u*eta2) + 7.04872294e-03*(u2*eta) + 1.43212631e-02*(eta*a12) + -1.00647316e-04*(u2*u) + -1.35969854e-03*(u2*a12) + -3.93712504e-02*(u2*eta2) + -7.78395207e-03*(u*eta*a12) + -1.13058512e-01*(eta2*a12) + -1.34352285e-02*(u2*eta*a1) + -1.25180209e-02*(eta*a13) + -1.01321063e-02*(u*eta2*eta) + 4.40550829e-04*(u3*a1) + 3.84892395e-04*(u3*u) + -7.80529026e-02*(eta3*a1) + -3.25433211e-02*(u*eta2*a1) + 1.79883592e-03*(u3*eta) + 6.29204610e-02*(u*eta2*a12) + 2.50379871e-01*(eta3*a12) + -8.32034716e-03*(u3*eta*a1) + 7.56521509e-02*(u*eta3*a1) + 7.28975324e-02*(u2*eta2*eta) + 5.37425071e-02*(u2*eta2*a1) + -1.12119169e-02*(u3*eta2) + -1.02025068e-03*(u4*a1) + -2.41141979e-04*(u3*a12) + 2.78877308e-03*(u*eta*a13) + 5.38601623e-04*(u2*a13) + 9.54037416e-02*(eta2*a13) + -6.12949349e-03*(u4*eta) + -3.65723636e-04*(u3*a13) + 1.04357606e-02*(u4*eta*a1) + -1.57067291e-01*(u*eta3*a12) + 5.73582747e-02*(u3*eta2*a1) + -2.07603889e-01*(eta3*a13) + 7.00338195e-02*(u2*eta2*a12) + 1.26623192e-03*(u4*a12) + -2.93338852e-02*(u*eta2*a13) + -8.34637581e-02*(u2*eta3*a1) + 2.41953141e-02*(u3*eta2*eta) + 3.49743145e-02*(u4*eta2) + 8.81404580e-03*(u2*eta*a13) + 4.85598379e-03*(u3*eta*a12) + -4.80782379e-03*(u4*eta*a12) + -1.05920066e-01*(u2*eta2*a13) + -5.10331930e-02*(u3*eta2*a12) + -8.75710297e-04*(u4*a13) + -6.64465638e-02*(u4*eta2*eta) + 6.66407126e-03*(u3*eta*a13) + -1.86904214e-01*(u2*eta3*a12) + -1.36060851e-01*(u3*eta3*a1) + -4.95540194e-02*(u4*eta2*a1) + 8.17319559e-02*(u*eta3*a13) + 2.44381585e-01*(u2*eta3*a13) + -2.06971940e-03*(u4*eta2*a12) + 9.77637782e-02*(u4*eta3*a1) + -2.42054471e-02*(u3*eta2*a13) + 1.64859194e-01*(u3*eta3*a12) + 3.71060241e-03*(u4*eta*a13);

	// Return answer
	return zeta1;

} // END of ZETA1 (2,2) fit implementation



// ZETA2 fit implementation 
// Header formatting for ZETA2 of (l,m)=(2,2) multipole
double IMRPhenomXCP_ZETA2_l2m2( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	UNUSED double eta3 = eta2*eta;
	double zeta2;

	// Evaluate fit for this parameter
	zeta2 = 2.38589050e+01*a1 + 1.49646293e+02*eta + -1.08077461e+01 + -1.98272513e+01*a12 + 5.19590590e+01*(u*eta) + -4.18121738e+01*(u*a1) + -6.12432443e+02*eta2 + 8.93335744e+01*(u*u) + 2.19662366e+02*(u*eta*a1) + -2.18132463e+03*(eta2*a1) + -3.33462448e+02*(u2*a1) + 1.32338893e+02*(u*a12) + -6.95947292e+02*(u*eta2) + -1.43194939e+03*(u2*eta) + 7.87290524e+02*(eta2*eta) + -8.44923064e+02*(eta*a12) + 5.38520565e+02*(u2*a12) + 7.79763026e+03*(u2*eta2) + -8.84567351e+02*(u*eta*a12) + 1.05694265e+04*(eta2*a12) + 4.25690943e+03*(u2*eta*a1) + 9.42785196e+02*(eta*a13) + -1.04016065e+02*(u*a13) + 2.11417669e+03*(u*eta2*eta) + 8.78513832e+01*(u3*a1) + -7.30195622e+01*(u3*u) + 7.08693433e+03*(eta3*a1) + 2.42379239e+03*(u*eta2*a1) + -4.99631088e+01*(u3*eta) + -4.45638531e+03*(u*eta2*a12) + -2.73666412e+04*(eta3*a12) + -5.37554390e+03*(u2*eta*a12) + -9.83736275e+02*(u3*eta*a1) + -1.17138325e+04*(u*eta3*a1) + -1.40709491e+04*(u2*eta2*eta) + -1.95616824e+04*(u2*eta2*a1) + 9.03100944e+02*(u3*eta2) + 2.03412767e+02*(u4*a1) + -3.26307623e+02*(u3*a12) + 6.82288963e+02*(u*eta*a13) + -3.15457440e+02*(u2*a13) + -9.79224584e+03*(eta2*a13) + 1.09441469e+03*(u4*eta) + 2.78500025e+02*(u3*a13) + -1.90234542e+03*(u4*eta*a1) + 2.60007757e+04*(u*eta3*a12) + 2.40446935e+04*(eta3*a13) + 1.79426884e+04*(u2*eta2*a12) + -2.68907843e+02*(u4*a12) + 3.13884198e+03*(u*eta2*a13) + 3.18823291e+04*(u2*eta3*a1) + -3.10751621e+03*(u3*eta2*eta) + -5.88480206e+03*(u4*eta2) + 2.58885777e+03*(u2*eta*a13) + 4.19394523e+03*(u3*eta*a12) + 1.14601064e+03*(u4*eta*a12) + -5.19325615e+03*(u2*eta2*a13) + -9.86224179e+03*(u3*eta2*a12) + 1.71473362e+02*(u4*a13) + 1.08077538e+04*(u4*eta2*eta) + -3.70449801e+03*(u3*eta*a13) + -2.10167333e+04*(u2*eta3*a12) + 1.21434706e+04*(u3*eta3*a1) + 7.75689930e+03*(u4*eta2*a1) + -1.83847722e+04*(u*eta3*a13) + -1.41709446e+04*(u4*eta3*a1) + 1.10008432e+04*(u3*eta2*a13) + -1.09074944e+04*(u3*eta3*a12) + -7.29320563e+02*(u4*eta*a13);

	// Return answer
	return zeta2;

} // END of ZETA2 (2,2) fit implementation



// NU0 fit implementation 
// Header formatting for NU0 of (l,m)=(2,2) multipole
double IMRPhenomXCP_NU0_l2m2( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	UNUSED double eta3 = eta2*eta;
	double nu0;

	// Evaluate fit for this parameter
	nu0 = 1.61301701e+01*a1 + -1.23389235e+01*u + -1.11294833e+00 + -6.22781979e+01*a12 + 3.70456454e+02*(u*eta) + 3.53831184e+01*(u*a1) + 2.21756515e+01*eta2 + -7.33746990e+01*(eta*a1) + -1.87675925e+03*(u*eta*a1) + 7.41953052e+01*(a12*a1) + -1.34539878e+01*(u2*a1) + 1.47162296e+02*(u*a12) + -2.78379431e+03*(u*eta2) + 1.61125032e+02*(u2*eta) + 4.67725211e+02*(eta*a12) + 1.54381118e+01*(u2*u) + 1.80517266e+02*(u2*a12) + -5.19821294e+02*(u2*eta2) + 4.43808889e+02*(u*eta*a12) + -1.04984634e+03*(eta2*a12) + -1.70607802e+03*(u2*eta*a1) + -7.72006796e+02*(eta*a13) + -3.03235796e+02*(u*a13) + 6.20416385e+03*(u*eta2*eta) + -3.59682021e+01*(u3*a1) + -7.99529174e+00*(u3*u) + 1.61065453e+04*(u*eta2*a1) + -4.99584410e+02*(u3*eta) + -1.55078216e+04*(u*eta2*a12) + 1.00755553e+03*(eta3*a12) + 2.92022364e+03*(u2*eta*a12) + 2.53159828e+03*(u3*eta*a1) + -3.83182653e+04*(u*eta3*a1) + -9.01641008e+02*(u2*eta2*eta) + 9.66198249e+03*(u2*eta2*a1) + 3.90747991e+03*(u3*eta2) + 6.82876439e+01*(u4*a1) + -2.62253985e+02*(u3*a12) + 3.11056625e+03*(u*eta*a13) + -3.41815489e+02*(u2*a13) + 2.81845427e+03*(eta2*a13) + 4.93507296e+02*(u3*a13) + 7.21326374e+02*(u4*eta*a1) + 4.94192687e+04*(u*eta3*a12) + -2.32518188e+04*(u3*eta2*a1) + -3.86649971e+03*(eta3*a13) + -2.43670600e+04*(u2*eta2*a12) + -3.21476677e+02*(u4*a12) + -8.27850070e+03*(u*eta2*a13) + -7.56985290e+03*(u2*eta3*a1) + -8.88093314e+03*(u3*eta2*eta) + -7.04396837e+02*(u4*eta2) + 4.63379687e+02*(u2*eta*a13) + -7.00543383e+02*(u4*eta*a12) + 9.30316297e+03*(u2*eta2*a13) + 2.15759460e+04*(u3*eta2*a12) + 4.63860352e+02*(u4*a13) + 3.88966384e+03*(u4*eta2*eta) + -5.29048698e+03*(u3*eta*a13) + 3.30773828e+04*(u2*eta3*a12) + 5.69448537e+04*(u3*eta3*a1) + -2.45824102e+03*(u4*eta2*a1) + -1.89687354e+04*(u2*eta3*a13) + 1.00262657e+04*(u4*eta2*a12) + -1.00783858e+04*(u4*eta3*a1) + 1.41835107e+04*(u3*eta2*a13) + -7.46773179e+04*(u3*eta3*a12) + -2.17628247e+03*(u4*eta*a13);

	// Return answer
	return nu0;

} // END of NU0 (2,2) fit implementation



// MU1 fit implementation 
// Header formatting for MU1 of (l,m)=(3,3) multipole
double IMRPhenomXCP_MU1_l3m3( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu1;

	// Evaluate fit for this parameter
	mu1 = -4.53857934e+00*a1 + 2.02622643e-01*u + -1.40834978e+00*eta + 1.21270021e-02 + 1.60740218e+01*a12 + -4.31596045e+00*(u*a1) + 5.10780490e+00*eta2 + 6.68631713e+00*(u*u) + 7.51884266e+01*(eta*a1) + 4.32392269e+01*(u*eta*a1) + -2.74700781e+02*(eta2*a1) + -1.39215540e+01*(a12*a1) + -2.18524983e+01*(u2*a1) + 1.35955932e+01*(u*a12) + 5.35168686e+00*(u*eta2) + -9.36113280e+01*(u2*eta) + -2.50248614e+02*(eta*a12) + -3.09988017e+00*(u2*u) + 7.34206489e+00*(u2*a12) + 2.99978196e+02*(u2*eta2) + -1.61206977e+02*(u*eta*a12) + 8.99163377e+02*(eta2*a12) + 2.72461923e+02*(u2*eta*a1) + 2.12593524e+02*(eta*a13) + -1.07885117e+01*(u*a13) + 2.84165415e+01*(u3*a1) + -9.94474532e+00*(u3*u) + -1.88779348e+02*(u*eta2*a1) + 4.13013673e+01*(u3*eta) + 6.49500354e+02*(u*eta2*a12) + -3.86292895e+02*(u3*eta*a1) + -6.99366534e+02*(u2*eta2*a1) + -1.46259030e+02*(u3*eta2) + 4.38391424e+01*(u4*a1) + -6.76149423e+01*(u3*a12) + 1.37136521e+02*(u*eta*a13) + 1.58093559e+01*(u2*a13) + -7.56122302e+02*(eta2*a13) + 1.38387621e+02*(u4*eta) + 4.72154247e+01*(u3*a13) + -5.76205002e+02*(u4*eta*a1) + 1.35939064e+03*(u3*eta2*a1) + -5.51818736e+02*(u2*eta2*a12) + -5.66576298e+01*(u4*a12) + -5.44106503e+02*(u*eta2*a13) + -4.33512288e+02*(u4*eta2) + -3.06070998e+02*(u2*eta*a13) + 9.38763144e+02*(u3*eta*a12) + 6.87408257e+02*(u4*eta*a12) + 1.42817412e+03*(u2*eta2*a13) + -3.32593998e+03*(u3*eta2*a12) + 1.79995581e+01*(u4*a13) + -6.68300907e+02*(u3*eta*a13) + 1.59468021e+03*(u4*eta2*a1) + -1.50814107e+03*(u4*eta2*a12) + 2.39316277e+03*(u3*eta2*a13) + -1.66756891e+02*(u4*eta*a13);

	// Return answer
	return mu1;

} // END of MU1 (3,3) fit implementation



// MU2 fit implementation 
// Header formatting for MU2 of (l,m)=(3,3) multipole
double IMRPhenomXCP_MU2_l3m3( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu2;

	// Evaluate fit for this parameter
	mu2 = 7.31346991e+00*a1 + -2.37789280e+00*u + 1.20431459e+01*eta + -6.94613712e-01 + -1.63335562e+01*a12 + 2.29540487e+01*(u*eta) + 2.22408234e+01*(u*a1) + -4.59010533e+01*eta2 + -1.32986121e+00*(u*u) + -1.21830549e+02*(eta*a1) + -2.23862577e+02*(u*eta*a1) + 4.45137692e+02*(eta2*a1) + 1.09570030e+01*(a12*a1) + -8.12894944e+00*(u2*a1) + -5.27817357e+01*(u*a12) + -3.85472111e+01*(u*eta2) + 2.98663571e+00*(u2*eta) + 2.77769469e+02*(eta*a12) + 2.72408145e+00*(u2*u) + 2.24690014e+01*(u2*a12) + 5.28931607e+02*(u*eta*a12) + -1.01674849e+03*(eta2*a12) + 2.20347674e+02*(u2*eta*a1) + -1.89288446e+02*(eta*a13) + 3.57962325e+01*(u*a13) + -2.81230892e+01*(u3*a1) + 2.82289045e+00*(u3*u) + 4.42583458e+02*(u*eta2*a1) + -2.01273818e+01*(u3*eta) + -1.05792918e+03*(u*eta2*a12) + -5.19799307e+02*(u2*eta*a12) + 2.52053535e+02*(u3*eta*a1) + -7.22122006e+02*(u2*eta2*a1) + 7.13335364e+01*(u3*a12) + -3.50135832e+02*(u*eta*a13) + -1.09934774e+01*(u2*a13) + 6.95980132e+02*(eta2*a13) + -2.52618341e+01*(u4*eta) + -5.09299431e+01*(u3*a13) + -9.69395602e+01*(u4*eta*a1) + -3.29798628e+02*(u3*eta2*a1) + 1.56064950e+03*(u2*eta2*a12) + -6.29989422e+00*(u4*a12) + 6.73289670e+02*(u*eta2*a13) + 7.15588580e+01*(u4*eta2) + 2.70940451e+02*(u2*eta*a13) + -6.77156015e+02*(u3*eta*a12) + 2.62979107e+02*(u4*eta*a12) + -6.94020410e+02*(u2*eta2*a13) + 1.10634560e+03*(u3*eta2*a12) + 4.90829716e+02*(u3*eta*a13) + 2.97641188e+02*(u4*eta2*a1) + -6.16841375e+02*(u4*eta2*a12) + -8.47195138e+02*(u3*eta2*a13) + -9.01572329e+01*(u4*eta*a13);

	// Return answer
	return mu2;

} // END of MU2 (3,3) fit implementation



// MU3 fit implementation 
// Header formatting for MU3 of (l,m)=(3,3) multipole
double IMRPhenomXCP_MU3_l3m3( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu3;

	// Evaluate fit for this parameter
	mu3 = 5.03435540e-01*a1 + 2.09371420e-02*u + 8.11064757e-01*eta + -5.66076153e-02 + -1.07466159e+00*a12 + -3.82277696e-01*(u*eta) + -7.65461333e-02*(u*a1) + -2.47074300e+00*eta2 + 2.52495292e-02*(u*u) + -7.17629592e+00*(eta*a1) + 1.64746902e+00*(u*eta*a1) + 2.14746185e+01*(eta2*a1) + 6.51955944e-01*(a12*a1) + -8.58755962e-01*(u2*a1) + 6.67928842e-02*(u*a12) + 1.93331873e+00*(u*eta2) + -4.02712239e-01*(u2*eta) + 1.54110191e+01*(eta*a12) + -2.92908806e-02*(u2*u) + 2.11265622e+00*(u2*a12) + 3.06317236e-01*(u2*eta2) + -2.09409193e+00*(u*eta*a12) + -4.58363438e+01*(eta2*a12) + 1.21444183e+01*(u2*eta*a1) + -9.41507028e+00*(eta*a13) + 8.87687821e-02*(u3*a1) + 9.75772078e-02*(u3*u) + -1.01519233e+01*(u*eta2*a1) + 6.11710914e-01*(u3*eta) + 1.66663401e+01*(u*eta2*a12) + -2.92705647e+01*(u2*eta*a12) + -2.71093191e+00*(u3*eta*a1) + -2.96683563e+01*(u2*eta2*a1) + -3.26578469e+00*(u3*eta2) + 3.78914976e-02*(u4*a1) + -7.35987033e-02*(u3*a12) + 6.97490870e-01*(u*eta*a13) + -1.35029546e+00*(u2*a13) + 2.77906167e+01*(eta2*a13) + -1.41813577e+00*(u4*eta) + 1.77262744e+01*(u3*eta2*a1) + 7.13101875e+01*(u2*eta2*a12) + -5.24684150e-01*(u4*a12) + -8.48528518e+00*(u*eta2*a13) + 5.33323533e+00*(u4*eta2) + 1.83144897e+01*(u2*eta*a13) + 4.01084556e+00*(u3*eta*a12) + 5.51909615e+00*(u4*eta*a12) + -4.24223784e+01*(u2*eta2*a13) + -3.14530962e+01*(u3*eta2*a12) + 4.60897509e-01*(u4*a13) + -1.84423596e+00*(u3*eta*a13) + -8.07573341e+00*(u4*eta2*a1) + 2.81635075e+00*(u4*eta2*a12) + 1.76241306e+01*(u3*eta2*a13) + -4.77244847e+00*(u4*eta*a13);

	// Return answer
	return mu3;

} // END of MU3 (3,3) fit implementation



// MU4 fit implementation 
// Header formatting for MU4 of (l,m)=(3,3) multipole
double IMRPhenomXCP_MU4_l3m3( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu4;

	// Evaluate fit for this parameter
	mu4 = 6.28363302e+01*a1 + -2.65751613e+00*u + 8.57117549e+01*eta + -5.85669025e+00 + -1.80686936e+02*a12 + -1.15837370e+02*(u*eta) + 7.79531542e+01*(u*a1) + -3.22243748e+02*eta2 + -3.09482897e+01*(u*u) + -8.65368123e+02*(eta*a1) + 2.96180380e+03*(eta2*a1) + 1.40811710e+02*(a12*a1) + -2.68013727e+01*(u2*a1) + -2.47104427e+02*(u*a12) + 8.87207761e+02*(u*eta2) + 4.76542519e+02*(u2*eta) + 2.42299624e+03*(eta*a12) + 5.16334940e+01*(u2*u) + 4.51994163e+02*(u2*a12) + -1.65543557e+03*(u2*eta2) + 1.14422477e+03*(u*eta*a12) + -7.97471130e+03*(eta2*a12) + -1.87643854e+03*(eta*a13) + 1.96299384e+02*(u*a13) + -4.39122214e+02*(u3*a1) + 7.29411925e+01*(u3*u) + -3.40581214e+03*(u*eta2*a1) + -4.95709085e+02*(u3*eta) + 3.04066375e+03*(u*eta2*a12) + -5.35856380e+03*(u2*eta*a12) + 4.43826807e+03*(u3*eta*a1) + 1.10730029e+03*(u2*eta2*a1) + 8.05330808e+02*(u3*eta2) + -2.22859025e+02*(u4*a1) + 9.77091874e+02*(u3*a12) + -1.26122868e+03*(u*eta*a13) + -4.74946053e+02*(u2*a13) + 6.07255646e+03*(eta2*a13) + -1.17066478e+03*(u4*eta) + -6.35147000e+02*(u3*a13) + 4.22543621e+03*(u4*eta*a1) + -8.95051985e+03*(u3*eta2*a1) + 1.41616807e+04*(u2*eta2*a12) + 4.17948642e+03*(u4*eta2) + 5.93904590e+03*(u2*eta*a13) + -1.00486166e+04*(u3*eta*a12) + -2.40795616e+03*(u4*eta*a12) + -1.66958985e+04*(u2*eta2*a13) + 2.17649983e+04*(u3*eta2*a12) + 2.29474989e+02*(u4*a13) + 6.57580939e+03*(u3*eta*a13) + -1.67844800e+04*(u4*eta2*a1) + 1.51950861e+04*(u4*eta2*a12) + -1.48164161e+04*(u3*eta2*a13) + -1.62424414e+03*(u4*eta*a13);

	// Return answer
	return mu4;

} // END of MU4 (3,3) fit implementation



// NU4 fit implementation 
// Header formatting for NU4 of (l,m)=(3,3) multipole
double IMRPhenomXCP_NU4_l3m3( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu4;

	// Evaluate fit for this parameter
	nu4 = 2.69907618e+02*a1 + -1.79681954e+01*u + 3.79298851e+02*eta + -3.45644795e+01 + -7.59114912e+02*a12 + 2.34022156e+02*(u*a1) + -1.11904423e+03*eta2 + 4.84356225e+01*(u*u) + -2.99382288e+03*(eta*a1) + -1.53398051e+03*(u*eta*a1) + 8.51576182e+03*(eta2*a1) + 6.05120131e+02*(a12*a1) + -6.44976847e+02*(u2*a1) + -5.02912425e+02*(u*a12) + 7.57798714e+02*(u*eta2) + 8.51841757e+03*(eta*a12) + 1.87091120e+02*(u2*u) + 2.42266466e+03*(u2*a12) + -1.47961451e+03*(u2*eta2) + 3.53929782e+03*(u*eta*a12) + -2.38636794e+04*(eta2*a12) + 3.59767194e+03*(u2*eta*a1) + -6.87220073e+03*(eta*a13) + 2.72608016e+02*(u*a13) + -1.44256120e+03*(u3*a1) + 5.27810947e+01*(u3*u) + -2.06595260e+03*(u3*eta) + -1.99099359e+03*(u*eta2*a12) + -2.07258869e+04*(u2*eta*a12) + 1.62328526e+04*(u3*eta*a1) + 5.13543555e+03*(u3*eta2) + 2.92617546e+03*(u3*a12) + -1.66810043e+03*(u*eta*a13) + -2.16293768e+03*(u2*a13) + 1.92291531e+04*(eta2*a13) + -1.69969190e+03*(u4*eta) + -1.72093889e+03*(u3*a13) + 7.07293824e+03*(u4*eta*a1) + -4.19078451e+04*(u3*eta2*a1) + 3.76005240e+04*(u2*eta2*a12) + -1.14898823e+03*(u4*a12) + 7.48772153e+03*(u4*eta2) + 2.08043128e+04*(u2*eta*a13) + -3.26500005e+04*(u3*eta*a12) + -4.54119084e+04*(u2*eta2*a13) + 8.44676872e+04*(u3*eta2*a12) + 1.33698431e+03*(u4*a13) + 1.88328303e+04*(u3*eta*a13) + -3.75105563e+04*(u4*eta2*a1) + 3.49049620e+04*(u4*eta2*a12) + -4.83163811e+04*(u3*eta2*a13) + -7.68973373e+03*(u4*eta*a13);

	// Return answer
	return nu4;

} // END of NU4 (3,3) fit implementation



// NU5 fit implementation 
// Header formatting for NU5 of (l,m)=(3,3) multipole
double IMRPhenomXCP_NU5_l3m3( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu5;

	// Evaluate fit for this parameter
	nu5 = 2.83207301e-01*a1 + -1.27541713e-01*u + 5.28993977e-01*eta + -3.71858343e-02 + -7.86314194e-01*a12 + 1.18317861e+00*(u*eta) + 1.09207541e+00*(u*a1) + -2.19774891e+00*eta2 + -1.16980723e-01*(u*u) + -4.57824565e+00*(eta*a1) + -1.07243889e+01*(u*eta*a1) + 1.78628695e+01*(eta2*a1) + 5.47322809e-01*(a12*a1) + 3.57181642e-01*(u2*a1) + -2.55566334e+00*(u*a12) + -2.20459650e+00*(u*eta2) + 1.89120274e+00*(u2*eta) + 1.08760274e+01*(eta*a12) + 2.27249473e-01*(u2*u) + -2.03399176e-01*(u2*a12) + -5.59095795e+00*(u2*eta2) + 2.59112590e+01*(u*eta*a12) + -4.06716592e+01*(eta2*a12) + -5.43966613e+00*(u2*eta*a1) + -7.59065908e+00*(eta*a13) + 1.65844584e+00*(u*a13) + -1.83355896e+00*(u3*a1) + 1.84519893e-01*(u3*u) + 2.34556229e+01*(u*eta2*a1) + -2.14078345e+00*(u3*eta) + -6.16800992e+01*(u*eta2*a12) + 5.72999213e+00*(u2*eta*a12) + 1.72511331e+01*(u3*eta*a1) + 8.84116946e+00*(u2*eta2*a1) + 4.03602944e+00*(u3*eta2) + -8.14802247e-01*(u4*a1) + 3.98592232e+00*(u3*a12) + -1.74967280e+01*(u*eta*a13) + 2.79555816e+01*(eta2*a13) + -2.73295481e+00*(u4*eta) + -2.44921014e+00*(u3*a13) + 1.01414841e+01*(u4*eta*a1) + -3.53412874e+01*(u3*eta2*a1) + 1.26844449e+00*(u4*a12) + 4.37284002e+01*(u*eta2*a13) + 7.31219991e+00*(u4*eta2) + -3.15342921e+00*(u2*eta*a13) + -3.75699743e+01*(u3*eta*a12) + -1.41720279e+01*(u4*eta*a12) + 8.02335311e+01*(u3*eta2*a12) + -6.86542622e-01*(u4*a13) + 2.29823644e+01*(u3*eta*a13) + -1.62074212e+01*(u4*eta2*a1) + 8.08765252e+00*(u4*eta2*a12) + -4.93368660e+01*(u3*eta2*a13) + 7.47061014e+00*(u4*eta*a13);

	// Return answer
	return nu5;

} // END of NU5 (3,3) fit implementation



// NU6 fit implementation 
// Header formatting for NU6 of (l,m)=(3,3) multipole
double IMRPhenomXCP_NU6_l3m3( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu6;

	// Evaluate fit for this parameter
	nu6 = 8.13036882e-02*eta + 7.30922921e-03 + 3.52496197e-01*a12 + 1.21666261e+00*(u*eta) + -2.87451551e-01*(u*a1) + 1.44533980e-01*(u*u) + -1.76856431e+00*(eta*a1) + -4.23520174e+00*(u*eta*a1) + 4.69416752e+00*(eta2*a1) + -4.58633993e-01*(a12*a1) + -3.99297622e-01*(u2*a1) + 9.49719509e-01*(u*a12) + -7.59253895e+00*(u*eta2) + -3.56665525e+00*(u2*eta) + -4.76432973e-01*(u2*u) + -1.76600556e+00*(u2*a12) + 1.29612292e+01*(u2*eta2) + 3.70029464e+00*(u*eta*a12) + 1.61994771e+01*(u2*eta*a1) + 2.75766036e+00*(eta*a13) + -7.76518709e-01*(u*a13) + 3.57506138e+00*(u3*a1) + -3.36205933e-01*(u3*u) + 3.64055899e+01*(u*eta2*a1) + 5.02452904e+00*(u3*eta) + -5.43955027e+01*(u*eta2*a12) + 3.54341538e+00*(u2*eta*a12) + -3.75744026e+01*(u3*eta*a1) + -5.92182610e+01*(u2*eta2*a1) + 1.44120361e+00*(u4*a1) + -1.11239789e+01*(u3*eta2) + -7.29414608e+00*(u3*a12) + 2.67164753e+00*(u2*a13) + -7.93195503e+00*(eta2*a13) + 6.95952605e+00*(u4*eta) + 4.40496893e+00*(u3*a13) + -3.55995150e+01*(u4*eta*a1) + 8.74569944e+01*(u3*eta2*a1) + 2.42142453e+01*(u*eta2*a13) + -2.66836491e+01*(u4*eta2) + -2.47634220e+01*(u2*eta*a13) + 7.54521342e+01*(u3*eta*a12) + 2.94422675e+01*(u4*eta*a12) + 7.48979362e+01*(u2*eta2*a13) + -1.78470858e+02*(u3*eta2*a12) + -1.78134919e+00*(u4*a13) + -4.44695743e+01*(u3*eta*a13) + 1.40602112e+02*(u4*eta2*a1) + -1.42687401e+02*(u4*eta2*a12) + 1.06031824e+02*(u3*eta2*a13) + 8.12662402e+00*(u4*eta*a13);

	// Return answer
	return nu6;

} // END of NU6 (3,3) fit implementation



// ZETA1 fit implementation 
// Header formatting for ZETA1 of (l,m)=(3,3) multipole
double IMRPhenomXCP_ZETA1_l3m3( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double zeta1;

	// Evaluate fit for this parameter
	zeta1 = -3.94085343e-02*a1 + 1.05305536e-02*u + -6.43818142e-02*eta + 4.59045525e-03 + 9.64671039e-02*a12 + -1.36047713e-01*(u*eta) + -8.29615890e-02*(u*a1) + 2.13194924e-01*eta2 + 5.76300568e-01*(eta*a1) + 1.09557147e+00*(u*eta*a1) + -1.91694439e+00*(eta2*a1) + -6.87608468e-02*(a12*a1) + 5.52951609e-02*(u2*a1) + 1.82981986e-01*(u*a12) + 3.99503304e-01*(u*eta2) + -1.39567406e+00*(eta*a12) + -2.78578749e-02*(u2*u) + -1.95596687e-01*(u2*a12) + -2.43592700e+00*(u*eta*a12) + 4.60713660e+00*(eta2*a12) + -7.90694221e-01*(u2*eta*a1) + 9.91555939e-01*(eta*a13) + -1.17023671e-01*(u*a13) + 2.11330472e-01*(u3*a1) + -1.77181662e-02*(u3*u) + -3.29740379e+00*(u*eta2*a1) + 3.65007990e-01*(u3*eta) + 7.46102754e+00*(u*eta2*a12) + 2.73794874e+00*(u2*eta*a12) + -2.78196927e+00*(u3*eta*a1) + 2.55446818e+00*(u2*eta2*a1) + -1.08759749e+00*(u3*eta2) + 6.76161670e-02*(u4*a1) + -4.51456154e-01*(u3*a12) + 1.58043426e+00*(u*eta*a13) + 1.67341520e-01*(u2*a13) + -3.26586082e+00*(eta2*a13) + 2.54812355e-01*(u4*eta) + 2.91970470e-01*(u3*a13) + -9.35045956e-01*(u4*eta*a1) + 8.36707448e+00*(u3*eta2*a1) + -8.62039875e+00*(u2*eta2*a12) + -5.89860144e-02*(u4*a12) + -4.93903010e+00*(u*eta2*a13) + -7.95088131e-01*(u4*eta2) + -2.30649257e+00*(u2*eta*a13) + 5.96471243e+00*(u3*eta*a12) + 8.01033738e-01*(u4*eta*a12) + 7.13445713e+00*(u2*eta2*a13) + -1.81013351e+01*(u3*eta2*a12) + -3.87193675e+00*(u3*eta*a13) + 2.77289850e+00*(u4*eta2*a1) + -2.31866940e+00*(u4*eta2*a12) + 1.18479971e+01*(u3*eta2*a13);

	// Return answer
	return zeta1;

} // END of ZETA1 (3,3) fit implementation



// ZETA2 fit implementation 
// Header formatting for ZETA2 of (l,m)=(3,3) multipole
double IMRPhenomXCP_ZETA2_l3m3( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double zeta2;

	// Evaluate fit for this parameter
	zeta2 = 9.29107126e+02*a1 + -2.56766535e+02*u + 1.24151063e+03*eta + -1.04558291e+02 + -2.60252727e+03*a12 + 2.75880428e+03*(u*eta) + 2.05275964e+03*(u*a1) + -3.72063332e+03*eta2 + 1.58671336e+02*(u*u) + -1.17669706e+04*(eta*a1) + -2.33687962e+04*(u*eta*a1) + 3.66881027e+04*(eta2*a1) + 2.04350440e+03*(a12*a1) + -2.89425040e+03*(u2*a1) + -4.35609237e+03*(u*a12) + -6.54916424e+03*(u*eta2) + 3.29382706e+04*(eta*a12) + 9.37516833e+02*(u2*u) + 9.59540450e+03*(u2*a12) + -4.38108860e+03*(u2*eta2) + 5.07119452e+04*(u*eta*a12) + -1.01747166e+05*(eta2*a12) + 2.28958462e+04*(u2*eta*a1) + -2.59999057e+04*(eta*a13) + 2.62055635e+03*(u*a13) + -7.11341067e+03*(u3*a1) + 3.05817124e+02*(u3*u) + 6.02213219e+04*(u*eta2*a1) + -1.14191467e+04*(u3*eta) + -1.35185555e+05*(u*eta2*a12) + -9.29441823e+04*(u2*eta*a12) + 8.79214604e+04*(u3*eta*a1) + -3.99813906e+04*(u2*eta2*a1) + 3.20416872e+04*(u3*eta2) + 1.47039016e+04*(u3*a12) + -3.07179878e+04*(u*eta*a13) + -7.97932887e+03*(u2*a13) + 8.01902629e+04*(eta2*a13) + -7.63570501e+03*(u4*eta) + -9.03780382e+03*(u3*a13) + 2.44025878e+04*(u4*eta*a1) + -2.51368617e+05*(u3*eta2*a1) + 2.15313265e+05*(u2*eta2*a12) + -3.91011900e+03*(u4*a12) + 8.35473695e+04*(u*eta2*a13) + 3.07460428e+04*(u4*eta2) + 8.31029618e+04*(u2*eta*a13) + -1.82839776e+05*(u3*eta*a12) + -2.06086010e+05*(u2*eta2*a13) + 5.28032481e+05*(u3*eta2*a12) + 4.28554074e+03*(u4*a13) + 1.12890274e+05*(u3*eta*a13) + -1.21745951e+05*(u4*eta2*a1) + 1.03390621e+05*(u4*eta2*a12) + -3.29145222e+05*(u3*eta2*a13) + -2.30271021e+04*(u4*eta*a13);

	// Return answer
	return zeta2;

} // END of ZETA2 (3,3) fit implementation



// NU0 fit implementation 
// Header formatting for NU0 of (l,m)=(3,3) multipole
double IMRPhenomXCP_NU0_l3m3( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	UNUSED double u2 = u*u;
	double u3 = u2*u;
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu0;

	// Evaluate fit for this parameter
	nu0 = 5.36687963e+03*eta + 1.70443090e+02 + 8.70547376e+03*a12 + 1.40137903e+04*(u*eta) + 4.09379466e+03*(u*a1) + -2.48129707e+04*eta2 + 2.30979121e+03*(u*u) + -4.77831622e+04*(eta*a1) + -1.40379010e+05*(u*eta*a1) + 2.01077912e+05*(eta2*a1) + -1.16503462e+04*(a12*a1) + -1.05248431e+04*(u2*a1) + -2.27099418e+04*(u*a12) + -7.52760213e+04*(u*eta2) + -8.10540199e+04*(u2*eta) + -2.26078101e+04*(u2*a12) + 3.71739788e+05*(u2*eta2) + 4.57485856e+05*(u*eta*a12) + -1.41683045e+05*(eta2*a12) + 5.27303092e+05*(u2*eta*a1) + 7.46041747e+04*(eta*a13) + 2.30125173e+04*(u*a13) + -2.74528821e+04*(u3*a1) + -5.79207070e+03*(u3*u) + 6.20868233e+05*(u*eta2*a1) + -2.65756855e+04*(u3*eta) + -1.72524321e+06*(u*eta2*a12) + -5.55957329e+05*(u2*eta*a12) + 5.22666265e+05*(u3*eta*a1) + -2.45621582e+06*(u2*eta2*a1) + 1.74592173e+04*(u4*a1) + 1.37832371e+05*(u3*eta2) + 1.23438577e+05*(u3*a12) + -3.95112574e+05*(u*eta*a13) + 4.28000178e+04*(u2*a13) + -1.25652475e+05*(eta2*a13) + 1.31685427e+05*(u4*eta) + -1.11193555e+05*(u3*a13) + -6.09859826e+05*(u4*eta*a1) + -1.91105176e+06*(u3*eta2*a1) + 3.53833660e+06*(u2*eta2*a12) + 4.10194256e+04*(u4*a12) + 1.37888198e+06*(u*eta2*a13) + -5.50421328e+05*(u4*eta2) + -1.92290349e+06*(u3*eta*a12) + 2.27362936e+05*(u4*eta*a12) + -1.23801700e+06*(u2*eta2*a13) + 6.37554388e+06*(u3*eta2*a12) + -6.73840896e+04*(u4*a13) + 1.64262461e+06*(u3*eta*a13) + 2.79154732e+06*(u4*eta2*a1) + -2.61257305e+06*(u4*eta2*a12) + -5.26928552e+06*(u3*eta2*a13) + 4.08345986e+05*(u4*eta*a13);

	// Return answer
	return nu0;

} // END of NU0 (3,3) fit implementation


#ifdef __cplusplus
}
#endif
