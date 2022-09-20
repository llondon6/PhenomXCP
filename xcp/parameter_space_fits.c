
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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double mu1;

	// Evaluate fit for this parameter
	mu1 = 2.15682581e+00*a1 + -1.28496412e+00*u + 5.83270450e-01*eta + -1.32517396e-01 + -1.51789734e+00*(u*u) + 2.45165736e+01*(u*eta) + -3.11460675e+01*(eta*a1) + -4.20022288e+00*a12 + 6.62848134e+00*(u*a1) + -6.74726032e+00*(u*a12) + 1.66988267e+02*(eta2*a1) + -1.43054987e+02*(u*eta*a1) + -3.22824879e-01*(u2*u) + -1.38670435e+02*(u*eta2) + 3.41320483e+01*(u2*eta) + 6.76804386e+01*(eta*a12) + 2.41453197e+00*(u3*a1) + 1.54354983e+02*(u*eta*a12) + -2.18591605e+02*(u2*eta2) + -2.75717314e+01*(u2*eta*a1) + 3.18893622e+00*(u2*a12) + 2.34122411e+02*(u*eta2*eta) + 8.79097605e+02*(u*eta2*a1) + -3.05513902e+02*(eta3*a1) + -3.65844597e+02*(eta2*a12) + -3.08756278e+01*(u2*eta*a12) + 2.44992156e+02*(u2*eta2*a1) + -9.73231053e+02*(u*eta2*a12) + -3.31359164e+00*(u3*a12) + -6.71251357e+00*(u3*eta*a1) + 4.44791686e+02*(u2*eta2*eta) + 6.49350493e+02*(eta3*a12) + -1.60164134e+03*(u*eta3*a1) + 4.47397009e+01*(u3*eta2*eta) + -3.75231921e+01*(u3*eta2*a1) + 1.81264655e+03*(u*eta3*a12) + -6.06558133e+02*(u2*eta3*a1) + 8.56777572e+01*(u2*eta2*a12) + 1.88756234e+01*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double mu2;

	// Evaluate fit for this parameter
	mu2 = -7.54164053e+00*u + -1.91516385e-01 + -4.18061332e+00*(u*u) + 1.55039994e+02*(u*eta) + 3.20399024e+00*(eta*a1) + -8.53142828e-01*a12 + 5.43122209e+01*(u*a1) + -6.36902248e+01*(u*a12) + -1.07378584e+03*(u*eta*a1) + -7.16400097e+00*(u2*u) + -9.60997353e+02*(u*eta2) + 7.96409152e+00*(u2*a1) + 6.68912601e+01*(u2*eta) + 1.20379699e+03*(u*eta*a12) + -3.96440207e+02*(u2*eta2) + -4.91323229e+01*(u2*eta*a1) + -9.16024335e+00*(u2*a12) + 1.83641227e+03*(u*eta2*eta) + 1.29899349e+02*(u3*eta) + 6.56170457e+03*(u*eta2*a1) + 7.50131222e+01*(u2*eta*a12) + -7.19880147e+03*(u*eta2*a12) + 1.67610628e+00*(u3*a12) + -7.73634153e+02*(u3*eta2) + 8.65119085e+02*(u2*eta2*eta) + -1.25208717e+04*(u*eta3*a1) + 1.53086671e+03*(u3*eta2*eta) + -3.72877742e+01*(u3*eta2*a1) + 1.35954995e+04*(u*eta3*a12) + -9.56493504e+01*(u2*eta2*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double mu3;

	// Evaluate fit for this parameter
	mu3 = 6.44905054e-02*a1 + 1.76863919e-02*u + 6.66706654e-03*eta + -3.19807241e-03 + -3.79911024e-03*(u*u) + -8.53135817e-01*(eta*a1) + -1.16891300e-01*(u*a1) + 1.52106345e-01*(u*a12) + 4.35962056e+00*(eta2*a1) + 5.72935240e-01*(u*eta*a1) + -7.60172185e-01*(u*eta2) + -5.50827313e-02*(u2*a1) + 8.71133458e-02*(u3*a1) + -9.01179170e-01*(u*eta*a12) + 1.05548067e+00*(u2*eta*a1) + 1.87056474e+00*(u*eta2*eta) + -2.07814905e-01*(u3*eta) + -7.79327670e+00*(eta3*a1) + -1.30087246e-01*(u2*eta*a12) + -5.50044904e+00*(u2*eta2*a1) + 1.08773010e+00*(u*eta2*a12) + -1.26442295e-01*(u3*a12) + 1.36581619e+00*(u3*eta2) + -2.52172811e-01*(u3*eta*a1) + -1.48820532e-01*(u2*eta2*eta) + -1.44914438e+00*(u*eta3*a1) + -2.17495557e+00*(u3*eta2*eta) + -4.31214126e-01*(u3*eta2*a1) + 9.76511974e+00*(u2*eta3*a1) + 3.05066950e-01*(u2*eta2*a12) + 5.27227564e-01*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu4;

	// Evaluate fit for this parameter
	nu4 = -3.09964085e-02*a1 + 1.57271240e-02*u + -4.62111242e-02*eta + 4.43403755e-03 + -7.24626894e-03*(u*u) + 1.07505792e-01*eta2 + -3.23948966e-01*(u*eta) + 3.71166548e-01*(eta*a1) + 3.19607415e-02*a12 + -7.69276357e-02*(u*a1) + 4.78116769e-02*(u*a12) + -1.30763107e+00*(eta2*a1) + 1.72556494e+00*(u*eta*a1) + 1.88646735e+00*(u*eta2) + 7.52390889e-02*(u2*a1) + -3.39932512e-01*(eta*a12) + -7.69252474e-03*(u3*a1) + -1.40806912e+00*(u*eta*a12) + 5.08362676e-01*(u2*eta2) + -7.33275104e-01*(u2*eta*a1) + -8.20112905e-02*(u2*a12) + -3.33109320e+00*(u*eta2*eta) + -1.06030817e+01*(u*eta2*a1) + 1.37500320e+00*(eta3*a1) + 8.20493332e-01*(eta2*a12) + 8.15309343e-01*(u2*eta*a12) + 1.77477691e+00*(u2*eta2*a1) + 9.38718692e+00*(u*eta2*a12) + 3.38509103e-02*(u3*a12) + 2.19073512e-01*(u3*eta2) + -9.15076848e-02*(u3*eta*a1) + -1.65464279e+00*(u2*eta2*eta) + 1.96269582e+01*(u*eta3*a1) + -1.00299232e+00*(u3*eta2*eta) + 5.84743480e-01*(u3*eta2*a1) + -1.80393010e+01*(u*eta3*a12) + -1.97128750e+00*(u2*eta2*a12) + -1.51449688e-01*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu5;

	// Evaluate fit for this parameter
	nu5 = -3.05660691e-02*u + 1.31352981e-01*eta + -8.35979819e-03 + -3.97631283e-02*(u*u) + -7.74561025e-01*eta2 + 4.67442854e-01*(u*eta) + -4.56835905e-01*(eta*a1) + -6.96499953e-02*a12 + 8.92504364e-02*(u*a1) + 1.34264675e+00*(eta2*eta) + -8.46986502e-02*(u*a12) + 3.51856126e+00*(eta2*a1) + -8.95585605e-01*(u*eta*a1) + 9.75793445e-03*(u2*u) + -2.38128939e+00*(u*eta2) + 1.30075523e-01*(u2*a1) + 6.68980911e-01*(u2*eta) + 1.34689385e+00*(eta*a12) + -4.06966260e-02*(u3*a1) + 3.42694043e-01*(u*eta*a12) + -4.23329268e+00*(u2*eta2) + -1.69354010e+00*(u2*eta*a1) + -5.38280413e-02*(u2*a12) + 3.92615903e+00*(u*eta2*eta) + -1.04810973e-01*(u3*eta) + 3.45498185e+00*(u*eta2*a1) + -6.33304481e+00*(eta3*a1) + -8.96445798e+00*(eta2*a12) + 2.94977443e-01*(u2*eta*a12) + 9.66621695e+00*(u2*eta2*a1) + 3.99054227e-02*(u3*a12) + 4.40717414e-01*(u3*eta*a1) + 9.16060017e+00*(u2*eta2*eta) + 1.83156510e+01*(eta3*a12) + -5.08260069e+00*(u*eta3*a1) + 1.46376607e+00*(u3*eta2*eta) + -1.36393119e+00*(u3*eta2*a1) + -2.13514623e+01*(u2*eta3*a1) + -1.28125539e-01*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu6;

	// Evaluate fit for this parameter
	nu6 = -7.91833101e-02*a1 + 3.65688008e-02*u + -3.11667273e-01*eta + 1.58221518e-02 + 1.90135838e+00*eta2 + -5.36918575e-01*(u*eta) + 1.62841708e+00*(eta*a1) + 1.54991878e-01*a12 + -3.54229177e-01*(u*a1) + -3.78788439e+00*(eta2*eta) + 4.50558046e-01*(u*a12) + -1.00492902e+01*(eta2*a1) + 5.78641907e+00*(u*eta*a1) + 4.74027862e-02*(u2*u) + 3.14680451e+00*(u*eta2) + -3.41885953e-02*(u2*a1) + -3.99761812e-02*(u2*eta) + -2.74630824e+00*(eta*a12) + -7.23995151e+00*(u*eta*a12) + 4.17428334e-01*(u2*eta*a1) + 4.96012269e-02*(u2*a12) + -6.31613054e+00*(u*eta2*eta) + -9.13937175e-01*(u3*eta) + -3.29305662e+01*(u*eta2*a1) + 1.96576021e+01*(eta3*a1) + 1.56527334e+01*(eta2*a12) + -7.02233835e-01*(u2*eta*a12) + 3.96779936e+01*(u*eta2*a12) + -7.15296358e-02*(u3*a12) + 4.66551632e+00*(u3*eta2) + 5.65996951e-01*(u3*eta*a1) + 2.82009876e-01*(u2*eta2*eta) + -2.88629987e+01*(eta3*a12) + 6.15962682e+01*(u*eta3*a1) + -7.16284104e+00*(u3*eta2*eta) + -2.13232667e+00*(u3*eta2*a1) + -7.13316197e+01*(u*eta3*a12) + -3.67931430e+00*(u2*eta3*a1) + 1.86112931e+00*(u2*eta2*a12) + 2.41635196e-01*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double zeta1;

	// Evaluate fit for this parameter
	zeta1 = 9.00838223e-05*a1 + -2.61569600e-05*u + 7.16021495e-05*eta + -6.53153489e-06 + 2.91218080e-05*(u*u) + -1.60054005e-04*eta2 + 6.07770789e-04*(u*eta) + -1.11726499e-03*(eta*a1) + -6.89971710e-05*a12 + 2.78523863e-04*(u*a1) + -3.28502043e-04*(u*a12) + 4.73841883e-03*(eta2*a1) + -5.75897175e-03*(u*eta*a1) + -4.10195362e-05*(u2*u) + -4.31141469e-03*(u*eta2) + -2.38007338e-04*(u2*a1) + -1.70136568e-04*(u2*eta) + 6.74949940e-04*(eta*a12) + 2.99354335e-05*(u3*a1) + 6.90036189e-03*(u*eta*a12) + 2.67585238e-03*(u2*eta*a1) + 1.67836942e-04*(u2*a12) + 9.22181273e-03*(u*eta2*eta) + 7.24517901e-04*(u3*eta) + 3.68125030e-02*(u*eta2*a1) + -7.09563751e-03*(eta3*a1) + -1.54635069e-03*(eta2*a12) + -1.30201508e-03*(u2*eta*a12) + -1.13454738e-02*(u2*eta2*a1) + -4.28832285e-02*(u*eta2*a12) + -3.61709857e-03*(u3*eta2) + -6.20621847e-04*(u3*eta*a1) + 1.13142583e-03*(u2*eta2*eta) + -7.28507021e-02*(u*eta3*a1) + 5.79768042e-03*(u3*eta2*eta) + 1.69019440e-03*(u3*eta2*a1) + 8.19251950e-02*(u*eta3*a12) + 1.73192253e-02*(u2*eta3*a1) + 2.60799878e-03*(u2*eta2*a12) + 8.49508560e-05*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double zeta2;

	// Evaluate fit for this parameter
	zeta2 = -1.78400964e+01*a1 + 3.50956181e+00*u + 9.60136924e+00*eta + 6.82931365e-01 + -4.03092002e+00*(u*u) + -1.25379881e+02*eta2 + -1.00242982e+02*(u*eta) + 1.80805286e+02*(eta*a1) + 1.82857576e+01*a12 + -3.30100396e+01*(u*a1) + 3.30699137e+02*(eta2*eta) + 3.68856564e+01*(u*a12) + -5.51384118e+02*(eta2*a1) + 7.98366884e+02*(u*eta*a1) + 6.66779765e+00*(u2*u) + 7.04004154e+02*(u*eta2) + 4.11392741e+01*(u2*a1) + -1.74965293e+02*(eta*a12) + -1.77442646e+01*(u3*a1) + -8.86161650e+02*(u*eta*a12) + 2.95946717e+02*(u2*eta2) + -4.28113243e+02*(u2*eta*a1) + -3.48216369e+01*(u2*a12) + -1.41609599e+03*(u*eta2*eta) + -7.39848186e+01*(u3*eta) + -5.20127369e+03*(u*eta2*a1) + 4.10553410e+02*(eta3*a1) + 4.13945565e+02*(eta2*a12) + 3.37924311e+02*(u2*eta*a12) + 1.26339786e+03*(u2*eta2*a1) + 5.67524228e+03*(u*eta2*a12) + 1.53612107e+01*(u3*a12) + 2.80681291e+02*(u3*eta2) + 1.22333354e+02*(u3*eta*a1) + -9.12632916e+02*(u2*eta2*eta) + 1.00216210e+04*(u*eta3*a1) + -4.38379022e+02*(u3*eta2*eta) + -1.22639277e+02*(u3*eta2*a1) + -1.07353624e+04*(u*eta3*a12) + -8.89972937e+02*(u2*eta3*a1) + -7.84231099e+02*(u2*eta2*a12) + -8.27258185e+01*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu0;

	// Evaluate fit for this parameter
	nu0 = 9.73061928e+02*a1 + 1.20543450e+02*u + 1.27484961e+03*eta + -1.39537396e+02 + 3.71950292e+02*(u*u) + -2.98191097e+03*eta2 + -8.58556058e+03*(eta*a1) + -3.84715641e+02*a12 + -2.64697029e+03*(u*a1) + 3.79489418e+03*(u*a12) + 1.90770015e+04*(eta2*a1) + 3.51129846e+04*(u*eta*a1) + 3.67130832e+02*(u2*u) + -2.79214444e+03*(u2*a1) + -3.25454984e+03*(u2*eta) + 1.11693629e+03*(u3*a1) + -5.45669002e+04*(u*eta*a12) + 2.96460120e+04*(u2*eta*a1) + 2.34027481e+03*(u2*a12) + -1.27304698e+04*(u*eta2*eta) + -9.07293457e+03*(u3*eta) + -1.79249757e+05*(u*eta2*a1) + 2.55795990e+04*(eta2*a12) + -2.40717405e+04*(u2*eta*a12) + -7.31512307e+04*(u2*eta2*a1) + 2.81123919e+05*(u*eta2*a12) + -1.08578690e+03*(u3*a12) + 5.02235441e+04*(u3*eta2) + -3.65979106e+03*(u3*eta*a1) + 2.63207469e+04*(u2*eta2*eta) + -7.73926053e+04*(eta3*a12) + 3.38873151e+05*(u*eta3*a1) + -7.19780980e+04*(u3*eta2*eta) + -8.16469730e+03*(u3*eta2*a1) + -5.02931114e+05*(u*eta3*a12) + 5.81371672e+04*(u2*eta2*a12) + 5.08606096e+03*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu1;

	// Evaluate fit for this parameter
	mu1 = 1.05577046e+00*a1 + -1.11593273e+00*u + 3.32603177e+00*eta + -3.19894268e-01 + 1.21806484e+00*(u*u) + -1.01150188e+01*eta2 + 2.43560932e+01*(u*eta) + -6.04555076e+00*(eta*a1) + 8.53210218e+00*(u*a1) + -2.03289935e+01*(u*a12) + -1.83407578e+02*(u*eta*a1) + -1.00859623e+00*(u2*u) + -8.66461824e+01*(u*eta2) + -7.80105454e+00*(u2*a1) + -1.71625487e+01*(u2*eta) + -1.71570781e+01*(eta*a12) + -1.48996563e+00*(a12*a1) + 8.00264943e+00*(u3*a1) + 4.17662826e+02*(u*eta*a12) + 5.91560558e+01*(u2*eta2) + 1.02257312e+02*(u2*eta*a1) + 1.44658864e+01*(u*a13) + 1.06087728e+01*(u2*a12) + 2.57446607e+00*(u3*eta) + 6.47054486e+02*(u*eta2*a1) + 3.20787039e+01*(eta*a13) + 1.09302100e+02*(eta2*a12) + -1.42990819e+02*(eta2*a13) + -1.22413830e+02*(u2*eta*a12) + -3.31994345e+02*(u2*eta2*a1) + -1.44653565e+03*(u*eta2*a12) + -2.84931604e+02*(u*eta*a13) + -1.38870216e+01*(u3*a12) + -2.72120215e+00*(u2*a13) + -2.60397958e+01*(u3*eta*a1) + 3.00615273e+01*(u3*eta2*a1) + 7.27169366e+00*(u3*a13) + 9.65521209e+02*(u*eta2*a13) + 1.60334161e+01*(u2*eta*a13) + 3.47648250e+02*(u2*eta2*a12) + 2.10940353e+01*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu2;

	// Evaluate fit for this parameter
	mu2 = 3.24429744e+00*a1 + 6.04813725e+00*eta + -3.91385818e-01 + -1.95790454e-01*(u*u) + -2.45282037e+01*eta2 + 3.83343922e+00*(u*eta) + -5.22065674e+01*(eta*a1) + -6.71558295e+00*a12 + 1.16235743e+00*(u*a1) + -1.76312177e+00*(u*a12) + 2.11893503e+02*(eta2*a1) + -4.19143811e+01*(u*eta*a1) + -9.20767234e-01*(u2*u) + -1.47544363e+01*(u*eta2) + 1.21186898e+00*(u2*a1) + 1.18086297e+02*(eta*a12) + 4.90606124e+00*(a12*a1) + 4.12943924e+00*(u3*a1) + 6.99588345e+01*(u*eta*a12) + -3.77755764e+00*(u2*a12) + 8.16132810e+00*(u3*eta) + 1.59224604e+02*(u*eta2*a1) + -8.98621006e+01*(eta*a13) + -4.98641940e+02*(eta2*a12) + 3.84316885e+02*(eta2*a13) + 1.37358401e+01*(u2*eta*a12) + -1.02088824e+01*(u2*eta2*a1) + -2.64142917e+02*(u*eta2*a12) + -2.03230396e+01*(u*eta*a13) + -6.63014940e+00*(u3*a12) + 2.92635353e+00*(u2*a13) + -3.05288531e+01*(u3*eta2) + -1.81646194e+01*(u3*eta*a1) + 5.77444255e+01*(u3*eta2*a1) + 3.73155877e+00*(u3*a13) + 7.15260594e+01*(u*eta2*a13) + -1.27628294e+01*(u2*eta*a13) + 5.22555176e+00*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu3;

	// Evaluate fit for this parameter
	mu3 = -1.68830509e-01*a1 + 7.96816079e-02*u + -3.37069996e-01*eta + 1.75233624e-02 + 2.01469739e-02*(u*u) + 1.41955608e+00*eta2 + -4.35453714e-01*(u*eta) + 3.00155196e+00*(eta*a1) + 4.79604415e-01*a12 + -4.30922283e-01*(u*a1) + 7.04699983e-01*(u*a12) + -1.27423755e+01*(eta2*a1) + 8.77446696e-01*(u*eta*a1) + -1.68453238e-01*(u2*u) + 9.04584094e-01*(u*eta2) + -1.57961048e-01*(u2*a1) + -3.31499286e-01*(u2*eta) + -8.01568734e+00*(eta*a12) + -4.35756215e-01*(a12*a1) + 9.37941824e-01*(u3*a1) + 8.93844216e-01*(u*eta*a12) + 6.61950014e-01*(u2*eta2) + 2.53046604e+00*(u2*eta*a1) + -4.56915773e-01*(u*a13) + 2.75737136e-01*(u2*a12) + 1.25202727e+00*(u3*eta) + 6.96133723e+00*(eta*a13) + 3.34179051e+01*(eta2*a12) + -2.81523151e+01*(eta2*a13) + -4.73626803e+00*(u2*eta*a12) + -4.67397189e+00*(u2*eta2*a1) + -4.92884506e+00*(u*eta2*a12) + -1.50656484e+00*(u3*a12) + -8.95161766e-02*(u2*a13) + -3.11071086e+00*(u3*eta2) + -4.44685065e+00*(u3*eta*a1) + 8.36227639e+00*(u3*eta2*a1) + 8.24068572e-01*(u3*a13) + -2.59195093e+00*(u*eta2*a13) + 1.92247325e+00*(u2*eta*a13) + 6.87539117e+00*(u2*eta2*a12) + 2.43103536e+00*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu4;

	// Evaluate fit for this parameter
	mu4 = -3.48916485e+01*a1 + 3.08122264e+00*u + -5.54001527e+01*eta + 3.04423956e+00 + -2.62600670e+00*(u*u) + 1.99228404e+02*eta2 + 5.51143151e+02*(eta*a1) + 8.58175073e+01*a12 + -3.57745024e+01*(u*a1) + 1.10810753e+02*(u*a12) + -2.00649903e+03*(eta2*a1) + 1.61824739e+02*(u*eta*a1) + 3.11902963e+01*(u*eta2) + 1.41803478e+01*(u2*a1) + 1.64573591e+01*(u2*eta) + -1.31089768e+03*(eta*a12) + -6.28989370e+01*(a12*a1) + 2.39587998e+01*(u3*a1) + -7.51669063e+02*(u*eta*a12) + -1.10517370e+02*(u2*eta2) + -8.30412416e+01*(u*a13) + -3.73627527e+01*(u2*a12) + -7.31434503e+01*(u3*eta) + -6.71139309e+02*(u*eta2*a1) + 9.30506228e+02*(eta*a13) + 4.73092348e+03*(eta2*a12) + -3.31029099e+03*(eta2*a13) + 2.17549524e+02*(u2*eta2*a1) + 2.57612636e+03*(u*eta2*a12) + 6.10058822e+02*(u*eta*a13) + -1.03459602e+02*(u3*a12) + 2.87820852e+01*(u2*a13) + 1.64636726e+02*(u3*eta2) + 2.99343308e+02*(u3*eta*a1) + -6.06315847e+02*(u3*eta2*a1) + 8.09197376e+01*(u3*a13) + -1.97560979e+03*(u*eta2*a13) + -2.04069455e+02*(u2*eta2*a12) + -1.14993489e+02*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu4;

	// Evaluate fit for this parameter
	nu4 = 1.23003585e+02*a1 + 3.07080122e+01*u + 2.05541220e+02*eta + -1.70416997e+01 + -9.97137773e+00*(u*u) + -6.16227852e+02*eta2 + -5.32378163e+02*(u*eta) + -1.53190984e+03*(eta*a1) + -3.33553614e+02*a12 + -1.79243525e+02*(u*a1) + 3.63855521e+02*(u*a12) + 4.37247978e+03*(eta2*a1) + 3.18743840e+03*(u*eta*a1) + 1.83311058e+01*(u2*u) + 1.82740448e+03*(u*eta2) + -2.14612702e+01*(u2*a1) + 2.17090432e+02*(u2*eta) + 3.99410204e+03*(eta*a12) + 2.66980118e+02*(a12*a1) + -1.20289751e+02*(u3*a1) + -6.19355947e+03*(u*eta*a12) + -1.29506516e+03*(u2*eta2) + -5.53294253e+02*(u2*eta*a1) + -2.39725989e+02*(u*a13) + 2.27506716e+02*(u2*a12) + -6.19975004e+01*(u3*eta) + -1.09190563e+04*(u*eta2*a1) + -3.13327694e+03*(eta*a13) + -1.10020071e+04*(eta2*a12) + 8.53008211e+03*(eta2*a13) + -8.02761969e+02*(u2*eta*a12) + 5.89918538e+03*(u2*eta2*a1) + 2.05277875e+04*(u*eta2*a12) + 3.84165497e+03*(u*eta*a13) + 2.02934481e+02*(u3*a12) + -2.80599196e+02*(u2*a13) + 3.90299030e+02*(u3*eta*a1) + -7.60010707e+01*(u3*eta2*a1) + -8.82883249e+01*(u3*a13) + -1.20499664e+04*(u*eta2*a13) + 2.03098124e+03*(u2*eta*a13) + -6.64162924e+03*(u2*eta2*a12) + -4.64835587e+02*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu5;

	// Evaluate fit for this parameter
	nu5 = -3.18083953e-01*a1 + -4.40013636e-01*eta + 2.82156221e-02 + -1.18955622e-02*(u*u) + 1.22838930e+00*eta2 + 1.88105095e-01*(u*eta) + 4.04840427e+00*(eta*a1) + 6.23463104e-01*a12 + -1.04171159e-01*(u*a1) + 3.65785694e-01*(u*a12) + -1.21436698e+01*(eta2*a1) + 3.06958212e-02*(u2*u) + -5.37338609e-01*(u*eta2) + 9.10638252e-02*(u2*a1) + -8.90377684e+00*(eta*a12) + -4.58883548e-01*(a12*a1) + -1.33974173e+00*(u*eta*a12) + -1.56845419e-01*(u2*eta*a1) + -3.02400490e-01*(u*a13) + -6.40565475e-02*(u2*a12) + -6.49637206e-01*(u3*eta) + 6.21454365e+00*(eta*a13) + 2.74780796e+01*(eta2*a12) + -1.91190619e+01*(eta2*a13) + 2.42153297e+00*(u*eta2*a12) + 8.20101732e-01*(u*eta*a13) + -3.20095609e-01*(u3*a12) + 1.72902660e+00*(u3*eta2) + 1.89560704e+00*(u3*eta*a1) + -3.62528767e+00*(u3*eta2*a1) + 3.31482514e-01*(u3*a13) + 1.99149235e-01*(u2*eta*a13) + -8.41843977e-01*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu6;

	// Evaluate fit for this parameter
	nu6 = -8.54066852e-02*a1 + -3.68102593e-02*u + -1.34539693e-01*eta + 1.48437188e-02 + 1.42066939e-02*(u*u) + 2.54621426e-01*eta2 + 6.88086421e-01*(u*eta) + 6.82061882e-01*(eta*a1) + 2.68882652e-01*a12 + 2.74424508e-01*(u*a1) + -5.83809242e-01*(u*a12) + -4.80962800e+00*(u*eta*a1) + -3.35409813e-02*(u2*u) + -2.62685143e+00*(u*eta2) + 9.40412756e-02*(u2*a1) + -4.54426512e-01*(u2*eta) + -1.94203426e+00*(eta*a12) + -2.21954980e-01*(a12*a1) + 1.44147362e-01*(u3*a1) + 9.58540541e+00*(u*eta*a12) + 2.81744332e+00*(u2*eta2) + 3.08636464e-01*(u*a13) + -4.31620356e-01*(u2*a12) + 3.49724954e-01*(u3*eta) + 1.77741156e+01*(u*eta2*a1) + 1.54291378e+00*(eta*a13) + -1.85539216e-01*(eta2*a13) + 3.09106890e+00*(u2*eta*a12) + -7.62862416e+00*(u2*eta2*a1) + -3.43508535e+01*(u*eta2*a12) + -5.19578062e+00*(u*eta*a13) + -1.46730230e-01*(u3*a12) + 3.81602476e-01*(u2*a13) + -1.00147722e+00*(u3*eta2) + -1.01020262e+00*(u3*eta*a1) + 2.00490657e+00*(u3*eta2*a1) + 5.69113996e-02*(u3*a13) + 1.83407486e+01*(u*eta2*a13) + -3.13081796e+00*(u2*eta*a13) + 5.32134301e+00*(u2*eta2*a12) + 4.68731310e-01*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double zeta1;

	// Evaluate fit for this parameter
	zeta1 = 2.12228543e-04 + 1.72765044e-03*(u*u) + 8.42796386e-03*(u*eta) + 4.29469866e-03*a12 + 5.87715482e-03*(u*a1) + -2.20492102e-02*(u*a12) + -2.69340783e-02*(eta2*a1) + -1.25502648e-01*(u*eta*a1) + -4.07456860e-03*(u2*u) + -3.77784160e-02*(u*eta2) + -6.71651789e-03*(u2*a1) + -2.90383421e-02*(u2*eta) + -5.77659664e-02*(eta*a12) + -6.09102904e-03*(a12*a1) + 1.82566085e-02*(u3*a1) + 3.76343406e-01*(u*eta*a12) + 1.40241278e-01*(u2*eta2) + 1.47343295e-01*(u2*eta*a1) + 2.08064136e-02*(u*a13) + 3.84612021e-02*(u3*eta) + 4.27172649e-01*(u*eta2*a1) + 8.70674721e-02*(eta*a13) + 2.56293164e-01*(eta2*a12) + -3.40743477e-01*(eta2*a13) + -1.25324880e-01*(u2*eta*a12) + -7.66102785e-01*(u2*eta2*a1) + -1.13210844e+00*(u*eta2*a12) + -3.04754122e-01*(u*eta*a13) + -1.94285312e-02*(u3*a12) + 1.44777143e-02*(u2*a13) + -8.78626802e-02*(u3*eta2) + -1.34626464e-01*(u3*eta*a1) + 2.54345881e-01*(u3*eta2*a1) + 8.03173880e-03*(u3*a13) + 8.16692351e-01*(u*eta2*a13) + -1.10497094e-01*(u2*eta*a13) + 9.52344098e-01*(u2*eta2*a12) + 5.59848397e-02*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double zeta2;

	// Evaluate fit for this parameter
	zeta2 = 1.17466045e+02*u + -4.71072239e-01 + -7.00444003e+01*(u*u) + -2.13529227e+03*(u*eta) + -4.81259441e+02*(eta*a1) + -2.54337833e+02*a12 + -8.25019324e+02*(u*a1) + 1.82289769e+03*(u*a12) + 2.66693164e+03*(eta2*a1) + 1.45360410e+04*(u*eta*a1) + 9.23710359e+01*(u2*u) + 7.45985689e+03*(u*eta2) + 2.50866542e+02*(u2*a1) + 1.03094794e+03*(u2*eta) + 4.10456314e+03*(eta*a12) + 3.40357143e+02*(a12*a1) + -6.16196181e+02*(u3*a1) + -3.04386531e+04*(u*eta*a12) + -4.55029015e+03*(u2*eta2) + -4.86230463e+03*(u2*eta*a1) + -1.28505090e+03*(u*a13) + -3.97035729e+02*(u3*eta) + -4.98921580e+04*(u*eta2*a1) + -4.92447931e+03*(eta*a13) + -1.51843861e+04*(eta2*a12) + 1.70234529e+04*(eta2*a13) + 3.72786250e+03*(u2*eta*a12) + 2.44729982e+04*(u2*eta2*a1) + 1.00754912e+05*(u*eta2*a12) + 2.02360239e+04*(u*eta*a13) + 8.85351065e+02*(u3*a12) + -4.43984070e+02*(u2*a13) + 2.91057271e+03*(u3*eta*a1) + -2.39608418e+03*(u3*eta2*a1) + -3.43969342e+02*(u3*a13) + -6.36647882e+04*(u*eta2*a13) + 3.55390358e+03*(u2*eta*a13) + -2.97822305e+04*(u2*eta2*a12) + -2.49152006e+03*(u3*eta*a12);

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
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu0;

	// Evaluate fit for this parameter
	nu0 = -2.15575051e+04*a1 + 1.38521730e+03*u + -3.94709851e+04*eta + 3.03679606e+03 + -4.30551177e+03*(u*u) + 1.21440151e+05*eta2 + -7.96183597e+03*(u*eta) + 2.90554449e+05*(eta*a1) + 4.88139438e+04*a12 + -2.17821844e+04*(u*a1) + 6.19225515e+04*(u*a12) + -8.97371851e+05*(eta2*a1) + 2.33882501e+05*(u*eta*a1) + -1.95896699e+03*(u2*u) + 2.85085615e+04*(u2*a1) + 4.36439849e+04*(u2*eta) + -6.42660293e+05*(eta*a12) + -3.36478632e+04*(a12*a1) + 1.11412359e+04*(u3*a1) + -7.49641009e+05*(u*eta*a12) + -8.76986804e+04*(u2*eta2) + -2.55139838e+05*(u2*eta*a1) + -4.71967496e+04*(u*a13) + -4.99654979e+04*(u2*a12) + 1.08725008e+04*(u3*eta) + -5.67272530e+05*(u*eta2*a1) + 4.35842624e+05*(eta*a13) + 1.96090218e+06*(eta2*a12) + -1.31875567e+06*(eta2*a13) + 3.54580633e+05*(u2*eta*a12) + 4.48084775e+05*(u2*eta2*a1) + 2.06672655e+06*(u*eta2*a12) + 6.03273858e+05*(u*eta*a13) + -7.93837442e+03*(u3*a12) + 2.88591348e+04*(u2*a13) + -8.27407566e+04*(u3*eta*a1) + 9.54969922e+04*(u3*eta2*a1) + -1.75651159e+06*(u*eta2*a13) + -1.61066419e+05*(u2*eta*a13) + -3.55507114e+05*(u2*eta2*a12) + 4.91595706e+04*(u3*eta*a12);

	// Return answer
	return nu0;

} // END of NU0 (3,3) fit implementation


#ifdef __cplusplus
}
#endif
