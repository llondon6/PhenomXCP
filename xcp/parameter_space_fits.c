
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
	mu1 = 1.67615435e+00*eta + 2.10903554e+00*a1 + -2.10405191e-01 + -2.15368369e+00*eta2 + -2.74600005e+01*(eta*a1) + -3.97281275e+00*a12 + -9.51780540e-01*(u*u) + -1.36313731e+00*(u*eta) + 1.15197580e+00*(u*a1) + -2.84091288e+01*(u*eta*a1) + -1.81265653e+00*(u*a12) + 1.84695048e+01*(u*eta2) + -4.93343590e-01*(u2*a1) + -5.40199062e-01*(u2*u) + 5.99421776e+01*(eta*a12) + 2.16438896e+01*(u2*eta) + 1.21874839e+02*(eta2*a1) + -1.90898093e+02*(eta3*a1) + 4.74309326e+01*(u*eta*a12) + 1.61109865e+02*(u*eta2*a1) + 3.24997973e+00*(u3*eta) + 3.13411227e+00*(u2*a12) + -1.49381152e+01*(u2*eta*a1) + -2.93651145e+02*(eta2*a12) + 2.83847008e+00*(u3*a1) + -6.01462558e+01*(u*eta2*eta) + -1.29130509e+02*(u2*eta2) + -2.83849908e+02*(u*eta2*a12) + -2.16456930e+02*(u*eta3*a1) + 1.42893120e+02*(u2*eta2*a1) + -2.94068529e+01*(u2*eta*a12) + -1.91357974e+01*(u3*eta*a1) + 2.42687147e+02*(u2*eta2*eta) + 4.76476395e+02*(eta3*a12) + -2.48164412e+00*(u3*a12) + -3.39815589e+02*(u2*eta3*a1) + 4.46931514e+02*(u*eta3*a12) + -1.00245853e+01*(u3*eta2*eta) + 7.50802592e+01*(u2*eta2*a12) + 1.38054997e+01*(u3*eta*a12) + 1.71111868e+01*(u3*eta2*a1);

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
	mu2 = 9.77476726e+00*eta + -5.67336980e-01 + -7.58223852e+01*eta2 + -1.03667717e+00*a12 + -1.99202361e+00*(u*u) + 1.95226932e+01*(u*a1) + -3.43529889e+02*(u*eta*a1) + -3.18223303e+01*(u*a12) + 2.04195924e+02*(eta2*eta) + 1.39907305e+01*(u2*a1) + -6.22564850e+00*(u2*u) + 5.65083871e+01*(eta2*a1) + -2.97831976e+02*(eta3*a1) + 5.24287797e+02*(u*eta*a12) + 1.96097697e+03*(u*eta2*a1) + 1.10986170e+02*(u3*eta) + -1.86685477e+01*(u2*a12) + -8.52685363e+01*(u2*eta*a1) + 1.83239900e+02*(u2*eta2) + -2.86499032e+03*(u*eta2*a12) + -3.57899444e+03*(u*eta3*a1) + -3.26888695e+02*(u2*eta2*a1) + 2.06584624e+02*(u2*eta*a12) + -4.01820346e+01*(u3*eta*a1) + -6.21272324e+02*(u2*eta2*eta) + 1.17498752e+02*(eta3*a12) + 6.60059213e+00*(u3*a12) + -5.80683611e+02*(u3*eta2) + 1.86367826e+03*(u2*eta3*a1) + 5.07681372e+03*(u*eta3*a12) + 9.14199594e+02*(u3*eta2*eta) + -5.43526012e+02*(u2*eta2*a12) + -2.93463287e+01*(u3*eta*a12) + 1.84920115e+02*(u3*eta2*a1);

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
	mu3 = 8.30075438e-03*eta + 1.95995198e-02*u + 6.05894863e-02*a1 + -4.05331175e-03 + -7.51367231e-01*(eta*a1) + -2.80601179e-02*(u*u) + -8.64252003e-02*(u*eta) + -1.03994068e-01*(u*a1) + 4.51527662e-01*(u*eta*a1) + 1.37289858e-01*(u*a12) + 4.76735647e-01*(u2*eta) + 3.69197556e+00*(eta2*a1) + -6.45957802e+00*(eta3*a1) + -6.62047594e-01*(u*eta*a12) + -1.46493908e-01*(u3*eta) + -1.07926903e-02*(u2*a12) + 6.69305062e-02*(u3*a1) + -3.03114304e+00*(u2*eta2) + 9.21763522e-01*(u2*eta2*a1) + 5.99635274e+00*(u2*eta2*eta) + -1.30278599e-01*(u3*a12) + 6.38421498e-01*(u3*eta2) + -2.86995492e+00*(u2*eta3*a1) + 1.21928810e+00*(u*eta3*a12) + 5.76558998e-01*(u3*eta*a12) + -1.30743768e+00*(u3*eta2*a1);

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
	nu4 = -5.28964673e-02*eta + 1.45587076e-02*u + -3.05353548e-02*a1 + 4.78315181e-03 + 1.14312673e-01*eta2 + 3.52904764e-01*(eta*a1) + 2.72274041e-02*a12 + -1.22313330e-02*(u*u) + -2.85261342e-01*(u*eta) + -6.70855046e-02*(u*a1) + 1.46807634e+00*(u*eta*a1) + 3.99952943e-02*(u*a12) + 1.59734801e+00*(u*eta2) + 8.44742867e-02*(u2*a1) + -2.33749854e-01*(eta*a12) + 1.10876247e-01*(u2*eta) + -1.04236543e+00*(eta2*a1) + 6.86147162e-01*(eta3*a1) + -1.19089706e+00*(u*eta*a12) + -8.83436203e+00*(u*eta2*a1) + -8.10956402e-02*(u2*a12) + -9.42365590e-01*(u2*eta*a1) + -1.18085994e-02*(u3*a1) + -2.69012513e+00*(u*eta2*eta) + -2.67952010e-01*(u2*eta2) + 7.81367953e+00*(u*eta2*a12) + 1.58381577e+01*(u*eta3*a1) + 3.20078217e+00*(u2*eta2*a1) + 8.10981705e-01*(u2*eta*a12) + 1.83440578e+00*(eta3*a12) + 3.26132644e-02*(u3*a12) + 1.17641949e-01*(u3*eta2) + -3.05783590e+00*(u2*eta3*a1) + -1.44375098e+01*(u*eta3*a12) + -5.32885613e-01*(u3*eta2*eta) + -1.94040129e+00*(u2*eta2*a12) + -1.55519991e-01*(u3*eta*a12) + 2.79719394e-01*(u3*eta2*a1);

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
	nu5 = -1.08097061e-01*eta + -1.63990385e-02*u + -6.56415474e-02*a1 + 4.32929530e-03 + 7.07571883e-01*eta2 + 8.76194707e-01*(eta*a1) + -2.31469507e-02*(u*u) + 1.03974979e-01*(u*eta) + -1.50210568e-02*(u*a1) + 1.64955198e+00*(u*eta*a1) + -1.33868961e+00*(eta2*eta) + 1.30881731e-01*(u2*a1) + 2.46966917e-02*(u2*u) + -1.15355150e-01*(eta*a12) + 2.48663844e-01*(u2*eta) + -4.91309559e+00*(eta2*a1) + 9.62322686e+00*(eta3*a1) + -1.90385071e+00*(u*eta*a12) + -1.34691594e+01*(u*eta2*a1) + -3.69651659e-01*(u3*eta) + -7.88460378e-02*(u2*a12) + -1.48309975e+00*(u2*eta*a1) + 4.64666973e-01*(eta2*a12) + -3.71962875e-02*(u3*a1) + -4.94424855e-01*(u*eta2*eta) + -6.40117543e-01*(u2*eta2) + 1.57495096e+01*(u*eta2*a12) + 2.80622277e+01*(u*eta3*a1) + 5.85490189e+00*(u2*eta2*a1) + 7.40440142e-01*(u2*eta*a12) + 8.72492431e-02*(u3*a12) + 2.24544666e+00*(u3*eta2) + -7.83723244e+00*(u2*eta3*a1) + -3.22487650e+01*(u*eta3*a12) + -5.10601301e+00*(u3*eta2*eta) + -1.71267097e+00*(u2*eta2*a12) + -4.21527750e-01*(u3*eta*a12) + 9.74610535e-01*(u3*eta2*a1);

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
	nu6 = -1.57810817e-02*eta + 8.03716030e-03*u + 2.09041142e-03 + 2.74606690e-02*(eta*a1) + 7.83568913e-02*a12 + -1.02540957e-02*(u*u) + -2.30355256e-01*(u*a1) + 3.39814551e+00*(u*eta*a1) + 3.29784680e-01*(u*a12) + 4.75695014e-02*(u2*u) + -1.22417837e+00*(eta*a12) + 1.07827854e-01*(u2*eta) + -4.88082410e+00*(u*eta*a12) + -1.85324929e+01*(u*eta2*a1) + -8.73335188e-01*(u3*eta) + 6.18584686e+00*(eta2*a12) + -6.01112365e-01*(u*eta2*eta) + -2.92836282e-01*(u2*eta2) + 2.53495758e+01*(u*eta2*a12) + 3.47151436e+01*(u*eta3*a1) + 5.73308318e-01*(u3*eta*a1) + -1.04904828e+01*(eta3*a12) + -8.09147762e-02*(u3*a12) + 4.20611185e+00*(u3*eta2) + -4.44865221e+01*(u*eta3*a12) + -5.78975494e+00*(u3*eta2*eta) + 3.28514924e-01*(u3*eta*a12) + -2.37522717e+00*(u3*eta2*a1);

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
	zeta1 = 1.84275409e-05*eta + 6.17393794e-05*u + 4.98413221e-05*a1 + -4.24529262e-06 + -3.31694593e-04*(eta*a1) + -3.35912027e-05*a12 + 5.14292355e-05*(u*u) + -1.02770422e-03*(u*eta) + -7.99122093e-05*(u*a1) + 8.08018763e-04*(u*eta*a1) + 4.87711411e-05*(u*a12) + 5.22803207e-03*(u*eta2) + -2.20894248e-04*(u2*a1) + -8.78997448e-05*(u2*u) + -7.03788114e-04*(u2*eta) + 5.20072456e-04*(eta2*a1) + -1.45011523e-03*(u*eta2*a1) + 1.37028084e-03*(u3*eta) + 1.42213003e-04*(u2*a12) + 2.68847892e-03*(u2*eta*a1) + 2.07966422e-03*(eta2*a12) + 1.08707225e-04*(u3*a1) + -8.79586277e-03*(u*eta2*eta) + 3.33220824e-03*(u2*eta2) + -2.72692533e-03*(u*eta2*a12) + -1.18866872e-02*(u2*eta2*a1) + -1.13259900e-03*(u2*eta*a12) + -1.12063142e-03*(u3*eta*a1) + -5.50511930e-03*(u2*eta2*eta) + -6.13465438e-03*(eta3*a12) + -7.82508550e-05*(u3*a12) + -6.55318946e-03*(u3*eta2) + 1.90196893e-02*(u2*eta3*a1) + 5.70727742e-03*(u*eta3*a12) + 1.08242192e-02*(u3*eta2*eta) + 2.19699824e-03*(u2*eta2*a12) + 5.88155169e-04*(u3*eta*a12) + 1.62657629e-03*(u3*eta2*a1);

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
	zeta2 = -3.62738517e+00*u + -9.65003656e+00*a1 + 1.54166524e+00 + -7.73663504e+01*eta2 + 5.80189578e+00*a12 + -1.32488537e+01*(u*u) + 2.77812175e+01*(u*eta) + 1.93168568e+02*(u*eta*a1) + 2.29984000e+02*(eta2*eta) + 4.88610142e+01*(u2*a1) + 1.31772989e+01*(u2*u) + 8.46790335e+01*(eta*a12) + 1.94437587e+02*(u2*eta) + 5.88748232e+02*(eta2*a1) + -1.77452309e+03*(eta3*a1) + -2.01105647e+02*(u*eta*a12) + -1.72901523e+03*(u*eta2*a1) + -1.64714625e+02*(u3*eta) + -3.01762981e+01*(u2*a12) + -6.37655308e+02*(u2*eta*a1) + -1.20235969e+03*(eta2*a12) + -2.72122271e+01*(u3*a1) + -1.66431248e+02*(u*eta2*eta) + -8.99544028e+02*(u2*eta2) + 1.67806297e+03*(u*eta2*a12) + 3.54849647e+03*(u*eta3*a1) + 2.65491187e+03*(u2*eta2*a1) + 2.97572457e+02*(u2*eta*a12) + 2.10865223e+02*(u3*eta*a1) + 1.41421637e+03*(u2*eta2*eta) + 3.11207827e+03*(eta3*a12) + 2.09830504e+01*(u3*a12) + 6.34284175e+02*(u3*eta2) + -3.79613803e+03*(u2*eta3*a1) + -3.18438279e+03*(u*eta3*a12) + -8.47960603e+02*(u3*eta2*eta) + -6.58418481e+02*(u2*eta2*a12) + -1.20328313e+02*(u3*eta*a12) + -2.69624494e+02*(u3*eta2*a1);

	// Return answer
	return zeta2;

} // END of ZETA2 (2,2) fit implementation



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
	mu1 = 5.63840601e+00*eta + -9.57588049e-01*u + 2.15027413e+00*a1 + -4.77487694e-01 + -1.73913871e+01*eta2 + -2.22034421e+01*(eta*a1) + -2.30942677e+00*a12 + 1.29176987e+00*(u*u) + 2.06807660e+01*(u*eta) + 6.94062613e+00*(u*a1) + -1.49026060e+02*(u*eta*a1) + -1.67566769e+01*(u*a12) + -7.29903118e+01*(u*eta2) + -8.24839091e+00*(u2*a1) + -8.99767735e-01*(u2*u) + 1.71958314e+01*(eta*a12) + -1.74680990e+01*(u2*eta) + 5.09236380e+01*(eta2*a1) + 1.22449051e+01*(u*a13) + 3.39169787e+02*(u*eta*a12) + 9.77848738e+00*(eta*a13) + 5.20982611e+02*(u*eta2*a1) + 2.62773317e+00*(u3*eta) + 1.21254166e+01*(u2*a12) + 1.02710429e+02*(u2*eta*a1) + 7.19410740e+00*(u3*a1) + 5.55209603e+01*(u2*eta2) + -1.15503022e+03*(u*eta2*a12) + -2.98139541e+02*(u2*eta2*a1) + -4.16158578e+00*(u2*a13) + -1.32143640e+02*(u2*eta*a12) + -2.34194647e+02*(u*eta*a13) + -7.13890926e+01*(eta2*a13) + -2.67202484e+01*(u3*eta*a1) + -1.20219639e+01*(u3*a12) + 2.53667955e+01*(u3*eta2*a1) + 7.73953276e+02*(u*eta2*a13) + 5.80508600e+00*(u3*a13) + 3.07254578e+02*(u2*eta2*a12) + 2.38077057e+01*(u3*eta*a12) + 2.95253768e+01*(u2*eta*a13);

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
	mu2 = 9.53727879e-01*eta + 5.09819658e-01*u + 4.28648836e-02 + -1.08940605e+01*eta2 + -1.46692505e+01*(eta*a1) + -4.89107432e-01*(u*u) + -3.22581302e+00*(u*eta) + -1.91408506e+00*(u*a1) + 4.21596066e+00*(u*a12) + 8.59812417e+00*(u*eta2) + 7.32105748e-01*(a12*a1) + 3.62380544e+00*(u2*a1) + -1.10051214e+00*(u2*u) + 4.16799162e+01*(eta*a12) + 1.92003509e+00*(u2*eta) + 1.12707964e+02*(eta2*a1) + -3.71512712e+00*(u*a13) + -1.08231907e+01*(u*eta*a12) + -4.31013110e+01*(eta*a13) + 2.16550830e+01*(u*eta2*a1) + 1.01975940e+01*(u3*eta) + -8.72026966e+00*(u2*a12) + -1.70375547e+01*(u2*eta*a1) + -3.00234981e+02*(eta2*a12) + 4.63186590e+00*(u3*a1) + 5.96414504e+00*(u2*a13) + 4.64380294e+01*(u2*eta*a12) + 2.96813968e+01*(u*eta*a13) + 2.64946471e+02*(eta2*a13) + -2.16561165e+01*(u3*eta*a1) + -7.16663037e+00*(u3*a12) + -3.67706066e+01*(u3*eta2) + 6.82879282e+01*(u3*eta2*a1) + -9.17422108e+01*(u*eta2*a13) + 4.06907691e+00*(u3*a13) + -1.08847218e+01*(u2*eta2*a12) + 5.30891258e+00*(u3*eta*a12) + -3.15130139e+01*(u2*eta*a13);

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
	mu3 = -2.79406833e-01*eta + 5.24013826e-02*u + -7.53461338e-02*a1 + 1.52384541e-02 + 1.08248201e+00*eta2 + 1.41973421e+00*(eta*a1) + 1.27783842e-01*a12 + 1.36594643e-02*(u*u) + -4.98696353e-01*(u*eta) + -2.52353915e-01*(u*a1) + 1.41801551e+00*(u*eta*a1) + 4.46087817e-01*(u*a12) + 1.75548774e+00*(u*eta2) + -1.03385649e-01*(a12*a1) + -1.45632318e-01*(u2*a1) + -7.02529639e-02*(u2*u) + -2.26330785e+00*(eta*a12) + -2.52285877e-01*(u2*eta) + -6.27157971e+00*(eta2*a1) + -2.49848045e-01*(u*a13) + -1.08170338e+00*(u*eta*a12) + 1.56684467e+00*(eta*a13) + -5.22898142e+00*(u*eta2*a1) + 5.85170979e-01*(u3*eta) + 2.46009732e-01*(u2*a12) + 2.39772548e+00*(u2*eta*a1) + 1.10421385e+01*(eta2*a12) + 3.06797938e-01*(u3*a1) + 6.44445443e-01*(u2*eta2) + 4.37725526e+00*(u*eta2*a12) + -5.24643652e+00*(u2*eta2*a1) + -4.56655277e-02*(u2*a13) + -4.33609110e+00*(u2*eta*a12) + -1.30859779e-01*(u*eta*a13) + -7.43807228e+00*(eta2*a13) + -7.66589188e-01*(u3*eta*a1) + -6.03081888e-01*(u3*a12) + -1.71914939e+00*(u3*eta2) + 2.07393027e+00*(u3*eta2*a1) + 4.12485642e-01*(u3*a13) + 7.54554403e+00*(u2*eta2*a12) + 7.83641980e-02*(u3*eta*a12) + 1.36351408e+00*(u2*eta*a13);

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
	mu4 = -1.95199610e+02*eta + 1.60038040e+01*u + -7.58496368e+01*a1 + 1.33586487e+01 + 6.41715564e+02*eta2 + 1.09946754e+03*(eta*a1) + 1.33701650e+02*a12 + -2.33286693e+02*(u*eta) + -3.12587555e+01*(u*a1) + 3.26299925e+02*(u*eta*a1) + 3.86629863e+01*(u*a12) + 9.19843801e+02*(u*eta2) + -7.89022643e+01*(a12*a1) + -1.93135649e+01*(u2*a1) + -2.86676588e+01*(u2*u) + -1.93070914e+03*(eta*a12) + -2.94697507e+01*(u2*eta) + -3.78503147e+03*(eta2*a1) + -3.51071055e+01*(u*a13) + 1.12080942e+03*(eta*a13) + -2.03703221e+03*(u*eta2*a1) + 3.86709138e+02*(u3*eta) + 2.46395330e+01*(u2*a12) + 4.66431652e+02*(u2*eta*a1) + 6.85148114e+03*(eta2*a12) + 6.16511865e+01*(u3*a1) + 1.64904738e+02*(u2*eta2) + 1.69437330e+03*(u*eta2*a12) + -1.48321485e+03*(u2*eta2*a1) + -6.82559489e+02*(u2*eta*a12) + -4.03414531e+03*(eta2*a13) + -4.96340753e+02*(u3*eta*a1) + -7.69225287e+01*(u3*a12) + -1.22522249e+03*(u3*eta2) + 1.76574925e+03*(u3*eta2*a1) + -8.79303618e+02*(u*eta2*a13) + 6.75508613e+01*(u3*a13) + 1.52988091e+03*(u2*eta2*a12) + -9.81075731e+01*(u3*eta*a12) + 1.83086276e+02*(u2*eta*a13);

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
	nu4 = 1.86528260e+01*u + 2.62648279e+00*a1 + -1.09839711e+00 + -1.45508472e+01*eta2 + -6.13494463e+01*a12 + -7.65844635e+00*(u*u) + -3.87232533e+02*(u*eta) + -6.82683804e+01*(u*a1) + 1.69182376e+03*(u*eta*a1) + 1.57020210e+02*(u*a12) + 1.42946152e+03*(u*eta2) + 8.07192221e+01*(a12*a1) + 9.00071095e+00*(u2*u) + 5.80094883e+02*(eta*a12) + 1.32601704e+02*(u2*eta) + -1.45914873e+02*(u*a13) + -3.19607644e+03*(u*eta*a12) + -8.27548367e+02*(eta*a13) + -6.35859899e+03*(u*eta2*a1) + 1.46882649e+02*(u3*eta) + 9.98919177e+01*(u2*a12) + -4.05420287e+02*(u2*eta*a1) + -1.44173294e+03*(eta2*a12) + -1.40266772e+02*(u3*a1) + -6.64369335e+02*(u2*eta2) + 1.08025466e+04*(u*eta2*a12) + 2.81825016e+03*(u2*eta2*a1) + -1.56064969e+02*(u2*a13) + 2.30842436e+03*(u*eta*a13) + 2.18385969e+03*(eta2*a13) + 2.93215052e+02*(u3*eta*a1) + 2.65065670e+02*(u3*a12) + -8.78593288e+02*(u3*eta2) + 1.25607404e+03*(u3*eta2*a1) + -6.61407075e+03*(u*eta2*a13) + -9.78586106e+01*(u3*a13) + -3.73715059e+03*(u2*eta2*a12) + -8.33488025e+02*(u3*eta*a12) + 9.15401671e+02*(u2*eta*a13);

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
	nu5 = -1.11946584e+00*eta + 8.57313802e-02*u + -4.67546981e-01*a1 + 6.90824265e-02 + 3.67409767e+00*eta2 + 7.00870315e+00*(eta*a1) + 6.67800273e-01*a12 + -9.80095695e-01*(u*eta) + -3.31679152e-01*(u*a1) + 2.33137168e+00*(u*eta*a1) + 4.52571772e-01*(u*a12) + 3.55060741e+00*(u*eta2) + -3.55191629e-01*(a12*a1) + -9.98040979e-02*(u2*a1) + -8.30752210e-02*(u2*u) + -1.14813182e+01*(eta*a12) + -2.40845954e+01*(eta2*a1) + -2.27940879e-01*(u*a13) + 6.06547416e+00*(eta*a13) + -8.81345943e+00*(u*eta2*a1) + 8.37635094e-01*(u3*eta) + 5.31919530e-01*(u2*a12) + 8.37146640e-01*(u2*eta*a1) + 4.17778528e+01*(eta2*a12) + 2.97418716e-01*(u3*a1) + -4.16195983e-01*(u2*a13) + -3.91334723e+00*(u2*eta*a12) + -2.12672258e+00*(u*eta*a13) + -2.29039859e+01*(eta2*a13) + -1.08417064e+00*(u3*eta*a1) + -5.41247310e-01*(u3*a12) + -2.62913059e+00*(u3*eta2) + 4.75458360e+00*(u3*eta2*a1) + 8.42369051e+00*(u*eta2*a13) + 4.64389741e-01*(u3*a13) + 1.83200547e+00*(u2*eta2*a12) + -5.73890116e-01*(u3*eta*a12) + 2.78650351e+00*(u2*eta*a13);

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
	nu6 = 1.20553806e+00*eta + -5.46938336e-02*u + 5.17835702e-01*a1 + -9.12366138e-02 + -3.63030647e+00*eta2 + -6.81874630e+00*(eta*a1) + -8.96304171e-01*a12 + 1.86030641e-02*(u*u) + 1.03972550e+00*(u*eta) + -3.61151882e-01*(u*a12) + -4.63076387e+00*(u*eta2) + 5.06849183e-01*(a12*a1) + 1.49259910e-01*(u2*u) + 1.23397318e+01*(eta*a12) + -2.02027961e-01*(u2*eta) + 2.15192126e+01*(eta2*a1) + 6.51836950e-01*(u*a13) + 1.58687937e+00*(u*eta*a12) + -7.28306706e+00*(eta*a13) + 5.59694838e+00*(u*eta2*a1) + -2.71061536e+00*(u3*eta) + 5.09497013e-02*(u2*a12) + -7.89431529e-01*(u2*eta*a1) + -4.07903762e+01*(eta2*a12) + -9.99896875e+00*(u*eta2*a12) + 5.84004833e+00*(u2*eta2*a1) + -5.34480863e+00*(u*eta*a13) + 2.48386637e+01*(eta2*a13) + 2.12990955e+00*(u3*eta*a1) + -2.08987874e-01*(u3*a12) + 1.00999928e+01*(u3*eta2) + -1.50260041e+01*(u3*eta2*a1) + 1.64888343e+01*(u*eta2*a13) + -1.54129354e-01*(u3*a13) + -3.59135921e+00*(u2*eta2*a12) + 2.87145177e+00*(u3*eta*a12) + 3.23188258e-01*(u2*eta*a13);

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
	zeta1 = 2.10733853e-02*eta + -1.39377765e-03*u + 9.01930059e-03*a1 + -1.36723901e-03 + -6.97658230e-02*eta2 + -1.12782033e-01*(eta*a1) + -1.29888960e-02*a12 + 3.11661250e-03*(u*u) + 3.10185675e-02*(u*eta) + 8.21695216e-03*(u*a1) + -1.52908747e-01*(u*eta*a1) + -3.39544995e-02*(u*a12) + -1.27196202e-01*(u*eta2) + 4.55262561e-03*(a12*a1) + -1.63864992e-02*(u2*a1) + -2.18860512e-03*(u2*u) + 1.47683674e-01*(eta*a12) + -4.62714328e-02*(u2*eta) + 3.24652803e-01*(eta2*a1) + 3.66046131e-02*(u*a13) + 4.93850855e-01*(u*eta*a12) + -3.45599588e-02*(eta*a13) + 5.35647106e-01*(u*eta2*a1) + -4.29723100e-03*(u3*eta) + 2.01496637e-02*(u2*a12) + 2.49704594e-01*(u2*eta*a1) + -3.49494935e-01*(eta2*a12) + 2.27376502e-02*(u3*a1) + 1.72288803e-01*(u2*eta2) + -1.41047083e+00*(u*eta2*a12) + -8.78738787e-01*(u2*eta2*a1) + 2.78934348e-03*(u2*a13) + -3.28229468e-01*(u2*eta*a12) + -4.74770015e-01*(u*eta*a13) + -1.41694332e-01*(u3*eta*a1) + -2.85694564e-02*(u3*a12) + 1.21787381e-01*(u3*eta2) + 1.22537244e+00*(u*eta2*a13) + 7.73652364e-03*(u3*a13) + 1.09717658e+00*(u2*eta2*a12) + 1.37537398e-01*(u3*eta*a12);

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
	zeta2 = -3.74703320e+02*eta + 6.81473948e+01*u + -2.15452850e+02*a1 + 2.54982207e+01 + 1.32507030e+03*eta2 + 2.46301023e+03*(eta*a1) + 2.73321595e+02*a12 + -6.79491469e+01*(u*u) + -1.51977893e+03*(u*eta) + -3.83084898e+02*(u*a1) + 8.30678662e+03*(u*eta*a1) + 1.10390875e+03*(u*a12) + 5.77111940e+03*(u*eta2) + -3.97641904e+01*(a12*a1) + 3.73336171e+02*(u2*a1) + 7.80576454e+01*(u2*u) + -2.87813732e+03*(eta*a12) + 9.36400027e+02*(u2*eta) + -6.89326273e+03*(eta2*a1) + -1.01995066e+03*(u*a13) + -1.93818274e+04*(u*eta*a12) + -2.97666331e+04*(u*eta2*a1) + 2.11818901e+02*(u3*eta) + -4.58392020e+02*(u2*a12) + -5.62796694e+03*(u2*eta*a1) + 6.38965545e+03*(eta2*a12) + -7.85943501e+02*(u3*a1) + -3.35994264e+03*(u2*eta2) + 6.15073786e+04*(u*eta2*a12) + 1.93574879e+04*(u2*eta2*a1) + -7.66603070e+01*(u2*a13) + 7.78509637e+03*(u2*eta*a12) + 1.53710256e+04*(u*eta*a13) + 2.29684725e+03*(eta2*a13) + 3.79092528e+03*(u3*eta*a1) + 1.13119324e+03*(u3*a12) + -3.76756581e+03*(u3*eta2) + -4.36806568e+04*(u*eta2*a13) + -3.43341707e+02*(u3*a13) + -2.57682849e+04*(u2*eta2*a12) + -4.55830519e+03*(u3*eta*a12) + 1.74109105e+03*(u3*eta2*a1);

	// Return answer
	return zeta2;

} // END of ZETA2 (3,3) fit implementation


#ifdef __cplusplus
}
#endif
