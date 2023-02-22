
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
	mu1 = -8.08951688e-01 + 1.37205020e+01*eta + 4.66663607e+00*a1 + -2.85808603e-01*u + -7.60464929e+01*eta2 + -7.91015855e+01*(eta*a1) + -6.20915494e+00*a12 + 2.52288263e+00*(u*eta) + -5.96662963e-01*(u*u) + 2.39469775e+00*(u*a1) + 1.06116952e+02*(eta*a12) + -3.63470027e+00*(u*a12) + -7.19664455e+00*(u*eta2) + -4.38731385e+01*(u*eta*a1) + 4.35463717e+02*(eta2*a1) + 1.39055965e+02*(eta2*eta) + 1.89226573e+01*(u2*eta) + -2.12094188e+00*(u2*a1) + 7.16417319e+01*(u*eta*a12) + 2.62549060e+02*(u*eta2*a1) + 2.08650484e+00*(u3*a1) + -7.75849371e+02*(eta3*a1) + -1.35457403e+02*(u2*eta2) + -4.74291075e+00*(u3*eta) + 4.85684871e+00*(u2*a12) + -5.76612202e+02*(eta2*a12) + 1.35770484e+02*(u2*eta2*a1) + -4.20986084e+02*(u*eta2*a12) + -1.51716458e+00*(u3*a12) + -4.67903984e+02*(u*eta3*a1) + 2.90383102e+02*(u2*eta2*eta) + -4.86846974e+01*(u2*eta*a12) + -1.71086546e+01*(u3*eta*a1) + 4.24675923e+01*(u3*eta2) + 1.00620363e+03*(eta3*a12) + 2.01811438e+01*(u3*eta2*a1) + -8.19247487e+01*(u3*eta2*eta) + 1.26487349e+02*(u2*eta2*a12) + -4.54388888e+02*(u2*eta3*a1) + 7.43176863e+02*(u*eta3*a12) + 9.37893124e+00*(u3*eta*a12);

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
	mu2 = -2.67700472e+00 + 5.96978178e+01*eta + 1.04284118e+01*a1 + -7.08227536e-01*u + -4.13874294e+02*eta2 + -2.67182780e+02*(eta*a1) + -1.20418365e+01*a12 + -5.49336799e+00*(u*u) + 3.87954349e+01*(u*a1) + 2.96129570e+02*(eta*a12) + -5.21863160e+01*(u*a12) + -6.44206173e+02*(u*eta*a1) + 1.97123813e+03*(eta2*a1) + 8.73668084e+02*(eta2*eta) + 7.96003736e+01*(u2*eta) + -1.06264513e+01*(u2*u) + 2.11330097e+01*(u2*a1) + 8.22762454e+02*(u*eta*a12) + 3.72164064e+03*(u*eta2*a1) + -5.37473165e+00*(u3*a1) + -4.29976147e+03*(eta3*a1) + 4.73861754e+01*(u*eta2*eta) + -3.97976493e+02*(u2*eta2) + -2.11714413e+02*(u2*eta*a1) + 2.14010943e+02*(u3*eta) + -2.44186972e+01*(u2*a12) + -2.17216620e+03*(eta2*a12) + 4.73848665e+02*(u2*eta2*a1) + -4.50819265e+03*(u*eta2*a12) + 1.38688856e+01*(u3*a12) + -7.09588771e+03*(u*eta3*a1) + 7.23776169e+02*(u2*eta2*eta) + 2.58332038e+02*(u2*eta*a12) + -2.23444932e+01*(u3*eta*a1) + -1.23974580e+03*(u3*eta2) + 4.74236050e+03*(eta3*a12) + 1.69293148e+02*(u3*eta2*a1) + 2.23960788e+03*(u3*eta2*eta) + -6.18925240e+02*(u2*eta2*a12) + 8.24494529e+03*(u*eta3*a12) + -5.55144959e+01*(u3*eta*a12);

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
	mu3 = -4.44622022e-02 + 6.53238149e-01*eta + 2.70599070e-01*a1 + 1.71367750e-02*u + -3.05012794e+00*eta2 + -4.09038151e+00*(eta*a1) + -2.41954824e-01*a12 + -4.90322793e-02*(u*eta) + -1.17949241e-01*(u*a1) + 3.89760528e+00*(eta*a12) + 2.37884771e-01*(u*a12) + 6.78513428e-01*(u*eta*a1) + 1.94129340e+01*(eta2*a1) + 4.42174673e+00*(eta2*eta) + -1.18368571e-01*(u2*a1) + -2.43035117e+00*(u*eta*a12) + -6.34873441e-01*(u*eta2*a1) + -2.90543071e+01*(eta3*a1) + -5.28602147e-01*(u*eta2*eta) + -7.80959120e-01*(u2*eta2) + 1.72528166e+00*(u2*eta*a1) + 1.12053546e-01*(u2*a12) + -1.87308048e+01*(eta2*a12) + -4.93109506e+00*(u2*eta2*a1) + 8.23519408e+00*(u*eta2*a12) + -1.04736381e-01*(u3*a12) + 2.91127498e+00*(u2*eta2*eta) + -1.70907120e+00*(u2*eta*a12) + 7.69581785e-01*(u3*eta*a1) + -1.26039917e+00*(u3*eta2) + 2.77095313e+01*(eta3*a12) + -3.43386592e+00*(u3*eta2*a1) + 5.40075316e+00*(u3*eta2*eta) + 5.01703159e+00*(u2*eta2*a12) + -9.88779375e+00*(u*eta3*a12) + 4.89717938e-01*(u3*eta*a12);

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
	nu4 = 3.80863025e-03 + -3.20366991e-02*eta + -1.15610971e-02*a1 + 8.32703313e-04*u + 7.36387742e-03*a12 + -4.54196144e-03*(u*eta) + -2.21691423e-02*(u*u) + -6.27129914e-02*(u*a1) + 1.27410244e-01*(eta*a12) + 5.77768096e-02*(u*a12) + 1.17717192e+00*(u*eta*a1) + 9.32982474e-01*(eta2*a1) + 2.21754703e-01*(eta2*eta) + 2.42511599e-01*(u2*eta) + 1.87251931e-02*(u2*u) + 1.09509418e-01*(u2*a1) + -1.15656628e+00*(u*eta*a12) + -6.74877137e+00*(u*eta2*a1) + -2.44003744e-02*(u3*a1) + -2.84903558e+00*(eta3*a1) + -8.24235314e-01*(u2*eta2) + -1.17367510e+00*(u2*eta*a1) + -3.26165542e-01*(u3*eta) + -1.04073876e-01*(u2*a12) + -1.94700813e+00*(eta2*a12) + 3.78642233e+00*(u2*eta2*a1) + 6.42475020e+00*(u*eta2*a12) + 1.70161915e-02*(u3*a12) + 1.22766545e+01*(u*eta3*a1) + 7.72559109e-01*(u2*eta2*eta) + 9.91792290e-01*(u2*eta*a12) + 2.37118837e-01*(u3*eta*a1) + 1.76261766e+00*(u3*eta2) + 5.16248152e+00*(eta3*a12) + -5.99281511e-01*(u3*eta2*a1) + -2.99197839e+00*(u3*eta2*eta) + -2.29782782e+00*(u2*eta2*a12) + -3.32337285e+00*(u2*eta3*a1) + -1.09601110e+01*(u*eta3*a12) + -6.19701167e-02*(u3*eta*a12);

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
	nu5 = 5.60386107e-03 + -1.09386310e-01*a1 + -9.76134307e-03*u + -6.71062874e-01*eta2 + 1.06043628e+00*(eta*a1) + 9.75096452e-02*a12 + -7.87088264e-03*(u*u) + -1.01758786e-01*(u*a1) + -1.25907546e+00*(eta*a12) + 3.34106089e+00*(u*eta*a1) + -2.51413694e+00*(eta2*a1) + 2.29944332e+00*(eta2*eta) + 6.04133041e-02*(u2*u) + 1.23865298e-01*(u2*a1) + -2.23410641e+00*(u*eta*a12) + -2.20146094e+01*(u*eta2*a1) + 8.62448438e-01*(u*eta2*eta) + 6.32924814e-01*(u2*eta2) + -1.21017505e+00*(u2*eta*a1) + -1.12755879e+00*(u3*eta) + -1.43412324e-01*(u2*a12) + 3.46069779e+00*(eta2*a12) + 2.94249774e+00*(u2*eta2*a1) + 1.84737693e+01*(u*eta2*a12) + 1.17803732e-01*(u3*a12) + 4.00894518e+01*(u*eta3*a1) + -2.11011862e+00*(u2*eta2*eta) + 1.60867202e+00*(u2*eta*a12) + -7.33265400e-01*(u3*eta*a1) + 7.58809705e+00*(u3*eta2) + 3.35254056e+00*(u3*eta2*a1) + -1.65637129e+01*(u3*eta2*eta) + -4.19186295e+00*(u2*eta2*a12) + -3.71958729e+01*(u*eta3*a12) + -5.70513147e-01*(u3*eta*a12);

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
	nu6 = 5.95383574e-02 + -1.17522261e+00*eta + -4.10668251e-01*a1 + 6.85791213e+00*eta2 + 8.33253217e+00*(eta*a1) + 5.73468307e-01*a12 + 6.37822009e-01*(u*eta) + 6.77667891e-02*(u*u) + -2.21965533e-01*(u*a1) + -1.12958001e+01*(eta*a12) + 1.83496379e-01*(u*a12) + -5.51679820e+00*(u*eta2) + 1.12119091e+00*(u*eta*a1) + -4.96322865e+01*(eta2*a1) + -1.24530137e+01*(eta2*eta) + -1.22806347e+00*(u2*eta) + 6.25101589e-02*(u2*u) + -1.70657057e-01*(u2*a1) + -1.10098577e+00*(u*eta2*a1) + 2.57990492e-01*(u3*a1) + 9.15265190e+01*(eta3*a1) + 1.21312978e+01*(u*eta2*eta) + 7.34805053e+00*(u2*eta2) + 2.23485134e+00*(u2*eta*a1) + -1.90056218e+00*(u3*eta) + 1.11490182e-01*(u2*a12) + 6.67632335e+01*(eta2*a12) + -1.07137695e+01*(u2*eta2*a1) + -5.42019511e+00*(u*eta2*a12) + -2.67677050e-01*(u3*a12) + -1.48533766e+01*(u2*eta2*eta) + -7.89221718e-01*(u2*eta*a12) + -1.38615096e+00*(u3*eta*a1) + 1.34637726e+01*(u3*eta2) + -1.22999151e+02*(eta3*a12) + 1.90584246e+00*(u3*eta2*a1) + -2.81630956e+01*(u3*eta2*eta) + 8.96504748e-01*(u2*eta2*a12) + 2.04007543e+01*(u2*eta3*a1) + 1.02279645e+01*(u*eta3*a12) + 9.95938528e-01*(u3*eta*a12);

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
	zeta1 = -6.98755515e-06 + 6.73068445e-05*u + 1.20938955e-03*eta2 + 7.40522845e-04*(eta*a1) + 9.93940487e-05*a12 + -1.91996667e-03*(u*eta) + 5.12072679e-05*(u*u) + 1.25088717e-04*(u*a1) + -2.60888702e-03*(eta*a12) + -5.64279088e-05*(u*a12) + 1.38952111e-02*(u*eta2) + -8.40923610e-03*(eta2*a1) + -4.34104186e-03*(eta2*eta) + -6.95810283e-04*(u2*eta) + -1.55100187e-04*(u2*u) + -3.10318275e-04*(u2*a1) + -2.77120842e-04*(u*eta*a12) + -8.22666841e-03*(u*eta2*a1) + -1.05758270e-04*(u3*a1) + 2.15732224e-02*(eta3*a1) + -2.93973437e-02*(u*eta2*eta) + 2.06871621e-03*(u2*eta2) + 5.37116081e-03*(u2*eta*a1) + 3.32545077e-03*(u3*eta) + 8.85373648e-05*(u2*a12) + 1.93406123e-02*(eta2*a12) + -2.71443207e-02*(u2*eta2*a1) + 6.46032374e-03*(u*eta2*a12) + 2.60740186e-02*(u*eta3*a1) + -1.18505133e-03*(u2*eta*a12) + 1.02084953e-03*(u3*eta*a1) + -2.14842153e-02*(u3*eta2) + -4.17935678e-02*(eta3*a12) + -2.78534671e-03*(u3*eta2*a1) + 4.30182627e-02*(u3*eta2*eta) + 3.39578939e-03*(u2*eta2*a12) + 4.19967865e-02*(u2*eta3*a1) + -1.90205252e-02*(u*eta3*a12) + 1.18748965e-04*(u3*eta*a12);

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
	zeta2 = 4.13491487e+00 + -3.84727611e+01*eta + -1.68532906e+01*a1 + 1.16023146e+02*(eta*a1) + 1.14013547e+01*a12 + 5.98778793e+01*(u*eta) + -1.50117382e+01*(u*u) + -5.43973840e+01*(u*a1) + 3.95041175e+01*(u*a12) + -6.04923546e+02*(u*eta2) + 8.10097285e+02*(u*eta*a1) + 3.00315041e+02*(eta2*a1) + 3.89675895e+02*(eta2*eta) + 1.84127385e+02*(u2*eta) + 2.28281853e+01*(u2*u) + 6.94689745e+01*(u2*a1) + -6.44768753e+02*(u*eta*a12) + -3.88754638e+03*(u*eta2*a1) + -4.94876917e+00*(u3*a1) + -2.04805112e+03*(eta3*a1) + 1.50728376e+03*(u*eta2*eta) + -5.01989298e+02*(u2*eta2) + -9.46751204e+02*(u2*eta*a1) + -4.16955342e+02*(u3*eta) + -4.78832575e+01*(u2*a12) + -1.02709937e+03*(eta2*a12) + 3.51253661e+03*(u2*eta2*a1) + 3.29645407e+03*(u*eta2*a12) + 1.47960097e+01*(u3*a12) + 5.82048564e+03*(u*eta3*a1) + 5.76705010e+02*(u2*eta*a12) + 2.47150678e+03*(u3*eta2) + 3.42312199e+03*(eta3*a12) + 1.80594510e+02*(u3*eta2*a1) + -4.75689263e+03*(u3*eta2*eta) + -1.54231811e+03*(u2*eta2*a12) + -3.32158560e+03*(u2*eta3*a1) + -5.11163668e+03*(u*eta3*a12) + -8.83908467e+01*(u3*eta*a12);

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
	nu0 = -8.34730807e+01 + 3.41633942e+03*eta + -9.79587101e+02*u + -1.04614442e+04*eta2 + 1.45504523e+04*(eta*a1) + 1.92217508e+03*a12 + 7.25112083e+03*(u*eta) + 2.90503462e+03*(u*u) + 7.17154065e+03*(u*a1) + -4.98627130e+04*(eta*a12) + -4.36915940e+03*(u*a12) + -9.75899438e+04*(u*eta*a1) + -1.61292414e+05*(eta2*a1) + -4.20163451e+04*(u2*eta) + -1.56022876e+03*(u2*u) + -9.41786186e+03*(u2*a1) + 6.33828928e+04*(u*eta*a12) + 4.22665864e+05*(u*eta2*a1) + -4.20614083e+02*(u3*a1) + 3.95172545e+05*(eta3*a1) + -5.93021268e+04*(u*eta2*eta) + 1.91569247e+05*(u2*eta2) + 1.33693135e+05*(u2*eta*a1) + 2.96399040e+04*(u3*eta) + 4.69422993e+03*(u2*a12) + 3.70231324e+05*(eta2*a12) + -5.82299440e+05*(u2*eta2*a1) + -2.80644357e+05*(u*eta2*a12) + -5.62902896e+05*(u*eta3*a1) + -2.72072542e+05*(u2*eta2*eta) + -5.40481819e+04*(u2*eta*a12) + -1.71817660e+05*(u3*eta2) + -7.92708486e+05*(eta3*a12) + 3.17306362e+05*(u3*eta2*eta) + 1.42115361e+05*(u2*eta2*a12) + 7.81516176e+05*(u2*eta3*a1) + 3.62517211e+05*(u*eta3*a12) + 1.95046844e+03*(u3*eta*a12);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu1;

	// Evaluate fit for this parameter
	mu1 = 8.50572933e-02 + -2.25863377e+00*eta + -5.07166639e+00*a1 + 3.05502445e-01*u + 7.40804108e+00*eta2 + 8.12913395e+01*(eta*a1) + 1.72494649e+01*a12 + -1.80861187e+00*(u*eta) + 6.18796135e+00*(u*u) + -4.72968492e+00*(u*a1) + -2.63742278e+02*(eta*a12) + 1.38695185e+01*(u*a12) + 1.26439108e+01*(u*eta2) + 5.09480761e+01*(u*eta*a1) + -2.90775577e+02*(eta2*a1) + -1.46974579e+01*(a12*a1) + -8.81247615e+01*(u2*eta) + -3.30853692e+00*(u2*u) + -1.82789325e+01*(u2*a1) + -1.67764145e+02*(u*eta*a12) + -2.20534177e+02*(u*eta2*a1) + 2.95452421e+01*(u3*a1) + -1.05799605e+01*(u*a13) + 2.21542534e+02*(eta*a13) + 2.87038121e+02*(u2*eta2) + 2.34212165e+02*(u2*eta*a1) + 4.49249118e+01*(u3*eta) + -9.59986064e+00*(u3*u) + 9.34599845e+02*(eta2*a12) + -6.13684517e+02*(u2*eta2*a1) + 6.78885478e+02*(u*eta2*a12) + 4.11198793e+01*(u4*a1) + -7.79615841e+02*(eta2*a13) + -6.93173506e+01*(u3*a12) + 1.35095819e+02*(u*eta*a13) + 7.60929360e+01*(u2*eta*a12) + -4.06365082e+02*(u3*eta*a1) + 2.01496443e+01*(u2*a13) + 1.35408447e+02*(u4*eta) + -1.60903259e+02*(u3*eta2) + 1.44055957e+03*(u3*eta2*a1) + -4.30013968e+02*(u4*eta2) + -5.10077050e+01*(u4*a12) + -5.52030001e+02*(u4*eta*a1) + -3.48393865e+02*(u2*eta*a13) + -7.09777053e+02*(u2*eta2*a12) + 9.70224678e+02*(u3*eta*a12) + -5.38819285e+02*(u*eta2*a13) + 4.78720686e+01*(u3*a13) + 1.48412158e+01*(u4*a13) + 1.56200615e+03*(u4*eta2*a1) + -6.81479585e+02*(u3*eta*a13) + -3.45455525e+03*(u3*eta2*a12) + 1.50239407e+03*(u2*eta2*a13) + 6.39351664e+02*(u4*eta*a12) + -1.44034422e+02*(u4*eta*a13) + 2.44851711e+03*(u3*eta2*a13) + -1.45599368e+03*(u4*eta2*a12);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu2;

	// Evaluate fit for this parameter
	mu2 = -2.91059924e-01 + 6.15801994e+00*eta + 4.36582211e+00*a1 + -3.93740384e+00*u + -2.68441571e+01*eta2 + -7.85989205e+01*(eta*a1) + -1.02163908e+01*a12 + 4.63121895e+01*(u*eta) + -3.84048896e+00*(u*u) + 3.46577083e+01*(u*a1) + 1.87568146e+02*(eta*a12) + -8.21545587e+01*(u*a12) + -1.20139633e+02*(u*eta2) + -4.10584712e+02*(u*eta*a1) + 3.04674381e+02*(eta2*a1) + 7.14069907e+00*(a12*a1) + 3.77550633e+01*(u2*eta) + 5.11552828e+00*(u2*u) + 8.00221890e+00*(u2*a1) + 9.73893926e+02*(u*eta*a12) + 1.09764264e+03*(u*eta2*a1) + -4.69545875e+01*(u3*a1) + 5.72304711e+01*(u*a13) + -1.32751771e+02*(eta*a13) + -1.01123188e+02*(u2*eta2) + -5.58735169e+01*(u3*eta) + -4.33242887e+00*(u2*a12) + 5.66152778e+00*(u3*u) + -7.22836757e+02*(eta2*a12) + -1.09146815e+02*(u2*eta2*a1) + -2.63263233e+03*(u*eta2*a12) + -1.74143094e+01*(u4*a1) + 5.11354263e+02*(eta2*a13) + 1.15348469e+02*(u3*a12) + -6.77764009e+02*(u*eta*a13) + -1.66979775e+02*(u2*eta*a12) + 5.34295936e+02*(u3*eta*a1) + -6.44501280e+01*(u4*eta) + 1.24590191e+02*(u3*eta2) + -1.31584725e+03*(u3*eta2*a1) + 1.84369778e+02*(u4*eta2) + 2.05849370e+01*(u4*a12) + 1.38088186e+02*(u4*eta*a1) + 1.42138845e+02*(u2*eta*a13) + 6.85620238e+02*(u2*eta2*a12) + -1.34112641e+03*(u3*eta*a12) + 1.84521665e+03*(u*eta2*a13) + -8.26129309e+01*(u3*a13) + -9.06176911e+00*(u4*a13) + -3.34526439e+02*(u4*eta2*a1) + 9.72768074e+02*(u3*eta*a13) + 3.44366411e+03*(u3*eta2*a12) + -4.98388064e+02*(u2*eta2*a13) + -7.91441564e+01*(u4*eta*a12) + -2.56072167e+03*(u3*eta2*a13) + 1.41302188e+02*(u4*eta2*a12);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu3;

	// Evaluate fit for this parameter
	mu3 = 1.07951856e-02 + -2.30911187e-01*eta + -3.20995807e-02*a1 + 4.91254343e-03*u + 8.43842378e-01*eta2 + 1.08745455e+00*(eta*a1) + 1.18287667e-02*a12 + -1.00242905e-01*(u*u) + -2.78050932e-02*(u*a1) + -1.41451939e+00*(eta*a12) + 2.84109605e-02*(u*a12) + 6.24662233e-02*(u*eta2) + -7.05612749e-02*(u*eta*a1) + -4.79463055e+00*(eta2*a1) + 1.85890599e+00*(u2*eta) + -4.94744484e-03*(u2*u) + 2.02904991e-01*(u2*a1) + 5.67490802e-01*(u*eta*a12) + -5.81447043e-01*(u*eta2*a1) + 6.46303345e-03*(u3*a1) + 1.70067062e-02*(u*a13) + 7.22842514e-01*(eta*a13) + -6.81817259e+00*(u2*eta2) + -6.66453808e+00*(u2*eta*a1) + 4.45184810e-02*(u3*eta) + 3.22107173e-01*(u2*a12) + 1.31408013e-01*(u3*u) + 7.70884651e+00*(eta2*a12) + 2.94433192e+01*(u2*eta2*a1) + -3.09041540e-01*(u4*a1) + -4.51658476e+00*(eta2*a13) + -9.21641490e-01*(u*eta*a13) + 4.00661619e+00*(u2*eta*a12) + -5.76938956e-01*(u2*a13) + -2.38966375e+00*(u4*eta) + -5.32312290e-01*(u3*eta2) + 3.20219557e+00*(u3*eta2*a1) + 8.34597545e+00*(u4*eta2) + -3.24942589e-01*(u4*a12) + 8.72849386e+00*(u4*eta*a1) + 2.48547651e+00*(u2*eta*a13) + -3.27000447e+01*(u2*eta2*a12) + -3.01235420e-01*(u3*eta*a12) + 1.95007224e+00*(u*eta2*a13) + -3.06950146e-02*(u3*a13) + 7.51190493e-01*(u4*a13) + -3.48513379e+01*(u4*eta2*a1) + 7.49300604e-01*(u3*eta*a13) + -5.78639402e+00*(u3*eta2*a12) + 6.27094210e+00*(u2*eta2*a13) + -5.29426728e+00*(u4*eta*a12) + -3.99363642e+00*(u4*eta*a13) + 1.61775302e+00*(u3*eta2*a13) + 3.44640913e+01*(u4*eta2*a12);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double mu4;

	// Evaluate fit for this parameter
	mu4 = -3.44711781e-01 + 2.37691177e+01*a1 + -6.68275884e+01*u + -2.50841961e+02*(eta*a1) + -9.23648456e+01*a12 + 9.05242487e+02*(u*eta) + -2.97795217e+01*(u*u) + 4.23474700e+02*(u*a1) + 1.02647410e+03*(eta*a12) + -8.33455241e+02*(u*a12) + -2.98867561e+03*(u*eta2) + -5.52216060e+03*(u*eta*a1) + 6.05552752e+02*(eta2*a1) + 7.95277112e+01*(a12*a1) + 4.73002054e+02*(u2*eta) + 1.48297121e+02*(u2*u) + 1.05392951e+04*(u*eta*a12) + 1.77218988e+04*(u*eta2*a1) + -9.27661661e+02*(u3*a1) + 5.16673635e+02*(u*a13) + -9.07557256e+02*(eta*a13) + -1.81487754e+03*(u2*eta2) + -5.92187985e+02*(u2*eta*a1) + -2.06525638e+03*(u3*eta) + 2.75481213e+02*(u2*a12) + 5.71463394e+01*(u3*u) + -2.61980490e+03*(eta2*a12) + 4.69714442e+03*(u2*eta2*a1) + -3.30910724e+04*(u*eta2*a12) + -1.10507602e+02*(u4*a1) + 2.38006082e+03*(eta2*a13) + 1.77137215e+03*(u3*a12) + -6.39468402e+03*(u*eta*a13) + -2.25044815e+03*(u2*eta*a12) + 1.25254804e+04*(u3*eta*a1) + -2.82733188e+02*(u2*a13) + -8.88114357e+02*(u4*eta) + 6.79518870e+03*(u3*eta2) + -4.05589375e+04*(u3*eta2*a1) + 3.23726484e+03*(u4*eta2) + -1.11643195e+02*(u4*a12) + 2.22242318e+03*(u4*eta*a1) + 2.75630650e+03*(u2*eta*a13) + -2.34107960e+04*(u3*eta*a12) + 1.97670908e+04*(u*eta2*a13) + -1.06001999e+03*(u3*a13) + 2.04263358e+02*(u4*a13) + -1.04222785e+04*(u4*eta2*a1) + 1.37965642e+04*(u3*eta*a13) + 7.49827137e+04*(u3*eta2*a12) + -3.54905612e+03*(u2*eta2*a13) + -1.77152088e+03*(u4*eta*a13) + -4.38924910e+04*(u3*eta2*a13) + 7.99691361e+03*(u4*eta2*a12);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu4;

	// Evaluate fit for this parameter
	nu4 = -6.22204388e+01 + 8.57935885e+02*eta + 5.08813003e+02*a1 + -7.49239291e+01*u + -2.74842720e+03*eta2 + -6.96789639e+03*(eta*a1) + -1.27361376e+03*a12 + 9.45758690e+02*(u*eta) + 3.10641649e+02*(u*a1) + 1.68150663e+04*(eta*a12) + -3.21417564e+02*(u*a12) + -2.89396750e+03*(u*eta2) + -3.48831301e+03*(u*eta*a1) + 2.16188178e+04*(eta2*a1) + 9.20199339e+02*(a12*a1) + 2.51069435e+02*(u2*u) + -6.95764238e+02*(u2*a1) + 2.77363129e+03*(u*eta*a12) + 1.05429794e+04*(u*eta2*a1) + -1.25882879e+03*(u3*a1) + 7.47436481e+01*(u*a13) + -1.18405899e+04*(eta*a13) + -5.86106720e+02*(u2*eta2) + 8.07525678e+03*(u2*eta*a1) + -3.44180351e+03*(u3*eta) + 2.55259302e+03*(u2*a12) + 2.49119858e+02*(u3*u) + -5.05501548e+04*(eta2*a12) + -1.65852069e+04*(u2*eta2*a1) + -8.31007173e+03*(u*eta2*a12) + -5.66559060e+02*(u4*a1) + 3.48690539e+04*(eta2*a13) + 1.97244431e+03*(u3*a12) + -2.90490447e+04*(u2*eta*a12) + 1.63205006e+04*(u3*eta*a1) + -2.09301975e+03*(u2*a13) + -3.60598597e+03*(u4*eta) + 1.10928853e+04*(u3*eta2) + -5.25630413e+04*(u3*eta2*a1) + 1.23490159e+04*(u4*eta2) + -3.47692983e+02*(u4*a12) + 9.41519710e+03*(u4*eta*a1) + 2.32728026e+04*(u2*eta*a13) + 6.43316199e+04*(u2*eta2*a12) + -2.40499344e+04*(u3*eta*a12) + -1.00440380e+03*(u3*a13) + 7.70181112e+02*(u4*a13) + -3.99949675e+04*(u4*eta2*a1) + 1.14798099e+04*(u3*eta*a13) + 7.73869745e+04*(u3*eta2*a12) + -5.14595699e+04*(u2*eta2*a13) + -6.58612673e+03*(u4*eta*a13) + -3.69303823e+04*(u3*eta2*a13) + 2.76748052e+04*(u4*eta2*a12);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu5;

	// Evaluate fit for this parameter
	nu5 = 5.47420530e-02 + -8.94959955e-01*eta + -4.48211948e-01*a1 + -2.21470356e-01*u + 2.42234914e+00*eta2 + 6.71731558e+00*(eta*a1) + 7.29212101e-01*a12 + 2.82963459e+00*(u*eta) + -2.49657050e-01*(u*u) + 1.60871578e+00*(u*a1) + -1.25507705e+01*(eta*a12) + -3.45244308e+00*(u*a12) + -8.82522769e+00*(u*eta2) + -1.99814549e+01*(u*eta*a1) + -1.87589724e+01*(eta2*a1) + -3.87338970e-01*(a12*a1) + 4.47021897e+00*(u2*eta) + 3.80340200e-01*(u2*u) + 1.57991428e+00*(u2*a1) + 4.23232699e+01*(u*eta*a12) + 6.08965476e+01*(u*eta2*a1) + -2.68760336e+00*(u3*a1) + 2.16596622e+00*(u*a13) + 6.86548870e+00*(eta*a13) + -1.42062196e+01*(u2*eta2) + -2.82748082e+01*(u2*eta*a1) + -4.76896655e+00*(u3*eta) + -2.40782694e+00*(u2*a12) + 1.62299907e-01*(u3*u) + 3.54074016e+01*(eta2*a12) + 8.50982105e+01*(u2*eta2*a1) + -1.28448797e+02*(u*eta2*a12) + -8.26178009e-01*(u4*a1) + -1.90815366e+01*(eta2*a13) + 5.47724096e+00*(u3*a12) + -2.69178341e+01*(u*eta*a13) + 4.87326264e+01*(u2*eta*a12) + 3.23318072e+01*(u3*eta*a1) + 1.06991991e+00*(u2*a13) + -3.15547322e+00*(u4*eta) + 1.44500261e+01*(u3*eta2) + -9.53715451e+01*(u3*eta2*a1) + 9.19093024e+00*(u4*eta2) + 8.02899512e-01*(u4*a12) + 1.60547765e+01*(u4*eta*a1) + -2.57247994e+01*(u2*eta*a13) + -1.45418579e+02*(u2*eta2*a12) + -6.46360726e+01*(u3*eta*a12) + 8.21871383e+01*(u*eta2*a13) + -3.29178151e+00*(u3*a13) + -3.91330215e+01*(u4*eta2*a1) + 3.86105027e+01*(u3*eta*a13) + 1.88828972e+02*(u3*eta2*a12) + 7.75676818e+01*(u2*eta2*a13) + -1.98839788e+01*(u4*eta*a12) + 5.52547807e+00*(u4*eta*a13) + -1.12419319e+02*(u3*eta2*a13) + 3.56323152e+01*(u4*eta2*a12);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu6;

	// Evaluate fit for this parameter
	nu6 = 1.54909480e-01 + -2.08714906e+00*eta + -1.25588088e+00*a1 + 5.07546962e-01*u + 6.33331173e+00*eta2 + 1.65217795e+01*(eta*a1) + 3.04124272e+00*a12 + -6.76162592e+00*(u*eta) + -1.64037365e-01*(u*u) + -2.81497262e+00*(u*a1) + -3.87853812e+01*(eta*a12) + 4.83386792e+00*(u*a12) + 2.25615893e+01*(u*eta2) + 3.55020645e+01*(u*eta*a1) + -4.86971576e+01*(eta2*a1) + -2.14294133e+00*(a12*a1) + 1.37416267e+00*(u2*eta) + -1.28339179e+00*(u2*u) + 2.12994363e+00*(u2*a1) + -5.68816054e+01*(u*eta*a12) + -1.15376038e+02*(u*eta2*a1) + 7.27743840e+00*(u3*a1) + -2.60911386e+00*(u*a13) + 2.68511417e+01*(eta*a13) + -2.22534951e+01*(u2*eta*a1) + 1.80292359e+01*(u3*eta) + -5.63482139e+00*(u2*a12) + -1.25979511e-01*(u3*u) + 1.12229992e+02*(eta2*a12) + 3.51554271e+01*(u2*eta2*a1) + 1.77777853e+02*(u*eta2*a12) + -7.16823689e-01*(u4*a1) + -7.70138974e+01*(eta2*a13) + -1.26466647e+01*(u3*a12) + 2.80344096e+01*(u*eta*a13) + 6.03483851e+01*(u2*eta*a12) + -9.83804616e+01*(u3*eta*a1) + 3.79888378e+00*(u2*a13) + 2.69848451e+00*(u4*eta) + -6.04160302e+01*(u3*eta2) + 3.25634558e+02*(u3*eta2*a1) + -1.37243243e+01*(u4*eta2) + 2.91765440e+00*(u4*a12) + 3.75739124e+00*(u4*eta*a1) + -3.96343117e+01*(u2*eta*a13) + -1.07181492e+02*(u2*eta2*a12) + 1.64595509e+02*(u3*eta*a12) + -8.24751805e+01*(u*eta2*a13) + 6.84390932e+00*(u3*a13) + -1.99009931e+00*(u4*a13) + 2.73884033e+01*(u4*eta2*a1) + -8.53249987e+01*(u3*eta*a13) + -5.37275220e+02*(u3*eta2*a12) + 6.60402931e+01*(u2*eta2*a13) + -2.78922681e+01*(u4*eta*a12) + 1.91644698e+01*(u4*eta*a13) + 2.73657124e+02*(u3*eta2*a13);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double zeta1;

	// Evaluate fit for this parameter
	zeta1 = 3.30946341e-05 + -1.62161945e-02*eta + 1.62971120e-02*u + 1.03827265e-01*eta2 + 1.32726762e-01*(eta*a1) + 7.09892037e-03*a12 + -1.68687468e-01*(u*eta) + 4.51264361e-02*(u*u) + -9.90657510e-02*(u*a1) + -3.42012253e-01*(eta*a12) + 2.04506535e-01*(u*a12) + 3.21465888e-01*(u*eta2) + 1.00354659e+00*(u*eta*a1) + -7.72528372e-01*(eta2*a1) + -7.84179495e-03*(a12*a1) + -4.85290336e-01*(u2*eta) + -2.60288829e-01*(u2*a1) + -2.07615806e+00*(u*eta*a12) + -2.02093225e+00*(u*eta2*a1) + -3.59832282e-02*(u3*a1) + -1.39864018e-01*(u*a13) + 2.45368998e-01*(eta*a13) + 1.17276200e+00*(u2*eta2) + 2.80086171e+00*(u2*eta*a1) + -1.53344692e-01*(u3*eta) + 4.83563087e-01*(u2*a12) + -5.41612504e-02*(u3*u) + 1.69778592e+00*(eta2*a12) + -7.12358090e+00*(u2*eta2*a1) + 4.47778938e+00*(u*eta2*a12) + 2.39033303e-01*(u4*a1) + -1.10544278e+00*(eta2*a13) + 6.36183705e-02*(u3*a12) + 1.47009000e+00*(u*eta*a13) + -5.26617393e+00*(u2*eta*a12) + 1.68311763e+00*(u3*eta*a1) + -2.98103807e-01*(u2*a13) + 5.31308951e-01*(u4*eta) + 1.00228016e+00*(u3*eta2) + -8.49897047e+00*(u3*eta2*a1) + -1.07620061e+00*(u4*eta2) + -3.87627923e-01*(u4*a12) + -1.74883077e+00*(u4*eta*a1) + 3.34572844e+00*(u2*eta*a13) + 1.41134786e+01*(u2*eta2*a12) + -3.54554972e+00*(u3*eta*a12) + -3.43159558e+00*(u*eta2*a13) + 2.46806714e-01*(u4*a13) + 2.09414054e+00*(u4*eta2*a1) + 1.80499409e+00*(u3*eta*a13) + 1.76720329e+01*(u3*eta2*a12) + -9.43276155e+00*(u2*eta2*a13) + 2.04986442e+00*(u4*eta*a12) + -1.26877754e+00*(u4*eta*a13) + -9.93436864e+00*(u3*eta2*a13);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double zeta2;

	// Evaluate fit for this parameter
	zeta2 = -1.38763654e+02 + 2.31329235e+03*eta + 1.18782378e+03*a1 + -5.91297318e+02*u + -8.35640207e+03*eta2 + -1.96785333e+04*(eta*a1) + -3.13958389e+03*a12 + 6.82576016e+03*(u*eta) + -8.46854527e+02*(u*u) + 3.62560714e+03*(u*a1) + 4.85576913e+04*(eta*a12) + -7.05976457e+03*(u*a12) + -1.69594385e+04*(u*eta2) + -4.10142862e+04*(u*eta*a1) + 6.92799573e+04*(eta2*a1) + 2.33064048e+03*(a12*a1) + 9.01738239e+03*(u2*eta) + 7.95842846e+02*(u2*u) + 2.38899116e+03*(u2*a1) + 7.97577367e+04*(u*eta*a12) + 1.03617626e+05*(u*eta2*a1) + -4.54257045e+03*(u3*a1) + 4.43479869e+03*(u*a13) + -3.45630035e+04*(eta*a13) + -2.31182354e+04*(u2*eta2) + -2.09031799e+04*(u2*eta*a1) + -8.74829214e+03*(u3*eta) + 1.94763541e+03*(u3*u) + -1.63531954e+05*(eta2*a12) + 5.21723868e+04*(u2*eta2*a1) + -2.06384787e+05*(u*eta2*a12) + -6.95531291e+03*(u4*a1) + 1.13143575e+05*(eta2*a13) + 8.82829787e+03*(u3*a12) + -5.06532268e+04*(u*eta*a13) + -1.52246916e+04*(u2*eta*a12) + 4.39153355e+04*(u3*eta*a1) + -2.01324004e+03*(u2*a13) + -2.42237327e+04*(u4*eta) + 2.00241085e+04*(u3*eta2) + -9.26313473e+04*(u3*eta2*a1) + 7.08091077e+04*(u4*eta2) + 6.49779672e+03*(u4*a12) + 7.85033425e+04*(u4*eta*a1) + 3.27073158e+04*(u2*eta*a13) + 4.14033836e+04*(u2*eta2*a12) + -7.95317751e+04*(u3*eta*a12) + 1.34954398e+05*(u*eta2*a13) + -5.83747483e+03*(u3*a13) + -1.39590277e+03*(u4*a13) + -2.27383521e+05*(u4*eta2*a1) + 5.26428846e+04*(u3*eta*a13) + 1.61992368e+05*(u3*eta2*a12) + -8.28407731e+04*(u2*eta2*a13) + -5.63508173e+04*(u4*eta*a12) + -1.11854150e+05*(u3*eta2*a13) + 1.61559087e+05*(u4*eta2*a12);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double nu0;

	// Evaluate fit for this parameter
	nu0 = 4.41851606e+02 + -6.68134748e+03*eta + 3.33257264e+04*eta2 + 2.37020668e+04*(eta*a1) + 2.31019185e+04*(u*eta) + 1.99964872e+04*(u*u) + -2.78992227e+04*(eta*a12) + -2.80590960e+03*(u*a12) + -1.80677505e+05*(u*eta2) + -1.32765158e+05*(u*eta*a1) + -1.39305533e+05*(eta2*a1) + 4.61027216e+02*(a12*a1) + -2.31194326e+05*(u2*eta) + 2.43779955e+04*(u2*u) + -1.25139795e+05*(u2*a1) + 2.58296201e+05*(u*eta*a12) + 9.80216629e+05*(u*eta2*a1) + -1.85654979e+05*(u3*a1) + 6.40768175e+05*(u2*eta2) + 1.51807603e+06*(u2*eta*a1) + -3.86323461e+05*(u3*eta) + 2.45394531e+05*(u2*a12) + -2.33037785e+04*(u3*u) + 1.48764385e+05*(eta2*a12) + -4.51948421e+06*(u2*eta2*a1) + -1.71249905e+06*(u*eta2*a12) + 1.05552917e+05*(u4*a1) + 4.09022566e+05*(u3*a12) + -1.12139527e+05*(u*eta*a13) + -3.09745338e+06*(u2*eta*a12) + 2.84438715e+06*(u3*eta*a1) + -1.52933937e+05*(u2*a13) + 2.41653152e+05*(u4*eta) + 1.44437044e+06*(u3*eta2) + -1.00993560e+07*(u3*eta2*a1) + -5.51290815e+05*(u4*eta2) + -1.45482110e+05*(u4*a12) + -9.48761885e+05*(u4*eta*a1) + 1.99686006e+06*(u2*eta*a13) + 9.67794515e+06*(u2*eta2*a12) + -6.14178070e+06*(u3*eta*a12) + 8.11896020e+05*(u*eta2*a13) + -2.62172718e+05*(u3*a13) + 6.43564849e+04*(u4*a13) + 1.96111791e+06*(u4*eta2*a1) + 3.90510933e+06*(u3*eta*a13) + 2.11661962e+07*(u3*eta2*a12) + -6.44367064e+06*(u2*eta2*a13) + 1.01281009e+06*(u4*eta*a12) + -2.86409142e+05*(u4*eta*a13) + -1.32632334e+07*(u3*eta2*a13) + -1.51976047e+06*(u4*eta2*a12);

	// Return answer
	return nu0;

} // END of NU0 (3,3) fit implementation


#ifdef __cplusplus
}
#endif
