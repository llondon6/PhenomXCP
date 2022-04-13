
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
// Header formatting for MU1
double IMRPhenomXCP_MU1( double theta, double eta, double a1 ){ 

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
	mu1 = -1.80179173e-01*u + 1.89554473e+00*a1 + -1.65265990e-01 + 1.01452835e+00*eta + 1.72461188e+00*(u*eta) + -2.49558718e+01*(eta*a1) + 1.65835935e+00*(u*a1) + -1.08263404e+00*(u*u) + -3.74459177e+00*a12 + 1.18464862e+02*(eta2*a1) + -2.43757329e+00*(u*a12) + 5.75928200e+01*(eta*a12) + -3.61965375e+01*(u*eta*a1) + 2.39534385e+01*(u2*eta) + -3.58107003e-01*(u2*u) + -2.07512457e+02*(eta3*a1) + 2.01953140e+02*(u*eta2*a1) + 2.65775905e+00*(u2*a12) + -1.40656182e+02*(u2*eta2) + -2.24269312e+01*(u2*eta*a1) + 2.78412852e+00*(u3*a1) + 5.75062627e+01*(u*eta*a12) + -2.35805812e+01*(u*eta2*eta) + -2.93051712e+02*(eta2*a12) + 5.00554091e+02*(eta3*a12) + -2.89384932e+02*(u*eta3*a1) + 2.10020046e+01*(u3*eta2) + 2.59579998e+02*(u2*eta2*eta) + -2.29356284e+01*(u2*eta*a12) + -2.36590205e+00*(u3*a12) + 1.71256132e+02*(u2*eta2*a1) + -1.96292483e+01*(u3*eta*a1) + -3.36308756e+02*(u*eta2*a12) + 5.56248818e+01*(u2*eta2*a12) + -5.44316335e+01*(u3*eta2*eta) + 5.37103619e+02*(u*eta3*a12) + -3.60413323e+02*(u2*eta3*a1) + 1.31943434e+01*(u3*eta*a12) + 2.05755395e+01*(u3*eta2*a1);

	// Return answer
	return mu1;

} // END of MU1 fit implementation



// MU3 fit implementation 
// Header formatting for MU3
double IMRPhenomXCP_MU3( double theta, double eta, double a1 ){ 

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
	mu3 = -3.20838086e-01*u + -2.23830563e+00*a1 + -3.64480327e-01 + 8.73539487e+00*eta + 2.33691645e+01*(eta*a1) + -8.25908591e+01*eta2 + 2.03876701e+01*(u*a1) + -2.32286486e+00*(u*u) + 2.27670817e+00*a12 + -3.36325032e+01*(u*a12) + -3.91022516e+01*(eta*a12) + -3.19487680e+02*(u*eta*a1) + 2.34609604e+02*(eta2*eta) + 1.74028620e+01*(u2*a1) + -5.63567183e+00*(u2*u) + -3.03129808e+02*(eta3*a1) + 1.76491589e+03*(u*eta2*a1) + 1.07060912e+02*(u3*eta) + -2.31925126e+01*(u2*a12) + 2.10735712e+02*(u2*eta2) + -1.15461378e+02*(u2*eta*a1) + -2.80374329e+00*(u3*a1) + 5.10194295e+02*(u*eta*a12) + 3.51311477e+01*(u*eta2*eta) + 1.33152724e+02*(eta2*a12) + -3.29056525e+03*(u*eta3*a1) + -5.53583396e+02*(u3*eta2) + -7.12886244e+02*(u2*eta2*eta) + 2.53678867e+02*(u2*eta*a12) + 1.11502140e+01*(u3*a12) + -2.93994125e+02*(u2*eta2*a1) + -3.93259233e+01*(u3*eta*a1) + -2.68316158e+03*(u*eta2*a12) + -6.61713767e+02*(u2*eta2*a12) + 8.09348456e+02*(u3*eta2*eta) + 4.73761943e+03*(u*eta3*a12) + 2.00636048e+03*(u2*eta3*a1) + -5.16081314e+01*(u3*eta*a12) + 2.46122082e+02*(u3*eta2*a1);

	// Return answer
	return mu3;

} // END of MU3 fit implementation



// MU4 fit implementation 
// Header formatting for MU4
double IMRPhenomXCP_MU4( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double delta = sqrt(1-4*eta);
	double delta2 = delta*delta;
	UNUSED double delta3 = delta2*delta;
	double mu4;

	// Evaluate fit for this parameter
	mu4 = 7.93754887e-04*u + 7.83328182e-03*a1 + -3.01041506e-03 + 6.87801146e-03*(u*delta) + -7.90856204e-03*(delta*a1) + -3.98495536e-03*(u*a1) + -3.53667106e-03*(u*u) + -5.45418403e-03*a12 + 1.86649581e-02*(delta*delta) + 4.41148885e-02*(delta*a12) + -7.10312941e-02*(delta2*a1) + -2.84195823e-02*(delta2*delta) + -1.22893575e-02*(u*delta*a1) + -1.38838211e-02*(u2*delta) + 5.73853805e-03*(u2*a1) + 2.98606646e-03*(u*a12) + 1.13058273e-02*(u*delta2*delta) + 1.43603356e-01*(delta3*a1) + -2.78864018e-03*(u2*a12) + -3.11857315e-02*(u3*delta) + 6.83390301e-02*(u2*delta*a1) + -9.98893224e-02*(u*delta2*a1) + -4.45900761e-02*(delta2*a12) + 1.36149051e-01*(u*delta2*a12) + -2.69163965e-02*(delta3*a12) + 1.18375757e-01*(u3*delta*a1) + 2.60356714e-02*(u3*delta*delta) + 2.58078380e-02*(u2*delta2*delta) + -1.08895656e-01*(u2*delta*a12) + 2.92872415e-02*(u*delta3*a1) + -9.06022044e-02*(u3*delta*a12) + 1.52551970e-01*(u2*delta2*a12) + -1.42455729e-01*(u2*delta3*a1) + -6.45524121e-02*(u3*delta2*a1) + -2.58512474e-02*(u*delta3*a12);

	// Return answer
	return mu4;

} // END of MU4 fit implementation



// NU4 fit implementation 
// Header formatting for NU4
double IMRPhenomXCP_NU4( double theta, double eta, double a1 ){ 

	/*
	Hola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.
	*/  

	// Preliminaries
	double u = cos(theta);
	double u2 = u*u;
	double u3 = u2*u;
	double a12 = a1*a1;
	double delta = sqrt(1-4*eta);
	double delta2 = delta*delta;
	UNUSED double delta3 = delta2*delta;
	double nu4;

	// Evaluate fit for this parameter
	nu4 = 3.78109197e-04*u + 3.76974105e-03*a1 + -1.56523110e-03 + -4.58086629e-03*(delta) + 5.14060530e-03*(u*delta) + 2.23157694e-02*(delta*a1) + -2.09233495e-03*(u*a1) + -1.85020469e-03*(u*u) + -2.71565589e-03*a12 + 1.31916954e-02*(delta*delta) + -2.80912431e-02*(delta*a12) + -5.07588060e-02*(delta2*a1) + -4.53644736e-03*(delta2*delta) + -3.91670402e-02*(u*delta*a1) + -2.63277785e-02*(u*delta*delta) + 7.85114257e-03*(u2*delta) + 2.84133658e-03*(u2*a1) + 2.80185978e-03*(u*a12) + 2.43287776e-02*(u*delta2*delta) + 3.69671245e-02*(u*delta*a12) + 8.24917382e-03*(delta3*a1) + -1.37984132e-02*(u2*delta*delta) + -1.26246289e-03*(u2*a12) + 5.74618982e-03*(u3*delta) + -3.92255186e-02*(u2*delta*a1) + 1.65595815e-03*(u3*a1) + 1.79467319e-01*(u*delta2*a1) + 5.20677398e-02*(delta2*a12) + -1.81287283e-01*(u*delta2*a12) + -1.74448542e-02*(u3*delta*a1) + -1.26950283e-02*(u3*delta*delta) + 7.14321678e-02*(u2*delta2*a1) + 5.32474395e-02*(u2*delta*a12) + -1.52724986e-01*(u*delta3*a1) + -2.89018251e-03*(u3*a12) + 6.22192610e-03*(u3*delta2*delta) + 1.86876291e-02*(u3*delta*a12) + -1.01105448e-01*(u2*delta2*a12) + 1.02732867e-02*(u2*delta3*a1) + 1.29823150e-02*(u3*delta2*a1) + 1.40586591e-01*(u*delta3*a12);

	// Return answer
	return nu4;

} // END of NU4 fit implementation



// NU5 fit implementation 
// Header formatting for NU5
double IMRPhenomXCP_NU5( double theta, double eta, double a1 ){ 

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
	nu5 = -1.22315452e-02*u + -1.39438234e-01*a1 + 1.84189403e-02 + -3.53876211e-01*eta + 2.13382834e+00*(eta*a1) + 2.06758387e+00*eta2 + -2.64841405e-02*(u*a1) + -2.73022262e-02*(u*u) + 7.94679758e-02*a12 + -1.17269816e+01*(eta2*a1) + -1.45037791e+00*(eta*a12) + 2.02246712e+00*(u*eta*a1) + 6.09265523e-01*(u*eta2) + 2.95911631e-01*(u2*eta) + -3.72973420e+00*(eta2*eta) + 1.57640988e-01*(u2*a1) + 2.68603153e-02*(u2*u) + 2.14172660e+01*(eta3*a1) + -1.59085311e+01*(u*eta2*a1) + -3.77327481e-01*(u3*eta) + -1.14945563e-01*(u2*a12) + -7.60681123e-01*(u2*eta2) + -1.76491871e+00*(u2*eta*a1) + -4.88347899e-02*(u3*a1) + -2.09462385e+00*(u*eta*a12) + -1.49441369e+00*(u*eta2*eta) + 7.59715228e+00*(eta2*a12) + -1.22054960e+01*(eta3*a12) + 3.24155208e+01*(u*eta3*a1) + 2.31795764e+00*(u3*eta2) + 1.15242438e+00*(u2*eta*a12) + 1.06512459e-01*(u3*a12) + 6.38379352e+00*(u2*eta2*a1) + 1.73359076e+01*(u*eta2*a12) + -2.80242653e+00*(u2*eta2*a12) + -5.49024499e+00*(u3*eta2*eta) + -3.53830406e+01*(u*eta3*a12) + -7.09928004e+00*(u2*eta3*a1) + -5.18019835e-01*(u3*eta*a12) + 1.24390880e+00*(u3*eta2*a1);

	// Return answer
	return nu5;

} // END of NU5 fit implementation



// NU6 fit implementation 
// Header formatting for NU6
double IMRPhenomXCP_NU6( double theta, double eta, double a1 ){ 

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
	nu6 = 1.39128045e-02*u + -1.72933273e-02*a1 + 2.46097164e-03 + -1.41820974e-02*eta + -5.84392465e-02*(u*eta) + 2.64351146e-01*(eta*a1) + -2.42691485e-01*(u*a1) + -7.90314669e-03*(u*u) + 9.59731614e-02*a12 + -1.05266799e+00*(eta2*a1) + 3.27150897e-01*(u*a12) + -1.48271224e+00*(eta*a12) + 3.43487281e+00*(u*eta*a1) + 9.46613106e-02*(u2*eta) + 4.26255512e-02*(u2*u) + 1.28251025e+00*(eta3*a1) + -1.77097865e+01*(u*eta2*a1) + -8.21529808e-01*(u3*eta) + -4.86141982e-01*(u2*eta2) + 1.55026950e-02*(u3*a1) + -4.66371009e+00*(u*eta*a12) + 7.30760801e+00*(eta2*a12) + -1.17247573e+01*(eta3*a12) + 3.15081286e+01*(u*eta3*a1) + 4.21822307e+00*(u3*eta2) + 8.52070157e-01*(u2*eta2*eta) + -1.54044517e-02*(u2*eta*a12) + -8.53548375e-02*(u3*a12) + 3.30227979e-01*(u2*eta2*a1) + 3.83864631e-01*(u3*eta*a1) + 2.33955972e+01*(u*eta2*a12) + -6.39658601e+00*(u3*eta2*eta) + -3.99023856e+01*(u*eta3*a12) + -1.15096556e+00*(u2*eta3*a1) + 3.43400055e-01*(u3*eta*a12) + -1.83929849e+00*(u3*eta2*a1);

	// Return answer
	return nu6;

} // END of NU6 fit implementation



// ZETA1 fit implementation 
// Header formatting for ZETA1
double IMRPhenomXCP_ZETA1( double theta, double eta, double a1 ){ 

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
	zeta1 = 1.67743416e-05*u + 6.71180503e-05*a1 + -4.44179970e-06 + 1.90620006e-05*eta + -1.59647404e-04*(u*eta) + -6.69108793e-04*(eta*a1) + 1.45927570e-04*(u*a1) + 2.70496670e-05*(u*u) + -4.94893751e-05*a12 + 2.66118826e-03*(eta2*a1) + -1.75545401e-04*(u*a12) + 3.17785805e-04*(eta*a12) + -3.58931535e-03*(u*eta*a1) + -2.24523040e-04*(u2*eta) + -1.87012999e-04*(u2*a1) + -8.87868729e-05*(u2*u) + -4.28534040e-03*(eta3*a1) + 2.51909156e-02*(u*eta2*a1) + 1.40416692e-03*(u3*eta) + 1.33676521e-04*(u2*a12) + 4.37032513e-04*(u2*eta2) + 2.05694774e-03*(u2*eta*a1) + 1.05676172e-04*(u3*a1) + 4.38376768e-03*(u*eta*a12) + 1.15222335e-03*(u*eta2*eta) + -1.86264443e-03*(eta3*a12) + -5.08952754e-02*(u*eta3*a1) + -6.81118451e-03*(u3*eta2) + -1.00984833e-03*(u2*eta*a12) + -7.67194398e-05*(u3*a12) + -8.35718441e-03*(u2*eta2*a1) + -1.09724008e-03*(u3*eta*a1) + -2.93696439e-02*(u*eta2*a12) + 1.83903177e-03*(u2*eta2*a12) + 1.13624784e-02*(u3*eta2*eta) + 5.67267029e-02*(u*eta3*a12) + 1.28144869e-02*(u2*eta3*a1) + 5.79148141e-04*(u3*eta*a12) + 1.59600797e-03*(u3*eta2*a1);

	// Return answer
	return zeta1;

} // END of ZETA1 fit implementation



// ZETA2 fit implementation 
// Header formatting for ZETA2
double IMRPhenomXCP_ZETA2( double theta, double eta, double a1 ){ 

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
	zeta2 = -3.10978687e+00*u + -9.99071747e+00*a1 + 1.57569134e+00 + 2.42202533e+01*(u*eta) + -7.52749200e+01*eta2 + -1.54054426e+01*(u*a1) + -1.32418745e+01*(u*u) + 6.07481534e+00*a12 + 5.94064420e+02*(eta2*a1) + 1.76707995e+01*(u*a12) + 8.74215049e+01*(eta*a12) + 4.71396662e+02*(u*eta*a1) + 1.94370276e+02*(u2*eta) + 2.18064853e+02*(eta2*eta) + 4.90475448e+01*(u2*a1) + 1.56056175e+01*(u2*u) + -1.77088226e+03*(eta3*a1) + -3.33987339e+03*(u*eta2*a1) + -2.19190190e+02*(u3*eta) + -2.97866125e+01*(u2*a12) + -9.10391228e+02*(u2*eta2) + -6.39677252e+02*(u2*eta*a1) + -2.56051271e+01*(u3*a1) + -5.26945210e+02*(u*eta*a12) + -1.46398474e+02*(u*eta2*eta) + -1.22982593e+03*(eta2*a12) + 3.15875668e+03*(eta3*a12) + 6.55501227e+03*(u*eta3*a1) + 9.79214949e+02*(u3*eta2) + 1.45848989e+03*(u2*eta2*eta) + 2.86235986e+02*(u2*eta*a12) + 1.94929788e+01*(u3*a12) + 2.70696734e+03*(u2*eta2*a1) + 2.03601441e+02*(u3*eta*a1) + 3.58851537e+03*(u*eta2*a12) + -6.18592204e+02*(u2*eta2*a12) + -1.50666827e+03*(u3*eta2*eta) + -6.77288453e+03*(u*eta3*a12) + -3.98799540e+03*(u2*eta3*a1) + -1.12422366e+02*(u3*eta*a12) + -2.75019418e+02*(u3*eta2*a1);

	// Return answer
	return zeta2;

} // END of ZETA2 fit implementation


#ifdef __cplusplus
}
#endif
