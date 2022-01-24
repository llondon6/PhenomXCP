
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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double mu1;

	// Evaluate fit for this parameter
	mu1 = -4.47904325e-02 + -3.92628579e+00*a1 + 6.13379166e-01*u + 8.24606449e+01*(eta*a1) + -6.78189639e+00*(u*eta) + -3.14013912e+00*(u*a1) + 1.60793743e+01*eta2 + 1.78871587e+01*a12 + 2.18562855e+00*(u*u) + -1.89587022e+01*(a12*a1) + -6.48906103e+01*(u2*eta) + 7.32350790e+00*(u*a12) + 5.12758004e+01*(u*eta2) + -4.98740967e+01*(eta2*eta) + -6.38199416e+02*(eta2*a1) + -3.49166131e+02*(eta*a12) + -9.43212213e+00*(u2*a1) + -8.20736840e-01*(u2*u) + 5.28975311e+00*(u3*a1) + 2.41325729e+03*(eta2*a12) + 3.19605372e+02*(u2*eta*a1) + -1.27833620e+02*(u*eta2*eta) + -6.17009553e+00*(u*a13) + 3.58029781e+02*(eta*a13) + 4.38158534e+02*(u2*eta2) + 1.42572397e+03*(eta3*a1) + 2.05352788e+00*(u3*eta) + 2.89280367e+01*(u*eta*a13) + 1.20239193e+02*(u*eta3*a1) + -9.42766418e+00*(u3*eta*a1) + -3.51570305e+02*(u2*eta*a12) + -1.95958327e+03*(u2*eta2*a1) + -2.33156258e+03*(eta2*a13) + 2.20397622e+01*(u4*eta) + 1.49010381e+01*(u2*a13) + -8.91720077e+02*(u2*eta2*eta) + -5.65774745e+01*(u*eta2*a12) + -5.08942948e+03*(eta3*a12) + -6.51961539e+00*(u4*a1) + -8.84834796e+00*(u3*a12) + 2.16179474e+03*(u2*eta2*a12) + -1.72427353e+02*(u4*eta2) + 4.99694255e+00*(u3*a13) + 4.73320684e+03*(eta3*a13) + -1.34288315e+01*(u2*eta*a13) + 3.70640753e+03*(u2*eta3*a1) + 1.62658781e+01*(u4*a12) + -1.11248301e+02*(u*eta2*a13) + 1.99216501e+00*(u4*eta*a12) + -2.94452655e+01*(u3*eta3*a1) + 3.83939435e+02*(u4*eta2*eta) + -4.01792957e+03*(u2*eta3*a12) + 1.84419670e+02*(u*eta3*a13) + 3.70355337e+01*(u3*eta2*a12) + -1.24289937e+01*(u4*a13);

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
	double a13 = a12*a1;
	double delta = sqrt(1-4*eta);
	double delta2 = delta*delta;
	UNUSED double delta3 = delta2*delta;
	double mu3;

	// Evaluate fit for this parameter
	mu3 = 1.80971808e-01 + 4.39118863e-01*(delta) + -4.54522952e-01*a1 + 2.98236052e-01*(u*a1) + -2.33349979e+00*(delta*a1) + 5.83772286e-01*(u*delta) + -7.04444952e-01*(delta*delta) + 2.22351979e-01*(a12*a1) + -3.34473856e+00*(u*delta*a1) + 2.84885302e-01*(u2*a1) + -2.36038249e-01*(u2*u) + 1.14313833e+01*(delta2*a12) + -5.27488061e-01*(u*delta2*delta) + -7.00534891e-01*(u*a13) + 2.35223681e+00*(delta*a13) + 9.27099881e+00*(u*delta2*a1) + -5.65759998e-01*(u2*delta*delta) + -1.12703817e+01*(delta2*a13) + 2.88883382e+00*(u2*delta2*a1) + -2.52040391e+00*(delta3*a12) + -2.15919799e+00*(u2*delta*a12) + 8.63575673e-01*(u3*a12) + -8.59781200e+00*(u*delta3*a1);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double delta = sqrt(1-4*eta);
	double delta2 = delta*delta;
	UNUSED double delta3 = delta2*delta;
	double mu4;

	// Evaluate fit for this parameter
	mu4 = -4.61569747e-03 + 2.13564486e-02*a1 + -6.80784050e-04*(u*a1) + -3.83295985e-02*a12 + -4.32850965e-03*(u*u) + 2.34576766e-02*(a12*a1) + 3.00925355e-02*(delta*a12) + 8.58018031e-03*(u2*a1) + -1.59607508e-02*(u2*delta*a1) + -1.22328575e-01*(delta2*a12) + 8.34271220e-03*(delta3*a1) + -1.77355689e-02*(u3*delta) + 2.91639938e-02*(delta2*a13) + 9.12180538e-02*(u2*delta2*a1) + 1.50329260e-01*(delta3*a12) + -1.49935265e-01*(u*delta2*a12) + -7.15261822e-03*(u2*a13) + -1.41989846e-02*(u*delta*a13) + 9.58795658e-02*(u3*delta*a1) + 7.34693632e-02*(u*delta3*a1) + -9.25397565e-02*(u2*delta3*a1) + -8.83665990e-02*(u*delta3*a12) + -1.21898775e-02*(u4*delta*a1) + -7.74940971e-02*(u3*delta*a12) + -5.07013739e-02*(u3*delta2*a1) + 6.27694277e-03*(u2*delta2*a12) + 1.28345540e-02*(u3*delta2*delta) + -6.82500858e-02*(delta3*a13) + 2.36973240e-01*(u*delta2*a13);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double delta = sqrt(1-4*eta);
	double delta2 = delta*delta;
	UNUSED double delta3 = delta2*delta;
	double nu4;

	// Evaluate fit for this parameter
	nu4 = -1.72164470e-03 + 3.48917607e-03*a1 + 1.32431394e-04*u + 2.14368467e-03*(delta*delta) + -2.81424339e-03*(a12*a1) + -4.76943879e-03*(delta*a12) + 8.78400971e-04*(u2*delta) + -2.27920402e-03*(u2*a1) + -5.02149264e-03*(u3*u) + -8.50426617e-03*(u*delta*a12) + -1.21391281e-03*(u*delta2*delta) + 9.02753824e-03*(u*delta2*a1) + -2.54559021e-02*(u2*delta2*a1) + -2.61796792e-02*(delta3*a12) + -8.87735224e-03*(u*delta2*a12) + 7.56937264e-03*(u4*delta) + 3.70496612e-03*(u2*a13) + 6.90147054e-03*(u2*delta*a12) + 7.47176963e-03*(u*delta*a13) + 2.89239474e-02*(u4*a1) + -4.64331583e-02*(u4*delta*a1) + -9.09602290e-03*(u3*delta2*a1) + 1.28459032e-01*(u2*delta2*a12) + 4.26366031e-02*(delta3*a13) + -5.36262829e-02*(u4*a12) + 5.42803165e-02*(u4*delta*a12) + -1.57937552e-01*(u2*delta2*a13) + 2.15082446e-02*(u3*delta2*a12) + 9.73587750e-03*(u2*delta3*a12) + -1.54482079e-02*(u*delta3*a13) + 3.08138593e-02*(u4*a13);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu5;

	// Evaluate fit for this parameter
	nu5 = 7.24466124e-02 + -1.72918658e+00*eta + -2.01172081e-01*a1 + -3.86925523e-01*u + 5.87117015e+00*(eta*a1) + 7.29716455e+00*(u*eta) + 2.49634444e+00*(u*a1) + 1.24864590e+01*eta2 + 1.39511942e-01*(u*u) + 1.07204154e-01*(a12*a1) + -4.89120234e+00*(u*a12) + -4.31728840e+01*(u*eta2) + -4.42342213e+01*(u*eta*a1) + -2.74324608e+01*(eta2*eta) + -4.90649920e+01*(eta2*a1) + -3.84861666e+00*(eta*a12) + -1.97163812e+00*(u2*a1) + 8.19093965e-01*(u3*a1) + 4.95512720e+01*(eta2*a12) + 1.70044660e+01*(u2*eta*a1) + 1.33838403e-02*(u3*u) + 8.32641328e+01*(u*eta*a12) + 8.29744792e+01*(u*eta2*eta) + 3.27406918e+00*(u*a13) + 2.41352000e+02*(u*eta2*a1) + -1.19999478e+01*(u2*eta2) + 1.18352560e+02*(eta3*a1) + 4.96951117e+00*(u2*a12) + -2.23089166e+00*(u3*eta) + -5.34309238e+01*(u*eta*a13) + -4.26348747e+02*(u*eta3*a1) + 2.08317810e+01*(u3*eta2) + -5.39209693e+01*(u2*eta*a12) + -1.45789010e+01*(eta2*a13) + -3.25967365e+00*(u2*a13) + 3.48513149e+01*(u2*eta2*eta) + -4.25870363e+02*(u*eta2*a12) + -1.42802470e+02*(eta3*a12) + -2.66233463e+00*(u3*a12) + 6.94074252e+02*(u*eta3*a12) + 1.68861786e+01*(u3*eta*a12) + 1.13377308e+02*(u2*eta2*a12) + 1.28324640e+00*(u3*a13) + 5.51506281e+01*(eta3*a13) + -4.12426408e+01*(u3*eta2*a1) + 3.69533800e+01*(u2*eta*a13) + -4.99636792e+01*(u3*eta2*eta) + -1.27428910e+02*(u2*eta3*a1) + 5.32273556e-01*(u4*a12) + -2.07511698e+00*(u4*eta*a1) + 2.58649283e+02*(u*eta2*a13) + -5.17536845e+00*(u3*eta*a13) + 1.21435229e+00*(u4*eta*a12) + 1.20428282e+02*(u3*eta3*a1) + -9.24138218e+01*(u2*eta2*a13) + 6.01824825e+00*(u4*eta2*eta) + 6.43165716e+01*(u2*eta3*a12) + -3.91952959e+02*(u*eta3*a13) + -2.63857273e+01*(u3*eta2*a12) + -4.43642707e-01*(u4*a13);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu6;

	// Evaluate fit for this parameter
	nu6 = -1.28192194e-02 + 2.49723423e-01*eta + 6.02295716e-02*a1 + -4.51681300e-02*u + -9.45039233e-01*(eta*a1) + 4.49620517e-01*(u*a1) + -1.68066378e+00*eta2 + 2.61369667e-02*(a12*a1) + -1.42499147e+00*(u*a12) + 2.36119331e+00*(u*eta2) + -1.07947522e+00*(u*eta*a1) + 3.51025512e+00*(eta2*eta) + 6.60052125e+00*(eta2*a1) + -3.10613562e-01*(eta*a12) + 4.53639235e-02*(u2*u) + -2.75992533e-01*(u3*a1) + -8.56547554e-03*(u3*u) + 8.45825758e+00*(u*eta*a12) + -6.71863947e+00*(u*eta2*eta) + 1.24467285e+00*(u*a13) + -1.33749703e+01*(u*eta2*a1) + -1.49957863e+01*(eta3*a1) + -1.80969264e-01*(u3*eta) + -9.72960648e+00*(u*eta*a13) + 4.25071122e+01*(u*eta3*a1) + 7.73753947e-02*(u4*eta) + -7.66506240e-02*(u2*eta2*eta) + 4.93601096e+00*(eta3*a12) + 8.17522669e-01*(u3*a12) + -4.48501841e+01*(u*eta3*a12) + -2.21130148e+00*(u3*eta*a12) + -1.74861030e-01*(u4*eta2) + -7.00087682e-01*(u3*a13) + -1.69344500e+00*(eta3*a13) + 6.10321478e+00*(u3*eta2*a1) + 1.91330060e+01*(u*eta2*a13) + 2.88026152e+00*(u3*eta*a13) + -6.01211521e+00*(u3*eta3*a1) + -4.76988290e+00*(u3*eta2*a12);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double zeta1;

	// Evaluate fit for this parameter
	zeta1 = 1.49586477e-06 + -5.51813953e-05*(u*eta) + 9.15449900e-05*a12 + -9.67201832e-06*(u*u) + -7.16860619e-05*(a12*a1) + 4.31228984e-04*(u*a12) + 3.72189254e-04*(u*eta*a1) + -1.61690535e-04*(eta2*eta) + -6.93084738e-04*(eta2*a1) + -9.89397701e-04*(eta*a12) + -3.73608484e-05*(u2*u) + 5.62372125e-03*(eta2*a12) + 4.85813230e-04*(u2*eta*a1) + 5.34808404e-05*(u3*u) + -8.82089526e-03*(u*eta*a12) + -5.31597449e-04*(u*a13) + 3.48231979e-04*(eta*a13) + 3.49549061e-03*(eta3*a1) + -1.90263693e-04*(u2*a12) + 5.83191892e-04*(u3*eta) + 1.07433465e-02*(u*eta*a13) + -3.00743376e-03*(u3*eta2) + 5.31879294e-04*(u3*eta*a1) + 1.82433210e-03*(u2*eta*a12) + -6.42404076e-04*(u4*eta) + 1.47380854e-04*(u2*a13) + 4.68629990e-02*(u*eta2*a12) + -1.41182169e-02*(eta3*a12) + -6.62820987e-05*(u4*a1) + -8.77585976e-05*(u3*a12) + -8.46094369e-02*(u*eta3*a12) + -1.40015074e-02*(u2*eta2*a12) + 2.56625301e-03*(u4*eta2) + 1.25133174e-04*(u3*a13) + -4.95753389e-03*(u3*eta2*a1) + -3.52771375e-04*(u2*eta*a13) + 6.92337615e-03*(u3*eta2*eta) + -4.81404049e-03*(u2*eta3*a1) + 4.74195119e-04*(u4*eta*a1) + -5.95212946e-02*(u*eta2*a13) + -1.07431362e-03*(u3*eta*a13) + -3.12736437e-03*(u4*eta2*eta) + 3.34592092e-02*(u2*eta3*a12) + -9.68145998e-04*(u4*eta2*a1) + 1.05854871e-01*(u*eta3*a13) + 6.11752972e-03*(u3*eta2*a12);

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
	UNUSED double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double zeta2;

	// Evaluate fit for this parameter
	zeta2 = -3.21970082e-01 + 7.85531992e+00*eta + -1.44916480e+01*u + -8.14133832e+00*(eta*a1) + 2.71313046e+02*(u*eta) + -1.89263361e+01*eta2 + -1.43019439e+01*a12 + 1.38411253e+00*(u*u) + 1.24362251e+01*(a12*a1) + 4.23738895e+01*(u*a12) + -1.50806641e+03*(u*eta2) + -4.39850071e+01*(u*eta*a1) + 1.05169145e+02*(eta*a12) + 1.29431216e+01*(u2*a1) + 3.76196721e+01*(u2*u) + -7.71320785e+01*(u3*a1) + -2.32712165e+02*(u2*eta*a1) + -1.38857203e+01*(u3*u) + -4.92628562e+02*(u*eta*a12) + 2.72327976e+03*(u*eta2*eta) + -1.88647113e+01*(u*a13) + -6.44830380e+02*(u3*eta) + 3.42647161e+03*(u3*eta2) + 1.20551451e+03*(u3*eta*a1) + 6.42770237e+02*(u2*eta2*a1) + -8.73549879e+02*(eta2*a13) + 1.84931931e+02*(u4*eta) + 2.56634075e+03*(u*eta2*a12) + -5.32609858e+02*(eta3*a12) + 1.72680092e+01*(u4*a1) + 3.21551549e+01*(u3*a12) + -3.77416197e+03*(u*eta3*a12) + -4.15978394e+02*(u3*eta*a12) + 1.24822494e+03*(u2*eta2*a12) + -7.65838482e+02*(u4*eta2) + 2.57140043e+03*(eta3*a13) + -5.51346442e+03*(u3*eta2*a1) + -9.78773914e+01*(u2*eta*a13) + -6.00728568e+03*(u3*eta2*eta) + -1.64122027e+02*(u4*eta*a1) + 7.51231031e+02*(u*eta2*a13) + 8.24579110e+01*(u3*eta*a13) + 6.31583087e+01*(u4*eta*a12) + 9.04979312e+03*(u3*eta3*a1) + 3.13895932e+02*(u2*eta2*a13) + 1.05475456e+03*(u4*eta2*eta) + -4.51957403e+03*(u2*eta3*a12) + 2.54365050e+02*(u4*eta2*a1) + -2.53090984e+03*(u*eta3*a13) + 5.40461479e+02*(u3*eta2*a12) + -8.39493042e+00*(u4*a13);

	// Return answer
	return zeta2;

} // END of ZETA2 fit implementation


#ifdef __cplusplus
}
#endif
