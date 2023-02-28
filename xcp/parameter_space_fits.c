
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
	mu1 = -2.28290120e+01*a1 + -4.33988384e+01*eta + 2.17935137e+00 + 1.04952699e+00*u + -9.34426006e+00*(u*a1) + 4.48553156e+02*(eta*a1) + 2.74802700e+02*eta2 + 5.93369552e+01*a12 + -2.80805870e+03*(eta2*a1) + -1.14301847e+03*(eta*a12) + 2.69511363e+01*(u2*a1) + -5.36480527e+02*(eta2*eta) + -4.56652563e+01*(a12*a1) + 3.05536771e+01*(u*a12) + -5.81977929e+01*(u*eta2) + -1.38212212e+00*(u2*u) + -2.74697380e+01*(u*a13) + -3.94092370e+01*(u2*eta2) + 5.12571687e+02*(u*eta2*a1) + -5.15057529e+02*(u2*eta*a1) + 5.47633023e+03*(eta3*a1) + -1.88580490e+02*(u*eta*a12) + 1.62413532e+02*(u*eta2*eta) + -1.14959026e+01*(u3*eta) + 1.48678953e+01*(u3*a1) + 8.64023501e+02*(eta*a13) + -9.31552639e+01*(u2*a12) + -6.72353600e+00*(u3*u) + 7.06865624e+03*(eta2*a12) + 1.14402675e+02*(u2*eta2*eta) + -1.44666649e+03*(u*eta3*a1) + -5.26472650e+03*(eta2*a13) + 7.76287119e+01*(u2*a13) + -4.21434542e+01*(u3*a12) + 3.46061528e+03*(u2*eta2*a1) + 1.88334070e+02*(u3*eta2) + 1.23191345e+02*(u4*eta) + -1.36923485e+04*(eta3*a12) + 1.91081043e+01*(u4*a1) + 2.75336134e+02*(u*eta*a13) + 1.64374319e+03*(u2*eta*a12) + -1.37594794e+01*(u4*a12) + -6.64259508e+02*(u4*eta2) + -1.01538923e+04*(u2*eta2*a12) + 1.00939925e+04*(eta3*a13) + -6.92519793e+03*(u2*eta3*a1) + -9.26621473e+02*(u*eta2*a13) + -8.96154990e+02*(u3*eta2*a1) + 1.10158889e+03*(u*eta3*a12) + 2.14301667e+02*(u3*eta*a12) + -3.02477199e+02*(u4*eta*a1) + -1.30417285e+03*(u2*eta*a13) + -4.75285404e+02*(u3*eta2*eta) + 3.49206427e+01*(u3*a13) + 5.70454324e+02*(u3*eta2*a12) + 1.92198841e+04*(u2*eta3*a12) + 2.38240303e+02*(u4*eta*a12) + 1.01026385e+03*(u*eta3*a13) + 1.12228457e+03*(u4*eta2*a1) + 1.17953537e+03*(u4*eta2*eta) + 7.66022017e+03*(u2*eta2*a13) + -3.02920012e+02*(u3*eta*a13) + 2.62089862e+03*(u3*eta3*a1) + -4.14721255e+02*(u4*eta2*a12) + -1.38365548e+03*(u4*eta3*a1) + -3.06562604e+03*(u3*eta3*a12) + -5.84026521e+01*(u4*eta*a13) + -1.39493710e+04*(u2*eta3*a13) + 6.71031388e+02*(u3*eta2*a13);

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
	mu2 = -1.78329318e+01*a1 + -2.02746273e+01*eta + 1.97740907e+00 + 9.62361325e+00*u + -3.49367711e+01*(u*a1) + 1.98991115e+02*(eta*a1) + -2.27194503e+02*(u*eta) + 6.18288194e+01*eta2 + 2.99287961e+01*a12 + -8.04605043e+00*(u*u) + -8.09419489e+02*(eta2*a1) + -2.62876472e+02*(eta*a12) + 4.65229521e+01*(u2*a1) + 9.69091703e+02*(u*eta*a1) + -4.24031416e+01*(eta2*eta) + -1.77363702e+01*(a12*a1) + 3.66781728e+01*(u*a12) + 1.69175665e+02*(u2*eta) + 1.75301998e+03*(u*eta2) + -8.77549866e+00*(u2*u) + -9.54442632e+02*(u2*eta2) + -8.75240367e+03*(u*eta2*a1) + -1.27547896e+03*(u2*eta*a1) + 1.13813484e+03*(eta3*a1) + -1.37193315e+03*(u*eta*a12) + -4.05627437e+03*(u*eta2*eta) + 2.43245867e+02*(u3*eta) + -2.53338868e+01*(u3*a1) + 1.21378972e+02*(eta*a13) + 3.30329085e+01*(u2*a12) + -8.24382725e+00*(u3*u) + 8.16117480e+02*(eta2*a12) + 1.57380128e+03*(u2*eta2*eta) + 2.21622313e+04*(u*eta3*a1) + 1.49313887e+04*(u*eta2*a12) + -2.27843536e+02*(eta2*a13) + -1.19182900e+02*(u2*a13) + 1.60862004e+02*(u3*a12) + 8.60666486e+03*(u2*eta2*a1) + -5.06547439e+01*(u3*eta*a1) + -2.18360227e+03*(u3*eta2) + 1.59050928e+02*(u4*eta) + -8.40841757e+02*(eta3*a12) + 1.22828179e+01*(u4*a1) + 4.29422797e+02*(u*eta*a13) + 6.18507604e+02*(u2*eta*a12) + -5.43424912e+01*(u4*a12) + -1.31709022e+03*(u4*eta2) + -8.30254719e+03*(u2*eta2*a12) + -1.65481080e+04*(u2*eta3*a1) + -7.16302555e+03*(u*eta2*a13) + 5.67943398e+03*(u3*eta2*a1) + -4.09855256e+04*(u*eta3*a12) + -1.97352337e+03*(u3*eta*a12) + 1.25076201e+02*(u4*eta*a1) + 1.21150665e+03*(u2*eta*a13) + 5.57189939e+03*(u3*eta2*eta) + -1.78744614e+02*(u3*a13) + 7.54192167e+01*(u4*a13) + 2.01115054e+04*(u2*eta3*a12) + -1.67154148e+02*(u4*eta*a12) + 2.23186221e+04*(u*eta3*a13) + 3.43736890e+03*(u4*eta2*eta) + -2.73405718e+03*(u2*eta2*a13) + 2.72600959e+03*(u3*eta*a13) + -2.13221083e+04*(u3*eta3*a1) + 2.10339752e+03*(u4*eta2*a12) + -3.93055059e+02*(u4*eta*a13) + -3.92673741e+03*(u4*eta3*a1) + 2.33216823e+04*(u3*eta3*a12) + -8.35564833e+03*(u3*eta2*a13);

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
	mu3 = -6.15504439e-02*a1 + 1.25508429e-02 + -3.54428996e-01*u + 2.85948012e+00*(u*a1) + -1.01562415e+00*(eta*a1) + 4.45659468e+00*(u*eta) + -1.21600106e+00*eta2 + 3.15197060e-01*a12 + -3.59220564e-01*(u*u) + 1.62068267e+01*(eta2*a1) + 2.10536912e-01*(eta*a12) + 1.74069984e+00*(u2*a1) + -3.41834016e+01*(u*eta*a1) + 3.79505800e+00*(eta2*eta) + -2.36164090e-01*(a12*a1) + -6.63566312e+00*(u*a12) + 4.77347912e+00*(u2*eta) + -1.65132584e+01*(u*eta2) + 3.96845632e-01*(u2*u) + 4.60278623e+00*(u*a13) + -2.26103731e+01*(u2*eta2) + 1.16758395e+02*(u*eta2*a1) + -1.87067862e+01*(u2*eta*a1) + -4.35493201e+01*(eta3*a1) + 7.73218783e+01*(u*eta*a12) + 1.69643517e+01*(u*eta2*eta) + -4.28988823e+00*(u3*eta) + -3.25856003e+00*(u3*a1) + -3.50535324e+00*(u2*a12) + 3.05309698e-01*(u3*u) + -2.81616998e+01*(eta2*a12) + 3.73336325e+01*(u2*eta2*eta) + -9.91660108e+01*(u*eta3*a1) + -2.52847578e+02*(u*eta2*a12) + 1.96641212e+01*(eta2*a13) + 2.45522353e+00*(u2*a13) + 7.81862211e+00*(u3*a12) + 7.25636838e+01*(u2*eta2*a1) + 3.43723268e+01*(u3*eta*a1) + 1.09894635e+01*(u3*eta2) + -3.98227889e+00*(u4*eta) + 8.73854209e+01*(eta3*a12) + -1.16149492e+00*(u4*a1) + -5.27325198e+01*(u*eta*a13) + 3.56377803e+01*(u2*eta*a12) + 1.76473881e+00*(u4*a12) + 1.94927424e+01*(u4*eta2) + -1.27506495e+02*(u2*eta2*a12) + -6.26224991e+01*(eta3*a13) + -1.04348709e+02*(u2*eta3*a1) + 1.67846185e+02*(u*eta2*a13) + -8.35250901e+01*(u3*eta2*a1) + 1.89475416e+02*(u*eta3*a12) + -8.23905645e+01*(u3*eta*a12) + 9.29044923e+00*(u4*eta*a1) + -2.55963630e+01*(u2*eta*a13) + -5.61845330e+00*(u3*a13) + 2.01832469e+02*(u3*eta2*a12) + -1.16136179e+00*(u4*a13) + 1.68218388e+02*(u2*eta3*a12) + -7.42381275e+00*(u4*eta*a12) + -1.15698702e+02*(u*eta3*a13) + -3.11087642e+01*(u4*eta2*a1) + -3.44682215e+01*(u4*eta2*eta) + 9.28936280e+01*(u2*eta2*a13) + 5.96346537e+01*(u3*eta*a13) + -1.35214512e+01*(u3*eta3*a1) + 4.83740230e+00*(u4*eta*a13) + 5.34352694e+01*(u4*eta3*a1) + 2.52700192e+01*(u3*eta3*a12) + -1.21431698e+02*(u2*eta3*a13) + -1.51185541e+02*(u3*eta2*a13);

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
	nu4 = -1.78377388e-01*a1 + -4.64888909e-01*eta + 1.80833391e-02 + -4.83609694e-02*u + 3.97235846e-01*(u*a1) + 4.37680784e+00*(eta*a1) + 8.63928826e-01*(u*eta) + 3.16247971e+00*eta2 + 4.11462902e-01*a12 + 7.33745267e-02*(u*u) + -2.93530138e+01*(eta2*a1) + -1.02891691e+01*(eta*a12) + -2.97941770e-01*(u2*a1) + -7.07705463e+00*(u*eta*a1) + -6.43597793e+00*(eta2*eta) + -2.79350656e-01*(a12*a1) + -9.63670444e-01*(u*a12) + -6.33133517e-01*(u2*eta) + -5.05777413e+00*(u*eta2) + 6.14000700e-03*(u2*u) + 6.64878258e-01*(u*a13) + 1.35093393e+00*(u2*eta2) + 4.16409446e+01*(u*eta2*a1) + 5.90890114e+01*(eta3*a1) + 1.73721986e+01*(u*eta*a12) + 9.60935198e+00*(u*eta2*eta) + -3.22510075e-02*(u3*eta) + -3.34068547e-02*(u3*a1) + 7.22654339e+00*(eta*a13) + 6.97149185e-01*(u2*a12) + -7.69966927e-02*(u3*u) + 6.94889214e+01*(eta2*a12) + -7.96533300e+01*(u*eta3*a1) + -1.03438298e+02*(u*eta2*a12) + -4.96111649e+01*(eta2*a13) + -5.35290234e-01*(u2*a13) + 1.62504317e+01*(u2*eta2*a1) + 8.73306970e-01*(u4*eta) + -1.40099718e+02*(eta3*a12) + 2.65413874e-01*(u4*a1) + -1.22702571e+01*(u*eta*a13) + 3.61322245e-01*(u2*eta*a12) + -5.34335284e-01*(u4*a12) + -3.65151454e+00*(u4*eta2) + -3.96511439e+01*(u2*eta2*a12) + 1.00926896e+02*(eta3*a13) + -4.59493345e+01*(u2*eta3*a1) + 8.13876572e-01*(u3*eta2*a1) + 7.39827527e+01*(u*eta2*a13) + 1.99651979e+02*(u*eta3*a12) + 9.75489388e-01*(u3*eta*a12) + -1.40584383e+00*(u4*eta*a1) + 5.27514418e-02*(u3*a13) + -4.59545524e+00*(u3*eta2*a12) + 3.92088649e-01*(u4*a13) + 1.08088573e+02*(u2*eta3*a12) + 2.11333163e+00*(u4*eta*a12) + -1.43667550e+02*(u*eta3*a13) + 1.45879759e+00*(u4*eta2*a1) + 5.49563540e+00*(u4*eta2*eta) + 2.83163602e+01*(u2*eta2*a13) + -1.22901922e+00*(u3*eta*a13) + -1.55571230e+00*(u4*eta*a13) + -7.89032782e+01*(u2*eta3*a13) + 4.55586865e+00*(u3*eta2*a13);

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
	nu5 = 5.87433014e-02*a1 + 1.45656015e-01*eta + -1.91411420e-02 + 1.40612163e-01*u + -1.27122004e+00*(u*a1) + -1.11961126e+00*(u*eta) + -3.01901284e-01*eta2 + -2.73472765e-01*a12 + 1.43641267e-01*(u*u) + -3.80906646e+00*(eta2*a1) + 1.68789746e+00*(eta*a12) + -3.91719696e-01*(u2*a1) + 1.05008309e+01*(u*eta*a1) + 1.92282621e-01*(a12*a1) + 3.15119946e+00*(u*a12) + -1.69695980e+00*(u2*eta) + -5.76944432e-01*(u*eta2) + -2.00533255e-01*(u2*u) + -2.31219430e+00*(u*a13) + 7.45723870e+00*(u2*eta2) + 2.08770777e+00*(u2*eta*a1) + 1.19688578e+01*(eta3*a1) + -2.65571830e+01*(u*eta*a12) + 1.15960042e+01*(u*eta2*eta) + 1.71730518e+00*(u3*eta) + 2.02561677e+00*(u3*a1) + -1.27587987e+00*(eta*a13) + 4.75353738e-01*(u2*a12) + -1.90206505e-01*(u3*u) + -1.10640551e+01*(u2*eta2*eta) + -8.97824357e+01*(u*eta3*a1) + 6.55950114e+00*(u*eta2*a12) + -3.17185692e-01*(u2*a13) + -5.29479571e+00*(u3*a12) + -3.92943630e+00*(u2*eta2*a1) + -2.01978248e+01*(u3*eta*a1) + 2.41520960e+00*(u4*eta) + -1.08788049e+01*(eta3*a12) + 6.42306433e-01*(u4*a1) + 1.98078523e+01*(u*eta*a13) + 7.40005371e-01*(u2*eta*a12) + -6.37821663e-01*(u4*a12) + -1.13354663e+01*(u4*eta2) + -1.42248571e+01*(u2*eta2*a12) + 9.05025681e+00*(eta3*a13) + -1.03811257e+01*(u*eta2*a13) + 3.26366880e+01*(u3*eta2*a1) + 2.04292763e+02*(u*eta3*a12) + 5.60879236e+01*(u3*eta*a12) + -5.83169938e+00*(u4*eta*a1) + -1.51774413e+01*(u3*eta2*eta) + 4.07672541e+00*(u3*a13) + -1.21227061e+02*(u3*eta2*a12) + 2.77929707e-01*(u4*a13) + 3.31233781e+01*(u2*eta3*a12) + 1.96204310e+00*(u4*eta*a12) + -1.32388357e+02*(u*eta3*a13) + 1.77157653e+01*(u4*eta2*eta) + 2.47709905e+01*(u4*eta2*a1) + 5.40520923e+00*(u2*eta2*a13) + -4.53481916e+01*(u3*eta*a13) + 6.74542966e+01*(u3*eta3*a1) + -3.81702003e+00*(u4*eta2*a12) + -3.67371847e+01*(u4*eta3*a1) + -8.41180490e+01*(u3*eta3*a12) + -1.48233092e+01*(u2*eta3*a13) + 1.17948116e+02*(u3*eta2*a13);

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
	nu6 = 1.74800255e-01*a1 + -9.00988999e-03 + -1.52561090e-01*u + 7.67164355e-01*(u*a1) + -1.67631916e+00*(eta*a1) + 3.23664046e+00*(u*eta) + -4.48626830e-01*a12 + 6.95814521e-02*(u*u) + 9.05290786e+00*(eta2*a1) + 4.85547257e+00*(eta*a12) + -1.22568944e+00*(u2*a1) + -1.67760486e+01*(u*eta*a1) + 4.06609923e-01*(eta2*eta) + 4.07279561e-01*(a12*a1) + -1.26597520e+00*(u*a12) + -1.14559089e+00*(u2*eta) + -2.27082932e+01*(u*eta2) + 7.90833918e-02*(u2*u) + 5.23379329e-01*(u*a13) + 6.16788198e+00*(u2*eta2) + 1.25153279e+02*(u*eta2*a1) + 2.27030927e+01*(u2*eta*a1) + -1.95058525e+01*(eta3*a1) + 2.88403216e+01*(u*eta*a12) + 4.93808836e+01*(u*eta2*eta) + -2.79216311e+00*(u3*eta) + 5.41803799e-01*(u3*a1) + -5.09379327e+00*(eta*a13) + 2.42947245e+00*(u2*a12) + 1.93430586e-01*(u3*u) + -2.63358359e+01*(eta2*a12) + -9.82965924e+00*(u2*eta2*eta) + -2.84913289e+02*(u*eta3*a1) + -2.29425009e+02*(u*eta2*a12) + 2.76675011e+01*(eta2*a13) + -1.06981723e+00*(u2*a13) + -2.74416220e+00*(u3*a12) + -1.36708104e+02*(u2*eta2*a1) + 2.55415627e+01*(u3*eta2) + -3.70938714e+00*(u4*eta) + 5.39754077e+01*(eta3*a12) + -2.44661425e-01*(u4*a1) + -1.32558395e+01*(u*eta*a13) + -4.65014079e+01*(u2*eta*a12) + 3.73436356e-01*(u4*a12) + 2.42912534e+01*(u4*eta2) + 2.92438457e+02*(u2*eta2*a12) + -5.34553059e+01*(eta3*a13) + 2.53278133e+02*(u2*eta3*a1) + 1.19598699e+02*(u*eta2*a13) + -6.89276262e+01*(u3*eta2*a1) + 5.44324166e+02*(u*eta3*a12) + 2.72766736e+01*(u3*eta*a12) + 2.54178482e+00*(u4*eta*a1) + 2.27154368e+01*(u2*eta*a13) + -6.41167115e+01*(u3*eta2*eta) + 2.92450668e+00*(u3*a13) + -7.39606675e-01*(u4*a13) + -5.59051290e+02*(u2*eta3*a12) + 1.39059039e+00*(u4*eta*a12) + -3.04007740e+02*(u*eta3*a13) + -1.91434429e+01*(u4*eta2*a1) + -5.26873689e+01*(u4*eta2*eta) + -1.56893839e+02*(u2*eta2*a13) + -3.73211448e+01*(u3*eta*a13) + 2.50402407e+02*(u3*eta3*a1) + -2.05597393e+01*(u4*eta2*a12) + 6.82385927e+01*(u4*eta3*a1) + -2.78597762e+02*(u3*eta3*a12) + 4.50144516e+00*(u4*eta*a13) + 3.15837214e+02*(u2*eta3*a13) + 1.05121504e+02*(u3*eta2*a13);

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
	zeta1 = 4.62354400e-04*a1 + 1.81566528e-05 + 1.15573027e-03*(u*a1) + -1.04600435e-02*(eta*a1) + -3.05523734e-03*(u*eta) + 9.16249734e-04*eta2 + -1.66412465e-03*a12 + -1.60900056e-03*(u*u) + 5.38699977e-02*(eta2*a1) + 3.45885363e-02*(eta*a12) + 6.89229329e-03*(u2*a1) + -4.88359685e-03*(eta2*eta) + 1.38427163e-03*(a12*a1) + -4.01631494e-03*(u*a12) + 2.40736900e-02*(u2*eta) + 3.40554983e-02*(u*eta2) + 3.68302713e-03*(u*a13) + -1.23908173e-01*(u2*eta2) + -1.17469273e-01*(u*eta2*a1) + -9.13986092e-02*(u2*eta*a1) + -7.72915331e-02*(eta3*a1) + 3.05162742e-02*(u*eta*a12) + -8.90168324e-02*(u*eta2*eta) + 3.71112964e-03*(u3*eta) + -2.05360752e-03*(u3*a1) + -2.81361118e-02*(eta*a13) + -1.03192113e-02*(u2*a12) + 1.54737015e-03*(u3*u) + -1.83761601e-01*(eta2*a12) + 2.15079016e-01*(u2*eta2*eta) + 4.08818852e-01*(u*eta3*a1) + 7.26694964e-02*(u*eta2*a12) + 1.52453963e-01*(eta2*a13) + 5.78711184e-03*(u2*a13) + 6.98656947e-03*(u3*a12) + 4.44481215e-01*(u2*eta2*a1) + 1.08720164e-02*(u3*eta*a1) + -4.35532343e-02*(u3*eta2) + -2.19482034e-02*(u4*eta) + 2.87244461e-01*(eta3*a12) + -5.36521564e-03*(u4*a1) + -4.02339989e-02*(u*eta*a13) + 1.24917223e-01*(u2*eta*a12) + 5.75213161e-03*(u4*a12) + 1.05522606e-01*(u4*eta2) + -6.01273011e-01*(u2*eta2*a12) + -2.47748985e-01*(eta3*a13) + -7.65408779e-01*(u2*eta3*a1) + 7.15336001e-02*(u*eta2*a13) + 1.00192400e-01*(u3*eta2*a1) + -5.49478634e-01*(u*eta3*a12) + -7.07743217e-02*(u3*eta*a12) + 5.90214568e-02*(u4*eta*a1) + -6.90765648e-02*(u2*eta*a13) + 1.17154344e-01*(u3*eta2*eta) + -6.34955859e-03*(u3*a13) + 5.15837488e-02*(u3*eta2*a12) + -2.82089388e-03*(u4*a13) + 1.08441725e+00*(u2*eta3*a12) + -3.50325843e-02*(u4*eta*a12) + 1.39326975e-01*(u*eta3*a13) + -2.29317515e-01*(u4*eta2*a1) + -1.72014646e-01*(u4*eta2*eta) + 3.44693270e-01*(u2*eta2*a13) + 7.91286294e-02*(u3*eta*a13) + -4.60698093e-01*(u3*eta3*a1) + 4.04516265e-02*(u4*eta2*a12) + 3.30549747e-01*(u4*eta3*a1) + 5.16187739e-01*(u3*eta3*a12) + 1.25054811e-02*(u4*eta*a13) + -6.54563344e-01*(u2*eta3*a13) + -2.20736127e-01*(u3*eta2*a13);

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
	zeta2 = -4.14484736e+01*a1 + -6.42097760e+00 + -2.40660756e+01*u + 1.58467742e+03*(eta*a1) + 8.72606521e+02*(u*eta) + 2.76391526e+02*eta2 + 1.72918084e+02*a12 + 2.80955509e+02*(u*u) + -1.10685433e+04*(eta2*a1) + -4.95107875e+03*(eta*a12) + -1.31107799e+03*(u2*a1) + -3.11151082e+03*(u*eta*a1) + -6.48258914e+02*(eta2*eta) + -1.43194695e+02*(a12*a1) + 1.74053201e+02*(u*a12) + -4.01249564e+03*(u2*eta) + -7.27018846e+03*(u*eta2) + -2.07080836e+02*(u*a13) + 1.93374883e+04*(u2*eta2) + 3.45116860e+04*(u*eta2*a1) + 1.67881084e+04*(u2*eta*a1) + 2.14652331e+04*(eta3*a1) + 3.47033304e+03*(u*eta*a12) + 1.69018666e+04*(u*eta2*eta) + -5.66982442e+02*(u3*eta) + 3.79736388e+02*(u3*a1) + 3.82280025e+03*(eta*a13) + 2.16716787e+03*(u2*a12) + -2.52412442e+02*(u3*u) + 3.25376492e+04*(eta2*a12) + -3.12644267e+04*(u2*eta2*eta) + -9.01180446e+04*(u*eta3*a1) + -5.63671903e+04*(u*eta2*a12) + -2.47636565e+04*(eta2*a13) + -1.24127102e+03*(u2*a13) + -1.27085855e+03*(u3*a12) + -7.47851889e+04*(u2*eta2*a1) + -2.57851694e+03*(u3*eta*a1) + 6.45761236e+03*(u3*eta2) + 3.40174229e+03*(u4*eta) + -6.18840800e+04*(eta3*a12) + 9.01040002e+02*(u4*a1) + -6.18319260e+02*(u*eta*a13) + -2.62782683e+04*(u2*eta*a12) + -1.03799834e+03*(u4*a12) + -1.52537427e+04*(u4*eta2) + 1.15654196e+05*(u2*eta2*a12) + 4.70249094e+04*(eta3*a13) + 1.15224420e+05*(u2*eta3*a1) + 2.78742823e+04*(u*eta2*a13) + -9.83677922e+03*(u3*eta2*a1) + 1.62865341e+05*(u*eta3*a12) + 1.34163435e+04*(u3*eta*a12) + -9.40065535e+03*(u4*eta*a1) + 1.50686903e+04*(u2*eta*a13) + -1.71460434e+04*(u3*eta2*eta) + 1.06439078e+03*(u3*a13) + -2.11088809e+04*(u3*eta2*a12) + 4.43128612e+02*(u4*a13) + -1.82818519e+05*(u2*eta3*a12) + 7.00652838e+03*(u4*eta*a12) + -9.09827408e+04*(u*eta3*a13) + 3.21463910e+04*(u4*eta2*a1) + 2.29667842e+04*(u4*eta2*eta) + -6.88261782e+04*(u2*eta2*a13) + -1.27567538e+04*(u3*eta*a13) + 5.91720457e+04*(u3*eta3*a1) + -1.04673082e+04*(u4*eta2*a12) + -1.93782327e+03*(u4*eta*a13) + -3.74163265e+04*(u4*eta3*a1) + -5.50958734e+04*(u3*eta3*a12) + 1.15022980e+05*(u2*eta3*a13) + 3.49903288e+04*(u3*eta2*a13);

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
	nu0 = 1.01998477e+03 + 5.63895091e+03*u + -2.81921363e+04*(u*a1) + -8.89833851e+04*(eta*a1) + -1.45136859e+05*(u*eta) + -3.71786756e+04*eta2 + -3.72081547e+03*a12 + -3.40664751e+04*(u*u) + 6.73185443e+05*(eta2*a1) + 2.40920834e+05*(eta*a12) + 1.78351367e+05*(u2*a1) + 7.92123671e+05*(u*eta*a1) + 9.19121950e+04*(eta2*eta) + 3.85402003e+03*(a12*a1) + 4.75590255e+04*(u*a12) + 4.74034671e+05*(u2*eta) + 1.03392575e+06*(u*eta2) + -2.29688973e+04*(u*a13) + -2.22335629e+06*(u2*eta2) + -5.92937011e+06*(u*eta2*a1) + -2.29264486e+06*(u2*eta*a1) + -1.29595261e+06*(eta3*a1) + -1.42664533e+06*(u*eta*a12) + -2.20029226e+06*(u*eta2*eta) + 6.71383185e+04*(u3*eta) + -4.15606353e+04*(u3*a1) + -1.74760915e+05*(eta*a13) + -3.27490509e+05*(u2*a12) + 2.96373435e+04*(u3*u) + -1.63272926e+06*(eta2*a12) + 3.52040171e+06*(u2*eta2*eta) + 1.30191382e+07*(u*eta3*a1) + 1.10769423e+07*(u*eta2*a12) + 1.13287584e+06*(eta2*a13) + 2.03811608e+05*(u2*a13) + 1.37273732e+05*(u3*a12) + 1.01174953e+07*(u2*eta2*a1) + 2.39231118e+05*(u3*eta*a1) + -7.16902198e+05*(u3*eta2) + -4.00603079e+05*(u4*eta) + 2.95052337e+06*(eta3*a12) + -1.03564984e+05*(u4*a1) + 7.79749488e+05*(u*eta*a13) + 4.14178416e+06*(u2*eta*a12) + 1.16727334e+05*(u4*a12) + 1.83996873e+06*(u4*eta2) + -1.84538893e+07*(u2*eta2*a12) + -1.99755919e+06*(eta3*a13) + -1.53607895e+07*(u2*eta3*a1) + -6.36598622e+06*(u*eta2*a13) + 1.20710378e+06*(u3*eta2*a1) + -2.48613279e+07*(u*eta3*a12) + -1.36874838e+06*(u3*eta*a12) + 1.07810118e+06*(u4*eta*a1) + -2.64678771e+06*(u2*eta*a13) + 1.82729007e+06*(u3*eta2*eta) + -1.10102938e+05*(u3*a13) + 2.08068918e+06*(u3*eta2*a12) + -5.26724143e+04*(u4*a13) + 2.87953576e+07*(u2*eta3*a12) + -7.40993721e+05*(u4*eta*a12) + 1.46582065e+07*(u*eta3*a13) + -3.83506533e+06*(u4*eta2*a1) + -2.88318638e+06*(u4*eta2*eta) + 1.23086688e+07*(u2*eta2*a13) + 1.25747028e+06*(u3*eta*a13) + -6.23711893e+06*(u3*eta3*a1) + 9.80657392e+05*(u4*eta2*a12) + 2.27901757e+05*(u4*eta*a13) + 4.95400413e+06*(u4*eta3*a1) + 5.27475605e+06*(u3*eta3*a12) + -2.00966651e+07*(u2*eta3*a13) + -3.34080403e+06*(u3*eta2*a13);

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
	mu1 = -4.01476264e+00*a1 + -6.39419682e-02 + 2.09686397e-01*u + -4.40488219e+00*(u*a1) + 6.51289796e+01*(eta*a1) + 1.50693987e+01*a12 + 6.34825864e+00*(u*u) + -2.37610387e+02*(eta2*a1) + -2.30152317e+02*(eta*a12) + -1.88133846e+01*(u2*a1) + 4.39635526e+01*(u*eta*a1) + -1.33377655e+01*(a12*a1) + 1.38006796e+01*(u*a12) + -9.18428883e+01*(u2*eta) + 4.98713766e+00*(u*eta2) + -2.97473019e+00*(u2*u) + -1.08909770e+01*(u*a13) + 3.01483113e+02*(u2*eta2) + -1.89099934e+02*(u*eta2*a1) + 2.53118207e+02*(u2*eta*a1) + -1.63088749e+02*(u*eta*a12) + 3.90668715e+01*(u3*eta) + 2.79007538e+01*(u3*a1) + 2.00458116e+02*(eta*a13) + -9.70839740e+00*(u3*u) + 8.23754453e+02*(eta2*a12) + 6.51599262e+02*(u*eta2*a12) + -7.09872503e+02*(eta2*a13) + 2.08533954e+01*(u2*a13) + -6.69565739e+01*(u3*a12) + -6.92727223e+02*(u2*eta2*a1) + -3.76793121e+02*(u3*eta*a1) + -1.36929882e+02*(u3*eta2) + 1.39044220e+02*(u4*eta) + 4.08862503e+01*(u4*a1) + 1.37893976e+02*(u*eta*a13) + 5.12751615e+01*(u2*eta*a12) + -4.88860938e+01*(u4*a12) + -4.47155510e+02*(u4*eta2) + -5.93064885e+02*(u2*eta2*a12) + 1.31836214e+03*(u3*eta2*a1) + -5.43740097e+02*(u*eta2*a13) + 9.26050225e+02*(u3*eta*a12) + -5.64185766e+02*(u4*eta*a1) + -3.42504592e+02*(u2*eta*a13) + 4.69216693e+01*(u3*a13) + -3.26840844e+03*(u3*eta2*a12) + 1.25383712e+01*(u4*a13) + 6.42048966e+02*(u4*eta*a12) + 1.63694748e+03*(u4*eta2*a1) + 1.46260063e+03*(u2*eta2*a13) + -6.62367134e+02*(u3*eta*a13) + -1.53370440e+03*(u4*eta2*a12) + -1.31938353e+02*(u4*eta*a13) + 2.36512196e+03*(u3*eta2*a13);

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
	mu2 = 7.11458013e+00*a1 + 1.16490092e+01*eta + -6.68054313e-01 + -2.28868880e+00*u + 2.14654851e+01*(u*a1) + -1.18912460e+02*(eta*a1) + 2.19460301e+01*(u*eta) + -4.45436048e+01*eta2 + -1.59025999e+01*a12 + -1.33560769e+00*(u*u) + 4.35196474e+02*(eta2*a1) + 2.71511634e+02*(eta*a12) + -7.98331642e+00*(u2*a1) + -2.14480915e+02*(u*eta*a1) + 1.06754877e+01*(a12*a1) + -5.10134125e+01*(u*a12) + 3.00817623e+00*(u2*eta) + -3.65926082e+01*(u*eta2) + 2.63442775e+00*(u2*u) + 3.46254895e+01*(u*a13) + 4.20329604e+02*(u*eta2*a1) + 2.18709748e+02*(u2*eta*a1) + 5.06981304e+02*(u*eta*a12) + -1.93551710e+01*(u3*eta) + -2.72111914e+01*(u3*a1) + -1.85233410e+02*(eta*a13) + 2.20096816e+01*(u2*a12) + 2.82340809e+00*(u3*u) + -9.95595548e+02*(eta2*a12) + -1.00258613e+03*(u*eta2*a12) + 6.82347863e+02*(eta2*a13) + -1.06417899e+01*(u2*a13) + 6.91653044e+01*(u3*a12) + -7.16810642e+02*(u2*eta2*a1) + 2.42145653e+02*(u3*eta*a1) + -2.53121889e+01*(u4*eta) + -3.35424939e+02*(u*eta*a13) + -5.14453543e+02*(u2*eta*a12) + -6.30583467e+00*(u4*a12) + 7.20675755e+01*(u4*eta2) + 1.54357513e+03*(u2*eta2*a12) + -3.12854103e+02*(u3*eta2*a1) + 6.35224644e+02*(u*eta2*a13) + -6.52109905e+02*(u3*eta*a12) + -9.68143847e+01*(u4*eta*a1) + 2.66893111e+02*(u2*eta*a13) + -4.94600474e+01*(u3*a13) + 1.05413120e+03*(u3*eta2*a12) + 2.63280599e+02*(u4*eta*a12) + 2.95107027e+02*(u4*eta2*a1) + -6.81477158e+02*(u2*eta2*a13) + 4.73342971e+02*(u3*eta*a13) + -6.14707679e+02*(u4*eta2*a12) + -9.05418955e+01*(u4*eta*a13) + -8.07860005e+02*(u3*eta2*a13);

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
	mu3 = -4.47722808e-02*a1 + -2.49245905e-01*eta + 1.26491268e-02 + 7.55211636e-03*u + 1.21273187e+00*(eta*a1) + -1.61469884e-01*(u*eta) + 8.56728844e-01*eta2 + 2.64670718e-02*a12 + -8.61634927e-02*(u*u) + -4.83297703e+00*(eta2*a1) + -1.51717695e+00*(eta*a12) + 9.67418101e-02*(u2*a1) + 3.59983292e-01*(u*eta*a1) + -5.47965581e-02*(u*a12) + 1.62290616e+00*(u2*eta) + 1.04463122e+00*(u*eta2) + -1.85694906e-02*(u2*u) + 5.24349424e-02*(u*a13) + -5.98652033e+00*(u2*eta2) + -5.03978515e+00*(u*eta2*a1) + -5.08953972e+00*(u2*eta*a1) + 4.21113352e-01*(u3*eta) + 3.67727936e-02*(u3*a1) + 6.64454101e-01*(eta*a13) + 5.35321536e-01*(u2*a12) + 1.12650942e-01*(u3*u) + 7.26847505e+00*(eta2*a12) + 8.34337798e+00*(u*eta2*a12) + -3.84608473e+00*(eta2*a13) + -6.72897102e-01*(u2*a13) + -1.70482402e-02*(u3*a12) + 2.45004295e+01*(u2*eta2*a1) + -1.68340178e+00*(u3*eta*a1) + -2.42341097e+00*(u3*eta2) + -2.08791389e+00*(u4*eta) + -1.44598402e-01*(u4*a1) + -2.44040799e-01*(u*eta*a13) + 1.27145694e+00*(u2*eta*a12) + -6.65037139e-01*(u4*a12) + 7.44407427e+00*(u4*eta2) + -2.58488763e+01*(u2*eta2*a12) + 1.30869577e+01*(u3*eta2*a1) + -4.69426783e+00*(u*eta2*a13) + 2.61590039e+00*(u3*eta*a12) + 6.38772487e+00*(u4*eta*a1) + 3.35287540e+00*(u2*eta*a13) + -2.46388242e+01*(u3*eta2*a12) + 9.10804177e-01*(u4*a13) + -9.90422817e-01*(u4*eta*a12) + -2.87204176e+01*(u4*eta2*a1) + 6.01779349e+00*(u2*eta2*a13) + -1.49952508e+00*(u3*eta*a13) + 2.57666090e+01*(u4*eta2*a12) + -5.53800143e+00*(u4*eta*a13) + 1.53420120e+01*(u3*eta2*a13);

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
	mu4 = -1.31119545e+01*a1 + -1.87597160e+01*eta + 2.08110833e+00 + -2.55570881e+01*u + 1.49978669e+02*(u*a1) + 1.43530341e+02*(eta*a1) + 2.15759641e+02*(u*eta) + -3.02908375e+01*(u*u) + -1.64212922e+02*(eta2*a1) + 9.15340981e+01*(u2*a1) + -1.10070663e+03*(u*eta*a1) + 1.76346922e+01*(a12*a1) + -2.72803311e+02*(u*a12) + 2.16768101e+02*(u2*eta) + -4.65386199e+02*(u*eta2) + 6.51619135e+01*(u2*u) + 1.55810964e+02*(u*a13) + -3.54897736e+02*(u2*eta2) + 1.88120653e+03*(u*eta2*a1) + -2.66719162e+02*(u2*eta*a1) + 1.66167744e+03*(u*eta*a12) + -7.37921592e+02*(u3*eta) + -3.74880261e+02*(u3*a1) + -2.14022998e+02*(eta*a13) + 6.45805541e+01*(u3*u) + -4.10747486e+02*(eta2*a12) + -1.74364238e+03*(u*eta2*a12) + 8.56460715e+02*(eta2*a13) + -7.71599096e+01*(u2*a13) + 6.44411192e+02*(u3*a12) + 3.91348573e+03*(u3*eta*a1) + 2.15168553e+03*(u3*eta2) + -6.40767905e+02*(u4*eta) + -2.12071831e+02*(u4*a1) + -7.52150451e+02*(u*eta*a13) + -1.36721911e+03*(u2*eta*a12) + 1.60310441e+02*(u4*a12) + 1.79245985e+03*(u4*eta2) + 3.54759024e+03*(u2*eta2*a12) + -1.09584964e+04*(u3*eta2*a1) + -6.12872603e+03*(u3*eta*a12) + 1.55544852e+03*(u4*eta*a1) + 1.58171536e+03*(u2*eta*a13) + -3.37284814e+02*(u3*a13) + 1.62072261e+04*(u3*eta2*a12) + -4.98132518e+03*(u4*eta2*a1) + -3.20203953e+03*(u2*eta2*a13) + 2.81250369e+03*(u3*eta*a13) + 2.91546821e+03*(u4*eta2*a12) + -1.05739791e+03*(u4*eta*a13) + -6.68759750e+03*(u3*eta2*a13);

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
	nu4 = 2.46020529e+02*a1 + 4.55311872e+02*eta + -3.08476483e+01 + -6.75160725e+01*u + 3.59504139e+02*(u*a1) + -3.60090887e+03*(eta*a1) + 5.15809727e+02*(u*eta) + -1.58140830e+03*eta2 + -7.01685600e+02*a12 + -9.19645264e+01*(u*u) + 1.18463025e+04*(eta2*a1) + 9.47988923e+03*(eta*a12) + 1.63564074e+02*(u2*a1) + -2.10126171e+03*(u*eta*a1) + 5.54758883e+02*(a12*a1) + -6.24140163e+02*(u*a12) + 6.69235288e+02*(u2*eta) + -8.30723071e+02*(u*eta2) + 1.39904831e+02*(u2*u) + 3.76180476e+02*(u*a13) + -1.22466275e+03*(u2*eta2) + 1.73565834e+03*(u*eta2*a1) + 2.89307920e+03*(u*eta*a12) + -1.37297643e+03*(u3*eta) + -6.58451851e+02*(u3*a1) + -7.14819304e+03*(eta*a13) + 6.08183709e+02*(u2*a12) + 2.59679383e+02*(u3*u) + -2.92244990e+04*(eta2*a12) + 2.12079966e+04*(eta2*a13) + -8.46344976e+02*(u2*a13) + 1.06549287e+03*(u3*a12) + 4.63177371e+03*(u3*eta*a1) + 3.69427785e+03*(u3*eta2) + -2.77568615e+03*(u4*eta) + -9.38034639e+02*(u4*a1) + -1.67432291e+03*(u*eta*a13) + -9.30817912e+03*(u2*eta*a12) + 7.36286483e+02*(u4*a12) + 7.33770733e+03*(u4*eta2) + 1.86938928e+04*(u2*eta2*a12) + -1.01644502e+04*(u3*eta2*a1) + -4.87517272e+03*(u3*eta*a12) + 8.67077698e+03*(u4*eta*a1) + 1.02223934e+04*(u2*eta*a13) + -6.02295824e+02*(u3*a13) + 6.03619667e+03*(u3*eta2*a12) + -4.17877197e+03*(u4*eta*a12) + -2.41937698e+04*(u4*eta2*a1) + -2.02102569e+04*(u2*eta2*a13) + 1.88138952e+03*(u3*eta*a13) + 1.67980895e+04*(u4*eta2*a12) + -2.19142017e+03*(u4*eta*a13);

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
	nu5 = -5.00656930e-01*a1 + -9.63416291e-01*eta + 6.24699447e-02 + -1.73351449e-01*u + 1.36045755e+00*(u*a1) + 7.15666591e+00*(eta*a1) + 1.89437767e+00*(u*eta) + 2.47709860e+00*eta2 + 8.10755815e-01*a12 + -2.72962580e-01*(u*u) + -1.88508010e+01*(eta2*a1) + -1.30771198e+01*(eta*a12) + 1.71690689e+00*(u2*a1) + -1.48969625e+01*(u*eta*a1) + -4.14671385e-01*(a12*a1) + -3.00610370e+00*(u*a12) + 4.48077989e+00*(u2*eta) + -4.89673834e+00*(u*eta2) + 2.69028403e-01*(u2*u) + 1.87274692e+00*(u*a13) + -1.33444416e+01*(u2*eta2) + 3.89559768e+01*(u*eta2*a1) + -2.80861786e+01*(u2*eta*a1) + 3.29026355e+01*(u*eta*a12) + -2.81406580e+00*(u3*eta) + -2.05200823e+00*(u3*a1) + 6.87002392e+00*(eta*a13) + -2.62397830e+00*(u2*a12) + 2.05584954e-01*(u3*u) + 3.42462123e+01*(eta2*a12) + -8.73075714e+01*(u*eta2*a12) + -1.72636879e+01*(eta2*a13) + 1.19666057e+00*(u2*a13) + 4.28252286e+00*(u3*a12) + 7.84694597e+01*(u2*eta2*a1) + 2.09403076e+01*(u3*eta*a1) + 6.72725281e+00*(u3*eta2) + -3.39952232e+00*(u4*eta) + -1.05887023e+00*(u4*a1) + -2.08137205e+01*(u*eta*a13) + 4.78213756e+01*(u2*eta*a12) + 1.18803553e+00*(u4*a12) + 9.18556788e+00*(u4*eta2) + -1.31507751e+02*(u2*eta2*a12) + -5.01603056e+01*(u3*eta2*a1) + 5.57147546e+01*(u*eta2*a13) + -4.29618656e+01*(u3*eta*a12) + 1.69406977e+01*(u4*eta*a1) + -2.53535250e+01*(u2*eta*a13) + -2.52858372e+00*(u3*a13) + 1.02515379e+02*(u3*eta2*a12) + -2.44512381e-01*(u4*a13) + -2.08623141e+01*(u4*eta*a12) + -3.70662526e+01*(u4*eta2*a1) + 7.07334137e+01*(u2*eta2*a13) + 2.47746454e+01*(u3*eta*a13) + 3.07194284e+01*(u4*eta2*a12) + 6.56983077e+00*(u4*eta*a13) + -5.74099880e+01*(u3*eta2*a13);

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
	nu6 = -3.89161859e-03 + 4.31056418e-01*u + -2.62874805e+00*(u*a1) + -4.76059710e+00*(u*eta) + 1.93050195e-01*eta2 + 3.58534628e-01*a12 + 1.65525616e-01*(u*u) + -3.46495511e+00*(eta*a12) + -5.50163766e-01*(u2*a1) + 2.84987942e+01*(u*eta*a1) + -4.39957187e-01*(a12*a1) + 5.01110518e+00*(u*a12) + 1.39553134e+01*(u*eta2) + -9.27041382e-01*(u2*u) + -3.03603787e+00*(u*a13) + -3.70647162e+00*(u2*eta2) + -8.31972005e+01*(u*eta2*a1) + -4.95305626e+00*(u2*eta*a1) + -5.29478700e+01*(u*eta*a12) + 1.16375034e+01*(u3*eta) + 5.37715839e+00*(u3*a1) + 4.42335843e+00*(eta*a13) + -5.10203712e-01*(u3*u) + 7.94516990e+00*(eta2*a12) + 1.52990081e+02*(u*eta2*a12) + -1.07203618e+01*(eta2*a13) + 4.37150141e-01*(u2*a13) + -9.62458687e+00*(u3*a12) + 2.86042065e+01*(u2*eta2*a1) + -6.50459911e+01*(u3*eta*a1) + -3.74366719e+01*(u3*eta2) + 3.63950719e+00*(u4*eta) + 2.12220772e+00*(u4*a1) + 3.11835037e+01*(u*eta*a13) + 2.03926918e+01*(u2*eta*a12) + -3.03577299e+00*(u4*a12) + -7.41516861e+00*(u4*eta2) + -7.30688194e+01*(u2*eta2*a12) + 2.08864324e+02*(u3*eta2*a1) + -8.88520523e+01*(u*eta2*a13) + 1.12398693e+02*(u3*eta*a12) + -8.75520950e+00*(u4*eta*a1) + -1.63383475e+01*(u2*eta*a13) + 5.36613602e+00*(u3*a13) + -3.59932494e+02*(u3*eta2*a12) + 1.62709067e+00*(u4*a13) + 3.45511729e+00*(u4*eta*a12) + 1.38088041e+01*(u4*eta2*a1) + 4.97658853e+01*(u2*eta2*a13) + -6.01320620e+01*(u3*eta*a13) + 1.91308756e+02*(u3*eta2*a13);

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
	zeta1 = -1.87315485e-02*eta + 1.89977505e-04 + 9.73150756e-03*u + -5.68962883e-02*(u*a1) + 1.35609469e-01*(eta*a1) + -7.03017805e-02*(u*eta) + 1.12086803e-01*eta2 + 5.58620515e-03*a12 + 4.51532216e-02*(u*u) + -7.85050597e-01*(eta2*a1) + -3.26496800e-01*(eta*a12) + -2.78892144e-01*(u2*a1) + 3.77294419e-01*(u*eta*a1) + -6.21828067e-03*(a12*a1) + 1.26117850e-01*(u*a12) + -4.65913559e-01*(u2*eta) + 1.42098293e-02*(u2*u) + -9.50196207e-02*(u*a13) + 1.07039174e+00*(u2*eta2) + 2.94197133e+00*(u2*eta*a1) + -9.11079431e-01*(u*eta*a12) + -3.65838844e-01*(u3*eta) + -1.25811824e-01*(u3*a1) + 2.25767790e-01*(eta*a13) + 5.51678175e-01*(u2*a12) + -4.69758598e-02*(u3*u) + 1.65532486e+00*(eta2*a12) + 7.30624517e-01*(u*eta2*a12) + -1.04689279e+00*(eta2*a13) + -3.54881299e-01*(u2*a13) + 2.29723628e-01*(u3*a12) + -7.35423232e+00*(u2*eta2*a1) + 3.02291350e+00*(u3*eta*a1) + 1.69286179e+00*(u3*eta2) + 3.91095595e-01*(u4*eta) + 2.20517629e-01*(u4*a1) + 8.01176031e-01*(u*eta*a13) + -5.99950936e+00*(u2*eta*a12) + -3.99039042e-01*(u4*a12) + -5.30870857e-01*(u4*eta2) + 1.60613259e+01*(u2*eta2*a12) + -1.28155536e+01*(u3*eta2*a1) + -1.27997460e+00*(u*eta2*a13) + -6.03069798e+00*(u3*eta*a12) + -1.25571931e+00*(u4*eta*a1) + 4.01344836e+00*(u2*eta*a13) + -9.53925671e-02*(u3*a13) + 2.56616910e+01*(u3*eta2*a12) + 2.76149010e-01*(u4*a13) + 1.78033773e+00*(u4*eta*a12) + -1.13359264e+01*(u2*eta2*a13) + 3.23762970e+00*(u3*eta*a13) + 1.77193933e+00*(u4*eta2*a12) + -1.43190521e+00*(u4*eta*a13) + -1.45423522e+01*(u3*eta2*a13);

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
	zeta2 = 8.71551143e+02*a1 + 1.91091575e+03*eta + -1.02366280e+02 + -4.24982182e+02*u + 2.59620555e+03*(u*a1) + -1.61050833e+04*(eta*a1) + 4.08465174e+03*(u*eta) + -7.28646829e+03*eta2 + -2.44606853e+03*a12 + -9.84188255e+02*(u*u) + 5.95890379e+04*(eta2*a1) + 4.06406335e+04*(eta*a12) + 4.07986701e+03*(u2*a1) + -2.44482713e+04*(u*eta*a1) + 1.88661576e+03*(a12*a1) + -5.24214701e+03*(u*a12) + 9.29362479e+03*(u2*eta) + -7.65624777e+03*(u*eta2) + 3.38855235e+02*(u2*u) + 3.45258851e+03*(u*a13) + -2.06327894e+04*(u2*eta2) + 4.86916115e+04*(u*eta2*a1) + -3.40376315e+04*(u2*eta*a1) + 5.05435527e+04*(u*eta*a12) + -1.59944928e+03*(u3*eta) + -1.71796436e+03*(u3*a1) + -2.94584715e+04*(eta*a13) + -4.43373044e+03*(u2*a12) + 1.76091981e+03*(u3*u) + -1.41855401e+05*(eta2*a12) + -1.10423069e+05*(u*eta2*a12) + 9.90755405e+04*(eta2*a13) + 1.16488042e+03*(u2*a13) + 3.72361930e+03*(u3*a12) + 7.41414614e+04*(u2*eta2*a1) + -3.26822850e+03*(u3*eta2) + -1.88301584e+04*(u4*eta) + -7.05791841e+03*(u4*a1) + -3.47659053e+04*(u*eta*a13) + 2.68206595e+04*(u2*eta*a12) + 8.64928832e+03*(u4*a12) + 4.68030179e+04*(u4*eta2) + -5.59425417e+04*(u2*eta2*a12) + 4.88186814e+04*(u3*eta2*a1) + 8.29932857e+04*(u*eta2*a13) + 6.30737137e+04*(u4*eta*a1) + -2.97543715e+03*(u3*a13) + -9.28832979e+04*(u3*eta2*a12) + -3.56388711e+03*(u4*a13) + -5.51443616e+04*(u4*eta*a12) + -1.40267938e+05*(u4*eta2*a1) + 7.88246491e+03*(u3*eta*a13) + 9.01533052e+04*(u4*eta2*a12) + 1.21235730e+04*(u4*eta*a13) + 3.12047115e+04*(u3*eta2*a13);

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
	nu0 = -1.25577184e+04*eta + 5.21457500e+02 + 7.14278538e+03*(u*a1) + 6.20947975e+04*(eta*a1) + 2.12885897e+04*(u*eta) + 6.14194093e+04*eta2 + 2.19645760e+04*(u*u) + -3.37904101e+05*(eta2*a1) + -1.14217123e+05*(eta*a12) + -1.44078342e+05*(u2*a1) + -2.21545043e+05*(u*eta*a1) + -2.61587425e+04*(u*a12) + -2.39136873e+05*(u2*eta) + -1.70381477e+05*(u*eta2) + 2.40398290e+04*(u2*u) + 1.82834203e+04*(u*a13) + 6.28412963e+05*(u2*eta2) + 1.23511156e+06*(u*eta2*a1) + 1.61160943e+06*(u2*eta*a1) + 5.65192464e+05*(u*eta*a12) + -3.79473974e+05*(u3*eta) + -1.93444623e+05*(u3*a1) + 6.37001323e+04*(eta*a13) + 2.88067990e+05*(u2*a12) + -2.52786009e+04*(u3*u) + 5.96014686e+05*(eta2*a12) + -2.64404379e+06*(u*eta2*a12) + -3.17266586e+05*(eta2*a13) + -1.79635712e+05*(u2*a13) + 4.37250467e+05*(u3*a12) + -4.47186222e+06*(u2*eta2*a1) + 2.94063266e+06*(u3*eta*a1) + 1.41890280e+06*(u3*eta2) + 2.51260566e+05*(u4*eta) + 1.26252378e+05*(u4*a1) + -3.57681803e+05*(u*eta*a13) + -3.29676191e+06*(u2*eta*a12) + -1.93621063e+05*(u4*a12) + -5.41577158e+05*(u4*eta2) + 9.50570781e+06*(u2*eta2*a12) + -1.03836350e+07*(u3*eta2*a1) + 1.57417465e+06*(u*eta2*a13) + -6.51516612e+06*(u3*eta*a12) + -1.07597971e+06*(u4*eta*a1) + 2.10401610e+06*(u2*eta*a13) + -2.84802438e+05*(u3*a13) + 2.23233226e+07*(u3*eta2*a12) + 9.52813432e+04*(u4*a13) + 1.30571565e+06*(u4*eta*a12) + 1.99921213e+06*(u4*eta2*a1) + -6.24939338e+06*(u2*eta2*a13) + 4.21114069e+06*(u3*eta*a13) + -1.60049166e+06*(u4*eta2*a12) + -4.63115672e+05*(u4*eta*a13) + -1.42280073e+07*(u3*eta2*a13);

	// Return answer
	return nu0;

} // END of NU0 (3,3) fit implementation


#ifdef __cplusplus
}
#endif
