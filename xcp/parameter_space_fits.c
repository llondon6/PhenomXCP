
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
	mu1 = -1.62091103e-01 + 1.50124214e+00*eta + 1.90706970e+00*a1 + -3.42699257e+00*a12 + -3.32248126e+00*eta2 + -2.68046241e+01*(eta*a1) + 1.37774445e+00*(u*eta) + -1.01504973e+01*(u*eta*a1) + 5.33187514e+01*(eta*a12) + -3.52308601e-01*(u*a12) + -1.82376393e-01*(u2*u) + 1.22665778e+02*(eta2*a1) + -1.98488852e+00*(u2*a1) + 3.78342203e+00*(u2*eta) + -1.07019788e+01*(u*eta2) + 1.49347104e+00*(u3*a1) + 1.64320952e+01*(u*eta2*eta) + -3.58726324e+01*(u2*eta2) + -1.82652708e+02*(eta3*a1) + 3.24180081e+00*(u2*a12) + 1.52590423e+01*(u2*eta*a1) + 1.91136227e+01*(u*eta*a12) + -2.62542215e+02*(eta2*a12) + 7.92417934e+01*(u*eta2*a1) + 8.61904241e+01*(u2*eta2*eta) + -3.53306344e+01*(u2*eta*a12) + -1.29041663e+02*(u*eta2*a12) + -1.30236306e+02*(u*eta3*a1) + 4.17819264e+02*(eta3*a12) + -1.21212306e+00*(u3*a12) + -1.16218902e+01*(u3*eta*a1) + 1.20929725e+01*(u3*eta2) + -3.12753539e+01*(u3*eta2*eta) + 9.18559260e+01*(u2*eta2*a12) + 1.39856046e+01*(u3*eta2*a1) + 2.09973607e+02*(u*eta3*a12) + 7.16516621e+00*(u3*eta*a12) + -1.28805048e+02*(u2*eta3*a1);

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
	mu3 = 3.68259524e-01 + -3.75893675e+00*eta + -1.66698038e+00*a1 + -5.17175349e-01*u + 1.55953680e+01*(u*a1) + -1.18589317e+00*(u*u) + 1.74715847e+01*(eta*a1) + -2.27365402e+02*(u*eta*a1) + -2.57325304e+01*(u*a12) + -2.19036044e+00*(u2*u) + 9.17204476e+00*(u2*a1) + 4.15009142e+01*(eta2*eta) + -5.31458900e+00*(u3*a1) + 3.96946254e+01*(u*eta2*eta) + 1.03511315e+02*(u2*eta2) + -1.91237172e+02*(eta3*a1) + -5.88934415e+01*(u2*eta*a1) + 5.57765929e+01*(u3*eta) + 3.77363119e+02*(u*eta*a12) + -7.80878355e+01*(eta2*a12) + 1.28021317e+03*(u*eta2*a1) + -1.24472719e+01*(u2*a12) + -3.45913475e+02*(u2*eta2*eta) + 1.36696757e+02*(u2*eta*a12) + -2.02730159e+03*(u*eta2*a12) + -2.50718388e+03*(u*eta3*a1) + 3.23842796e+02*(eta3*a12) + -1.73188952e+02*(u2*eta2*a1) + 9.23171718e+00*(u3*a12) + -3.26412671e+02*(u3*eta2) + 5.37808839e+02*(u3*eta2*eta) + -3.53367947e+02*(u2*eta2*a12) + 1.00416299e+02*(u3*eta2*a1) + 3.73671042e+03*(u*eta3*a12) + -3.92414588e+01*(u3*eta*a12) + 1.08775034e+03*(u2*eta3*a1);

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
	mu4 = -3.06325774e-04 + 1.50029814e-04*u + 2.16410229e-04*a12 + -7.07189505e-04*(u*a1) + -6.13221967e-04*(u*u) + 1.58067476e-02*(u*delta) + 1.03763559e-02*(delta*delta) + -1.72769515e-02*(delta*a1) + -1.54278336e-02*(delta2*delta) + -4.99852933e-02*(u*delta*a1) + 4.86137066e-02*(delta*a12) + -2.37995922e-02*(u*delta*delta) + 3.60373808e-02*(delta3*a1) + 2.18604378e-02*(u*delta2*a1) + 1.90506007e-02*(u*delta2*delta) + 1.55919299e-02*(u2*delta*a1) + -9.81164557e-02*(delta2*a12) + 3.84549792e-02*(u*delta*a12) + -2.78856088e-02*(u3*delta) + -2.04083085e-02*(u2*delta*delta) + 6.41979786e-02*(delta3*a12) + 9.47622491e-02*(u3*delta*a1) + 5.85421103e-02*(u2*delta2*a1) + 2.96727076e-02*(u2*delta2*delta) + 3.61549306e-02*(u3*delta*delta) + -3.29874225e-02*(u*delta3*a1) + -5.31585194e-02*(u2*delta*a12) + 5.30972381e-02*(u*delta3*a12) + -7.21832050e-02*(u3*delta*a12) + -5.46538561e-02*(u3*delta2*a1) + 7.24508198e-02*(u2*delta2*a12) + -1.19900615e-01*(u2*delta3*a1) + -1.34243683e-02*(u3*delta2*delta);

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
	nu4 = -1.20466596e-04 + -2.39270601e-04*a1 + 5.99136063e-04*u + -3.34252335e-03*(delta) + 1.58092356e-04*a12 + -3.08649396e-03*(u*a1) + -3.25184573e-04*(u*u) + -2.02670998e-03*(u*delta) + 6.48928228e-03*(delta*delta) + 1.95156468e-02*(delta*a1) + 3.56027095e-03*(u*a12) + -7.18682292e-04*(u2*u) + -1.20750940e-03*(delta2*delta) + -2.51679651e-02*(delta*a12) + -3.44053295e-02*(delta2*a1) + 5.40421301e-03*(u2*delta) + -3.06269848e-03*(u*delta*delta) + 3.41781207e-03*(delta3*a1) + 3.93427934e-03*(u3*a1) + 4.90788533e-02*(u*delta2*a1) + 4.37523178e-03*(u*delta2*delta) + -3.24199734e-02*(u2*delta*a1) + 4.05103565e-02*(delta2*a12) + 6.66543153e-03*(u3*delta) + -8.13534283e-03*(u2*delta*delta) + -2.25748365e-02*(u3*delta*a1) + -4.44761476e-03*(u3*a12) + -5.41957304e-02*(u*delta2*a12) + 5.34779824e-02*(u2*delta2*a1) + -2.18052835e-03*(u2*delta2*delta) + -9.61532884e-03*(u3*delta*delta) + -4.17217170e-02*(u*delta3*a1) + 4.49297906e-02*(u2*delta*a12) + 3.43916360e-02*(u*delta3*a12) + 2.12660975e-02*(u3*delta*a12) + 1.17844056e-02*(u3*delta2*a1) + -7.91128555e-02*(u2*delta2*a12) + 9.58931682e-03*(u2*delta3*a1) + 4.11837351e-03*(u3*delta2*delta);

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
	nu5 = 1.53731487e-02 + -2.70626633e-01*eta + -6.41506741e-02*a1 + -1.24706185e-02*u + 1.45141329e+00*eta2 + 6.26371907e-02*(u*a1) + -2.56898624e-02*(u*u) + 1.12579161e+00*(eta*a1) + 5.44015520e-02*(u*eta) + -3.86449529e-01*(eta*a12) + -8.80535527e-02*(u*a12) + 3.12945490e-03*(u2*u) + -5.78645575e+00*(eta2*a1) + 1.06465817e-01*(u2*a1) + 3.92396314e-01*(u2*eta) + -2.46328654e+00*(eta2*eta) + -4.23705636e-02*(u3*a1) + -1.95398023e+00*(u2*eta2) + 9.29302876e+00*(eta3*a1) + -1.51192062e+00*(u2*eta*a1) + 5.16009013e-02*(u3*eta) + 1.49286841e+00*(eta2*a12) + -2.18174318e+00*(u*eta2*a1) + -5.36885180e-02*(u2*a12) + 3.19638655e+00*(u2*eta2*eta) + 6.60967802e-01*(u2*eta*a12) + 3.10642469e+00*(u*eta2*a12) + 4.39125021e+00*(u*eta3*a1) + 7.06095137e+00*(u2*eta2*a1) + 8.95002040e-02*(u3*a12) + -2.99268353e-01*(u3*eta2) + -1.74012403e+00*(u2*eta2*a12) + 8.52442641e-01*(u3*eta2*a1) + -6.51220730e+00*(u*eta3*a12) + -3.93199115e-01*(u3*eta*a12) + -1.09380549e+01*(u2*eta3*a1);

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
	nu6 = 9.45580851e-03 + -1.46245926e-01*eta + -5.85459891e-02*a1 + 2.04903702e-02*u + 1.18886306e-01*a12 + 7.18434981e-01*eta2 + -1.96469661e-01*(u*a1) + -3.47437496e-03*(u*u) + 9.27264136e-01*(eta*a1) + -2.44473255e-01*(u*eta) + 2.83226750e+00*(u*eta*a1) + -1.84117341e+00*(eta*a12) + 2.61118689e-01*(u*a12) + 1.35631746e-02*(u2*u) + -4.75451162e+00*(eta2*a1) + 4.32286154e-02*(u2*eta) + -1.14703929e+00*(eta2*eta) + 1.23616837e+00*(u*eta2) + 3.45136737e-02*(u3*a1) + -2.39178442e+00*(u*eta2*eta) + -1.21766072e-01*(u2*eta2) + 7.94957547e+00*(eta3*a1) + -1.28503219e-03*(u2*a12) + -3.52776796e-01*(u3*eta) + -3.73918285e+00*(u*eta*a12) + 9.51536332e+00*(eta2*a12) + -1.49608747e+01*(u*eta2*a1) + 1.90696827e+01*(u*eta2*a12) + 2.72868988e+01*(u*eta3*a1) + -1.61850047e+01*(eta3*a12) + -7.83961156e-02*(u3*a12) + 1.15765718e-01*(u3*eta*a1) + 1.92766626e+00*(u3*eta2) + -2.86470589e+00*(u3*eta2*eta) + -1.08735842e+00*(u3*eta2*a1) + -3.33036428e+01*(u*eta3*a12) + 3.23402582e-01*(u3*eta*a12);

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
	zeta1 = -2.00997903e-05 + 3.17149324e-04*eta + 8.73154180e-05*a1 + -5.26263677e-05*a12 + -1.63554003e-03*eta2 + 7.95231216e-05*(u*a1) + 4.65020619e-05*(u*u) + -1.24570299e-03*(eta*a1) + 6.45701379e-05*(u*eta) + -2.07423871e-03*(u*eta*a1) + 6.03341340e-04*(eta*a12) + -8.42009040e-05*(u*a12) + -2.42011859e-05*(u2*u) + 5.78645573e-03*(eta2*a1) + -1.81558573e-04*(u2*a1) + -7.45708659e-04*(u2*eta) + 2.75384100e-03*(eta2*eta) + -5.62500370e-04*(u*eta2) + 2.83807292e-05*(u3*a1) + 1.07742853e-03*(u*eta2*eta) + 3.99286548e-03*(u2*eta2) + -8.79197946e-03*(eta3*a1) + 9.48075872e-05*(u2*a12) + 2.61207901e-03*(u2*eta*a1) + 4.26155729e-04*(u3*eta) + 2.47565717e-03*(u*eta*a12) + -1.57239637e-03*(eta2*a12) + 1.39176744e-02*(u*eta2*a1) + -7.02759247e-03*(u2*eta2*eta) + -8.65336745e-04*(u2*eta*a12) + -1.67284457e-02*(u*eta2*a12) + -2.68619322e-02*(u*eta3*a1) + -1.29768466e-02*(u2*eta2*a1) + -4.53942788e-05*(u3*a12) + -2.20054330e-04*(u3*eta*a1) + -2.39753072e-03*(u3*eta2) + 4.66248796e-03*(u3*eta2*eta) + 1.97060407e-03*(u2*eta2*a12) + 3.19845960e-02*(u*eta3*a12) + 2.86548707e-04*(u3*eta*a12) + 2.16785100e-02*(u2*eta3*a1);

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
	zeta2 = 2.07487222e+00 + -2.04839349e+01*eta + -1.23852521e+01*a1 + 1.06555891e+01*a12 + 4.94725341e+01*eta2 + -6.87323846e+00*(u*a1) + -5.20752874e+00*(u*u) + 1.24919758e+02*(eta*a1) + -1.36278711e+01*(u*eta) + 2.25697919e+02*(u*eta*a1) + -8.49005568e+01*(eta*a12) + 6.09566406e+00*(u*a12) + 3.21945636e+00*(u2*u) + -3.00999981e+02*(eta2*a1) + 2.41046932e+01*(u2*a1) + 7.44609997e+01*(u2*eta) + 1.18044287e+02*(u*eta2) + -7.83228475e+00*(u3*a1) + -2.35129467e+02*(u*eta2*eta) + -3.40842472e+02*(u2*eta2) + -3.20377640e+02*(u2*eta*a1) + -4.27601356e+01*(u3*eta) + -2.38744185e+02*(u*eta*a12) + -1.57895187e+03*(u*eta2*a1) + -1.57078212e+01*(u2*a12) + 5.19008491e+02*(u2*eta2*eta) + 1.63131912e+02*(u2*eta*a12) + 1.71086173e+03*(u*eta2*a12) + 3.03009484e+03*(u*eta3*a1) + 6.78102010e+02*(eta3*a12) + 1.37324163e+03*(u2*eta2*a1) + 9.84203127e+00*(u3*a12) + 4.48017416e+01*(u3*eta*a1) + 1.97954003e+02*(u3*eta2) + -3.51863611e+02*(u3*eta2*eta) + -3.92652812e+02*(u2*eta2*a12) + -3.28141202e+03*(u*eta3*a12) + -5.41981873e+01*(u3*eta*a12) + -1.95219863e+03*(u2*eta3*a1);

	// Return answer
	return zeta2;

} // END of ZETA2 fit implementation


#ifdef __cplusplus
}
#endif
