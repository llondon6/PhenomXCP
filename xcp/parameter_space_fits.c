
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
	double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double mu1;

	// Evaluate fit for this parameter
	mu1 = 1.09687476e+01*u + 2.90462866e+00 + -5.85895387e+01*eta + -2.75251411e+01*a1 + -2.03910735e+02*(u*eta) + 3.74239671e+02*eta2 + 5.45850193e+02*(eta*a1) + -8.84044874e+01*(u*a1) + 7.02394289e+01*a12 + -6.16219937e+00*(u*u) + -7.32932130e+02*(eta2*eta) + 2.08680975e+02*(u*a12) + 1.63703671e+03*(u*eta*a1) + 1.14616615e+02*(u2*eta) + 6.42125325e+01*(u2*a1) + -1.87532123e+01*(u2*u) + -3.43854812e+03*(eta2*a1) + -1.36527113e+03*(eta*a12) + 1.23659757e+03*(u*eta2) + -5.34519038e+01*(a12*a1) + -1.46749634e+02*(u*a13) + 1.02052552e+03*(eta*a13) + 3.43182035e+02*(u3*eta) + -1.64822628e+02*(u2*a12) + -6.91744639e+02*(u2*eta2) + 1.51880971e+02*(u3*a1) + 8.48184756e+03*(eta2*a12) + -3.89161189e+03*(u*eta*a12) + -9.91326061e+03*(u*eta2*a1) + -1.19503076e+03*(u2*eta*a1) + 6.71005708e+03*(eta3*a1) + -2.39502241e+03*(u*eta2*eta) + 1.22356951e+02*(u2*a13) + 2.93084759e+03*(u2*eta*a12) + 1.25804298e+03*(u2*eta2*eta) + -6.24608488e+03*(eta2*a13) + -2.81858669e+03*(u3*eta*a1) + -3.53015677e+02*(u3*a12) + 7.21492348e+03*(u2*eta2*a1) + 2.36348029e+04*(u*eta2*a12) + 2.75932025e+03*(u*eta*a13) + 1.91827454e+04*(u*eta3*a1) + -2.04788077e+03*(u3*eta2) + -1.64163676e+04*(eta3*a12) + -1.73214055e+01*(u4*a1) + -1.67976371e+04*(u*eta2*a13) + 2.44445941e+02*(u3*a13) + -1.70290231e+04*(u2*eta2*a12) + 3.20556493e+02*(u4*eta*a1) + -4.57337100e+04*(u*eta3*a12) + 4.49187926e+01*(u4*a12) + 2.30086808e+01*(u4*eta2) + 1.19664248e+04*(eta3*a13) + 6.63457102e+03*(u3*eta*a12) + 3.90956502e+03*(u3*eta2*eta) + -2.09626113e+03*(u2*eta*a13) + -1.31798863e+04*(u2*eta3*a1) + 1.69678747e+04*(u3*eta2*a1) + -4.02283938e+04*(u3*eta2*a12) + -6.73693937e+02*(u4*eta*a12) + -3.12554608e+01*(u4*a13) + -4.63614541e+03*(u3*eta*a13) + 1.17809396e+04*(u2*eta2*a13) + 3.00330033e+04*(u2*eta3*a12) + -3.25733488e+04*(u3*eta3*a1) + 3.24849281e+04*(u*eta3*a13) + -2.02582011e+03*(u4*eta2*a1) + 7.74982458e+04*(u3*eta3*a12) + 2.82340294e+04*(u3*eta2*a13) + 3.44239226e+03*(u4*eta2*a12) + 3.26054008e+03*(u4*eta3*a1) + 3.75138117e+02*(u4*eta*a13) + -2.01302395e+04*(u2*eta3*a13) + -3.77580817e+03*(u4*eta3*a12) + -1.38379156e+03*(u4*eta2*a13) + -5.44765615e+04*(u3*eta3*a13);

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
	double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double mu3;

	// Evaluate fit for this parameter
	mu3 = -6.37699948e+00*u + 1.09050724e+01 + -2.33100626e+02*eta + -8.55286929e+01*a1 + 1.27607115e+02*(u*eta) + 1.43935350e+03*eta2 + 1.81134110e+03*(eta*a1) + 5.54679247e+01*(u*a1) + 1.78936543e+02*a12 + -1.32029368e+02*(u*u) + -2.68513151e+03*(eta2*eta) + -1.04925839e+02*(u*a12) + -1.02804401e+03*(u*eta*a1) + 2.72740556e+03*(u2*eta) + 9.85284492e+02*(u2*a1) + 8.67605739e+00*(u2*u) + -1.12504141e+04*(eta2*a1) + -3.82088078e+03*(eta*a12) + -6.22406515e+02*(u*eta2) + -1.15497175e+02*(a12*a1) + 1.93956120e+02*(u3*u) + 4.77259167e+01*(u*a13) + 2.45685059e+03*(eta*a13) + -1.29368491e+02*(u3*eta) + -2.06371853e+03*(u2*a12) + -1.65187577e+04*(u2*eta2) + -1.03333699e+02*(u3*a1) + 2.39154155e+04*(eta2*a12) + 1.84778942e+03*(u*eta*a12) + 4.66469822e+03*(u*eta2*a1) + -2.05135200e+04*(u2*eta*a1) + 2.12441494e+04*(eta3*a1) + 8.80039628e+02*(u*eta2*eta) + 1.30934772e+03*(u2*a13) + 4.34068412e+04*(u2*eta*a12) + 3.07322274e+04*(u2*eta2*eta) + -1.54285125e+04*(eta2*a13) + 1.49319782e+03*(u3*eta*a1) + 2.44010033e+02*(u3*a12) + 1.24827554e+05*(u2*eta2*a1) + -3.95045830e+03*(u4*eta) + -7.29803072e+03*(u*eta2*a12) + -8.38564926e+02*(u*eta*a13) + -5.85624671e+03*(u*eta3*a1) + 3.64185600e+02*(u3*eta2) + -4.54935086e+04*(eta3*a12) + -1.45147101e+03*(u4*a1) + 2.57205533e+03*(u*eta2*a13) + -1.66334905e+02*(u3*a13) + -2.65810818e+05*(u2*eta2*a12) + 2.96359003e+04*(u4*eta*a1) + 6.47513174e+03*(u*eta3*a12) + 3.06698935e+03*(u4*a12) + 2.37819570e+04*(u4*eta2) + 2.94899777e+04*(eta3*a13) + -3.42019354e+03*(u3*eta*a12) + -2.77505384e+04*(u2*eta*a13) + -2.32715731e+05*(u2*eta3*a1) + -4.61032113e+03*(u3*eta2*a1) + 9.73665287e+03*(u3*eta2*a12) + -6.28820794e+04*(u4*eta*a12) + -1.96498518e+03*(u4*a13) + 2.34650800e+03*(u3*eta*a13) + 1.70882004e+05*(u2*eta2*a13) + 4.97325477e+05*(u2*eta3*a12) + 1.36454049e+03*(u3*eta3*a1) + -1.78440245e+05*(u4*eta2*a1) + -4.41856578e+04*(u4*eta2*eta) + -6.71648215e+03*(u3*eta2*a13) + 3.79362924e+05*(u4*eta2*a12) + 3.30996058e+05*(u4*eta3*a1) + 4.04317843e+04*(u4*eta*a13) + -3.20847364e+05*(u2*eta3*a13) + -7.03807605e+05*(u4*eta3*a12) + -2.44570541e+05*(u4*eta2*a13) + 4.54355429e+05*(u4*eta3*a13);

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
	double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double mu4;

	// Evaluate fit for this parameter
	mu4 = -6.33079589e-01*u + 7.45031542e-05 + 1.78358365e-01*eta + 1.04348406e+01*(u*eta) + -2.02773604e+00*eta2 + -1.55953801e+00*(eta*a1) + 5.09333318e+00*(u*a1) + 2.12220280e-01*a12 + -1.15241701e-01*(u*u) + 4.87065539e+00*(eta2*eta) + -1.16824166e+01*(u*a12) + -8.27415128e+01*(u*eta*a1) + 7.10928781e-01*(u2*eta) + 5.24144886e-01*(u2*a1) + 8.41376366e-01*(u2*u) + 1.67063819e+01*(eta2*a1) + 6.02094248e-01*(eta*a12) + -5.53828143e+01*(u*eta2) + -1.75620863e-01*(a12*a1) + 8.00439802e+00*(u*a13) + -1.40852562e+01*(u3*eta) + -1.83342086e+00*(u2*a12) + -2.13633620e+00*(u2*eta2) + -6.69557493e+00*(u3*a1) + -2.34910136e+01*(eta2*a12) + 1.87730597e+02*(u*eta*a12) + 4.33941225e+02*(u*eta2*a1) + -3.99108820e+01*(eta3*a1) + 9.49810560e+01*(u*eta2*eta) + 1.73043540e+00*(u2*a13) + 1.32686643e+01*(u2*eta*a12) + 5.63183050e+00*(u2*eta2*eta) + 1.45180076e+01*(eta2*a13) + 1.11486175e+02*(u3*eta*a1) + 1.53374633e+01*(u3*a12) + 1.46968021e+00*(u4*eta) + -1.07303214e+01*(u2*eta2*a1) + -9.75934907e+02*(u*eta2*a12) + -1.27274786e+02*(u*eta*a13) + -7.36932571e+02*(u*eta3*a1) + 7.55161350e+01*(u3*eta2) + 6.70458096e+01*(eta3*a12) + 3.57116967e-01*(u4*a1) + 6.56259118e+02*(u*eta2*a13) + -1.05093237e+01*(u3*a13) + -5.55796319e+01*(u2*eta2*a12) + -1.67320694e+01*(u4*eta*a1) + 1.64517310e+03*(u*eta3*a12) + -9.71202751e+00*(u4*eta2) + -4.45643476e+01*(eta3*a13) + -2.53824966e+02*(u3*eta*a12) + -1.30344558e+02*(u3*eta2*eta) + -1.82563379e+01*(u2*eta*a13) + -5.96019092e+02*(u3*eta2*a1) + 1.35036000e+03*(u3*eta2*a12) + 2.06700193e+01*(u4*eta*a12) + -7.30613897e-01*(u4*a13) + 1.72637690e+02*(u3*eta*a13) + 9.08175714e+01*(u2*eta2*a13) + 1.49469329e+02*(u2*eta3*a12) + 1.02726991e+03*(u3*eta3*a1) + 9.79801572e+01*(u4*eta2*a1) + -1.09901833e+03*(u*eta3*a13) + 1.37675357e+01*(u4*eta2*eta) + -2.31818687e+03*(u3*eta3*a12) + -9.13227494e+02*(u3*eta2*a13) + -1.12510483e+02*(u4*eta2*a12) + -1.35408788e+02*(u4*eta3*a1) + -1.96957262e+02*(u2*eta3*a13) + 9.27100915e+01*(u4*eta3*a12) + -6.62517782e+00*(u4*eta2*a13) + 1.56070105e+03*(u3*eta3*a13) + 9.11801036e+01*(u4*eta3*a13);

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
	double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu4;

	// Evaluate fit for this parameter
	nu4 = 4.93062055e-03*u + 3.23498481e-02 + -6.94332931e-01*eta + -2.75859165e-01*a1 + -2.40809794e-01*(u*eta) + 4.26345464e+00*eta2 + 5.94470182e+00*(eta*a1) + -2.80212300e-02*(u*a1) + 5.78736353e-01*a12 + -8.11333563e+00*(eta2*eta) + 1.60774012e+00*(u*eta*a1) + 4.89134876e-01*(u2*eta) + 2.41186799e-01*(u2*a1) + -1.64878769e-02*(u2*u) + -3.68643757e+01*(eta2*a1) + -1.29233535e+01*(eta*a12) + 1.90655533e+00*(u*eta2) + -3.62842073e-01*(a12*a1) + 8.46743870e+00*(eta*a13) + -6.23103649e-02*(u3*u) + 5.44987571e-01*(u3*eta) + -3.01199179e-01*(u2*a12) + -3.89265593e+00*(u2*eta2) + 1.77003576e-01*(u3*a1) + 8.16883393e+01*(eta2*a12) + -1.99758856e+00*(u*eta*a12) + -1.29647628e+01*(u*eta2*a1) + -8.68540623e+00*(u2*eta*a1) + 7.02971673e+01*(eta3*a1) + -4.07418837e+00*(u*eta2*eta) + 1.67511186e+01*(u2*eta*a12) + 7.67653705e+00*(u2*eta2*eta) + -5.48947317e+01*(eta2*a13) + -5.07268439e+00*(u3*eta*a1) + 7.56598461e-01*(u4*eta) + -4.65241625e-01*(u3*a12) + 5.90837937e+01*(u2*eta2*a1) + 1.76268151e+01*(u*eta2*a12) + 8.37457158e-01*(u*eta*a13) + 2.79101192e+01*(u*eta3*a1) + -4.08215498e+00*(u3*eta2) + -1.57396139e+02*(eta3*a12) + 1.73054505e-01*(u4*a1) + -7.11337771e+00*(u*eta2*a13) + 3.73734021e-01*(u3*a13) + -1.22687207e+02*(u2*eta2*a12) + -3.85884536e+01*(u*eta3*a12) + -5.76302395e-01*(u4*a12) + -3.79500748e+00*(u4*eta2) + 1.07487961e+02*(eta3*a13) + 1.19570301e+01*(u3*eta*a12) + 8.65551264e+00*(u3*eta2*eta) + -8.73881466e+00*(u2*eta*a13) + -1.12479878e+02*(u2*eta3*a1) + 3.60261918e+01*(u3*eta2*a1) + -8.06540684e+01*(u3*eta2*a12) + 2.35625891e+00*(u4*eta*a12) + 5.89952064e-01*(u4*a13) + -8.52602512e+00*(u3*eta*a13) + 7.27210634e+01*(u2*eta2*a13) + 2.41382375e+02*(u2*eta3*a12) + -7.41231922e+01*(u3*eta3*a1) + 1.51567699e+01*(u*eta3*a13) + -4.28467469e+00*(u4*eta2*a1) + 7.04875997e+00*(u4*eta2*eta) + 1.60719278e+02*(u3*eta3*a12) + 5.44658835e+01*(u3*eta2*a13) + 5.76349404e+00*(u4*eta3*a1) + -4.44247062e+00*(u4*eta*a13) + -1.50739335e+02*(u2*eta3*a13) + 1.31228702e+01*(u4*eta2*a13) + -1.05165885e+02*(u3*eta3*a13) + -1.96198340e+01*(u4*eta3*a13);

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
	double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu5;

	// Evaluate fit for this parameter
	nu5 = -5.46663782e-01*u + -6.24533439e-05 + 2.87268336e-01*eta + 1.14277640e+01*(u*eta) + -2.72873890e+00*eta2 + -2.26422782e+00*(eta*a1) + 3.07986460e+00*(u*a1) + -1.52124025e-01*a12 + 1.25752822e+00*(u*u) + 6.16852102e+00*(eta2*eta) + -5.36625251e+00*(u*a12) + -6.39270278e+01*(u*eta*a1) + -2.77572295e+01*(u2*eta) + -6.97398003e+00*(u2*a1) + 2.44379632e-01*(u2*u) + 2.10230658e+01*(eta2*a1) + 7.51709212e+00*(eta*a12) + -7.51066825e+01*(u*eta2) + 1.80744791e-01*(a12*a1) + -1.05702384e+00*(u3*u) + 3.23207819e+00*(u*a13) + -6.47200945e+00*(eta*a13) + -8.69094437e+00*(u3*eta) + 1.19814365e+01*(u2*a12) + 1.84422573e+02*(u2*eta2) + -6.02857516e+01*(eta2*a12) + 1.11475216e+02*(u*eta*a12) + 4.16489885e+02*(u*eta2*a1) + 1.56398350e+02*(u2*eta*a1) + -4.68333820e+01*(eta3*a1) + 1.56136698e+02*(u*eta2*eta) + -6.77914431e+00*(u2*a13) + -2.74912221e+02*(u2*eta*a12) + -3.81589361e+02*(u2*eta2*eta) + 4.88915926e+01*(eta2*a13) + 2.98532928e+01*(u3*eta*a1) + -2.19362806e+00*(u3*a12) + -1.03916128e+03*(u2*eta2*a1) + 2.80160022e+01*(u4*eta) + -7.23274520e+02*(u*eta2*a12) + -6.54209931e+01*(u*eta*a13) + -8.58804776e+02*(u*eta3*a1) + 7.12807202e+01*(u3*eta2) + 1.28236962e+02*(eta3*a12) + 3.83119173e+00*(u4*a1) + 4.16397159e+02*(u*eta2*a13) + 1.55490046e+00*(u3*a13) + 1.83461053e+03*(u2*eta2*a12) + -1.26957286e+02*(u4*eta*a1) + 1.48359142e+03*(u*eta3*a12) + -2.62773599e+00*(u4*a12) + -2.03250136e+02*(u4*eta2) + -1.02051234e+02*(eta3*a13) + -2.40388728e+01*(u3*eta*a12) + -1.66147494e+02*(u3*eta2*eta) + 1.56920400e+02*(u2*eta*a13) + 2.13871192e+03*(u2*eta3*a1) + -3.13572455e+02*(u3*eta2*a1) + 4.29362855e+02*(u3*eta2*a12) + 1.64046077e+02*(u4*eta*a12) + 1.10917425e+01*(u3*eta*a13) + -1.04649225e+03*(u2*eta2*a13) + -3.76791822e+03*(u2*eta3*a12) + 8.00218042e+02*(u3*eta3*a1) + -8.42019622e+02*(u*eta3*a13) + 9.86734543e+02*(u4*eta2*a1) + 4.42912785e+02*(u4*eta2*eta) + -1.23214261e+03*(u3*eta3*a12) + -2.38653641e+02*(u3*eta2*a13) + -1.43810975e+03*(u4*eta2*a12) + -2.21475587e+03*(u4*eta3*a1) + -7.05142024e+01*(u4*eta*a13) + 2.14097010e+03*(u2*eta3*a13) + 3.38387515e+03*(u4*eta3*a12) + 6.95753877e+02*(u4*eta2*a13) + 6.99059646e+02*(u3*eta3*a13) + -1.70109241e+03*(u4*eta3*a13);

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
	double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double nu6;

	// Evaluate fit for this parameter
	nu6 = -8.25429752e-02*u + -3.91729253e-03 + 1.16825375e-01*eta + 7.11647171e-02*a1 + 7.02943505e-01*(u*eta) + -4.98807800e-01*eta2 + -1.51057791e+00*(eta*a1) + 6.39305397e-01*(u*a1) + 5.76621804e-01*(u*u) + -1.70750517e+00*(u*a12) + -4.93453027e+00*(u*eta*a1) + -1.16730752e+01*(u2*eta) + -5.58039686e+00*(u2*a1) + 6.32353615e-02*(u2*u) + 7.84794562e+00*(eta2*a1) + 7.59396306e-01*(eta*a12) + -1.53631199e+00*(u*eta2) + 1.39631365e+00*(u*a13) + -3.75147239e-01*(eta*a13) + -4.41424867e-01*(u3*eta) + -9.96941878e-01*(u3*u) + 1.28481989e+01*(u2*a12) + 6.63384721e+01*(u2*eta2) + -3.60755788e-01*(u3*a1) + -4.18229516e+00*(eta2*a12) + 1.50133723e+01*(u*eta*a12) + 9.59623217e+00*(u*eta2*a1) + 1.10675020e+02*(u2*eta*a1) + -9.42331216e+00*(eta3*a1) + -8.56670561e+00*(u2*a13) + -2.54664029e+02*(u2*eta*a12) + -1.13401853e+02*(u2*eta2*eta) + 2.22591438e+00*(eta2*a13) + 1.45283624e+00*(u3*eta*a1) + 1.97122437e+01*(u4*eta) + 8.45529745e-01*(u3*a12) + -6.32869029e+02*(u2*eta2*a1) + -4.25747863e+01*(u*eta2*a12) + -1.37113566e+01*(u*eta*a13) + 7.90121855e-01*(u3*eta2) + 9.21762256e+00*(u4*a1) + 4.66645973e+01*(u*eta2*a13) + -6.54343648e-01*(u3*a13) + 1.46669336e+03*(u2*eta2*a12) + -1.80198708e+02*(u4*eta*a1) + 3.93908845e+01*(u*eta3*a12) + -1.12687769e+02*(u4*eta2) + -3.18414297e+00*(u3*eta*a12) + 1.70665903e+02*(u2*eta*a13) + 1.10270200e+03*(u2*eta3*a1) + -2.15603868e+01*(u4*a12) + 4.20304291e+02*(u4*eta*a12) + 1.46508489e+01*(u4*a13) + 2.68349810e+00*(u3*eta*a13) + -9.91229312e+02*(u2*eta2*a13) + -2.58221473e+03*(u2*eta3*a12) + 1.03334119e+03*(u4*eta2*a1) + -5.67755677e+01*(u*eta3*a13) + 1.96708968e+02*(u4*eta2*eta) + -4.61368755e+00*(u3*eta3*a12) + -2.41814471e+03*(u4*eta2*a12) + -1.82320552e+03*(u4*eta3*a1) + -2.86194289e+02*(u4*eta*a13) + 1.76139821e+03*(u2*eta3*a13) + 4.29164796e+03*(u4*eta3*a12) + 1.65484056e+03*(u4*eta2*a13) + -2.95444187e+03*(u4*eta3*a13);

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
	double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double zeta1;

	// Evaluate fit for this parameter
	zeta1 = 2.69377655e-04*u + 2.30095590e-05 + -2.37317912e-04*eta + -5.25352758e-03*(u*eta) + 5.60763922e-04*eta2 + -1.49057504e-03*(eta*a1) + -9.96727320e-04*(u*a1) + -1.44485866e-04*a12 + -2.17421479e-04*(u*u) + 8.80204894e-04*(u*a12) + 1.96711998e-02*(u*eta*a1) + 3.47341288e-03*(u2*eta) + 5.90491363e-04*(u2*a1) + -4.91755106e-04*(u2*u) + 1.55154011e-02*(eta2*a1) + 6.93721694e-03*(eta*a12) + 3.12920106e-02*(u*eta2) + 1.76748080e-04*(a12*a1) + -6.67723445e-03*(eta*a13) + 9.45881267e-03*(u3*eta) + 4.27253564e-04*(u3*u) + -5.78898118e-04*(u2*a12) + -2.04430782e-02*(u2*eta2) + 1.92543100e-03*(u3*a1) + -5.83357042e-02*(eta2*a12) + -1.76506232e-02*(u*eta*a12) + -1.16955749e-01*(u*eta2*a1) + -6.51123483e-03*(u2*eta*a1) + -3.75126727e-02*(eta3*a1) + -5.97760990e-02*(u*eta2*eta) + 2.41767568e-04*(u2*a13) + 4.11626389e-02*(u2*eta2*eta) + 5.29782666e-02*(eta2*a13) + -3.75523933e-02*(u3*eta*a1) + -7.62828710e-03*(u4*eta) + -1.75394813e-03*(u3*a12) + 3.86937241e-02*(u2*eta2*a1) + 1.02444697e-01*(u*eta2*a12) + 5.70691431e-04*(u*eta*a13) + 2.26268759e-01*(u*eta3*a1) + -5.67839754e-02*(u3*eta2) + 1.30178244e-01*(eta3*a12) + -2.03500476e-03*(u4*a1) + 3.66092213e-02*(u4*eta*a1) + -2.00673802e-01*(u*eta3*a12) + 3.81574116e-03*(u4*a12) + 4.64608484e-02*(u4*eta2) + -1.15553657e-01*(eta3*a13) + 3.42324286e-02*(u3*eta*a12) + 1.09910260e-01*(u3*eta2*eta) + 3.55954319e-03*(u2*eta*a13) + -9.83684361e-02*(u2*eta3*a1) + 2.28156740e-01*(u3*eta2*a1) + -2.09468388e-01*(u3*eta2*a12) + -7.02795851e-02*(u4*eta*a12) + -2.60097017e-03*(u4*a13) + -2.22130925e-02*(u2*eta2*a13) + 6.10382644e-02*(u2*eta3*a12) + -4.51541660e-01*(u3*eta3*a1) + -2.34092419e-01*(u4*eta2*a1) + -9.34034199e-02*(u4*eta2*eta) + 4.30805666e-01*(u3*eta3*a12) + 4.64522539e-01*(u4*eta2*a12) + 5.00416394e-01*(u4*eta3*a1) + 4.85850454e-02*(u4*eta*a13) + -1.02391509e+00*(u4*eta3*a12) + -3.24413057e-01*(u4*eta2*a13) + -1.45137082e-02*(u3*eta3*a13) + 7.18325640e-01*(u4*eta3*a13);

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
	double u4 = u3*u;
	double a12 = a1*a1;
	double a13 = a12*a1;
	double eta2 = eta*eta;
	double eta3 = eta2*eta;
	double zeta2;

	// Evaluate fit for this parameter
	zeta2 = -3.69589539e+01*u + 1.93713228e+00 + -6.55043435e+01*eta + -5.28133717e+01*a1 + 7.32880744e+02*(u*eta) + 5.35055604e+02*eta2 + 1.25464178e+03*(eta*a1) + 1.08122299e+02*(u*a1) + 1.32819093e+02*a12 + 1.50631848e+01*(u*u) + -1.18219387e+03*(eta2*eta) + -4.96034010e+01*(u*a12) + -2.35859551e+03*(u*eta*a1) + -1.73886677e+02*(u2*eta) + 1.40663102e+02*(u2*a1) + 7.12386229e+01*(u2*u) + -8.50579346e+03*(eta2*a1) + -3.27633464e+03*(eta*a12) + -4.54544140e+03*(u*eta2) + -9.16745350e+01*(a12*a1) + -5.57880952e+01*(u3*u) + -3.66957005e+01*(u*a13) + 2.36372749e+03*(eta*a13) + -1.36149312e+03*(u3*eta) + -4.32015453e+02*(u2*a12) + 9.98174962e+02*(u2*eta2) + -2.33756487e+02*(u3*a1) + 2.23154517e+04*(eta2*a12) + 1.80123692e+03*(u*eta*a12) + 1.60343420e+04*(u*eta2*a1) + -3.42835554e+03*(u2*eta*a1) + 1.71931049e+04*(eta3*a1) + 9.00458746e+03*(u*eta2*eta) + 2.92504030e+02*(u2*a13) + 9.90548340e+03*(u2*eta*a12) + -2.26054629e+03*(u2*eta2*eta) + -1.63522067e+04*(eta2*a13) + 4.76949896e+03*(u3*eta*a1) + 1.48375862e+02*(u3*a12) + 2.00151495e+04*(u2*eta2*a1) + 9.70972860e+02*(u4*eta) + -1.60491475e+04*(u*eta2*a12) + -3.45979148e+04*(u*eta3*a1) + 8.28992510e+03*(u3*eta2) + -4.51018410e+04*(eta3*a12) + 9.21753667e+01*(u4*a1) + 3.97885152e+03*(u*eta2*a13) + 4.57886903e+01*(u3*a13) + -5.76892458e+04*(u2*eta2*a12) + -1.41584847e+03*(u4*eta*a1) + 4.10486155e+04*(u*eta3*a12) + -1.48234138e+01*(u4*a12) + -5.88237724e+03*(u4*eta2) + 3.33207121e+04*(eta3*a13) + -3.87233695e+03*(u3*eta*a12) + -1.62854582e+04*(u3*eta2*eta) + -6.72535766e+03*(u2*eta*a13) + -3.32420282e+04*(u2*eta3*a1) + -3.13731094e+04*(u3*eta2*a1) + 3.13183133e+04*(u3*eta2*a12) + 3.91134203e+04*(u2*eta2*a13) + 9.73621220e+04*(u2*eta3*a12) + 6.65177903e+04*(u3*eta3*a1) + -1.50234948e+04*(u*eta3*a13) + 1.11939608e+04*(u4*eta2*a1) + 1.19138075e+04*(u4*eta2*eta) + -7.73776804e+04*(u3*eta3*a12) + -6.16590359e+03*(u3*eta2*a13) + -8.51133294e+03*(u4*eta2*a12) + -3.07112735e+04*(u4*eta3*a1) + -6.57208215e+04*(u2*eta3*a13) + 4.10894112e+04*(u4*eta3*a12) + 6.70223204e+03*(u4*eta2*a13) + 2.43758699e+04*(u3*eta3*a13) + -3.09059530e+04*(u4*eta3*a13);

	// Return answer
	return zeta2;

} // END of ZETA2 fit implementation


#ifdef __cplusplus
}
#endif
