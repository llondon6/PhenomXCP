

// Import usefuls
#include <math.h>



// MU1 fit implementation 
// Header formatting for MU1
static double IMRPhenomXCP_MU1( double theta, double eta, double a1 ){ 

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
	mu1 = -4.20618435e+00*a1 + 1.49177246e+00*u + 6.14648534e-02 + 1.63760994e+01*a12 + 1.45341069e+00*eta2 + 7.15666896e+01*(eta*a1) + -1.10330150e+01*(u*eta) + -1.02823271e+01*(u*a1) + -1.74771159e+00*(u2*u) + 4.45935726e+01*(u*eta*a1) + -1.68391260e+01*(a12*a1) + -2.16859678e+01*(u2*eta) + 2.21740921e+01*(u*a12) + -2.81478116e+02*(eta*a12) + 3.61439607e+01*(u*eta2) + -4.61790140e+02*(eta2*a1) + 3.90831711e+00*(u2*a1) + 1.92402755e+02*(u2*eta2) + 2.90639764e+02*(eta*a13) + -2.06058728e+01*(u2*a12) + 1.75385779e+03*(eta2*a12) + 1.59844231e+01*(u3*a1) + 9.17937044e+02*(eta3*a1) + -9.29295122e+01*(u*eta*a12) + -1.37300249e+01*(u*a13) + 7.50511732e+01*(u2*eta*a1) + -7.28499377e+01*(u*eta2*eta) + -1.75373704e+03*(eta2*a13) + -7.14470339e+02*(u2*eta2*a1) + 1.80064108e+01*(u4*eta) + -4.55949916e+02*(u2*eta2*eta) + -3.44368872e+03*(eta3*a12) + 6.95429997e+01*(u3*eta2) + -5.27128858e+00*(u4*a1) + 2.32963072e+01*(u2*a13) + 5.60176234e+01*(u*eta*a13) + -3.32611166e+01*(u3*a12) + -6.62343206e+01*(u3*eta*a1) + 1.52174113e+02*(u3*eta*a12) + -1.44616319e+02*(u2*eta*a13) + 1.95560265e+01*(u3*a13) + -1.52041823e+02*(u3*eta2*eta) + 1.32840766e+01*(u4*a12) + -1.42550037e+02*(u4*eta2) + -7.17019839e+01*(u3*eta2*a1) + 1.72190949e+03*(u2*eta3*a1) + 6.54102508e+02*(u2*eta2*a12) + 3.36581803e+03*(eta3*a13) + 2.36359528e+02*(u3*eta3*a1) + 3.23439391e+02*(u4*eta2*eta) + -4.78577957e+01*(u3*eta2*a12) + 3.59725652e+02*(u2*eta2*a13) + -2.12498513e+03*(u2*eta3*a12) + -9.84157144e+00*(u4*a13) + -8.02985078e+01*(u3*eta*a13);

	// Return answer
	return mu1;

} // END of MU1 fit implementation



// MU2 fit implementation 
// Header formatting for MU2
static double IMRPhenomXCP_MU2( double theta, double eta, double a1 ){ 

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
	double mu2;

	// Evaluate fit for this parameter
	mu2 = -1.70409578e+00*a1 + 1.60243369e+00*u + 1.84761916e-01 + 1.84964222e+00*(u*u) + -2.99104692e+01*eta2 + 2.09616550e+01*(eta*a1) + -3.55083307e+01*(u*eta) + 7.92033422e-01*(u2*u) + 9.10681555e-01*(a12*a1) + 9.65990084e+01*(eta2*eta) + -1.65032761e+01*(u2*eta) + 2.34415629e+00*(u*a12) + 2.10670680e+02*(u*eta2) + -2.11340748e+00*(u2*a1) + -1.18768540e+02*(eta2*a12) + -2.97180916e+00*(u3*a1) + -6.94119989e-01*(u3*u) + -1.74649516e+02*(eta3*a1) + -2.02308427e+00*(u*a13) + -6.42473415e+00*(u3*eta) + -3.65406917e+02*(u*eta2*eta) + 1.70366906e+02*(u2*eta2*a1) + 9.17471626e+00*(u4*eta) + 1.21484714e+02*(u2*eta2*eta) + 3.81665670e+02*(eta3*a12) + -2.74234838e+01*(u*eta2*a12) + 1.64521085e+01*(u3*eta2) + 2.63917471e+01*(u3*eta*a1) + -1.10223690e+02*(u*eta3*a1) + -8.37834886e+00*(u4*eta*a1) + -8.75724956e+00*(u3*eta*a12) + -4.89420075e+00*(u2*eta*a13) + 1.60248477e+00*(u3*a13) + 1.85319097e+02*(u*eta3*a12) + 1.94890802e+00*(u4*a12) + -1.69846511e+01*(u4*eta2) + -5.23935753e+01*(u3*eta2*a1) + -4.80159135e+02*(u2*eta3*a1);

	// Return answer
	return mu2;

} // END of MU2 fit implementation



// MU3 fit implementation 
// Header formatting for MU3
static double IMRPhenomXCP_MU3( double theta, double eta, double a1 ){ 

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
	mu3 = -8.45957335e+01*a1 + 1.26418699e+01*u + -2.00057814e+02*eta + 1.16889641e+01 + 1.56015413e+02*a12 + -8.72775398e+00*(u*u) + 1.06234462e+03*eta2 + 1.44676006e+03*(eta*a1) + -1.46178216e+02*(u*eta) + -7.60455014e+01*(u*a1) + -1.57812510e+01*(u2*u) + 6.67228485e+02*(u*eta*a1) + -8.87357838e+01*(a12*a1) + -1.79087097e+03*(eta2*eta) + 5.40054138e+02*(u*eta2) + 1.53289966e+02*(u2*eta) + 1.90541214e+02*(u*a12) + -2.65210596e+03*(eta*a12) + -7.82182612e+03*(eta2*a1) + 4.05305467e+01*(u2*a1) + -8.54854329e+02*(u2*eta2) + 1.48067473e+03*(eta*a13) + -1.45432864e+03*(u*eta2*a1) + 1.43375639e+04*(eta2*a12) + 9.66533014e+01*(u3*a1) + 1.35345037e+04*(eta3*a1) + -1.79652205e+03*(u*eta*a12) + -1.51442811e+02*(u*a13) + 1.46886949e+02*(u3*eta) + -7.04444357e+02*(u2*eta*a1) + -6.23273500e+02*(u*eta2*eta) + -7.94229828e+03*(eta2*a13) + 3.82686153e+03*(u2*eta2*a1) + 1.52030535e+03*(u2*eta2*eta) + -2.48920641e+04*(eta3*a12) + 5.78366609e-01*(u4*a1) + -3.38396769e+02*(u3*eta2) + -5.54448614e+01*(u2*a13) + 1.56807621e+03*(u*eta*a13) + -2.30063200e+02*(u3*a12) + -6.42826205e+02*(u3*eta*a1) + 4.92998740e+03*(u*eta2*a12) + 1.50237404e+03*(u3*eta*a12) + 9.80727029e+02*(u2*eta*a13) + 1.65543168e+02*(u3*a13) + -5.21222929e+03*(u*eta2*a13) + -3.21316681e+03*(u*eta3*a12) + -6.63489618e+03*(u2*eta3*a1) + 1.37516909e+04*(eta3*a13) + 4.14847297e+03*(u3*eta3*a1) + -5.39919574e+03*(u2*eta2*a13) + 5.49787535e+03*(u*eta3*a13) + -1.07766812e+03*(u3*eta*a13) + -5.47422593e+00*(u4*eta*a13) + -9.32142605e+03*(u3*eta3*a12) + 9.51221993e+03*(u2*eta3*a13) + 6.62675043e+03*(u3*eta3*a13);

	// Return answer
	return mu3;

} // END of MU3 fit implementation



// MU4 fit implementation 
// Header formatting for MU4
static double IMRPhenomXCP_MU4( double theta, double eta, double a1 ){ 

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
	mu4 = 9.93058882e-02*a1 + -8.67304804e-02*u + 1.14264485e-01*eta + -1.46752867e-02 + -2.93153971e-02*a12 + -3.14593593e-02*(u*u) + -3.12268295e-01*eta2 + -1.01346298e+00*(eta*a1) + 3.94313967e-01*(u*eta) + 7.80565839e-01*(u*a1) + 1.12853454e-01*(u2*u) + -3.52795716e+00*(u*eta*a1) + 5.22212311e-01*(u2*eta) + -1.90570685e+00*(u*a12) + 4.40968433e+00*(eta2*a1) + -1.99704827e-02*(u2*a1) + -2.86992240e+00*(u2*eta2) + 6.37979623e-02*(u2*a12) + -5.58283849e-01*(eta2*a12) + -9.95680887e-01*(u3*a1) + -5.78129910e+00*(eta3*a1) + 8.58049010e+00*(u*eta*a12) + 1.38617147e+00*(u*a13) + -4.19512789e-01*(u3*eta) + 1.72330366e+00*(eta2*a13) + 5.27735563e+00*(u2*eta2*eta) + -1.25767388e+00*(u3*eta2) + -6.37327421e+00*(u*eta*a13) + 2.31966865e+00*(u3*a12) + 4.80109841e+00*(u3*eta*a1) + -1.05576831e+01*(u3*eta*a12) + -2.15918317e-01*(u2*eta*a13) + -1.65464733e+00*(u3*a13) + 6.49394261e-01*(u*eta2*a13) + 3.46849893e+00*(u3*eta2*eta) + -1.10009417e+00*(u3*eta2*a1) + -4.15861523e+00*(eta3*a13) + -1.97535326e-01*(u4*eta2*eta) + -6.42400716e-02*(u4*a13) + 7.50775242e+00*(u3*eta*a13) + 2.75973944e-01*(u4*eta*a13);

	// Return answer
	return mu4;

} // END of MU4 fit implementation



// NU4 fit implementation 
// Header formatting for NU4
static double IMRPhenomXCP_NU4( double theta, double eta, double a1 ){ 

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
	nu4 = -2.38475790e-01*a1 + -6.73712192e-03*u + -6.30444504e-01*eta + 2.95331732e-02 + 4.83182187e-01*a12 + 3.84227872e+00*eta2 + 5.13928941e+00*(eta*a1) + 3.24371442e-02*(u*eta) + 4.91819161e-02*(u*a1) + 1.82163989e-02*(u2*u) + -2.31447178e-01*(u*eta*a1) + -2.95654967e-01*(a12*a1) + -7.27754529e+00*(eta2*eta) + 3.79294627e-01*(u2*eta) + -1.64947515e-01*(u*a12) + -1.08885374e+01*(eta*a12) + -3.16869008e+01*(eta2*a1) + 1.93391101e-01*(u2*a1) + -2.78538835e+00*(u2*eta2) + 7.04684079e+00*(eta*a13) + -2.42598719e-01*(u2*a12) + 6.87114677e+01*(eta2*a12) + -7.71862729e-02*(u3*a1) + -5.33772583e-02*(u3*u) + 6.01247773e+01*(eta3*a1) + 1.91840058e+00*(u*eta*a12) + 1.08890890e-01*(u*a13) + -2.02193630e-01*(u3*eta) + -6.82305924e+00*(u2*eta*a1) + -4.58883426e+01*(eta2*a13) + 4.39933188e+01*(u2*eta2*a1) + 7.04585349e-01*(u4*eta) + 5.09176601e+00*(u2*eta2*eta) + -1.32034286e+02*(eta3*a12) + 8.61084119e-01*(u3*eta2) + 1.20504823e-01*(u4*a1) + -1.72799881e+00*(u*eta*a13) + 9.64355323e-02*(u3*a12) + 3.82614205e-01*(u3*eta*a1) + -9.86990092e+00*(u*eta2*a12) + 1.35205363e+01*(u2*eta*a12) + -7.39079779e+00*(u2*eta*a13) + 1.07857869e+01*(u*eta2*a13) + 2.04348184e+01*(u*eta3*a12) + -1.52234416e+00*(u3*eta2*eta) + -3.39605145e-01*(u4*a12) + -3.98176008e+00*(u4*eta2) + -7.95463847e+01*(u2*eta3*a1) + -9.47649252e+01*(u2*eta2*a12) + 8.99705473e+01*(eta3*a13) + 7.94486925e+00*(u4*eta2*eta) + -2.18614989e+00*(u3*eta2*a12) + 5.93499185e+01*(u2*eta2*a13) + 1.79413233e+02*(u2*eta3*a12) + 3.57616914e-01*(u4*a13) + -2.30861823e+01*(u*eta3*a13) + -6.25243878e-01*(u3*eta*a13) + -7.09041523e+00*(u4*eta3*a1) + -1.34426226e+00*(u4*eta*a13) + 2.90666856e+00*(u3*eta2*a13) + 5.03651952e+00*(u4*eta2*a12) + -1.20177486e+02*(u2*eta3*a13);

	// Return answer
	return nu4;

} // END of NU4 fit implementation



// NU5 fit implementation 
// Header formatting for NU5
static double IMRPhenomXCP_NU5( double theta, double eta, double a1 ){ 

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
	nu5 = -6.41734659e-01*a1 + -5.25926227e-01*u + -1.32184645e+00*eta + 8.36855008e-02 + 1.23309799e+00*a12 + 6.51385836e-01*(u*u) + 6.89795594e+00*eta2 + 1.00694735e+01*(eta*a1) + 1.11209865e+01*(u*eta) + 2.95722155e+00*(u*a1) + -6.22943073e+01*(u*eta*a1) + -7.16354209e-01*(a12*a1) + -1.20666107e+01*(eta2*eta) + -7.37837999e+01*(u*eta2) + -1.59372498e+01*(u2*eta) + -5.13194446e+00*(u*a12) + 2.40388927e-01*(u2*u) + -1.91396627e+01*(eta*a12) + -5.28418001e+01*(eta2*a1) + -2.24657105e+00*(u2*a1) + 1.12667160e+02*(u2*eta2) + 1.08147951e+01*(eta*a13) + 1.61522808e+00*(u2*a12) + 4.10530523e+02*(u*eta2*a1) + 9.96148908e+01*(eta2*a12) + -4.15940373e-01*(u3*u) + 9.31774130e+01*(eta3*a1) + 1.08627460e+02*(u*eta*a12) + 3.09371839e+00*(u*a13) + -8.72736078e+00*(u3*eta) + 6.42100384e+01*(u2*eta*a1) + 1.54421670e+02*(u*eta2*eta) + -5.49413048e+01*(eta2*a13) + -4.79540875e+02*(u2*eta2*a1) + 1.54095905e+01*(u4*eta) + -2.44097046e+02*(u2*eta2*eta) + -1.75261962e+02*(eta3*a12) + 7.22074252e+01*(u3*eta2) + -1.17575384e+00*(u4*a1) + -6.38622728e+01*(u*eta*a13) + -2.12186430e+00*(u3*a12) + 3.06768126e+01*(u3*eta*a1) + -7.14693153e+02*(u*eta2*a12) + -8.53218144e+02*(u*eta3*a1) + -7.25267144e+01*(u2*eta*a12) + -2.82704062e+01*(u4*eta*a1) + -2.71731630e+01*(u3*eta*a12) + 2.43726540e+01*(u2*eta*a13) + 1.47115628e+00*(u3*a13) + 4.12629145e+02*(u*eta2*a13) + 1.47975911e+03*(u*eta3*a12) + -1.68984832e+02*(u3*eta2*eta) + 8.36997088e+00*(u4*a12) + -1.26176140e+02*(u4*eta2) + -3.23643722e+02*(u3*eta2*a1) + 1.06641094e+03*(u2*eta3*a1) + 6.04501449e+02*(u2*eta2*a12) + 9.52399801e+01*(eta3*a13) + 8.27200543e+02*(u3*eta3*a1) + 2.94419569e+02*(u4*eta2*eta) + 4.59082605e+02*(u3*eta2*a12) + -2.39773587e+02*(u2*eta2*a13) + -5.34015992e+01*(u4*eta*a12) + -1.40812658e+03*(u2*eta3*a12) + -7.20139921e+00*(u4*a13) + 3.82405298e+02*(u4*eta2*a1) + -8.42956944e+02*(u*eta3*a13) + 1.38055377e+01*(u3*eta*a13) + -1.04902648e+03*(u4*eta3*a1) + 7.22875200e+01*(u4*eta*a13) + -2.61854741e+02*(u3*eta2*a13) + -1.03192484e+02*(u4*eta2*a12) + -1.30528125e+03*(u3*eta3*a12) + 5.91563946e+02*(u2*eta3*a13) + -1.82836107e+02*(u4*eta2*a13) + 7.53691969e+02*(u3*eta3*a13) + 8.03707654e+02*(u4*eta3*a12);

	// Return answer
	return nu5;

} // END of NU5 fit implementation



// NU6 fit implementation 
// Header formatting for NU6
static double IMRPhenomXCP_NU6( double theta, double eta, double a1 ){ 

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
	nu6 = 7.54346094e-02*a1 + -6.90363459e-02*u + 3.45459387e-01*eta + -1.70858413e-02 + -2.24217344e+00*eta2 + -1.36422701e+00*(eta*a1) + 7.28306767e-01*(u*eta) + 4.62077882e-01*(u*a1) + 5.48410204e-02*(u2*u) + -3.33464294e+00*(u*eta*a1) + 1.65589304e-02*(a12*a1) + 4.37279019e+00*(eta2*eta) + -2.78844955e+00*(u*eta2) + -1.38448860e+00*(u*a12) + 9.12947099e+00*(eta2*a1) + 6.46293033e+00*(u*eta2*a1) + -1.97952260e+00*(eta2*a12) + -1.55375785e-01*(u3*a1) + -1.09887306e-02*(u3*u) + -1.85640424e+01*(eta3*a1) + 1.21293126e+01*(u*eta*a12) + 1.25090548e+00*(u*a13) + -6.90791646e-01*(u3*eta) + 3.55517898e+00*(u*eta2*eta) + 1.02676081e-01*(u4*eta) + -6.91999364e-02*(u2*eta2*eta) + 6.20258819e+00*(eta3*a12) + 3.87434660e+00*(u3*eta2) + -1.28077145e+01*(u*eta*a13) + 4.54069650e-01*(u3*a12) + -3.59878116e+01*(u*eta2*a12) + -4.83302392e-01*(u3*a13) + 4.66591592e+01*(u*eta2*a13) + 3.44130214e+01*(u*eta3*a12) + -7.54439972e+00*(u3*eta2*eta) + -2.44565199e-01*(u4*eta2) + 7.46415899e+00*(u3*eta3*a1) + -6.01793115e+00*(u3*eta2*a12) + -5.91077149e+01*(u*eta3*a13) + 1.72140101e+00*(u3*eta*a13);

	// Return answer
	return nu6;

} // END of NU6 fit implementation



// ZETA1 fit implementation 
// Header formatting for ZETA1
static double IMRPhenomXCP_ZETA1( double theta, double eta, double a1 ){ 

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
	zeta1 = -1.78667812e-04*a1 + -5.89149196e-04*eta + 3.46911473e-05 + 2.76239000e-04*a12 + -5.73601853e-06*(u*u) + 3.02206667e-03*eta2 + 2.93623309e-03*(eta*a1) + 3.27445083e-05*(u*eta) + -1.12735024e-05*(u2*u) + -9.28763206e-05*(a12*a1) + -4.90887846e-03*(eta2*eta) + 1.52840388e-05*(u2*eta) + -3.38090318e-03*(eta*a12) + -4.26768790e-04*(u*eta2) + -1.35843277e-02*(eta2*a1) + 9.12681311e-03*(eta2*a12) + 9.84069008e-05*(u3*u) + 1.89344571e-02*(eta3*a1) + 3.73633734e-05*(u*a13) + 9.23692445e-03*(eta2*a13) + -1.89902239e-03*(u4*eta) + -7.48446187e-04*(u4*a1) + 3.85206528e-05*(u2*a13) + -4.55185919e-04*(u*eta*a13) + 3.11884513e-04*(u3*a12) + 6.92757414e-04*(u3*eta*a1) + 8.41387501e-03*(u*eta3*a1) + 1.66410724e-02*(u4*eta*a1) + -7.61574137e-03*(u3*eta*a12) + -1.31238637e-04*(u2*eta*a13) + -4.63709069e-04*(u3*a13) + 3.42632684e-03*(u*eta2*a13) + -1.58463667e-02*(u*eta3*a12) + 3.01466555e-03*(u3*eta2*eta) + 1.92632248e-03*(u4*a12) + 1.19003299e-02*(u4*eta2) + -3.13021186e-02*(eta3*a13) + -2.63430067e-02*(u3*eta3*a1) + -2.42816716e-02*(u4*eta2*eta) + 3.71395861e-02*(u3*eta2*a12) + -4.50708506e-02*(u4*eta*a12) + -1.62098903e-03*(u4*a13) + -1.10716531e-01*(u4*eta2*a1) + 1.00720443e-02*(u3*eta*a13) + 2.31459511e-01*(u4*eta3*a1) + 3.74094800e-02*(u4*eta*a13) + -5.68342362e-02*(u3*eta2*a13) + 3.04626564e-01*(u4*eta2*a12) + -1.93617752e-02*(u3*eta3*a12) + -2.51793208e-01*(u4*eta2*a13) + 8.13686981e-02*(u3*eta3*a13) + -6.36818085e-01*(u4*eta3*a12) + 5.23150421e-01*(u4*eta3*a13);

	// Return answer
	return zeta1;

} // END of ZETA1 fit implementation



// ZETA2 fit implementation 
// Header formatting for ZETA2
static double IMRPhenomXCP_ZETA2( double theta, double eta, double a1 ){ 

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
	zeta2 = -2.09192114e+01*u + 2.64975006e+01*eta + -1.04106211e+00 + -2.45599248e+00*a12 + -7.20890468e+00*(u*u) + -1.60843235e+02*eta2 + 3.47579850e+02*(u*eta) + 3.87328647e+01*(u*a1) + 4.30581508e+01*(u2*u) + -4.47359210e+02*(u*eta*a1) + 3.17605685e+02*(eta2*eta) + -1.75639294e+03*(u*eta2) + 1.07200362e+02*(u2*eta) + -1.91704545e+01*(u*a12) + -1.78007566e+02*(eta*a12) + 1.90883157e+01*(u2*a1) + -3.87724315e+02*(u2*eta2) + 2.85206821e+02*(eta*a13) + 1.26904433e+01*(u2*a12) + 9.84669777e+02*(u*eta2*a1) + 1.69631578e+03*(eta2*a12) + -1.10687440e+02*(u3*a1) + -1.23950896e+02*(eta3*a1) + 1.06081234e+01*(u*a13) + -6.99282272e+02*(u3*eta) + -2.98310187e+02*(u2*eta*a1) + 2.91003323e+03*(u*eta2*eta) + -2.59605631e+03*(eta2*a13) + 7.00396575e+02*(u2*eta2*a1) + 4.37201837e+02*(u2*eta2*eta) + -3.57157720e+03*(eta3*a12) + 3.57588800e+03*(u3*eta2) + -1.91942319e+01*(u2*a13) + -1.43826079e+02*(u*eta*a13) + 8.54502161e+01*(u3*a12) + 2.37096613e+03*(u*eta2*a12) + 1.50500503e+03*(u3*eta*a1) + -7.78192223e+02*(u3*eta*a12) + -2.68873207e+01*(u3*a13) + -6.81601618e+03*(u*eta3*a12) + -6.10308006e+03*(u3*eta2*eta) + -6.23854971e+03*(u3*eta2*a1) + 9.61558534e+02*(u2*eta2*a12) + 5.71368386e+03*(eta3*a13) + 9.51247996e+03*(u3*eta3*a1) + 1.00347735e+03*(u3*eta2*a12) + 1.00998648e+02*(u2*eta2*a13) + -3.32455506e+03*(u2*eta3*a12) + 1.06553224e+00*(u4*a13) + 7.87714979e+02*(u*eta3*a13) + 2.14635462e+02*(u3*eta*a13);

	// Return answer
	return zeta2;

} // END of ZETA2 fit implementation

