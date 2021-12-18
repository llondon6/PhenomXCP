

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
	mu1 = 9.51107662e-01 + 1.00075474e+00*u + -7.74201014e+00*a1 + -2.06264435e+01*eta + -3.33555289e+00*(u*a1) + -1.66951302e+00*(u*u) + -1.49540922e+01*(u*eta) + 1.91212764e+01*a12 + 1.66209959e+02*(eta*a1) + 1.40156104e+02*eta2 + 1.57851783e+01*(u*eta*a1) + -1.56511774e+01*(a12*a1) + -8.56254432e-01*(u2*u) + 1.91354444e+00*(u2*a1) + 5.53066651e+00*(u*a12) + 2.87287969e+01*(u2*eta) + -2.90115390e+02*(eta2*eta) + 9.32767138e+01*(u*eta2) + -1.14241108e+03*(eta2*a1) + -3.89932807e+02*(eta*a12) + -3.23332437e+00*(u*a13) + -1.90005385e+02*(u*eta2*eta) + 2.81469721e+00*(u3*eta) + 2.41300458e+03*(eta3*a1) + -4.20089800e+00*(u2*a12) + -2.73514537e+01*(u2*eta*a1) + -1.47177665e+02*(u2*eta2) + 2.64500684e+03*(eta2*a12) + -8.30578175e+01*(u*eta2*a1) + 5.29729130e+00*(u3*a1) + 3.03441812e+02*(eta*a13) + 5.40396332e-01*(u3*u) + 2.05151191e+02*(u*eta3*a1) + -1.41101137e+01*(u3*eta*a1) + 2.54394318e+02*(u2*eta2*eta) + -8.94047716e+00*(u4*eta) + 5.96139582e+00*(u2*a13) + -5.57946745e+03*(eta3*a12) + -2.25492236e+01*(u*eta2*a12) + 1.17382688e+00*(u4*a1) + -1.98515651e+03*(eta2*a13) + 1.61407692e+02*(u2*eta2*a1) + -8.91929273e+00*(u3*a12) + 7.94789239e+01*(u2*eta2*a12) + 2.37357338e+01*(u4*eta2) + 1.36771212e+01*(u3*eta*a12) + 4.09599706e+03*(eta3*a13) + -2.20796082e+01*(u2*eta*a13) + -4.06974656e+02*(u2*eta3*a1) + -1.27053330e+00*(u4*a12) + 4.44605380e+00*(u3*a13);

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
	mu2 = 6.28985741e-01 + 2.16418274e+00*u + -1.55285059e+00*a1 + -7.22291053e+00*eta + -1.57699578e-01*(u*u) + -4.36341565e+01*(u*eta) + 8.94083778e-01*a12 + 2.61249004e+01*(eta*a1) + 3.53858282e+01*(u*eta*a1) + -1.34905694e-01*(u2*a1) + -6.48440808e+00*(u*a12) + 5.63886907e+01*(eta2*eta) + 2.40487241e+02*(u*eta2) + -5.89001179e+01*(eta2*a1) + -1.15228674e+01*(eta*a12) + 5.44595979e+00*(u*a13) + -3.92580889e+02*(u*eta2*eta) + 2.42844026e+01*(u*eta*a12) + -2.32903880e+02*(u*eta2*a1) + 2.73576987e+02*(u*eta3*a1) + -1.76908649e+00*(u3*eta*a1) + 5.04913323e+01*(u*eta2*a12) + 9.31334740e-02*(u4*a1) + 1.59201726e+01*(eta2*a13) + 3.12248315e+00*(u2*eta2*a1) + 4.50771370e+00*(u3*eta2) + -2.77153711e+01*(u*eta*a13);

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
	mu3 = -2.87256313e-01 + -7.52514926e+00*a1 + 1.97617975e+01*(u*a1) + 3.83665868e+01*a12 + 1.56047730e+02*(eta*a1) + -3.04537559e+02*(u*eta*a1) + -4.43556395e+01*(a12*a1) + -1.97638841e-01*(u2*u) + -1.71984289e+00*(u2*a1) + -4.94259080e+01*(u*a12) + -9.81594206e+02*(eta2*a1) + -6.62149243e+02*(eta*a12) + 2.05184040e+03*(eta3*a1) + 7.66078145e+02*(u*eta*a12) + 4.99525950e+00*(u2*a12) + 3.76481835e+03*(eta2*a12) + 1.49403536e+03*(u*eta2*a1) + -5.34462667e+00*(u3*a1) + 7.23169321e+02*(eta*a13) + -1.38146382e-01*(u3*u) + -2.32839616e+03*(u*eta3*a1) + -7.24193434e+03*(eta3*a12) + -4.03437868e+03*(u*eta2*a12) + 9.04937746e-01*(u4*a1) + -3.88237900e+03*(eta2*a13) + 3.79974093e+01*(u*eta*a13) + 2.78065746e+01*(u3*a12) + 6.95541906e+03*(u*eta3*a12) + -1.47727655e+02*(u2*eta2*a12) + -3.51226242e+02*(u3*eta*a12) + 6.99844955e+03*(eta3*a13) + 6.39263579e+02*(u3*eta2*a1) + -5.76423849e+02*(u*eta3*a13) + -2.20366549e+03*(u3*eta3*a1) + 4.89933453e+00*(u4*eta*a12) + 4.81658147e+02*(u2*eta3*a12) + -4.13470149e+00*(u4*a13) + 9.67212619e+02*(u3*eta2*a12);

	// Return answer
	return mu3;

} // END of MU3 fit implementation



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
	nu4 = -4.28479112e-03 + -3.37230238e-03*u + 1.32324965e-02*eta + 2.60736697e-02*(u*a1) + 8.18004031e-02*(u*u) + 6.14360552e-01*(eta*a1) + 1.71667392e-02*(u2*u) + -2.94468850e-01*(u2*a1) + -8.06622007e-02*(u*a12) + -1.15034008e+00*(u2*eta) + -4.72430034e+00*(eta2*a1) + -1.74027035e+00*(eta*a12) + 2.73838185e-01*(u*eta2*eta) + -2.24850038e-01*(u3*eta) + 9.09006960e+00*(eta3*a1) + 7.82846860e-01*(u*eta*a12) + 6.82004537e-01*(u2*a12) + 2.06710870e+00*(u2*eta*a1) + 6.37090777e+00*(u2*eta2) + 1.42351342e+01*(eta2*a12) + -7.67956625e-02*(u3*a1) + 1.45810104e+00*(eta*a13) + -8.20536003e-02*(u3*u) + -1.98000483e+00*(u*eta3*a1) + 6.65311266e-01*(u3*eta*a1) + -1.25590544e+01*(u2*eta2*eta) + 1.22807637e+00*(u4*eta) + -6.00081077e-01*(u2*a13) + -2.88813293e+01*(eta3*a12) + -6.59510161e+00*(u*eta2*a12) + 1.98887521e-01*(u4*a1) + -1.26071006e+01*(eta2*a13) + -8.31878716e+00*(u2*eta2*a1) + 1.23714121e+00*(u3*eta2) + -2.92244385e+00*(u2*eta*a12) + 6.07956740e-02*(u3*a12) + 1.96004740e+01*(u*eta3*a12) + -7.25958506e+00*(u4*eta2) + 2.68878104e+01*(eta3*a13) + 3.12363172e+00*(u2*eta*a13) + -1.10366396e+00*(u4*eta*a1) + 2.06376786e+01*(u2*eta3*a1) + -4.02066153e-01*(u4*a12) + -2.54429990e+00*(u3*eta2*eta) + 8.87787641e-02*(u3*a13) + 3.36704958e+00*(u*eta2*a13) + -3.43269395e+00*(u3*eta2*a1) + -1.37803670e+01*(u*eta3*a13) + 8.76242146e+00*(u3*eta3*a1) + -1.80559221e+00*(u3*eta*a13) + 4.34261545e-01*(u4*a13) + 6.60087367e+00*(u4*eta2*a1) + 1.48011949e+01*(u4*eta2*eta) + 2.39798288e+00*(u3*eta2*a12) + -9.78761288e+00*(u2*eta3*a13) + 7.20031463e+00*(u4*eta2*a12) + -1.51855567e+01*(u3*eta3*a12) + -1.84409806e+00*(u4*eta*a13) + -2.29717911e+01*(u4*eta3*a1) + 6.07566779e+00*(u3*eta2*a13);

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
	nu5 = 1.18064679e-02 + -4.84505107e-01*u + -1.29000653e-01*a1 + 2.53121977e+00*(u*a1) + 8.98822868e-01*(u*u) + 1.03182098e+01*(u*eta) + 2.49346199e-01*a12 + 6.85817399e-01*(eta*a1) + -5.91743254e-01*eta2 + -5.41246372e+01*(u*eta*a1) + -1.59977292e-01*(a12*a1) + 2.00615192e-01*(u2*u) + -3.99577498e+00*(u2*a1) + -3.88815714e+00*(u*a12) + -2.00484921e+01*(u2*eta) + 1.33781657e+00*(eta2*eta) + -6.92881307e+01*(u*eta2) + -1.27053985e+00*(eta*a12) + 2.00095764e+00*(u*a13) + 1.46636090e+02*(u*eta2*eta) + -7.90707202e+00*(u3*eta) + -7.99795659e-01*(eta3*a1) + 8.49340677e+01*(u*eta*a12) + 5.11875630e+00*(u2*a12) + 9.27374116e+01*(u2*eta*a1) + 1.32870602e+02*(u2*eta2) + 3.64787970e+02*(u*eta2*a1) + 4.47042887e-01*(u3*a1) + 8.21065640e-01*(eta*a13) + -5.83439762e-01*(u3*u) + -7.73743497e+02*(u*eta3*a1) + 2.17222159e+01*(u3*eta*a1) + -2.74644492e+02*(u2*eta2*eta) + 1.78428245e+01*(u4*eta) + -2.27447328e+00*(u2*a13) + -5.82022428e+02*(u*eta2*a12) + -6.14836983e+02*(u2*eta2*a1) + 6.75448152e+01*(u3*eta2) + -4.31777582e+01*(u*eta*a13) + -1.28216526e+02*(u2*eta*a12) + -3.50222173e+00*(u3*a12) + 1.24869728e+03*(u*eta3*a12) + 8.55242991e+02*(u2*eta2*a12) + -1.35531627e+02*(u4*eta2) + 5.99714345e+01*(u2*eta*a13) + -4.46106753e+01*(u4*eta*a1) + 1.26038211e+03*(u2*eta3*a1) + 5.95301239e+00*(u4*a12) + -1.60920314e+02*(u3*eta2*eta) + 2.74057950e+00*(u3*a13) + 2.96891067e+02*(u*eta2*a13) + -2.72899440e+02*(u3*eta2*a1) + -6.41118941e+02*(u*eta3*a13) + 7.38954740e+02*(u3*eta3*a1) + -1.08088185e+01*(u3*eta*a13) + -2.14303614e+01*(u4*eta*a12) + -1.73764563e+03*(u2*eta3*a12) + -5.47588107e+00*(u4*a13) + 4.37869531e+02*(u4*eta2*a1) + -3.94246752e+02*(u2*eta2*a13) + 3.03165513e+02*(u4*eta2*eta) + 3.05346895e+02*(u3*eta2*a12) + 7.80687226e+02*(u2*eta3*a13) + -1.94685885e+02*(u4*eta2*a12) + -1.03696087e+03*(u3*eta3*a12) + 4.92987513e+01*(u4*eta*a13) + -1.08086002e+03*(u4*eta3*a1) + -1.22988705e+02*(u3*eta2*a13) + 8.03745639e+02*(u4*eta3*a12) + -1.16895404e+02*(u4*eta2*a13) + 5.10963614e+02*(u3*eta3*a13);

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
	nu6 = -5.63268247e-03 + -2.44264325e-02*u + 5.01441815e-02*a1 + 2.95284195e-01*(u*a1) + -1.00238243e-01*a12 + 5.64218668e-01*eta2 + -1.16322582e+00*(u*eta*a1) + 2.09096053e-01*(a12*a1) + 5.32038207e-03*(u2*a1) + -9.83609927e-01*(u*a12) + -2.11652701e+00*(eta2*eta) + 2.07926500e+00*(u*eta2) + -4.00525565e+00*(eta2*a1) + 9.27680511e-01*(u*a13) + -6.96265681e+00*(u*eta2*eta) + 1.37957306e+01*(eta3*a1) + 7.84521737e+00*(u*eta*a12) + -4.78969881e-01*(u2*eta2) + 9.06993644e+00*(eta2*a12) + -1.22582736e+01*(u*eta2*a1) + -2.38669966e+00*(eta*a13) + -1.43958108e-02*(u3*u) + 5.00958156e+01*(u*eta3*a1) + 2.40011354e+00*(u2*eta2*eta) + 2.35376119e-01*(u4*eta) + -3.12892057e+01*(eta3*a12) + -2.19386912e-02*(u4*a1) + 6.40740888e+00*(eta2*a13) + -9.85877416e+00*(u*eta*a13) + -6.54690316e+01*(u*eta3*a12) + -8.75368007e-01*(u4*eta2) + 1.59999005e-01*(u4*eta*a1) + -1.38010328e+00*(u2*eta3*a1) + 2.50633047e+01*(u*eta2*a13);

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
	zeta1 = 2.81426565e-05 + 2.93827232e-04*u + -1.53059805e-04*a1 + -4.85373582e-04*eta + -1.22197974e-03*(u*a1) + -3.17541676e-05*(u*u) + -5.75290725e-03*(u*eta) + 3.83083033e-04*a12 + 2.52558978e-03*(eta*a1) + 2.53977801e-03*eta2 + 2.41577262e-02*(u*eta*a1) + -2.52579162e-04*(a12*a1) + -4.68055198e-04*(u2*u) + 1.49214749e-03*(u*a12) + 6.56402043e-04*(u2*eta) + -4.22830271e-03*(eta2*eta) + 3.47435227e-02*(u*eta2) + -1.17888263e-02*(eta2*a1) + -5.37716840e-03*(eta*a12) + -5.13827682e-04*(u*a13) + -6.71604656e-02*(u*eta2*eta) + 8.92212558e-03*(u3*eta) + 1.68442974e-02*(eta3*a1) + -2.96119984e-02*(u*eta*a12) + -2.29972816e-04*(u2*a12) + -4.43738677e-03*(u2*eta2) + 2.12752478e-02*(eta2*a12) + -1.47190056e-01*(u*eta2*a1) + 1.78090960e-03*(u3*a1) + 2.89922255e-03*(eta*a13) + 4.88202286e-05*(u3*u) + 2.89845339e-01*(u*eta3*a1) + -3.40234385e-02*(u3*eta*a1) + 9.09162168e-03*(u2*eta2*eta) + -9.18091471e-04*(u4*eta) + 1.83017626e-04*(u2*a13) + -2.40522023e-02*(eta3*a12) + 1.80610873e-01*(u*eta2*a12) + 6.32201036e-05*(u4*a1) + -7.51498007e-03*(eta2*a13) + -5.37512513e-02*(u3*eta2) + 1.04783375e-02*(u*eta*a13) + 1.95332471e-03*(u2*eta*a12) + -1.55551723e-03*(u3*a12) + -3.61115672e-01*(u*eta3*a12) + 5.44924738e-03*(u4*eta2) + 2.85158012e-02*(u3*eta*a12) + -5.29260195e-04*(u4*eta*a1) + -8.40587671e-05*(u4*a12) + 1.04816201e-01*(u3*eta2*eta) + -6.27817362e-02*(u*eta2*a13) + 2.08002006e-01*(u3*eta2*a1) + 1.25479115e-01*(u*eta3*a13) + -4.17808708e-01*(u3*eta3*a1) + 1.41398961e-03*(u3*eta*a13) + -1.57341671e-02*(u2*eta3*a12) + 2.28567905e-03*(u4*eta2*a1) + -1.43834554e-02*(u2*eta2*a13) + -1.06617064e-02*(u4*eta2*eta) + -1.76050055e-01*(u3*eta2*a12) + 4.54266838e-02*(u2*eta3*a13) + -1.58539452e-03*(u4*eta2*a12) + 3.75075082e-01*(u3*eta3*a12) + 4.93017239e-04*(u4*eta*a13) + -9.13890494e-03*(u3*eta2*a13);

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
	zeta2 = -9.97328890e+00 + -1.33992899e+01*u + 4.14628907e+01*a1 + 1.61580632e+02*eta + 1.15203728e+01*(u*a1) + 3.95212363e+01*(u*u) + 2.25251022e+02*(u*eta) + -8.21468689e+01*a12 + -5.40629810e+02*(eta*a1) + -8.07314221e+02*eta2 + -9.65252328e+01*(u*eta*a1) + 5.41068438e+01*(a12*a1) + 2.82669341e+01*(u2*u) + -7.96946979e+01*(u2*a1) + -6.39068515e+02*(u2*eta) + 1.32491298e+03*(eta2*eta) + -1.16791680e+03*(u*eta2) + 2.10842064e+03*(eta2*a1) + 8.11516465e+02*(eta*a12) + 1.71353870e+01*(u*a13) + 2.03789994e+03*(u*eta2*eta) + -4.57100279e+02*(u3*eta) + -2.63845882e+03*(eta3*a1) + 1.01222235e+02*(u2*a12) + 6.98337959e+02*(u2*eta*a1) + 3.74725375e+03*(u2*eta2) + -1.83605185e+03*(eta2*a12) + -2.21671629e+01*(u3*a1) + -4.05917402e+02*(eta*a13) + -6.38360022e+01*(u3*u) + 1.29037259e+02*(u3*eta*a1) + -7.39199116e+03*(u2*eta2*eta) + 1.12203899e+03*(u4*eta) + -7.10752348e+01*(u2*a13) + 1.90011863e+02*(u4*a1) + -4.13723624e+03*(u2*eta2*a1) + 2.38557422e+03*(u3*eta2) + -4.99090319e+02*(u*eta*a13) + -9.07236987e+01*(u3*a12) + 1.60900949e+03*(u*eta3*a12) + -6.76990552e+03*(u4*eta2) + 1.91012375e+03*(u3*eta*a12) + 2.80138542e+03*(eta3*a13) + -3.19144815e+03*(u4*eta*a1) + 1.15353653e+04*(u2*eta3*a1) + -2.77803370e+02*(u4*a12) + -4.28204154e+03*(u3*eta2*eta) + 8.31805059e+01*(u3*a13) + 3.57826509e+03*(u*eta2*a13) + -8.37226720e+03*(u*eta3*a13) + 1.01485804e+03*(u3*eta3*a1) + -1.46268178e+03*(u3*eta*a13) + 4.72128616e+03*(u4*eta*a12) + -9.23333092e+03*(u2*eta3*a12) + 1.80968602e+02*(u4*a13) + 2.13977199e+04*(u4*eta2*a1) + 1.35612722e+04*(u4*eta2*eta) + -1.08369385e+04*(u3*eta2*a12) + 6.52795988e+03*(u2*eta3*a13) + -3.54686913e+04*(u4*eta2*a12) + 1.53976735e+04*(u3*eta3*a12) + -3.22681275e+03*(u4*eta*a13) + -4.93675042e+04*(u4*eta3*a1) + 7.41235808e+03*(u3*eta2*a13) + 9.01842121e+04*(u4*eta3*a12) + 2.51161166e+04*(u4*eta2*a13) + -9.72451815e+03*(u3*eta3*a13) + -6.45004999e+04*(u4*eta3*a13);

	// Return answer
	return zeta2;

} // END of ZETA2 fit implementation

