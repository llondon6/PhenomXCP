

#
def generate_model_params(theta,eta,a1):

	'''
	Hola, soy un codigo escribido por "4b_document_fits.py". 
	~lionel.london@kcl.ac.uk/pilondon2@gmail.com
	'''  

	# Import usefuls
	from numpy import cos, sqrt

	# Preliminaries
	u = cos(theta)
	u2 = u*u
	u3 = u2*u
	u4 = u3*u
	eta2 = eta*eta
	eta3 = eta2*eta
	delta = sqrt(1-4*eta)
	delta2 = delta*delta
	delta3 = delta2*delta

	# mu1, 39 terms
	mu1_l2m2 = -1.28496412e+00*(u) + 5.83270450e-01*(eta) + 2.15682581e+00*(a1) + -1.32517396e-01 + 2.45165736e+01*(u*eta) + -1.51789734e+00*(u*u) + -3.11460675e+01*(eta*a1) + -4.20022288e+00*(a1*a1) + 6.62848134e+00*(u*a1) + -3.22824879e-01*(u2*u) + 6.76804386e+01*(eta*a1*a1) + 1.66988267e+02*(eta2*a1) + -1.43054987e+02*(u*eta*a1) + -6.74726032e+00*(u*a1*a1) + -1.38670435e+02*(u*eta*eta) + 3.41320483e+01*(u2*eta) + 3.18893622e+00*(u2*a1*a1) + -3.65844597e+02*(eta2*a1*a1) + 2.34122411e+02*(u*eta2*eta) + 1.54354983e+02*(u*eta*a1*a1) + 8.79097605e+02*(u*eta2*a1) + -2.75717314e+01*(u2*eta*a1) + 2.41453197e+00*(u3*a1) + -3.05513902e+02*(eta3*a1) + -2.18591605e+02*(u2*eta*eta) + 4.44791686e+02*(u2*eta2*eta) + -3.08756278e+01*(u2*eta*a1*a1) + 2.44992156e+02*(u2*eta2*a1) + -9.73231053e+02*(u*eta2*a1*a1) + -6.71251357e+00*(u3*eta*a1) + -1.60164134e+03*(u*eta3*a1) + 6.49350493e+02*(eta3*a1*a1) + -3.31359164e+00*(u3*a1*a1) + -6.06558133e+02*(u2*eta3*a1) + 1.88756234e+01*(u3*eta*a1*a1) + 8.56777572e+01*(u2*eta2*a1*a1) + -3.75231921e+01*(u3*eta2*a1) + 4.47397009e+01*(u3*eta2*eta) + 1.81264655e+03*(u*eta3*a1*a1)

	# mu2, 30 terms
	mu2_l2m2 = -7.54164053e+00*(u) + -1.91516385e-01 + 1.55039994e+02*(u*eta) + -4.18061332e+00*(u*u) + 3.20399024e+00*(eta*a1) + -8.53142828e-01*(a1*a1) + 5.43122209e+01*(u*a1) + -7.16400097e+00*(u2*u) + -1.07378584e+03*(u*eta*a1) + -6.36902248e+01*(u*a1*a1) + -9.60997353e+02*(u*eta*eta) + 7.96409152e+00*(u2*a1) + 6.68912601e+01*(u2*eta) + 1.29899349e+02*(u3*eta) + -9.16024335e+00*(u2*a1*a1) + 1.83641227e+03*(u*eta2*eta) + 1.20379699e+03*(u*eta*a1*a1) + 6.56170457e+03*(u*eta2*a1) + -4.91323229e+01*(u2*eta*a1) + -3.96440207e+02*(u2*eta*eta) + 8.65119085e+02*(u2*eta2*eta) + -7.73634153e+02*(u3*eta*eta) + 7.50131222e+01*(u2*eta*a1*a1) + -7.19880147e+03*(u*eta2*a1*a1) + -1.25208717e+04*(u*eta3*a1) + 1.67610628e+00*(u3*a1*a1) + -9.56493504e+01*(u2*eta2*a1*a1) + -3.72877742e+01*(u3*eta2*a1) + 1.53086671e+03*(u3*eta2*eta) + 1.35954995e+04*(u*eta3*a1*a1)

	# mu3, 31 terms
	mu3_l2m2 = 1.76863919e-02*(u) + 6.66706654e-03*(eta) + 6.44905054e-02*(a1) + -3.19807241e-03 + -3.79911024e-03*(u*u) + -8.53135817e-01*(eta*a1) + -1.16891300e-01*(u*a1) + 4.35962056e+00*(eta2*a1) + 5.72935240e-01*(u*eta*a1) + 1.52106345e-01*(u*a1*a1) + -7.60172185e-01*(u*eta*eta) + -5.50827313e-02*(u2*a1) + -2.07814905e-01*(u3*eta) + 1.87056474e+00*(u*eta2*eta) + -9.01179170e-01*(u*eta*a1*a1) + 1.05548067e+00*(u2*eta*a1) + 8.71133458e-02*(u3*a1) + -7.79327670e+00*(eta3*a1) + -1.48820532e-01*(u2*eta2*eta) + 1.36581619e+00*(u3*eta*eta) + -1.30087246e-01*(u2*eta*a1*a1) + -5.50044904e+00*(u2*eta2*a1) + 1.08773010e+00*(u*eta2*a1*a1) + -1.44914438e+00*(u*eta3*a1) + -1.26442295e-01*(u3*a1*a1) + -2.52172811e-01*(u3*eta*a1) + 9.76511974e+00*(u2*eta3*a1) + 5.27227564e-01*(u3*eta*a1*a1) + 3.05066950e-01*(u2*eta2*a1*a1) + -4.31214126e-01*(u3*eta2*a1) + -2.17495557e+00*(u3*eta2*eta)

	# nu4, 38 terms
	nu4_l2m2 = 1.57271240e-02*(u) + -4.62111242e-02*(eta) + -3.09964085e-02*(a1) + 4.43403755e-03 + -3.23948966e-01*(u*eta) + 1.07505792e-01*(eta*eta) + -7.24626894e-03*(u*u) + 3.71166548e-01*(eta*a1) + 3.19607415e-02*(a1*a1) + -7.69276357e-02*(u*a1) + -3.39932512e-01*(eta*a1*a1) + -1.30763107e+00*(eta2*a1) + 1.72556494e+00*(u*eta*a1) + 4.78116769e-02*(u*a1*a1) + 1.88646735e+00*(u*eta*eta) + 7.52390889e-02*(u2*a1) + -8.20112905e-02*(u2*a1*a1) + 8.20493332e-01*(eta2*a1*a1) + -3.33109320e+00*(u*eta2*eta) + -1.40806912e+00*(u*eta*a1*a1) + -1.06030817e+01*(u*eta2*a1) + -7.33275104e-01*(u2*eta*a1) + -7.69252474e-03*(u3*a1) + 1.37500320e+00*(eta3*a1) + 5.08362676e-01*(u2*eta*eta) + -1.65464279e+00*(u2*eta2*eta) + 2.19073512e-01*(u3*eta*eta) + 8.15309343e-01*(u2*eta*a1*a1) + 1.77477691e+00*(u2*eta2*a1) + 9.38718692e+00*(u*eta2*a1*a1) + -9.15076848e-02*(u3*eta*a1) + 1.96269582e+01*(u*eta3*a1) + 3.38509103e-02*(u3*a1*a1) + -1.51449688e-01*(u3*eta*a1*a1) + -1.97128750e+00*(u2*eta2*a1*a1) + 5.84743480e-01*(u3*eta2*a1) + -1.00299232e+00*(u3*eta2*eta) + -1.80393010e+01*(u*eta3*a1*a1)

	# nu5, 39 terms
	nu5_l2m2 = -3.05660691e-02*(u) + 1.31352981e-01*(eta) + -8.35979819e-03 + 4.67442854e-01*(u*eta) + -7.74561025e-01*(eta*eta) + -3.97631283e-02*(u*u) + -4.56835905e-01*(eta*a1) + -6.96499953e-02*(a1*a1) + 8.92504364e-02*(u*a1) + 9.75793445e-03*(u2*u) + 1.34689385e+00*(eta*a1*a1) + 3.51856126e+00*(eta2*a1) + -8.95585605e-01*(u*eta*a1) + -8.46986502e-02*(u*a1*a1) + 1.34264675e+00*(eta2*eta) + -2.38128939e+00*(u*eta*eta) + 1.30075523e-01*(u2*a1) + 6.68980911e-01*(u2*eta) + -1.04810973e-01*(u3*eta) + -5.38280413e-02*(u2*a1*a1) + -8.96445798e+00*(eta2*a1*a1) + 3.92615903e+00*(u*eta2*eta) + 3.42694043e-01*(u*eta*a1*a1) + 3.45498185e+00*(u*eta2*a1) + -1.69354010e+00*(u2*eta*a1) + -4.06966260e-02*(u3*a1) + -6.33304481e+00*(eta3*a1) + -4.23329268e+00*(u2*eta*eta) + 9.16060017e+00*(u2*eta2*eta) + 2.94977443e-01*(u2*eta*a1*a1) + 9.66621695e+00*(u2*eta2*a1) + 4.40717414e-01*(u3*eta*a1) + -5.08260069e+00*(u*eta3*a1) + 1.83156510e+01*(eta3*a1*a1) + 3.99054227e-02*(u3*a1*a1) + -2.13514623e+01*(u2*eta3*a1) + -1.28125539e-01*(u3*eta*a1*a1) + -1.36393119e+00*(u3*eta2*a1) + 1.46376607e+00*(u3*eta2*eta)

	# nu6, 40 terms
	nu6_l2m2 = 3.65688008e-02*(u) + -3.11667273e-01*(eta) + -7.91833101e-02*(a1) + 1.58221518e-02 + -5.36918575e-01*(u*eta) + 1.90135838e+00*(eta*eta) + 1.62841708e+00*(eta*a1) + 1.54991878e-01*(a1*a1) + -3.54229177e-01*(u*a1) + 4.74027862e-02*(u2*u) + -2.74630824e+00*(eta*a1*a1) + -1.00492902e+01*(eta2*a1) + 5.78641907e+00*(u*eta*a1) + 4.50558046e-01*(u*a1*a1) + -3.78788439e+00*(eta2*eta) + 3.14680451e+00*(u*eta*eta) + -3.41885953e-02*(u2*a1) + -3.99761812e-02*(u2*eta) + -9.13937175e-01*(u3*eta) + 4.96012269e-02*(u2*a1*a1) + 1.56527334e+01*(eta2*a1*a1) + -6.31613054e+00*(u*eta2*eta) + -7.23995151e+00*(u*eta*a1*a1) + -3.29305662e+01*(u*eta2*a1) + 4.17428334e-01*(u2*eta*a1) + 1.96576021e+01*(eta3*a1) + 2.82009876e-01*(u2*eta2*eta) + 4.66551632e+00*(u3*eta*eta) + -7.02233835e-01*(u2*eta*a1*a1) + 3.96779936e+01*(u*eta2*a1*a1) + 6.15962682e+01*(u*eta3*a1) + -2.88629987e+01*(eta3*a1*a1) + -7.15296358e-02*(u3*a1*a1) + 5.65996951e-01*(u3*eta*a1) + -3.67931430e+00*(u2*eta3*a1) + 2.41635196e-01*(u3*eta*a1*a1) + 1.86112931e+00*(u2*eta2*a1*a1) + -2.13232667e+00*(u3*eta2*a1) + -7.16284104e+00*(u3*eta2*eta) + -7.13316197e+01*(u*eta3*a1*a1)

	# zeta1, 40 terms
	zeta1_l2m2 = -2.61569600e-05*(u) + 7.16021495e-05*(eta) + 9.00838223e-05*(a1) + -6.53153489e-06 + 6.07770789e-04*(u*eta) + -1.60054005e-04*(eta*eta) + 2.91218080e-05*(u*u) + -1.11726499e-03*(eta*a1) + -6.89971710e-05*(a1*a1) + 2.78523863e-04*(u*a1) + -4.10195362e-05*(u2*u) + 6.74949940e-04*(eta*a1*a1) + 4.73841883e-03*(eta2*a1) + -5.75897175e-03*(u*eta*a1) + -3.28502043e-04*(u*a1*a1) + -4.31141469e-03*(u*eta*eta) + -2.38007338e-04*(u2*a1) + -1.70136568e-04*(u2*eta) + 7.24517901e-04*(u3*eta) + 1.67836942e-04*(u2*a1*a1) + -1.54635069e-03*(eta2*a1*a1) + 9.22181273e-03*(u*eta2*eta) + 6.90036189e-03*(u*eta*a1*a1) + 3.68125030e-02*(u*eta2*a1) + 2.67585238e-03*(u2*eta*a1) + 2.99354335e-05*(u3*a1) + -7.09563751e-03*(eta3*a1) + 1.13142583e-03*(u2*eta2*eta) + -3.61709857e-03*(u3*eta*eta) + -1.30201508e-03*(u2*eta*a1*a1) + -1.13454738e-02*(u2*eta2*a1) + -4.28832285e-02*(u*eta2*a1*a1) + -7.28507021e-02*(u*eta3*a1) + -6.20621847e-04*(u3*eta*a1) + 1.73192253e-02*(u2*eta3*a1) + 8.49508560e-05*(u3*eta*a1*a1) + 2.60799878e-03*(u2*eta2*a1*a1) + 1.69019440e-03*(u3*eta2*a1) + 5.79768042e-03*(u3*eta2*eta) + 8.19251950e-02*(u*eta3*a1*a1)

	# zeta2, 42 terms
	zeta2_l2m2 = 3.50956181e+00*(u) + 9.60136924e+00*(eta) + -1.78400964e+01*(a1) + 6.82931365e-01 + -1.00242982e+02*(u*eta) + -1.25379881e+02*(eta*eta) + -4.03092002e+00*(u*u) + 1.80805286e+02*(eta*a1) + 1.82857576e+01*(a1*a1) + -3.30100396e+01*(u*a1) + 6.66779765e+00*(u2*u) + -1.74965293e+02*(eta*a1*a1) + -5.51384118e+02*(eta2*a1) + 7.98366884e+02*(u*eta*a1) + 3.68856564e+01*(u*a1*a1) + 3.30699137e+02*(eta2*eta) + 7.04004154e+02*(u*eta*eta) + 4.11392741e+01*(u2*a1) + -7.39848186e+01*(u3*eta) + -3.48216369e+01*(u2*a1*a1) + 4.13945565e+02*(eta2*a1*a1) + -1.41609599e+03*(u*eta2*eta) + -8.86161650e+02*(u*eta*a1*a1) + -5.20127369e+03*(u*eta2*a1) + -4.28113243e+02*(u2*eta*a1) + -1.77442646e+01*(u3*a1) + 4.10553410e+02*(eta3*a1) + 2.95946717e+02*(u2*eta*eta) + -9.12632916e+02*(u2*eta2*eta) + 2.80681291e+02*(u3*eta*eta) + 3.37924311e+02*(u2*eta*a1*a1) + 1.26339786e+03*(u2*eta2*a1) + 5.67524228e+03*(u*eta2*a1*a1) + 1.22333354e+02*(u3*eta*a1) + 1.00216210e+04*(u*eta3*a1) + 1.53612107e+01*(u3*a1*a1) + -8.89972937e+02*(u2*eta3*a1) + -8.27258185e+01*(u3*eta*a1*a1) + -7.84231099e+02*(u2*eta2*a1*a1) + -1.22639277e+02*(u3*eta2*a1) + -4.38379022e+02*(u3*eta2*eta) + -1.07353624e+04*(u*eta3*a1*a1)

	# nu0, 37 terms
	nu0_l2m2 = 1.20543450e+02*(u) + 1.27484961e+03*(eta) + 9.73061928e+02*(a1) + -1.39537396e+02 + -2.98191097e+03*(eta*eta) + 3.71950292e+02*(u*u) + -8.58556058e+03*(eta*a1) + -3.84715641e+02*(a1*a1) + -2.64697029e+03*(u*a1) + 3.67130832e+02*(u2*u) + 1.90770015e+04*(eta2*a1) + 3.51129846e+04*(u*eta*a1) + 3.79489418e+03*(u*a1*a1) + -2.79214444e+03*(u2*a1) + -3.25454984e+03*(u2*eta) + -9.07293457e+03*(u3*eta) + 2.34027481e+03*(u2*a1*a1) + 2.55795990e+04*(eta2*a1*a1) + -1.27304698e+04*(u*eta2*eta) + -5.45669002e+04*(u*eta*a1*a1) + -1.79249757e+05*(u*eta2*a1) + 2.96460120e+04*(u2*eta*a1) + 1.11693629e+03*(u3*a1) + 2.63207469e+04*(u2*eta2*eta) + 5.02235441e+04*(u3*eta*eta) + -2.40717405e+04*(u2*eta*a1*a1) + -7.31512307e+04*(u2*eta2*a1) + 2.81123919e+05*(u*eta2*a1*a1) + 3.38873151e+05*(u*eta3*a1) + -7.73926053e+04*(eta3*a1*a1) + -1.08578690e+03*(u3*a1*a1) + -3.65979106e+03*(u3*eta*a1) + 5.08606096e+03*(u3*eta*a1*a1) + 5.81371672e+04*(u2*eta2*a1*a1) + -8.16469730e+03*(u3*eta2*a1) + -7.19780980e+04*(u3*eta2*eta) + -5.02931114e+05*(u*eta3*a1*a1)

	# mu1, 53 terms
	mu1_l3m3 = -3.97201083e+00*(a1) + -4.74585442e-02 + 6.15840277e+00*(u*u) + 6.26778106e+01*(eta*a1) + 1.47669704e+01*(a1*a1) + -2.47370680e+00*(u*a1) + -1.29961099e+01*(a1*a1*a1) + -2.82496042e+00*(u2*u) + -2.21682691e+02*(eta*a1*a1) + -2.28924058e+02*(eta2*a1) + 3.90144916e+01*(u*eta*a1) + 9.14387768e+00*(u*a1*a1) + 1.08839756e+01*(u*eta*eta) + -1.84899028e+01*(u2*a1) + -9.16326585e+01*(u2*eta) + 4.16671139e+01*(u3*eta) + -9.03780786e+00*(u3*u) + 7.93698872e+02*(eta2*a1*a1) + 1.92861740e+02*(eta*a1*a1*a1) + 2.98936067e+02*(u2*eta*eta) + -1.46405043e+02*(u*eta*a1*a1) + -7.78226294e+00*(u*a1*a1*a1) + -2.16592109e+02*(u*eta2*a1) + 2.70998766e+02*(u2*eta*a1) + 2.58430899e+01*(u3*a1) + 3.75918579e+01*(u4*a1) + 1.26471019e+02*(u*eta*a1*a1*a1) + -1.55823155e+02*(u3*eta*eta) + -7.49200945e+02*(u2*eta2*a1) + 1.32575882e+02*(u4*eta) + 6.91776335e+02*(u*eta2*a1*a1) + -6.83034770e+02*(eta2*a1*a1*a1) + -3.82806192e+02*(u3*eta*a1) + -6.15356015e+01*(u3*a1*a1) + 2.08507594e+01*(u2*a1*a1*a1) + 4.31808360e+01*(u3*a1*a1*a1) + 9.25412533e+02*(u3*eta*a1*a1) + -4.20559859e+02*(u2*eta2*a1*a1) + -4.23255825e+02*(u4*eta*eta) + 1.41001232e+03*(u3*eta2*a1) + -5.46478144e+02*(u4*eta*a1) + -5.68991980e+02*(u*eta2*a1*a1*a1) + -4.29912123e+01*(u4*a1*a1) + -3.08979423e+02*(u2*eta*a1*a1*a1) + -6.59305466e+02*(u3*eta*a1*a1*a1) + 1.34885404e+03*(u2*eta2*a1*a1*a1) + 6.26793281e+02*(u4*eta*a1*a1) + 1.57899898e+03*(u4*eta2*a1) + 8.92961928e+00*(u4*a1*a1*a1) + -3.41677564e+03*(u3*eta2*a1*a1) + -1.50277733e+03*(u4*eta2*a1*a1) + -1.26595079e+02*(u4*eta*a1*a1*a1) + 2.45185338e+03*(u3*eta2*a1*a1*a1)

	# mu2, 54 terms
	mu2_l3m3 = -2.28752621e+00*(u) + 1.19221294e+01*(eta) + 7.57940532e+00*(a1) + -7.33906513e-01 + 2.20133668e+01*(u*eta) + -4.30180611e+01*(eta*eta) + -9.67668912e-01*(u*u) + -1.21006823e+02*(eta*a1) + -1.70268577e+01*(a1*a1) + 2.15062714e+01*(u*a1) + 1.15538782e+01*(a1*a1*a1) + 2.60459086e+00*(u2*u) + 2.78558570e+02*(eta*a1*a1) + 4.24418916e+02*(eta2*a1) + -2.15718092e+02*(u*eta*a1) + -5.12425996e+01*(u*a1*a1) + -3.72364324e+01*(u*eta*eta) + -9.67166213e+00*(u2*a1) + -1.91255093e+01*(u3*eta) + 2.56895044e+01*(u2*a1*a1) + 2.61389395e+00*(u3*u) + -9.82381697e+02*(eta2*a1*a1) + -1.92236955e+02*(eta*a1*a1*a1) + 5.11880431e+02*(u*eta*a1*a1) + 3.48590420e+01*(u*a1*a1*a1) + 4.27957718e+02*(u*eta2*a1) + 2.30738821e+02*(u2*eta*a1) + -2.71823217e+01*(u3*a1) + -3.40414405e+02*(u*eta*a1*a1*a1) + -5.50747683e+02*(u2*eta*a1*a1) + -6.88831297e+02*(u2*eta2*a1) + -2.35861836e+01*(u4*eta) + -1.02678743e+03*(u*eta2*a1*a1) + 6.82651504e+02*(eta2*a1*a1*a1) + 2.43367441e+02*(u3*eta*a1) + 6.95032043e+01*(u3*a1*a1) + -1.37668487e+01*(u2*a1*a1*a1) + -4.99061541e+01*(u3*a1*a1*a1) + -6.61324255e+02*(u3*eta*a1*a1) + 1.53878249e+03*(u2*eta2*a1*a1) + 7.41693783e+01*(u4*eta*eta) + -3.22206427e+02*(u3*eta2*a1) + -9.32343574e+01*(u4*eta*a1) + 6.58595136e+02*(u*eta2*a1*a1*a1) + -5.73104298e+00*(u4*a1*a1) + 3.06731952e+02*(u2*eta*a1*a1*a1) + 4.83616845e+02*(u3*eta*a1*a1*a1) + -7.47535646e+02*(u2*eta2*a1*a1*a1) + 2.59546920e+02*(u4*eta*a1*a1) + 2.33962894e+02*(u4*eta2*a1) + 1.09785323e+03*(u3*eta2*a1*a1) + -5.23096151e+02*(u4*eta2*a1*a1) + -1.00500234e+02*(u4*eta*a1*a1*a1) + -8.53235749e+02*(u3*eta2*a1*a1*a1)

	# mu3, 58 terms
	mu3_l3m3 = 1.03794087e-01*(u) + -7.25068914e-01*(eta) + -3.17207943e-01*(a1) + 4.32830894e-02 + -1.40182949e+00*(u*eta) + 2.53176202e+00*(eta*eta) + -2.38823518e-01*(u*u) + 5.31041760e+00*(eta*a1) + 6.92167653e-01*(a1*a1) + -9.18777931e-01*(u*a1) + -4.85579123e-01*(a1*a1*a1) + -1.38026848e-01*(u2*u) + -1.13594813e+01*(eta*a1*a1) + -1.91860161e+01*(eta2*a1) + 1.26537428e+01*(u*eta*a1) + 2.42553938e+00*(u*a1*a1) + 5.54040120e+00*(u*eta*eta) + 1.33493452e+00*(u2*a1) + 3.56203076e+00*(u2*eta) + 1.85865525e+00*(u3*eta) + -1.86541718e+00*(u2*a1*a1) + 3.11596849e-01*(u3*u) + 4.14252019e+01*(eta2*a1*a1) + 7.73582324e+00*(eta*a1*a1*a1) + -1.04915834e+01*(u2*eta*eta) + -3.41942575e+01*(u*eta*a1*a1) + -1.93970817e+00*(u*a1*a1*a1) + -4.96097179e+01*(u*eta2*a1) + -2.06386905e+01*(u2*eta*a1) + 1.16938901e+00*(u3*a1) + -1.79737292e+00*(u4*a1) + 2.80327709e+01*(u*eta*a1*a1*a1) + -7.63978162e+00*(u3*eta*eta) + 2.88122419e+01*(u2*eta*a1*a1) + 5.99259044e+01*(u2*eta2*a1) + -4.67993978e+00*(u4*eta) + 1.33408097e+02*(u*eta2*a1*a1) + -2.81540220e+01*(eta2*a1*a1*a1) + -1.58124265e+01*(u3*eta*a1) + -3.02653408e+00*(u3*a1*a1) + 4.16192190e-01*(u2*a1*a1*a1) + 2.37455081e+00*(u3*a1*a1*a1) + 4.17810669e+01*(u3*eta*a1*a1) + -7.33493169e+01*(u2*eta2*a1*a1) + 1.34008098e+01*(u4*eta*eta) + 6.47569418e+01*(u3*eta2*a1) + 2.78528278e+01*(u4*eta*a1) + -1.09286935e+02*(u*eta2*a1*a1*a1) + 2.58004619e+00*(u4*a1*a1) + -5.85659435e+00*(u2*eta*a1*a1*a1) + -3.36641582e+01*(u3*eta*a1*a1*a1) + -4.03331462e+01*(u4*eta*a1*a1) + -7.76138066e+01*(u4*eta2*a1) + -6.12408157e-01*(u4*a1*a1*a1) + -1.69737750e+02*(u3*eta2*a1*a1) + 9.63660324e+01*(u4*eta2*a1*a1) + 9.36339176e+00*(u4*eta*a1*a1*a1) + 1.36338394e+02*(u3*eta2*a1*a1*a1)

	# mu4, 58 terms
	mu4_l3m3 = -5.17317824e+00*(u) + 6.90771831e+01*(eta) + 6.09021855e+01*(a1) + -6.64914977e+00 + 8.45897256e+01*(u*eta) + -1.67304672e+02*(eta*eta) + 1.22477360e+01*(u*u) + -6.76446739e+02*(eta*a1) + -1.53665348e+02*(a1*a1) + 3.88092780e+01*(u*a1) + 1.06715594e+02*(a1*a1*a1) + 1.65770568e+01*(u2*u) + 1.76124733e+03*(eta*a1*a1) + 1.56590417e+03*(eta2*a1) + -6.43185463e+02*(u*eta*a1) + -8.16802574e+01*(u*a1*a1) + -1.66480863e+02*(u*eta*eta) + -2.60159360e+02*(u2*a1) + -9.02514480e+01*(u2*eta) + -2.54527157e+02*(u3*eta) + 7.62559143e+02*(u2*a1*a1) + 1.50267473e+01*(u3*u) + -4.18400934e+03*(eta2*a1*a1) + -1.24534299e+03*(eta*a1*a1*a1) + 1.41548780e+03*(u*eta*a1*a1) + 5.99609819e+01*(u*a1*a1*a1) + 1.33121718e+03*(u*eta2*a1) + 2.80310696e+03*(u2*eta*a1) + -1.20848054e+02*(u3*a1) + 5.53447485e+01*(u4*a1) + -1.04402340e+03*(u*eta*a1*a1*a1) + 6.22890585e+02*(u3*eta*eta) + -8.63486419e+03*(u2*eta*a1*a1) + -5.93768982e+03*(u2*eta2*a1) + -3.27786849e+02*(u4*eta) + -3.06891230e+03*(u*eta2*a1*a1) + 2.98645348e+03*(eta2*a1*a1*a1) + 1.93548221e+03*(u3*eta*a1) + 2.59641090e+02*(u3*a1*a1) + -5.73662444e+02*(u2*a1*a1*a1) + -1.83662949e+02*(u3*a1*a1*a1) + -4.33586092e+03*(u3*eta*a1*a1) + 2.00402634e+04*(u2*eta2*a1*a1) + 1.29099903e+03*(u4*eta*eta) + -4.88402674e+03*(u3*eta2*a1) + 4.30407584e+02*(u4*eta*a1) + 2.45350577e+03*(u*eta2*a1*a1*a1) + -3.20930458e+02*(u4*a1*a1) + 6.61978212e+03*(u2*eta*a1*a1*a1) + 3.13574702e+03*(u3*eta*a1*a1*a1) + -1.57564112e+04*(u2*eta2*a1*a1*a1) + 1.51518904e+03*(u4*eta*a1*a1) + -4.45498268e+03*(u4*eta2*a1) + 2.84738269e+02*(u4*a1*a1*a1) + 1.13569736e+04*(u3*eta2*a1*a1) + 3.30169005e+03*(u4*eta2*a1*a1) + -1.87975981e+03*(u4*eta*a1*a1*a1) + -8.59667277e+03*(u3*eta2*a1*a1*a1)

	# nu4, 56 terms
	nu4_l3m3 = -2.40128981e+01*(u) + 5.12723958e+02*(eta) + 3.42050294e+02*(a1) + -4.20171044e+01 + 1.89948030e+02*(u*eta) + -1.50696434e+03*(eta*eta) + 3.31192346e+01*(u*u) + -4.28980386e+03*(eta*a1) + -8.70271555e+02*(a1*a1) + 1.99048000e+02*(u*a1) + 6.46280447e+02*(a1*a1*a1) + 1.06016214e+02*(u2*u) + 1.07852692e+04*(eta*a1*a1) + 1.24559106e+04*(eta2*a1) + -1.77245591e+03*(u*eta*a1) + -3.78892091e+02*(u*a1*a1) + -3.34909916e+02*(u*eta*eta) + -5.69349485e+02*(u2*a1) + -1.21979684e+03*(u3*eta) + 1.90786409e+03*(u2*a1*a1) + 3.10062562e+01*(u3*u) + -3.09175658e+04*(eta2*a1*a1) + -7.93646729e+03*(eta*a1*a1*a1) + -1.27505932e+03*(u2*eta*eta) + 3.47417286e+03*(u*eta*a1*a1) + 2.01300150e+02*(u*a1*a1*a1) + 3.83670261e+03*(u*eta2*a1) + 4.33607659e+03*(u2*eta*a1) + -7.26503919e+02*(u3*a1) + -1.85811193e+03*(u*eta*a1*a1*a1) + 3.46658421e+03*(u3*eta*eta) + -1.77866928e+04*(u2*eta*a1*a1) + -4.29777703e+03*(u2*eta2*a1) + -8.81904607e+02*(u4*eta) + -8.03807277e+03*(u*eta2*a1*a1) + 2.26200448e+04*(eta2*a1*a1*a1) + 8.34167298e+03*(u3*eta*a1) + 1.39303183e+03*(u3*a1*a1) + -1.61653052e+03*(u2*a1*a1*a1) + -7.94848150e+02*(u3*a1*a1*a1) + -1.59636372e+04*(u3*eta*a1*a1) + 3.39641515e+04*(u2*eta2*a1*a1) + 3.53219437e+03*(u4*eta*eta) + -2.37376419e+04*(u3*eta2*a1) + 3.21323384e+03*(u4*eta*a1) + 4.67709650e+03*(u*eta2*a1*a1*a1) + -5.40781097e+02*(u4*a1*a1) + 1.61219831e+04*(u2*eta*a1*a1*a1) + 9.13773557e+03*(u3*eta*a1*a1*a1) + -3.52329090e+04*(u2*eta2*a1*a1*a1) + -1.60301022e+04*(u4*eta2*a1) + 6.19735931e+02*(u4*a1*a1*a1) + 4.58117724e+04*(u3*eta2*a1*a1) + 1.47161050e+04*(u4*eta2*a1*a1) + -3.32726118e+03*(u4*eta*a1*a1*a1) + -2.68286745e+04*(u3*eta2*a1*a1*a1)

	# nu5, 57 terms
	nu5_l3m3 = -8.24878309e-02*(u) + -5.07644770e-01*(eta) + -2.66827485e-01*(a1) + 1.89838094e-02 + 1.46045005e+00*(u*eta) + 1.63293906e+00*(eta*eta) + 1.33321919e-01*(u*u) + 4.86048986e+00*(eta*a1) + 4.84924109e-01*(a1*a1) + 5.10633539e-01*(u*a1) + -3.57511537e-01*(a1*a1*a1) + 1.62965822e-01*(u2*u) + -1.03284910e+01*(eta*a1*a1) + -1.65210855e+01*(eta2*a1) + -9.46604547e+00*(u*eta*a1) + -9.46218752e-01*(u*a1*a1) + -4.77199303e+00*(u*eta*eta) + -6.38748746e-01*(u2*a1) + -8.40180128e-01*(u2*eta) + -2.64759847e+00*(u3*eta) + 1.35512613e+00*(u2*a1*a1) + -1.95132677e-01*(u3*u) + 3.63654748e+01*(eta2*a1*a1) + 6.99973693e+00*(eta*a1*a1*a1) + 2.60179287e+00*(u2*eta*eta) + 1.88066756e+01*(u*eta*a1*a1) + 5.41529909e-01*(u*a1*a1*a1) + 3.14435546e+01*(u*eta2*a1) + 9.12187971e-01*(u2*eta*a1) + -1.01369556e+00*(u3*a1) + 9.60281195e-01*(u4*a1) + -1.21027329e+01*(u*eta*a1*a1*a1) + 8.27827701e+00*(u3*eta*eta) + 1.61092537e+00*(u4*eta) + -6.43391156e+01*(u*eta2*a1*a1) + -2.46634264e+01*(eta2*a1*a1*a1) + 1.71321296e+01*(u3*eta*a1) + 1.89967469e+00*(u3*a1*a1) + -8.02700176e-01*(u2*a1*a1*a1) + -1.12630475e+00*(u3*a1*a1*a1) + -3.40307902e+01*(u3*eta*a1*a1) + -1.07682438e+01*(u2*eta2*a1*a1) + -5.76871627e+00*(u4*eta*eta) + -5.33777238e+01*(u3*eta2*a1) + -4.92695576e+00*(u4*eta*a1) + 4.27259375e+01*(u*eta2*a1*a1*a1) + -1.71078938e+00*(u4*a1*a1) + -1.44126372e+00*(u2*eta*a1*a1*a1) + 2.16955099e+01*(u3*eta*a1*a1*a1) + 1.37598246e+01*(u2*eta2*a1*a1*a1) + 5.62308390e+00*(u4*eta*a1*a1) + 1.74106458e+01*(u4*eta2*a1) + 8.74777467e-01*(u4*a1*a1*a1) + 1.07945784e+02*(u3*eta2*a1*a1) + -1.80877359e+01*(u4*eta2*a1*a1) + -4.62989152e-01*(u4*eta*a1*a1*a1) + -7.03264243e+01*(u3*eta2*a1*a1*a1)

	# nu6, 52 terms
	nu6_l3m3 = -5.79213752e-02*(u) + -1.04694533e+00*(eta) + -6.49775256e-01*(a1) + 8.12733140e-02 + 9.08626191e-01*(u*eta) + 3.08908199e+00*(eta*eta) + -6.48875872e-02*(u*u) + 8.46763806e+00*(eta*a1) + 1.54526166e+00*(a1*a1) + 4.97274067e-01*(u*a1) + -1.07256879e+00*(a1*a1*a1) + -1.96522626e+01*(eta*a1*a1) + -2.41815472e+01*(eta2*a1) + -7.32186633e+00*(u*eta*a1) + -1.20421778e+00*(u*a1*a1) + -3.25869784e+00*(u*eta*eta) + 1.08921545e+00*(u2*a1) + 5.22169761e-01*(u2*eta) + -2.77944499e+00*(u2*a1*a1) + -5.92919148e-02*(u3*u) + 5.52316222e+01*(eta2*a1*a1) + 1.33883003e+01*(eta*a1*a1*a1) + 1.68846788e+01*(u*eta*a1*a1) + 7.87647562e-01*(u*a1*a1*a1) + 2.51066222e+01*(u*eta2*a1) + -1.31483397e+01*(u2*eta*a1) + -2.10290460e-01*(u3*a1) + 4.86715744e-02*(u4*a1) + -1.09644279e+01*(u*eta*a1*a1*a1) + 3.52228520e+01*(u2*eta*a1*a1) + 3.14781281e+01*(u2*eta2*a1) + 9.83692433e-01*(u4*eta) + -5.58984280e+01*(u*eta2*a1*a1) + -3.73053886e+01*(eta2*a1*a1*a1) + 2.98152733e+00*(u3*eta*a1) + 8.40320165e-01*(u3*a1*a1) + 1.95716611e+00*(u2*a1*a1*a1) + -7.05386547e-01*(u3*a1*a1*a1) + -1.11372815e+01*(u3*eta*a1*a1) + -9.28972715e+01*(u2*eta2*a1*a1) + -3.39456920e+00*(u4*eta*eta) + -9.64370934e+00*(u3*eta2*a1) + -1.58869029e+00*(u4*eta*a1) + 3.54781969e+01*(u*eta2*a1*a1*a1) + -2.51338683e+01*(u2*eta*a1*a1*a1) + 9.17626547e+00*(u3*eta*a1*a1*a1) + 6.91056737e+01*(u2*eta2*a1*a1*a1) + 7.96315024e-01*(u4*eta*a1*a1) + 7.64022860e+00*(u4*eta2*a1) + 3.42607269e+01*(u3*eta2*a1*a1) + -5.76206910e+00*(u4*eta2*a1*a1) + -2.72628565e+01*(u3*eta2*a1*a1*a1)

	# zeta1, 49 terms
	zeta1_l3m3 = -3.31763878e-02*(eta) + -2.06546504e-02*(a1) + 2.35785912e-03 + 1.16610041e-01*(eta*eta) + 5.96858115e-03*(u*u) + 3.14764494e-01*(eta*a1) + 5.59003772e-02*(a1*a1) + -4.28701133e-02*(a1*a1*a1) + -7.90037091e-03*(u2*u) + -8.31623079e-01*(eta*a1*a1) + -1.10955973e+00*(eta2*a1) + 2.05794355e-02*(u*eta*a1) + 1.54574085e-03*(u*a1*a1) + -8.41316299e-02*(u2*eta) + 1.03343751e-01*(u3*eta) + -7.69087251e-02*(u2*a1*a1) + -1.36764324e-02*(u3*u) + 2.86678618e+00*(eta2*a1*a1) + 6.32625982e-01*(eta*a1*a1*a1) + 2.73605865e-01*(u2*eta*eta) + -9.12999785e-02*(u*eta*a1*a1) + -1.17597787e-01*(u*eta2*a1) + 5.69293840e-02*(u3*a1) + 5.87361288e-02*(u4*a1) + 7.50929121e-02*(u*eta*a1*a1*a1) + -3.06383738e-01*(u3*eta*eta) + 1.03132139e+00*(u2*eta*a1*a1) + 1.89904524e-01*(u4*eta) + 5.21787118e-01*(u*eta2*a1*a1) + -2.15822006e+00*(eta2*a1*a1*a1) + -7.61514469e-01*(u3*eta*a1) + -1.20771926e-01*(u3*a1*a1) + 9.46083103e-02*(u2*a1*a1*a1) + 8.26540002e-02*(u3*a1*a1*a1) + 1.65836583e+00*(u3*eta*a1*a1) + -3.10833306e+00*(u2*eta2*a1*a1) + -5.66789427e-01*(u4*eta*eta) + 2.33491889e+00*(u3*eta2*a1) + -7.83833838e-01*(u4*eta*a1) + -4.95620167e-01*(u*eta2*a1*a1*a1) + -5.45118632e-02*(u4*a1*a1) + -1.25967508e+00*(u2*eta*a1*a1*a1) + -1.16291296e+00*(u3*eta*a1*a1*a1) + 3.76245794e+00*(u2*eta2*a1*a1*a1) + 7.19797395e-01*(u4*eta*a1*a1) + 2.23186462e+00*(u4*eta2*a1) + -5.27263959e+00*(u3*eta2*a1*a1) + -2.03341856e+00*(u4*eta2*a1*a1) + 3.80673478e+00*(u3*eta2*a1*a1*a1)

	# zeta2, 54 terms
	zeta2_l3m3 = -2.83814989e+01*(u) + 1.44124042e+03*(eta) + 9.68068557e+02*(a1) + -1.10729482e+02 + 8.68915897e+01*(u*eta) + -4.34995948e+03*(eta*eta) + -8.39172020e+01*(u*u) + -1.33529576e+04*(eta*a1) + -2.58161270e+03*(a1*a1) + 2.05601781e+02*(u*a1) + 1.96327026e+03*(a1*a1*a1) + 3.52742572e+02*(u2*u) + 3.51042652e+04*(eta*a1*a1) + 4.20391938e+04*(eta2*a1) + -1.35332098e+03*(u*eta*a1) + -2.67808817e+02*(u*a1*a1) + -9.93489962e+02*(u2*a1) + 2.20575949e+03*(u2*eta) + -4.32328256e+03*(u3*eta) + 4.98799385e+03*(u2*a1*a1) + 3.74373961e+02*(u3*u) + -1.09882550e+05*(eta2*a1*a1) + -2.65523993e+04*(eta*a1*a1*a1) + -9.90824931e+03*(u2*eta*eta) + 2.16244904e+03*(u*eta*a1*a1) + 3.88803409e+03*(u*eta2*a1) + 6.06023655e+03*(u2*eta*a1) + -2.48029505e+03*(u3*a1) + -1.24497445e+03*(u4*a1) + 1.30813716e+04*(u3*eta*eta) + -5.09709545e+04*(u2*eta*a1*a1) + -6.34249408e+03*(u4*eta) + -1.04477890e+04*(u*eta2*a1*a1) + 8.30296680e+04*(eta2*a1*a1*a1) + 3.12861431e+04*(u3*eta*a1) + 4.71290738e+03*(u3*a1*a1) + -4.78198498e+03*(u2*a1*a1*a1) + -2.72128542e+03*(u3*a1*a1*a1) + -6.12779307e+04*(u3*eta*a1*a1) + 1.16931529e+05*(u2*eta2*a1*a1) + 2.17607416e+04*(u4*eta*eta) + -9.75684475e+04*(u3*eta2*a1) + 2.46398616e+04*(u4*eta*a1) + 5.81577721e+03*(u*eta2*a1*a1*a1) + 5.38140113e+04*(u2*eta*a1*a1*a1) + 3.67899472e+04*(u3*eta*a1*a1*a1) + -1.39670185e+05*(u2*eta2*a1*a1*a1) + -1.68546392e+04*(u4*eta*a1*a1) + -8.79309854e+04*(u4*eta2*a1) + 1.29029033e+03*(u4*a1*a1*a1) + 1.98401960e+05*(u3*eta2*a1*a1) + 8.02493792e+04*(u4*eta2*a1*a1) + -6.36130670e+03*(u4*eta*a1*a1*a1) + -1.24809928e+05*(u3*eta2*a1*a1*a1)

	# nu0, 58 terms
	nu0_l3m3 = -2.18165286e+03*(u) + -2.31665259e+04*(eta) + -1.38861147e+04*(a1) + 1.87231600e+03 + 3.66359141e+04*(u*eta) + 6.56144338e+04*(eta*eta) + 2.58775506e+03*(u*u) + 1.85827683e+05*(eta*a1) + 3.53853017e+04*(a1*a1) + 2.11568803e+04*(u*a1) + -2.66195486e+04*(a1*a1*a1) + 7.37656576e+03*(u2*u) + -4.60902900e+05*(eta*a1*a1) + -5.44921063e+05*(eta2*a1) + -3.22110807e+05*(u*eta*a1) + -5.88593695e+04*(u*a1*a1) + -1.16794257e+05*(u*eta*eta) + -1.02916880e+04*(u2*a1) + -5.39551682e+04*(u2*eta) + -1.10621290e+05*(u3*eta) + -1.07670319e+03*(u2*a1*a1) + -6.32311176e+03*(u3*u) + 1.34274732e+06*(eta2*a1*a1) + 3.41144260e+05*(eta*a1*a1*a1) + 2.67288975e+05*(u2*eta*eta) + 8.35303379e+05*(u*eta*a1*a1) + 4.50112679e+04*(u*a1*a1*a1) + 9.88354536e+05*(u*eta2*a1) + 2.69940368e+05*(u2*eta*a1) + -8.49671969e+04*(u3*a1) + 3.10003330e+04*(u4*a1) + -6.14738474e+05*(u*eta*a1*a1*a1) + 3.42978149e+05*(u3*eta*eta) + -2.86902443e+05*(u2*eta*a1*a1) + -1.45470905e+06*(u2*eta2*a1) + 9.02943783e+04*(u4*eta) + -2.47372416e+06*(u*eta2*a1*a1) + -9.91251471e+05*(eta2*a1*a1*a1) + 1.18552862e+06*(u3*eta*a1) + 2.44558921e+05*(u3*a1*a1) + 1.56471744e+04*(u2*a1*a1*a1) + -1.85960178e+05*(u3*a1*a1*a1) + -3.30534611e+06*(u3*eta*a1*a1) + 2.17392710e+06*(u2*eta2*a1*a1) + -3.42722828e+05*(u4*eta*eta) + -3.57255945e+06*(u3*eta2*a1) + -4.15952807e+05*(u4*eta*a1) + 1.77277837e+06*(u*eta2*a1*a1*a1) + -2.35832781e+04*(u4*a1*a1) + 2.47973907e+06*(u3*eta*a1*a1*a1) + -8.34743680e+05*(u2*eta2*a1*a1*a1) + 2.87026687e+05*(u4*eta*a1*a1) + 1.61879747e+06*(u4*eta2*a1) + -7.61359211e+03*(u4*a1*a1*a1) + 9.81339241e+06*(u3*eta2*a1*a1) + -1.53743639e+06*(u4*eta2*a1*a1) + 1.29599075e+05*(u4*eta*a1*a1*a1) + -7.29982038e+06*(u3*eta2*a1*a1*a1)

	#
	return (mu1_l2m2,mu2_l2m2,mu3_l2m2,nu4_l2m2,nu5_l2m2,nu6_l2m2,zeta1_l2m2,zeta2_l2m2,nu0_l2m2,mu1_l3m3,mu2_l3m3,mu3_l3m3,mu4_l3m3,nu4_l3m3,nu5_l3m3,nu6_l3m3,zeta1_l3m3,zeta2_l3m3,nu0_l3m3)
