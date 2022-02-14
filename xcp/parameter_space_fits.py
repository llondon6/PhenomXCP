

#
def generate_model_params(theta,eta,a1):

	'''
	Hola, soy un codigo escribido por "4b_document_fits.py". 
	~londonl@mit.edu/pilondon2@gmail.com 2020
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

	# mu1, 36 terms
	mu1 = -1.62202638e-01 + 5.47728895e-02*(u) + 1.49833799e+00*(eta) + 2.02086038e+00*(a1) + -3.32381754e+00*(eta*eta) + -2.89494869e+01*(eta*a1) + -3.60925797e+00*(a1*a1) + 5.67372703e+01*(eta*a1*a1) + 3.84559112e+00*(u2*eta) + -2.77866603e-01*(u*a1*a1) + -2.06041199e+00*(u2*a1) + -2.09673795e-01*(u2*u) + -9.82040922e+00*(u*eta*a1) + 1.35393041e+02*(eta2*a1) + -2.82683957e+02*(eta2*a1*a1) + 1.21764885e+00*(u3*eta) + -3.67710398e+01*(u2*eta*eta) + 1.30761014e+00*(u3*a1) + 7.28137724e+01*(u*eta2*a1) + 3.39563536e+00*(u2*a1*a1) + -8.39282254e+00*(u*eta2*eta) + -2.06153971e+02*(eta3*a1) + 1.80402643e+01*(u*eta*a1*a1) + 1.61004480e+01*(u2*eta*a1) + -1.22136863e+02*(u*eta2*a1*a1) + 8.93482503e+01*(u2*eta2*eta) + 4.54994772e+02*(eta3*a1*a1) + -3.75632057e+01*(u2*eta*a1*a1) + -9.08587635e+00*(u3*eta*a1) + -1.20928743e+00*(u3*a1*a1) + -1.08986793e+02*(u*eta3*a1) + 7.16054469e+00*(u3*eta2*a1) + 9.94016574e+01*(u2*eta2*a1*a1) + 1.94647457e+02*(u*eta3*a1*a1) + -1.41495571e+02*(u2*eta3*a1) + 6.97503312e+00*(u3*eta*a1*a1)

	# mu3, 36 terms
	mu3 = 3.98713523e-01 + -4.68010864e-01*(u) + -4.07981135e+00*(eta) + -1.75001833e+00*(a1) + -1.22155365e+00*(u*u) + 1.51643408e+01*(u*a1) + 1.85019419e+01*(eta*a1) + -2.47808895e+01*(u*a1*a1) + 9.37385395e+00*(u2*a1) + -2.33459774e+00*(u2*u) + 4.48692814e+01*(eta2*eta) + -2.24536280e+02*(u*eta*a1) + -8.19768170e+01*(eta2*a1*a1) + 5.70746704e+01*(u3*eta) + 1.07454157e+02*(u2*eta*eta) + -4.84515782e+00*(u3*a1) + 1.26416925e+03*(u*eta2*a1) + -1.26123659e+01*(u2*a1*a1) + 3.51615208e+01*(u*eta2*eta) + -2.02841322e+02*(eta3*a1) + 3.67594399e+02*(u*eta*a1*a1) + -6.08970925e+01*(u2*eta*a1) + -1.97991924e+03*(u*eta2*a1*a1) + -1.72144528e+02*(u2*eta2*a1) + -3.59422780e+02*(u2*eta2*eta) + 3.39698028e+02*(eta3*a1*a1) + 1.38675327e+02*(u2*eta*a1*a1) + 8.33303751e+00*(u3*a1*a1) + -3.33374673e+02*(u3*eta*eta) + -2.45363456e+03*(u*eta3*a1) + 5.55469975e+02*(u3*eta2*eta) + 9.16851566e+01*(u3*eta2*a1) + -3.57734963e+02*(u2*eta2*a1*a1) + 3.63525379e+03*(u*eta3*a1*a1) + 1.09987703e+03*(u2*eta3*a1) + -3.53473985e+01*(u3*eta*a1*a1)

	# mu4, 34 terms
	mu4 = -3.25079745e-04 + 3.27205799e-03*(delta) + -5.70163483e-04*(u*u) + -2.67404657e-02*(delta*a1) + 1.31134879e-02*(u*delta) + 2.26372484e-04*(a1*a1) + -3.38335239e-03*(u2*delta) + -4.34054603e-04*(u*a1*a1) + -7.22215933e-03*(delta2*delta) + -4.96719194e-02*(u*delta*a1) + -1.43174832e-02*(u*delta*delta) + 2.96306456e-02*(delta2*a1) + 5.55304129e-02*(delta*a1*a1) + -1.06622649e-02*(u2*delta*delta) + -2.21773351e-02*(u3*delta) + -2.92531370e-04*(u3*a1) + -1.19733281e-01*(delta2*a1*a1) + 2.20921166e-02*(u2*delta*a1) + 1.21663040e-02*(delta3*a1) + 2.26273634e-02*(u*delta2*a1) + 3.80427706e-02*(u*delta*a1*a1) + 1.11179462e-02*(u*delta2*delta) + -5.57439513e-02*(u2*delta*a1*a1) + -3.33680406e-02*(u*delta3*a1) + 2.25287374e-02*(u2*delta2*delta) + 8.22484389e-02*(delta3*a1*a1) + 4.30344812e-02*(u2*delta2*a1) + 1.88340289e-02*(u3*delta*delta) + 8.95861381e-02*(u3*delta*a1) + -6.87235434e-02*(u3*delta*a1*a1) + 7.50311377e-02*(u2*delta2*a1*a1) + -1.09921953e-01*(u2*delta3*a1) + -4.98462473e-02*(u3*delta2*a1) + 5.00435515e-02*(u*delta3*a1*a1)

	# nu4, 36 terms
	nu4 = -1.53680301e-04 + 5.22321996e-04*(u) + -2.81221022e-03*(delta) + -7.90634008e-05*(a1) + -3.25508078e-04*(u*u) + -2.71247654e-03*(u*a1) + 1.81341084e-02*(delta*a1) + -2.60949090e-03*(u*delta) + 4.93022264e-03*(delta*delta) + 4.78423147e-03*(u2*delta) + 3.21524321e-03*(u*a1*a1) + -3.13219416e-02*(delta2*a1) + -6.20738293e-04*(u2*u) + -2.44022939e-02*(delta*a1*a1) + -5.85515935e-03*(u2*delta*delta) + 7.15904170e-03*(u3*delta) + 3.46834653e-03*(u3*a1) + 3.97145823e-02*(delta2*a1*a1) + -3.13107141e-02*(u2*delta*a1) + 1.43877003e-03*(delta3*a1) + 4.47528552e-02*(u*delta2*a1) + 1.60630616e-03*(u*delta2*delta) + 4.48453250e-02*(u2*delta*a1*a1) + -5.04101375e-02*(u*delta2*a1*a1) + -3.72216612e-02*(u*delta3*a1) + -4.11422878e-03*(u2*delta2*delta) + 4.95322219e-02*(u2*delta2*a1) + -1.23509758e-02*(u3*delta*delta) + -4.02705031e-03*(u3*a1*a1) + -2.13706940e-02*(u3*delta*a1) + 6.49277532e-03*(u3*delta2*delta) + 2.03367650e-02*(u3*delta*a1*a1) + -7.89520445e-02*(u2*delta2*a1*a1) + 1.28385345e-02*(u2*delta3*a1) + 1.15171188e-02*(u3*delta2*a1) + 3.05261921e-02*(u*delta3*a1*a1)

	# nu5, 23 terms
	nu5 = 3.02281676e-04 + -1.09824909e-02*(u) + -1.72457987e-03*(u*u) + 5.90172392e-02*(u*a1) + 4.85970139e-02*(u*eta) + -4.66607802e-02*(a1*a1) + 3.79216120e-01*(eta*a1*a1) + -8.37699428e-02*(u*a1*a1) + 1.05240003e-02*(u2*a1) + -1.99917654e+00*(eta2*a1*a1) + 7.19893618e-02*(u3*eta) + -3.60036715e-02*(u3*a1) + -2.21513137e+00*(u*eta2*a1) + 3.06530620e+00*(u*eta2*a1*a1) + 5.48090801e-02*(u2*eta2*eta) + 4.71516170e+00*(eta3*a1*a1) + 8.10350408e-02*(u3*a1*a1) + -3.33798743e-01*(u3*eta*eta) + 4.74535619e+00*(u*eta3*a1) + 7.67315351e-01*(u3*eta2*a1) + -6.59573852e+00*(u*eta3*a1*a1) + -4.84473497e-01*(u2*eta3*a1) + -3.64883290e-01*(u3*eta*a1*a1)

	# nu6, 37 terms
	nu6 = 9.21007521e-03 + 2.06823624e-02*(u) + -1.41419639e-01*(eta) + -5.66950864e-02*(a1) + -3.53897470e-03*(u*u) + -1.97222508e-01*(u*a1) + 6.88221776e-01*(eta*eta) + -2.50551155e-01*(u*eta) + 8.92612970e-01*(eta*a1) + 1.16783944e-01*(a1*a1) + -1.80113085e+00*(eta*a1*a1) + 4.39443063e-02*(u2*eta) + 2.62056061e-01*(u*a1*a1) + 1.36320110e-02*(u2*u) + -1.08689390e+00*(eta2*eta) + 2.86800573e+00*(u*eta*a1) + -4.54537267e+00*(eta2*a1) + 1.27743489e+00*(u*eta*eta) + 9.27137473e+00*(eta2*a1*a1) + -3.51819191e-01*(u3*eta) + -1.23716424e-01*(u2*eta*eta) + 3.31288530e-02*(u3*a1) + -1.52346367e+01*(u*eta2*a1) + -1.29203218e-03*(u2*a1*a1) + -2.47103743e+00*(u*eta2*eta) + 7.54487568e+00*(eta3*a1) + -3.78266517e+00*(u*eta*a1*a1) + 1.94061986e+01*(u*eta2*a1*a1) + -1.57105931e+01*(eta3*a1*a1) + 1.20343591e-01*(u3*eta*a1) + -7.65646672e-02*(u3*a1*a1) + 1.92159870e+00*(u3*eta*eta) + 2.78539767e+01*(u*eta3*a1) + -2.86274788e+00*(u3*eta2*eta) + -1.07883930e+00*(u3*eta2*a1) + 3.14879845e-01*(u3*eta*a1*a1) + -3.40103705e+01*(u*eta3*a1*a1)

	# zeta1, 42 terms
	zeta1 = -1.73573456e-05 + 2.64805206e-04*(eta) + 7.33932118e-05*(a1) + 4.58257478e-05*(u*u) + 8.81447882e-05*(u*a1) + -1.32385197e-03*(eta*eta) + 1.01473363e-04*(u*eta) + -9.81760595e-04*(eta*a1) + -3.85356831e-05*(a1*a1) + 3.37264221e-04*(eta*a1*a1) + -7.34773057e-04*(u2*eta) + -1.04072214e-04*(u*a1*a1) + -1.79090376e-04*(u2*a1) + -3.13702399e-05*(u2*u) + 2.16738021e-03*(eta2*eta) + -2.31520423e-03*(u*eta*a1) + 4.22251647e-03*(eta2*a1) + -1.02945101e-03*(u*eta*eta) + 5.06489481e-04*(u3*eta) + 3.93965250e-03*(u2*eta*eta) + 4.71255683e-05*(u3*a1) + 1.58748174e-02*(u*eta2*a1) + 9.25805601e-05*(u2*a1*a1) + 2.39437711e-03*(u*eta2*eta) + -5.86086680e-03*(eta3*a1) + 2.85564226e-03*(u*eta*a1*a1) + 2.57808774e-03*(u2*eta*a1) + -1.90254425e-02*(u*eta2*a1*a1) + -1.28533140e-02*(u2*eta2*a1) + -6.94945581e-03*(u2*eta2*eta) + -2.94121402e-03*(eta3*a1*a1) + -8.37071097e-04*(u2*eta*a1*a1) + -4.69555935e-04*(u3*eta*a1) + -4.46377370e-05*(u3*a1*a1) + -2.49281356e-03*(u3*eta*eta) + -3.14778281e-02*(u*eta3*a1) + 4.16181276e-03*(u3*eta2*eta) + 7.23384477e-04*(u3*eta2*a1) + 1.88675597e-03*(u2*eta2*a1*a1) + 2.83855648e-04*(u3*eta*a1*a1) + 3.63891528e-02*(u*eta3*a1*a1) + 2.15937747e-02*(u2*eta3*a1)

	# zeta2, 42 terms
	zeta2 = 2.39832164e+00 + -2.84503972e+01*(eta) + -1.12041788e+01*(a1) + -6.25973611e+00*(u*u) + -7.37746757e+00*(u*a1) + 1.07401661e+02*(eta*eta) + -1.53629108e+01*(u*eta) + 1.10191473e+02*(eta*a1) + 7.53323965e+00*(a1*a1) + -3.24208055e+01*(eta*a1*a1) + 9.75919314e+01*(u2*eta) + 7.44943979e+00*(u*a1*a1) + 2.51619857e+01*(u2*a1) + 3.54967055e+00*(u2*u) + -1.26618116e+02*(eta2*eta) + 2.38223363e+02*(u*eta*a1) + -2.59544875e+02*(eta2*a1) + 1.41581365e+02*(u*eta*eta) + -2.70279098e+02*(eta2*a1*a1) + -4.63037636e+01*(u3*eta) + -4.96721284e+02*(u2*eta*eta) + -8.79437134e+00*(u3*a1) + -1.68252013e+03*(u*eta2*a1) + -1.38074536e+01*(u2*a1*a1) + -3.03781923e+02*(u*eta2*eta) + -2.62098096e+02*(u*eta*a1*a1) + -3.53945001e+02*(u2*eta*a1) + 1.84857874e+03*(u*eta2*a1*a1) + 1.65494555e+03*(u2*eta2*a1) + 8.43322928e+02*(u2*eta2*eta) + 1.11544909e+03*(eta3*a1*a1) + 1.38884132e+02*(u2*eta*a1*a1) + 5.88679098e+01*(u3*eta*a1) + 9.58716588e+00*(u3*a1*a1) + 1.98387317e+02*(u3*eta*eta) + 3.28197032e+03*(u*eta3*a1) + -3.15084499e+02*(u3*eta2*eta) + -4.19833276e+01*(u3*eta2*a1) + -3.22987256e+02*(u2*eta2*a1*a1) + -5.32551278e+01*(u3*eta*a1*a1) + -3.54766448e+03*(u*eta3*a1*a1) + -2.62454194e+03*(u2*eta3*a1)

	#
	return mu1,mu3,mu4,nu4,nu5,nu6,zeta1,zeta2
