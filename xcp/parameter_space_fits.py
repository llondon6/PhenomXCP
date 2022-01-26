

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

	# mu1, 39 terms
	mu1 = -1.53435960e-01 + 1.64884144e+00*(a1) + 9.53337422e-01*(eta) + -1.46903558e-01*(u) + -2.11016525e+01*(eta*a1) + -1.05098767e+00*(u*u) + 1.45397894e+00*(u*a1) + -3.45568913e+00*(a1*a1) + 1.50931043e+00*(u*eta) + 2.30385835e+01*(u2*eta) + -3.40623097e+01*(u*eta*a1) + -2.28328785e+00*(u*a1*a1) + -3.76401056e-01*(u2*u) + 9.70634760e+01*(eta2*a1) + 5.29422956e+01*(eta*a1*a1) + 2.90439519e+00*(u3*a1) + -1.34160650e+02*(u2*eta*eta) + 1.95165232e+02*(u*eta2*a1) + -2.66787474e+02*(eta2*a1*a1) + -1.66312045e+02*(eta3*a1) + 5.56315501e+01*(u*eta*a1*a1) + -2.12125055e+01*(u2*eta*a1) + 2.62047852e+00*(u2*a1*a1) + -2.21896429e+01*(u*eta2*eta) + -2.42273569e+00*(u3*a1*a1) + 2.13939425e+01*(u3*eta*eta) + -3.28633843e+02*(u*eta2*a1*a1) + 4.49840254e+02*(eta3*a1*a1) + -2.03538288e+01*(u3*eta*a1) + 1.61727922e+02*(u2*eta2*a1) + -2.83785475e+02*(u*eta3*a1) + -2.31576629e+01*(u2*eta*a1*a1) + 2.46536821e+02*(u2*eta2*eta) + -5.49964867e+01*(u3*eta2*eta) + 5.26907628e+02*(u*eta3*a1*a1) + -3.44063401e+02*(u2*eta3*a1) + 5.77539964e+01*(u2*eta2*a1*a1) + 1.33757040e+01*(u3*eta*a1*a1) + 2.17774586e+01*(u3*eta2*a1)

	# mu3, 40 terms
	mu3 = 2.74868402e-02 + -2.93437590e+00*(a1) + -9.74895539e-01*(u) + -2.66593769e+01*(eta*eta) + 4.43671476e+01*(eta*a1) + -2.10967680e+00*(u*u) + 2.12584077e+01*(u*a1) + 2.41286839e+00*(a1*a1) + 1.64925320e+01*(u*eta) + -3.60182114e+02*(u*eta*a1) + -3.39517307e+01*(u*a1*a1) + 1.25779511e+02*(eta2*eta) + 1.57746635e+01*(u2*a1) + -5.67058542e+00*(u2*u) + -1.49061184e+02*(eta2*a1) + -1.03094881e+02*(u*eta*eta) + -4.89434594e+01*(eta*a1*a1) + 1.95735623e+02*(u2*eta*eta) + 2.04922197e+03*(u*eta2*a1) + 2.13231524e+02*(eta2*a1*a1) + 1.00659336e+02*(u3*eta) + 5.40890861e+02*(u*eta*a1*a1) + -1.04047271e+02*(u2*eta*a1) + -2.18684929e+01*(u2*a1*a1) + 2.24085759e+02*(u*eta2*eta) + 8.72924364e+00*(u3*a1*a1) + -5.02792105e+02*(u3*eta*eta) + -2.91265767e+03*(u*eta2*a1*a1) + -1.69333511e+02*(eta3*a1*a1) + -5.33588390e+01*(u3*eta*a1) + -2.96407962e+02*(u2*eta2*a1) + -3.82515190e+03*(u*eta3*a1) + 2.42962860e+02*(u2*eta*a1*a1) + -6.67390774e+02*(u2*eta2*eta) + 7.14233188e+02*(u3*eta2*eta) + 5.17826212e+03*(u*eta3*a1*a1) + 1.94265031e+03*(u2*eta3*a1) + -6.41501169e+02*(u2*eta2*a1*a1) + -4.09489419e+01*(u3*eta*a1*a1) + 2.53217305e+02*(u3*eta2*a1)

	# mu4, 35 terms
	mu4 = -2.97929651e-03 + 7.73227806e-03*(a1) + 6.70211745e-04*(u) + 1.67688422e-02*(delta*delta) + -3.59477801e-03*(u*u) + -3.39513834e-03*(u*a1) + 1.13058547e-02*(u*delta) + -5.38127516e-03*(a1*a1) + -1.08025711e-02*(delta*a1) + -2.51772906e-02*(delta2*delta) + 2.40476407e-03*(u*a1*a1) + -3.21178169e-02*(u*delta*a1) + -1.28311303e-02*(u2*delta) + -2.23620450e-02*(u*delta*delta) + 5.94823662e-03*(u2*a1) + 4.78343675e-02*(delta*a1*a1) + -5.15514204e-02*(delta2*a1) + -6.62199728e-02*(delta2*a1*a1) + 2.98858901e-02*(u*delta2*delta) + -2.74793014e-02*(u3*delta) + 1.18251183e-01*(delta3*a1) + 6.41510885e-02*(u2*delta*a1) + -2.95182518e-03*(u2*a1*a1) + 1.53328629e-02*(u*delta*a1*a1) + -4.79500036e-02*(u*delta3*a1) + 2.50893773e-02*(u3*delta*delta) + 5.26730187e-02*(u*delta2*a1*a1) + 2.29331281e-02*(u2*delta2*delta) + -1.02621346e-01*(u2*delta*a1*a1) + 1.02706104e-01*(u3*delta*a1) + -1.29248359e-01*(u2*delta3*a1) + -6.76597958e-02*(u3*delta2*a1) + 1.40458256e-01*(u2*delta2*a1*a1) + -7.31499850e-02*(u3*delta*a1*a1) + 3.72761277e-02*(u*delta3*a1*a1)

	# nu4, 41 terms
	nu4 = -1.65988643e-03 + 4.24557131e-03*(a1) + 3.75036246e-04*(u) + -5.33305472e-03*(delta) + 1.71360962e-02*(delta*delta) + -1.60047429e-03*(u*u) + -2.07487426e-03*(u*a1) + 5.25879771e-03*(u*delta) + -3.19230576e-03*(a1*a1) + 2.58282006e-02*(delta*a1) + -8.29621736e-03*(delta2*delta) + -6.97311018e-02*(delta2*a1) + 2.79135668e-03*(u*a1*a1) + -3.97229924e-02*(u*delta*a1) + 6.90363372e-03*(u2*delta) + -2.67052006e-02*(u*delta*delta) + 1.58141436e-03*(u2*a1) + -3.15251865e-02*(delta*a1*a1) + 1.64259884e-03*(u3*a1) + 1.81449891e-01*(u*delta2*a1) + 7.09238072e-02*(delta2*a1*a1) + 2.46217444e-02*(u*delta2*delta) + 5.60041633e-03*(u3*delta) + 2.65194856e-02*(delta3*a1) + -3.41494039e-02*(u2*delta*a1) + 3.73639367e-02*(u*delta*a1*a1) + -1.29541876e-02*(u2*delta*delta) + -2.88054908e-03*(u3*a1*a1) + -1.54472246e-01*(u*delta3*a1) + -1.23211499e-02*(u3*delta*delta) + -1.83013493e-01*(u*delta2*a1*a1) + -1.82543489e-02*(delta3*a1*a1) + 4.81973054e-02*(u2*delta*a1*a1) + 6.63218583e-02*(u2*delta2*a1) + -1.71253367e-02*(u3*delta*a1) + 1.09030444e-02*(u2*delta3*a1) + 6.03569236e-03*(u3*delta2*delta) + 1.24775840e-02*(u3*delta2*a1) + -9.64801213e-02*(u2*delta2*a1*a1) + 1.86330633e-02*(u3*delta*a1*a1) + 1.42370752e-01*(u*delta3*a1*a1)

	# nu5, 39 terms
	nu5 = 1.25596929e-02 + -7.05753317e-02*(a1) + -3.69010705e-01*(eta) + -5.80986794e-02*(u) + 3.15798836e+00*(eta*eta) + 1.71667746e+00*(eta*a1) + 4.27184678e-01*(u*a1) + 1.13811027e+00*(u*eta) + 1.47442550e-01*(u2*eta) + -7.92212871e+00*(u*eta*a1) + -4.14668096e-01*(u*a1*a1) + -8.06478160e+00*(eta2*eta) + 2.22906890e-02*(u2*u) + -1.37002132e+01*(eta2*a1) + -7.37840406e+00*(u*eta*eta) + -6.88585980e-01*(eta*a1*a1) + -2.29332954e-01*(u3*a1) + -2.16729587e+00*(u2*eta*eta) + 4.61009275e+01*(u*eta2*a1) + 8.47252776e+00*(eta2*a1*a1) + 3.35893092e+01*(eta3*a1) + 8.80719863e+00*(u*eta*a1*a1) + -8.00973937e-01*(u2*eta*a1) + 2.84723613e-01*(u2*a1*a1) + 1.63833131e+01*(u*eta2*eta) + -4.06126211e-02*(u3*a1*a1) + -1.64336505e+00*(u3*eta*eta) + -5.43368367e+01*(u*eta2*a1*a1) + -2.43214187e+01*(eta3*a1*a1) + 2.49817041e+00*(u3*eta*a1) + 1.09390338e+01*(u2*eta2*a1) + -8.91357197e+01*(u*eta3*a1) + -3.42978838e+00*(u2*eta*a1*a1) + 6.82192092e+00*(u2*eta2*eta) + 3.85542725e+00*(u3*eta2*eta) + 1.06079823e+02*(u*eta3*a1*a1) + -3.27704788e+01*(u2*eta3*a1) + 9.55365840e+00*(u2*eta2*a1*a1) + -5.42603018e+00*(u3*eta2*a1)

	# nu6, 41 terms
	nu6 = -1.66788123e-03 + 1.40629568e-02*(a1) + 2.48809582e-02*(u) + 2.92427756e-01*(eta*eta) + 6.22172886e-03*(u*u) + -2.86648408e-01*(u*a1) + 6.34284487e-02*(a1*a1) + -2.67072140e-01*(u*eta) + 4.12357727e+00*(u*eta*a1) + 3.76688401e-01*(u*a1*a1) + -1.14865980e+00*(eta2*eta) + -1.06045443e-01*(u2*a1) + 3.77041179e-02*(u2*u) + -1.34625990e+00*(eta2*a1) + 1.41646717e+00*(u*eta*eta) + -1.14910222e+00*(eta*a1*a1) + 3.06435734e-02*(u3*a1) + -9.16152565e-01*(u2*eta*eta) + -2.18564544e+01*(u*eta2*a1) + 7.15849800e+00*(eta2*a1*a1) + -7.44647182e-01*(u3*eta) + 4.82125499e+00*(eta3*a1) + -5.39580663e+00*(u*eta*a1*a1) + 1.15629010e+00*(u2*eta*a1) + 1.20632204e-01*(u2*a1*a1) + -3.19467079e+00*(u*eta2*eta) + -1.04795458e-01*(u3*a1*a1) + 3.65658702e+00*(u3*eta*eta) + 2.75952933e+01*(u*eta2*a1*a1) + -1.45196742e+01*(eta3*a1*a1) + 3.42599834e-01*(u3*eta*a1) + -1.70183272e+00*(u2*eta2*a1) + 4.04166442e+01*(u*eta3*a1) + -1.48837813e+00*(u2*eta*a1*a1) + 3.26555398e+00*(u2*eta2*eta) + -4.86795213e+00*(u3*eta2*eta) + -4.85543635e+01*(u*eta3*a1*a1) + -5.18006571e+00*(u2*eta3*a1) + 4.06923477e+00*(u2*eta2*a1*a1) + 4.48845343e-01*(u3*eta*a1*a1) + -2.06703148e+00*(u3*eta2*a1)

	# zeta1, 39 terms
	zeta1 = -4.38426256e-06 + 6.74227037e-05*(a1) + 1.88743945e-05*(eta) + 1.68182782e-05*(u) + -6.75107731e-04*(eta*a1) + 2.77624114e-05*(u*u) + 1.51068748e-04*(u*a1) + -4.97782144e-05*(a1*a1) + -1.62714917e-04*(u*eta) + -2.34380212e-04*(u2*eta) + -3.64527401e-03*(u*eta*a1) + -1.82043737e-04*(u*a1*a1) + -1.90222526e-04*(u2*a1) + -8.85169284e-05*(u2*u) + 2.68211602e-03*(eta2*a1) + 3.21820332e-04*(eta*a1*a1) + 1.00276352e-04*(u3*a1) + 4.66100088e-04*(u2*eta*eta) + 2.54983158e-02*(u*eta2*a1) + 1.40671137e-03*(u3*eta) + -4.29741279e-03*(eta3*a1) + 4.45939919e-03*(u*eta*a1*a1) + 2.09661305e-03*(u2*eta*a1) + 1.37058897e-04*(u2*a1*a1) + 1.20670227e-03*(u*eta2*eta) + -7.08454938e-05*(u3*a1*a1) + -6.82956729e-03*(u3*eta*eta) + -2.97667734e-02*(u*eta2*a1*a1) + -1.90370141e-03*(eta3*a1*a1) + -1.07427239e-03*(u3*eta*a1) + -8.45591249e-03*(u2*eta2*a1) + -5.15995028e-02*(u*eta3*a1) + -1.05221928e-03*(u2*eta*a1*a1) + 1.13715053e-02*(u3*eta2*eta) + 5.75607032e-02*(u*eta3*a1*a1) + 1.27699256e-02*(u2*eta3*a1) + 1.95596446e-03*(u2*eta2*a1*a1) + 5.53004132e-04*(u3*eta*a1*a1) + 1.60059078e-03*(u3*eta2*a1)

	# zeta2, 41 terms
	zeta2 = 1.55516896e+00 + -9.91154235e+00*(a1) + -3.09475897e+00*(u) + -7.50730740e+01*(eta*eta) + -1.35318132e+01*(u*u) + -1.68055477e+01*(u*a1) + 5.96009842e+00*(a1*a1) + 2.47880962e+01*(u*eta) + 1.99756820e+02*(u2*eta) + 4.85759899e+02*(u*eta*a1) + 1.94803997e+01*(u*a1*a1) + 2.18660320e+02*(eta2*eta) + 4.97703662e+01*(u2*a1) + 1.56163919e+01*(u2*u) + 5.95522900e+02*(eta2*a1) + 8.83470132e+01*(eta*a1*a1) + -2.43440466e+01*(u3*a1) + -9.37956721e+02*(u2*eta*eta) + -3.41152563e+03*(u*eta2*a1) + -1.23810043e+03*(eta2*a1*a1) + -2.21096664e+02*(u3*eta) + -1.78197068e+03*(eta3*a1) + -5.47132285e+02*(u*eta*a1*a1) + -6.52316875e+02*(u2*eta*a1) + -3.02929530e+01*(u2*a1*a1) + -1.58313136e+02*(u*eta2*eta) + 1.80176708e+01*(u3*a1*a1) + 9.89514502e+02*(u3*eta*eta) + 3.68592782e+03*(u*eta2*a1*a1) + 3.18451760e+03*(eta3*a1*a1) + 1.99040792e+02*(u3*eta*a1) + 2.76070583e+03*(u2*eta2*a1) + 6.71051548e+03*(u*eta3*a1) + 2.94252098e+02*(u2*eta*a1*a1) + 1.50030416e+03*(u2*eta2*eta) + -1.51616085e+03*(u3*eta2*eta) + -6.96405576e+03*(u*eta3*a1*a1) + -4.04399795e+03*(u2*eta3*a1) + -6.43216031e+02*(u2*eta2*a1*a1) + -1.05867696e+02*(u3*eta*a1*a1) + -2.79546751e+02*(u3*eta2*a1)

	#
	return mu1,mu3,mu4,nu4,nu5,nu6,zeta1,zeta2
