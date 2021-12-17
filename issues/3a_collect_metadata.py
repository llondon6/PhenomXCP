#!/usr/bin/env python3

# Setup python environment
from numpy import *
from positive import *
from nrutils import scsearch
from numpy.linalg import norm
from xcp import determine_data_fitting_region,calibration_catalog
import pickle
from glob import glob
import xcp

#
l = 2

#
package_dir = parent( xcp.__path__[0] )
datadir = package_dir + 'data/version2/'
files = glob( datadir+'q*l%i*.txt'%l )
files.sort()

#
alert('We have %s files to consider.'%red(str(len(files))))

#
metadata = []
simnames = []

#
alert('Using %s and txt files at %s to gather metadata...'%(magenta('calibration_catalog'),magenta(datadir)))
for f in files:

    #
    file_name = f.split('/')[-1].split('.')[0].split('_l')[0]

    #
    A = scsearch( keyword=file_name, verbose=not True, catalog=calibration_catalog )

    #
    a = A[0]
    
    #
    l = a.L/norm(a.L)
    test_print = ('q1a' in file_name) and (norm(a.X1)<norm(a.X2))
    if test_print:
        warning(magenta('Flipping 1-2 labels in euqal mass case'),fname=file_name)
        a.X1,a.X2 = [ array(k) for k in (a.X2,a.X1) ]
        a.m1,a.m2 = [ float(k) for k in (a.m2,a.m1) ]
        print( '\t * |chi1| = ',norm(a.X1))
        print( '\t * chi1_l = ',dot(a.X1,l))
        print( '\t * chi2_l = ',dot(a.X2,l))

    #
    m1,m2 = [ k/(a.m1+a.m2) for k in (a.m1,a.m2) ] 
    eta = m1*m2/(m1+m2)
    
    #
    X1,X2,L,S = a.X1,a.X2,a.L,a.S
    
    #
    a1,a2 = norm(a.X1),norm(a.X2)
    
    #
    l = L/norm(L)
    s = S/norm(S)
    
    # Determin angles to rotate spins into frame where L || z
    # ---
    lx,ly,lz = l
    alpha = 0
    beta  = -arccos( lz )
    gamma = -arctan2( ly, lx )
    chi1_vec = rotate3( a.X1, alpha, beta, gamma )
    chi2_vec = rotate3( a.X2, alpha, beta, gamma )
    if test_print:
        print( '\t * chi1_vec = ',chi1_vec)
        print( '\t * chi2_vec = ',chi2_vec)
    
    
    # NOTE that theta is in radians
    theta = arccos( dot( l, s ) ) 
    
    #
    chi1 = dot(X1,l)
    chi2 = dot(X2,l)
    
    #
    chi_p   = calc_chi_p(   m1,X1, m2,X2, L )
    chi_eff = calc_chi_eff( m1,X1, m2,X2, L )
    
    #
    delta = (m1-m2)/(m1+m2)
    
    #
    simnames.append(file_name)
    metadata.append( [ theta,
                       m1,
                       m2,
                       eta,
                       delta,
                       chi_eff,
                       chi_p,
                       chi1,
                       chi2,
                       a1,
                       a2,
                       chi1_vec[0],
                       chi1_vec[1],
                       chi1_vec[2],
                       chi2_vec[0],
                       chi2_vec[1],
                       chi2_vec[2] ] )

#
print( 'Done.')

#
metadata_array = array(metadata)

#
keys = [ 'theta', 'm1', 'm2', 'eta', 'delta', 'chi_eff', 'chi_p', 'chi1', 'chi2', 'a1', 'a2', 'chi1_vec_x', 'chi1_vec_y', 'chi1_vec_z', 'chi2_vec_x', 'chi2_vec_y', 'chi2_vec_z' ]

#
metadata_dict = {}
metadata_dict = { keys[k]: metadata_array[:,k] for k in range(len(keys)) }
metadata_dict['simname'] = simnames

#
metadata_dict['array_data_columns'] = keys
metadata_dict['array_data'] = metadata_array

#
metadata_path = '/Users/book/KOALA/PhenomXCP/data/metadata_dict.pickle'
alert('Saving metadata dictionary to %s'%magenta(metadata_path))
pickle.dump( metadata_dict, open( metadata_path, "wb" ) )

#
alert('Saving complete.')
    
        

