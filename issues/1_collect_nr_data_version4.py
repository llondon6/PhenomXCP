#!/usr/bin/env python3

'''
#Load process and store NR data for coprecessing frame modeling.
londonl@mit.edu 2020

## Outline

This is the same as version 2, but will be used as a sanity check for the back-end code. Small changes listed below:
* unwrap_dphi input is now used in data formatting routine

'''

# Setup python environment
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import h5py
from os import path
import pickle
from xcp import *

# Preliminaries 
# --

# Define path for file IO
package_dir = parent( xcp.__path__[0] )
data_dir = package_dir + 'data/version4/'
alert(data_dir)
mkdir(data_dir,verbose=True)

# Find and load data
# --

# Define simulations to load
A = scsearch(keyword='pnr-catalog',verbose=True)
# Sort according to mass ration and spin
A.sort( key = lambda a: round( linalg.norm( a.m1 / a.m2 ), 1 ) + round( linalg.norm( a.S1 / a.m1**2 ), 2 ) + round( linalg.norm( a.S2 / a.m2**2 ), 2 ) )

# Reverse the sorting so that higher mass ratios are first
A = A[::-1]

#
catalog_path = package_dir + 'data/calibration_catalog.pickle'
alert('Saving scentry catalog list to %s'%magenta(catalog_path))
pickle.dump( A, open( catalog_path, "wb" ) )

# Let the people know.
alert('We have found %i simulations.'%len(A))

# Define loading parameters 
lmax = 3 # 4
pad = 1000
clean = True
dt = 0.5
kind = 'psi4'
gwf_smooth_width = 30 # samples

# Load and process simulations 
# --

# For all sims 
for a in [a for a in A if ('q2a02t150' in a.simname) ]:
# for a in [ b for b in A if 'q8a06t150dP25' in b.simname ]:#a in A:
# for a in A:
    
    # #
    # txt_file_path = data_dir+'%s_l%im%i.txt'%(a.simname,2,2)
    # if path.exists(txt_file_path):
    #     warning('It seems that %s already exists, so we\'re moving on ...'%magenta(txt_file_path),header=True)
    #     continue
    
    #
    alert('Processing: %s'%magenta(a.simname),header=True)
    
    # Load
    y_raw = gwylm(a,lmax=lmax,dt=dt,pad=pad,clean=clean,verbose=False)
    
    # Manage frames using dict defined below
    frame = {}
    frame['raw'] = y_raw

    # Put in initial J frame
    frame['init-j'] = y_raw.__calc_initial_j_frame__()
    
    # # Symmetrize the psi4 time domain coprecessing frame waveform, and return to the init-j frame
    # frame['star-init-j'] = gwylmo_cpclean( frame['init-j'], cp_domain='td' )
    
    # NOTE that although the angles model uses the j(t) frame, 
    # we do NOT use this here as the coprecessing frame is uniquely 
    # defined and the j(t) frame only adds problematic noise
    
    alert('Computing coprecessing frames',header=True)
    # Solve optimal emission problem for l subsets
    subframe       = {}
    cp_subframe_fd = {}
    s_cp_subframe_fd = {}
    # cp_subframe_td = {}
    subangles      = {}
    # subangles_td   = {}
    for ll in range(2,lmax+1):
            
        #
        mm = ll
        
        #
        txt_file_path = data_dir+'%s_l%im%i.txt'%(a.simname,ll,mm)
        the_case_has_been_processed = path.isfile(txt_file_path)
        
        #
        if the_case_has_been_processed:
            
            #
            alert('Skipping (%i,%i) for %s because the case has already been processed'%(ll,mm,cyan(a.simname)))
            
        else:
            
            #
            alert('Calculating coprecessing frame for l=%i subset ...'%ll)
            
            # Select only multipoles with l = ll
            raw_subframe = frame['init-j'].selectlm( [(ll,m_) for m_ in range(-ll,ll+1)] )
            
            # ** Go to coprecessing frame, symmetrize, revert back to original frame **
            subframe[ll] = gwylmo_cpclean( raw_subframe )
            
            '''NOTE that scrubbing is turned off below because it causes some angles to be artificially flipped'''
            # # Scrubbing away long duration noise
            # subframe[ll].scrub(apply=True)
            
            # Smoothing remaining noise at tolerance
            # NOTE that this object will be used for output
            subframe[ll] = subframe[ll].__smooth__(gwf_smooth_width)
            
            # Solve the optimal emission problem and rotate multipoles
            cp_subframe_fd[ll] =  subframe[ll].__calc_coprecessing_frame__( kind=kind, 
                                                                        transform_domain='fd' )
            
            # # Smoothing remaining noise at tolerance
            # # NOTE that this object will be used for output
            # s_cp_subframe_fd[ll] = cp_subframe_fd[ll].__smooth__(gwf_smooth_width)
            
            # Store angles for this ll
            foo = cp_subframe_fd[ll].previous_radiation_axis_info
            subangles[ll] = [ foo.radiation_axis[k] for k in ('fd_alpha','fd_beta','fd_gamma') ]
            # bar = cp_subframe_td[ll].previous_radiation_axis_info
            # subangles_td[ll] = [ foo.radiation_axis[k] for k in ('td_alpha','td_beta','td_gamma') ]
            
            # **
            alert('Saving diagnostic plot for l=m=%i ...'%ll)
            fig, ax, output_data, format_tags = collect_nr_data_plotting_helper( ll, cp_subframe_fd[ll], subangles[ll], unwrap_dphi=False )
            #
            png_file_path = data_dir+'%s_l%im%i.png'%(cp_subframe_fd[ll].simname,ll,mm)
            savefig( png_file_path, pad_inches=0 )
            alert('Diagnostic plot saved to "%s"'%yellow(png_file_path))
            close('all')
            
            #
            heaber = header = ',        '.join(format_tags)
            alert('Saving related data to "%s"'%yellow(txt_file_path))
            savetxt( txt_file_path, output_data, header=header )


alert('All done.')