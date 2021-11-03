#!/usr/bin/env python3

'''
#Load process and store NR data for coprecessing frame modeling.
londonl@mit.edu 2020

## Outline

1. Load simulation 
2. Compute TD Psi4 Co-Precessing Frame 
3. Symmetrize the co-precessing frame data 
4. Output the symmetrized FD amplitude and phase along with diagnostic plots
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
from mpa import gwylmo_cpclean

# Preliminaries 
# --

# Define path for file IO
package_dir = parent( xcp.__path__[0] )
data_dir = package_dir + 'data/version0/'
alert(data_dir)
mkdir(data_dir,verbose=True)

# Find and load data
# --

# Define simulations to load
A = scsearch(keyword='pnr-catalog',verbose=True)

#
catalog_path = package_dir + 'data/calibration_catalog.pickle'
alert('Saving scentry catalog list to %s'%magenta(catalog_path))
pickle.dump( A, open( catalog_path, "wb" ) )

# Let the people know.
alert('We have found %i simulations.'%len(A))

# Define loading parameters 
lmax = 4 # NOTE --- if we want moments up to ell, then lmax=ell+1 is needed for minimal accuracy
pad = 1800
clean = True
dt = 0.5

# Define indices of interest (ie ones that will be output after processing)
lm_of_interest = [ (2,1), (2,2), (3,2), (3,3), (4,4) ]

# Select multipole moments to use when calculating coprecessing frame
cpframe_select_lm = [ (2,1), (2,2), (2,-1), (2,-2), (2,0) ]
        
#
save_keys = ['cp-y-td-scrubbed-sym-z','cp-y-td-scrubbed-sym']

# Load and process simulations 
# --

# For all sims 
for a in A:
    
    #
    alert('Processing: %s'%magenta(a.simname),header=True)
    
    # Load
    y_raw = gwylm(a,lmax=lmax,dt=dt,pad=pad,clean=clean,verbose=False)
    
    # Manage frames using dict defined below
    frame = {}
    frame['raw'] = y_raw

    # Put in initial J frame
    frame['init-j'] = y_raw.__calc_initial_j_frame__()
    
    # NOTE that although the angles model uses the j(t) frame, 
    # we do NOT use this here as the coprecessing frame is uniquely 
    # defined and the j(t) frame only adds problematic noise
    
    # Compute TD adn FD coprecessing psi4 frames using new nrutils functoins
    frame['cp-y-fd'] = frame['init-j'].__calc_coprecessing_frame__( transform_domain='fd', kind='psi4',select_lm_list=cpframe_select_lm )
    frame['cp-y-td'] = frame['init-j'].__calc_coprecessing_frame__( transform_domain='td', kind='psi4',select_lm_list=cpframe_select_lm )
    
    #
    frame['cp-y-td-scrubbed'] = frame['cp-y-td'].scrub(lm=[(2,2),(2,1),(2,-1),(2,-2),(3,3),(3,-3),(4,4),(4,-4)])
    
    # Symmetrize both
    frame['cp-y-fd-sym'] = frame['cp-y-fd'].__symmetrize__()
    frame['cp-y-td-sym'] = frame['cp-y-td'].__symmetrize__()
    frame['cp-y-td-scrubbed-sym'] = frame['cp-y-td-scrubbed'].__symmetrize__()
    # --- #
    frame['cp-y-fd-sym-z'] = frame['cp-y-fd'].__symmetrize__(zparity=True,__heuristic__=False)
    frame['cp-y-td-sym-z'] = frame['cp-y-td'].__symmetrize__(zparity=True,__heuristic__=False)
    frame['cp-y-td-scrubbed-sym-z'] = frame['cp-y-td-scrubbed'].__symmetrize__(zparity=True,__heuristic__=False)
    
    # # -- Do the same using old workflow for comparison -- #
    
    # # Symmetrize the psi4 time domain coprecessing frame waveform, and return to the init-j frame
    # frame['star-init-j'] = gwylmo_cpclean( frame['init-j'], cp_domain='td' )
    # 
    # # Calculate the coprecessing frame for the case above
    # # Compute TD adn FD coprecessing psi4 frames
    # frame['star-cp-y-fd-sym'] = frame['star-init-j'].__calc_coprecessing_frame__( transform_domain='fd', kind='psi4' )
    # frame['star-cp-y-td-sym'] = frame['star-init-j'].__calc_coprecessing_frame__( transform_domain='td', kind='psi4' )
    
    # Produce diagnostic plots 
    # Produce diagnostic plots 
    def plot_amp_dphi(frame,l=2,m=2):
        
        alert((l,m),header=True)

        #
        fig = figure( figsize=4*figaspect(0.8) )

        #
        y0,y1 = -inf,1e4
        kind = 'psi4' # kind used to measure smoothness of phase derivate
        D0 = mean( frame[save_keys[0]][l,m][kind].fd_dphi )
        smoothness_measure = {}
        # case_test = lambda k: ('cp' in k) and ( not ('init' in k) )
        case_test = lambda k: ('sym' in k) and ('cp' in k) and ( not ('star' in k) ) # and ('fd' in k)
        
        #
        key_list = list(filter( case_test, sort(list(frame.keys())) ))

        #
        for k in key_list:

            #
            f = frame[k].f
            mask = ((f)>0.02*max(m,2)/2) & ((f)<0.06*max(m,2)/2)
            this_smoothness_measure = abs( std( frame[k][l,m][kind].fd_dphi[mask]-smooth(frame[k][l,m][kind].fd_dphi[mask],width=60).answer ) )
            smoothness_measure[k] = this_smoothness_measure 

        #
        ymin = inf
        ymax = -inf
        qnm_lines_applied = False
        
        #
        for k in key_list:
            
            #
            is_best = smoothness_measure[k]==min(list(smoothness_measure.values()))
            if is_best: alert(k)

            #
            alp = 1 if ('td' in k) else 0.5
            ls = '--' if not ('td' in k) else '-'
            #ls = '-.' if ('td' in k) and not ('star' in k) else ls
            lw = 7 if not ('td' in k) else 2
            lw = 1 if 'star' in k else lw


            # --- #
            ax1 = subplot(2,1,1)
            # --- #


            #
            kind = 'strain'
            f = frame[k].f
            mask = abs(f)<0.1*max(m,2)/2
            #ln = plot( f, frame[k][l,m][kind].fd_amp, label=k, alpha=1 if is_best else 0.2, ls=ls if not is_best else '-', lw=lw if not is_best else 2 )
            ln = plot( f, frame[k][l,m][kind].fd_amp, label= ('*' if is_best else '') + k, alpha=alp, ls=ls, lw=lw )
            ylabel(r'$|\tilde{h}_{%i%i}|$'%(l,m))
            yscale('log')
            xlabel('$fM$')
            xscale('log')
            xlim( 0.008*m/2, 0.115*max(m,2)/2 )

            mask = (f>min(xlim())) & (f<max(xlim()))
            ylim( lim( frame[k][l,m][kind].fd_amp[mask], dilate=1.1, dilate_with_multiply=True ) )
            ymin = min( min(ylim()), ymin )
            ymax = max( max(ylim()), ymax )
            ylim( ymin, ymax )

            if k == key_list[-1]:
                axvline(frame[k][l,m][kind].qnm_prograde_fring,color='dodgerblue',lw=1,alpha=0.5,label=r'$f_\mathrm{Ring}^\mathrm{Prograde}$')
                axvline(frame[k][l,m][kind].qnm_retrograde_fring*-1,color='red',lw=1,alpha=0.5,label=r'$f_\mathrm{Ring}^\mathrm{Retrograde}$')
                if (abs(m)<=(l-1)) and (abs(m)>=2):
                    axvline(frame[k][l-1,m][kind].qnm_prograde_fring,ls='--',color='k',alpha=0.12,lw=4)
                    axvline(frame[k][l-1,m][kind].qnm_retrograde_fring,ls=':',color='k',alpha=0.12,lw=4)
            #
            legend(ncol=3)


            # --- #
            ax2 = subplot(2,1,2)
            # --- #


            f = frame[k].f
            mask = ((f)>0.02*m/2) & ((f)<0.06*max(m,2)/2)
            kind = 'psi4'
            dphi_for_plotting = frame[k][l,m][kind].fd_dphi-D0
            #smoothness_measure = abs( std( frame[k][2,2][kind].fd_dphi[mask]-smooth(frame[k][2,2][kind].fd_dphi[mask]).answer ) )
            if True: #smoothness_measure[k]<10:
                plot( f, dphi_for_plotting, label=('*' if is_best else '') + k, alpha=alp, ls=ls, lw=lw, color=ln[0].get_color() )
                #plot( f, dphi_for_plotting, label=k, alpha=1 if is_best else 0.2, ls=ls if not is_best else '-', lw=lw if not is_best else 2, color=ln[0].get_color() )
                xlabel('$fM$')
                xlim( 0.02*m/2, 0.11*max(m,2)/2 )
                ya,yb = lim( frame[k][2,2][kind].fd_dphi[mask]-D0 )
                y0 = max( ya,y0 )
                y1 = min( yb,y1 )
                b = 0.25*abs(y1-y0)
                ylim( y0-b,y1+b )
            #
            xlabel('$fM$')
            # xscale('log')
            ylabel(r'$\frac{d}{df}\arg(\tilde{\psi}_{%i%i})$'%(l,m))

            if is_best:
                mask = (f>min(xlim())) & (f<max(xlim()))
                ylim( lim( dphi_for_plotting[mask], dilate=0.1 ) )

            # yscale('log')
            mask = (f>min(xlim())) & (f<max(xlim()))
            ylim( lim( dphi_for_plotting[mask], dilate=1.1, dilate_with_multiply=True ) )

            if k == key_list[-1]:
                axvline(frame[k][l,m][kind].qnm_prograde_fring,color='dodgerblue',lw=1,alpha=0.5,label=r'$f_\mathrm{Ring}^\mathrm{Prograde}$')
                axvline(frame[k][l,m][kind].qnm_retrograde_fring*-1,color='red',lw=1,alpha=0.5,label=r'$f_\mathrm{Ring}^\mathrm{Retrograde}$')
                # if (abs(m)<=(l-1)) and (abs(m)>=2):
                #     axvline(frame[k][l-1,m][kind].qnm_prograde_fring,ls='--',color='k',alpha=0.12,lw=4)
                #     axvline(frame[k][l-1,m][kind].qnm_retrograde_fring,ls=':',color='k',alpha=0.12,lw=4)

                #
                qnm_lines_applied = True
            #
            legend(ncol=3)


        #
        subplot(2,1,1)
        title(y_raw.simname)
        
        #
        return fig,[ax1,ax2]
    
    # Save comparison plots for multipole moments of interest
    for l,m in lm_of_interest:
        
        #
        fig,ax = plot_amp_dphi(frame,l,m)
        file_path = data_dir+'%s_l%im%i.png'%(a.simname,l,m)
        alert('Saving diagnostic plot to "%s"'%yellow(file_path))
        savefig( file_path )
        close('all')
    
        # Select and output amplitude and phase data
        f = frame['raw'].f
        mask = (f>0.03) & (f<0.12) 
        
        #
        for save_key in save_keys:

            #
            td_amp  = frame[save_key][l,m]['strain'].fd_amp
            td_dphi = frame[save_key][l,m]['psi4'].fd_dphi
            td_phi = frame[save_key][l,m]['psi4'].fd_phi
            
            #
            shift = min(smooth(td_dphi[mask]).answer)
            td_dphi -= shift
            td_phi  -= f * shift
            td_phi -= td_phi[ sum(f<0.03)-1+argmin(smooth(td_dphi[mask]).answer) ]

            data_array = array([ f, td_amp, td_dphi ]).T

            #
            txt_file_path = data_dir+'%s_%s_l%im%i.txt'%(a.simname,save_key.replace('-','_'),l,m)
            alert('Saving waveform data to "%s"'%yellow(txt_file_path))
            # pickle.dump( data_array, open( file_path, "wb" ) )
            savetxt( txt_file_path, data_array, header='[ f, td_amp, td_dphi ], here td and fd refer to the frame type used; frequencies are positive, waveform info are symmetrized in the psi4 TD/FD coprecessing frame from NR simulation at %s'%frame['raw'].simdir )
    
alert('All done.')