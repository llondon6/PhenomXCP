#!/usr/bin/env python3

# Setup python environment
from matplotlib.pyplot import *
from numpy import *
from positive import *
from glob import glob
from xcp import determine_data_fitting_region
import xcp


# Preliminaries
# --

# Define data location
package_dir = parent( xcp.__path__[0] )
datadir = package_dir + 'data/version2/'

#
# l,m = 2,2
l,m = 3,3

#
files = glob( datadir+'q*l%im%i.txt'%(l,m) )
files.sort()

#
if l==3:
    files = [ f for f in files if 'q1' not in f  ]

#
data = []
for f in files:
    alert('Loading %s'%red(f))
    data.append(loadtxt(f).T)

#
n = 4
fig,ax = subplots( int(ceil(len(data)/n)), n, figsize=12*figaspect(  0.618 * ceil(len(data)/n) / n  ), sharex=False )
ax = ax.flatten()

#
tight_layout(0,2,4)

#
alert('Plotting')
for k in range(len(data)):
    
    #
    # print( '.',end='')
    print(files[k].split('/')[-1].split('.')[0].split('_l')[0])

    # SELECT AND UNPACK DATA
    # ---
    
    #
    if k==len(ax): break
    sca( ax[k] )

    # Select and unpack
    f,amp_fd,dphi_fd,alpha,beta,gamma = data[k]
    
    # DETERMINE DATA FITTING REGION
    # ---
    calibration_data, dphi_lorentzian_min, f_min, f_max, f_lorentzian_min = determine_data_fitting_region( data[k] )
    
    #
    calibration_f, calibration_amp_fd, calibration_dphi_fd, calibration_alpha, calibration_beta, calibration_gamma = calibration_data.T
    
    # Plot RAW data
    # ---
    domain  = f
    codomain = dphi_fd - dphi_lorentzian_min
    plot(domain,codomain,lw=1,color='k',alpha=0.5,label='raw data')
    
    # Note location of min
    # ---
    
    plot( f_lorentzian_min, dphi_lorentzian_min - dphi_lorentzian_min, marker='o', ms='12', mfc='none', mec='k', mew=2,label='min(dphi_fd)' )
    
    # Plot calibration region
    domain4 = calibration_f
    codomain4 = calibration_dphi_fd - dphi_lorentzian_min
    plot( domain4, codomain4, lw=4, color='r', zorder=-10, alpha=1, ls='--', label='calibration data' )
    axvline( min(domain4), color='k', lw=6, ls=':' , alpha=0.3, label='fit region bounds')
    axvline( max(domain4), color='k', lw=6, ls=':' , alpha=0.3)
    
    # Other annotations
    # ---
    
    #
    xlim( lim(domain4,dilate=0.2) )
    ylim( limy(domain4,codomain4,dilate=0.1) )
    xlabel('$f M$')
    ylabel('$d\phi_{%i%i}/df$'%(l,m))
    
    #
    title(files[k].split('/')[-1].split('.')[0].split('_l')[0],loc='left',size=12)
    legend()

    
#
file_path = datadir+'dphi_fitting_region_diagnostic_l%im%i.pdf'%(l,m)
print('')
alert('Saving batch plot to %s'%magenta(file_path))
savefig(file_path,pad_inches=2, bbox_inches = "tight")
alert('\nDone.')