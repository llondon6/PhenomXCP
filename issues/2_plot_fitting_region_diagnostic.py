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

#
package_dir = parent( xcp.__path__[0] )
datadir = package_dir + 'data/version2/'

#
l,m = 2,2

#
files = glob( datadir+'q*l%im%i.txt'%(l,m) )
files.sort()

#
data = []
for f in files:
    alert('Loading %s'%red(f))
    data.append(loadtxt(f).T)

#
fig,ax = subplots( int(len(data)), 2, figsize=12*figaspect(len(data)*0.618/2), sharex=True )
ax = ax.flatten()

#
tight_layout(1,2,10)

#
alert('Plotting ...')
for k in range(len(data)):
    
    #
    print( '.',end='')

    #
    if k==len(ax): break
    sca( ax[k] )

    #
    f,amp_fd,dphi_fd,alpha,beta,gamma = data[k]
    
    #
    new_data,new_knot,new_fmin,new_fmax,_ = determine_data_fitting_region(data[k])
    plot( f[new_knot], amp_fd[new_knot], color='k', mfc='none', marker='o', ms=20, mew=4, alpha = 0.15  )
    axvline( new_fmin, color='k', ls='-',lw=8,alpha = 0.15 )
    axvline( new_fmax, color='k', ls='-',lw=8,alpha = 0.15 )

    #
    fmin,fmax = 0.03, 0.12
    mask = (f>=fmin) & (f<=fmax)
    
    #
    x = log(f[mask])
    y = smooth(dphi_fd[mask]).answer
    knot = argmin(y)
    plot( f[mask][knot], amp_fd[mask][knot], color='b', mfc='none', marker='o', ms=10, mew=2  )

    #
    plot( f, amp_fd, color='k', ls='-', label=r'cp-$\psi_4$-fd',lw=2 )
        
    #
    f_new,amp_fd_new,dphi_fd_new,alpha_new,beta_new,gamma_new = new_data.T
    
    #
    plot( f_new, amp_fd_new, alpha=0.8, color='orange', ls='-', label=r'calibration data',lw=2 )

    #
    xlim(0.002,0.2)
    ylim(1e-4,1e2)

    #
    axvline( fmin, color='k', ls=':',lw=1 )
    axvline( fmax, color='k', ls=':',lw=1 )
    
    fmin_new = exp(x[knot]) * 0.5
    fmax_new = exp(x[knot]) + 0.020 # * 1.315
    axvline( fmin_new, color='b', ls='--',lw=2 )
    axvline( fmax_new, color='b', ls='--',lw=2 )

    #
    xscale('log')
    yscale('log')
    
    #
    legend(ncol=2,loc=3)
    ylabel(r'$|\tilde{h}_{22}(f)|$')
    # if (k+1==len(data)) or (k+1==int(len(data)-1)):
    #     xlabel('$fM$')
    # else:
    #     xticks([])
    title(files[k].split('/')[-1].split('.')[0],loc='left',size=12)

#
alert('\nDone.')
file_path = datadir+'amp_fitting_region_diagnostic_l%im%i.pdf'%(l,m)
alert('Saving batch plot to %s'%magenta(file_path))
savefig(file_path,pad_inches=2, bbox_inches = "tight")