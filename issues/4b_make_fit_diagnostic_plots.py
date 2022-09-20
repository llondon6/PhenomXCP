#!/usr/bin/env python3

# Setup python environment
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import xcp
from xcp import determine_data_fitting_region,calibration_catalog,metadata_dict,advanced_gmvx_plot,template_amp_phase, gc

# # Load and unpuack physical parameter space
# raw_domain = loadtxt(datadir+'fit_initial_binary_parameters.txt')
# theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z = raw_domain.T
    
# For all pairs of l and m in the global config file
for ll,mm in gc.lmlist:

    # --------------------------------------- #
    # Preliminaries
    # --------------------------------------- #

    #Load parameter space fit data
    alert('Loading parameter space fit data.')
    package_dir = parent( xcp.__path__[0] )
    datadir = package_dir + 'data/version4/'
    foo_path = datadir+'parameter_space_fits_l%im%i.pickle'%(ll,mm)
    foo = pickle.load( open( foo_path, "rb" ) )

    # --------------------------------------- #
    # Plot ans save fits 
    # --------------------------------------- #

    #
    scarecrow = template_amp_phase(0.5, 0.5,zeros(3),zeros(3),lm=(ll,mm),include_nu0=True,floor_dphi=False)
    parameter_names_in_order = scarecrow.__code__.co_varnames[1:scarecrow.__code__.co_argcount]
    # if ll==2:
    #     parameter_names_in_order = tuple([ parameter_names_in_order[j] for j,k in enumerate(parameter_names_in_order) if k!='mu4' ])
    fit_object = { k:foo[k] for k in parameter_names_in_order }

    #
    for key in fit_object:
        
        # Summary figure for internal diagnostics 
        fit_object[key].labels={'python':[key,('u', 'eta', 'delta', 'a1'),''],'latex':['\\'+key,(r'\cos(\theta)', r'\eta', r'\delta', r'a_1'),'']}
        
        # Generate diagnostic figures
        summary_fig,eta_set_fig,theta_set_fig,a1_set_fig = advanced_gmvx_plot( fit_object[key], ll, mm )
                
        # Save summary figure
        figure_path = datadir + key+'_fit_diagnostic_1_summary_l%im%i.pdf'%(ll,mm)
        alert('Saving diagnostic plot to %s'%magenta(figure_path))
        summary_fig.savefig( figure_path, pad_inches=0, bbox_inches = 'tight' )
        
        # Save eta space figure
        figure_path = datadir + key+'_fit_diagnostic_2_eta_sets_l%im%i.pdf'%(ll,mm)
        alert('Saving diagnostic plot to %s'%magenta(figure_path))
        eta_set_fig.savefig( figure_path, pad_inches=0, bbox_inches = 'tight' )
        
        # Save theta space figure
        figure_path = datadir + key+'_fit_diagnostic_3_theta_sets_l%im%i.pdf'%(ll,mm)
        alert('Saving diagnostic plot to %s'%magenta(figure_path))
        theta_set_fig.savefig( figure_path, pad_inches=0, bbox_inches = 'tight' )
        
        # Save a1 space figure
        figure_path = datadir + key+'_fit_diagnostic_4_a1_sets_l%im%i.pdf'%(ll,mm)
        alert('Saving diagnostic plot to %s'%magenta(figure_path))
        a1_set_fig.savefig( figure_path, pad_inches=0, bbox_inches = 'tight' )
        
        close('all')

#
alert('All done.')
