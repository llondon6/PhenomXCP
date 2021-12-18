#!/usr/bin/env python3

# Setup python environment
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import xcp
from xcp import determine_data_fitting_region,calibration_catalog,metadata_dict,advanced_gmvx_plot,template_amp_phase

# --------------------------------------- #
# Preliminaries
# --------------------------------------- #

#Load parameter space fit data
alert('Loading parameter space fit data.')
package_dir = parent( xcp.__path__[0] )
datadir = package_dir + 'data/version2/'
foo_path = datadir+'parameter_space_fits.pickle'
foo = pickle.load( open( foo_path, "rb" ) )

# Load and unpuack physical parameter space
raw_domain = loadtxt(datadir+'fit_initial_binary_parameters.txt')
theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z = raw_domain.T

# --------------------------------------- #
# Plot ans save fits 
# --------------------------------------- #

#
scarecrow = template_amp_phase(0.5, 0.5,zeros(3),zeros(3),ell=2)
parameter_names_in_order = scarecrow.__code__.co_varnames[1:scarecrow.__code__.co_argcount]
fit_object = { k:foo[k] for k in parameter_names_in_order }

#
for key in fit_object:
    
    # Summary figure for internal diagnostics 
    fit_object[key].labels={'python':[key,('u', 'eta', 'delta', 'a1'),''],'latex':['\\'+key,(r'\cos(\theta)', r'\eta', r'\delta', r'a_1'),'']}
    
    # Generate diagnostic figures
    summary_fig,eta_set_fig,theta_set_fig = advanced_gmvx_plot( fit_object[key] )
            
    # Save summary figure
    figure_path = datadir + key+'_fit_diagnostic_1_summary.pdf'
    alert('Saving diagnostic plot to %s'%magenta(figure_path))
    summary_fig.savefig( figure_path, pad_inches=0, bbox_inches = 'tight' )
    
    # Save eta space figure
    figure_path = datadir + key+'_fit_diagnostic_2_eta_sets.pdf'
    alert('Saving diagnostic plot to %s'%magenta(figure_path))
    eta_set_fig.savefig( figure_path, pad_inches=0, bbox_inches = 'tight' )
    
    # Save theta space figure
    figure_path = datadir + key+'_fit_diagnostic_3_theta_sets.pdf'
    alert('Saving diagnostic plot to %s'%magenta(figure_path))
    theta_set_fig.savefig( figure_path, pad_inches=0, bbox_inches = 'tight' )
    
    close('all')
    
   
# --------------------------------------- #
# Generate fit python code 
# --------------------------------------- #

#
code_string = ['\n\n#\ndef generate_model_params(theta,eta,a1):\n\n',
               '\t\'\'\'\n\tHola, soy un codigo escribido por "4b_document_fits.py". \n\t~londonl@mit.edu/pilondon2@gmail.com 2020\n\t\'\'\'  \n\n',
               '\t# Import usefuls\n',
               '\tfrom numpy import cos, sqrt\n\n',
               '\t# Preliminaries\n',
               '\tu = cos(theta)\n',
               '\tu2 = u*u\n', 
               '\tu3 = u2*u\n', 
               '\tu4 = u3*u\n', 
               '\teta2 = eta*eta\n', 
               '\teta3 = eta2*eta\n' 
               '\tdelta = sqrt(1-4*eta)\n',
               '\tdelta2 = delta*delta\n', 
               '\tdelta3 = delta2*delta\n\n' 
              ]

# determine list of fitted variables and sort
fit_var = sort( [ k for k in foo.keys() if ('mu' in k)or('nu' in k)or('zeta' in k) ] )

#
for k in fit_var:
    
    # Store python code for fit
    code_string.append( '\t# %s\n'%k )

    #
    this_code_string = foo[k].__str_python__()
    this_code_string = this_code_string.replace('lambda u,eta,delta,a1: ','')
    this_code_string = this_code_string.replace('u*u*','u2*')
    this_code_string = this_code_string.replace('u2*u*','u3*')
    this_code_string = this_code_string.replace('u3*u*','u4*')
    this_code_string = this_code_string.replace('u2*u2*','u4*')
    this_code_string = this_code_string.replace('eta*eta*','eta2*')
    this_code_string = this_code_string.replace('eta2*eta*','eta3*')
    this_code_string = this_code_string.replace('delta*delta*','delta2*')
    this_code_string = this_code_string.replace('delta2*delta*','delta3*')

    #
    code_string.append( '\t'+this_code_string+'\n\n' )
        
#
code_string.append( '\t#\n' )
code_string.append( '\treturn %s\n'%(','.join(fit_var)) )

# Write fit equations to file 
codedir = package_dir+'xcp/'
code_path = codedir+'parameter_space_fits.py'
alert('Write fit equations to file at %s'%magenta(code_path))
f = open(code_path,'w+')
f.writelines(code_string)
f.close()

#
alert('All done.')
