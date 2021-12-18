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

# #
# for key in fit_object:
    
#     # Summary figure for internal diagnostics 
#     fit_object[key].labels={'python':[key,('u', 'eta', 'delta', 'a1'),''],'latex':['\\'+key,(r'\cos(\theta)', r'\eta', r'\delta', r'a_1'),'']}
    
#     # Generate diagnostic figures
#     summary_fig,eta_set_fig,theta_set_fig = advanced_gmvx_plot( fit_object[key] )
            
#     # Save summary figure
#     figure_path = datadir + key+'_fit_diagnostic_1_summary.pdf'
#     alert('Saving diagnostic plot to %s'%magenta(figure_path))
#     summary_fig.savefig( figure_path, pad_inches=0, bbox_inches = 'tight' )
    
#     # Save eta space figure
#     figure_path = datadir + key+'_fit_diagnostic_2_eta_sets.pdf'
#     alert('Saving diagnostic plot to %s'%magenta(figure_path))
#     eta_set_fig.savefig( figure_path, pad_inches=0, bbox_inches = 'tight' )
    
#     # Save theta space figure
#     figure_path = datadir + key+'_fit_diagnostic_3_theta_sets.pdf'
#     alert('Saving diagnostic plot to %s'%magenta(figure_path))
#     theta_set_fig.savefig( figure_path, pad_inches=0, bbox_inches = 'tight' )
    
#     close('all')
    
   
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

#
for k in parameter_names_in_order:
    
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
code_string.append( '\treturn %s\n'%(','.join(parameter_names_in_order)) )

# Write fit equations to file 
codedir = package_dir+'xcp/'
code_path = codedir+'parameter_space_fits.py'
alert('Write fit equations to file at %s'%magenta(code_path))
f = open(code_path,'w+')
f.writelines(code_string)
f.close() 
   
# --------------------------------------- #
# Generate fit C code 
# --------------------------------------- #


#
header_string = '\n\n/* Header file for IMRPhenomXCP\'s tuned parameters */\n\n' 

#
code_string = ['\n\n// Import usefuls\n', '#include <math.h>\n\n']

#
for k in parameter_names_in_order:
    
    #
    this_model_string = foo[k].__str_python__()+';' 
    this_model_string = this_model_string.replace('lambda u,eta,delta,a1: ','')
    
    #
    this_header_string = '// Header formatting for %s\nstatic double IMRPhenomXCP_%s( double theta, double eta, double a1 )'%(k.upper(),k.upper())
    header_string += this_header_string + ';\n'
    
    #
    this_code_string = ''.join(['\n\n// %s fit implementation \n%s{ \n\n'%(k.upper(),this_header_string),
                        '\t/*\n\tHola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.\n\t*/  \n\n',
                        '\t// Preliminaries\n',
                        '\tdouble u = cos(theta);\n' if 'u' in this_model_string else '',
                        '\tdouble u2 = u*u;\n' if 'u*u' in this_model_string else '', 
                        '\tdouble u3 = u2*u;\n' if 'u*u*u' in this_model_string else '', 
                        '\tdouble u4 = u3*u;\n' if 'u*u*u*u' in this_model_string else '', 
                        '\tdouble a12 = a1*a1;\n' if 'a1*a1' in this_model_string else '',
                        '\tdouble a13 = a12*a1;\n' if 'a1*a1*a1' in this_model_string else '',
                        '\tdouble eta2 = eta*eta;\n' if 'eta*eta' in this_model_string else '', 
                        '\tdouble eta3 = eta2*eta;\n' if 'eta*eta*eta' in this_model_string else '', 
                        '\tdouble delta = sqrt(1-4*eta);\n' if 'delta' in this_model_string else '',
                        '\tdouble delta2 = delta*delta;\n' if 'delta*delta' in this_model_string else '', 
                        '\tdouble delta3 = delta2*delta;\n' if 'delta2*delta' in this_model_string else '' ,
                        '\tdouble %s;\n'%k,
                        '\n\t// Evaluate fit for this parameter\n\t'
                    ])
    
    # # Store python code for fit
    # this_code_string.append( '\t// %s\n'%k )

    #
    this_model_string = this_model_string.replace('(eta)','eta')
    this_model_string = this_model_string.replace('(u)','u')
    this_model_string = this_model_string.replace('(a1)','a1')
    this_model_string = this_model_string.replace('(a1*a1)','a12')
    this_model_string = this_model_string.replace('a1*a1*','a12*')
    this_model_string = this_model_string.replace('*a1*a1','*a12')
    this_model_string = this_model_string.replace('*a12*a1','*a13')
    this_model_string = this_model_string.replace('u*u*','u2*')
    this_model_string = this_model_string.replace('u2*u*','u3*')
    this_model_string = this_model_string.replace('u3*u*','u4*')
    this_model_string = this_model_string.replace('u2*u2*','u4*')
    this_model_string = this_model_string.replace('(eta*eta)','eta2')
    this_model_string = this_model_string.replace('eta*eta*','eta2*')
    this_model_string = this_model_string.replace('*eta*eta','*eta2')
    this_model_string = this_model_string.replace('eta2*eta*','eta3*')
    this_model_string = this_model_string.replace('delta*delta*','delta2*')
    this_model_string = this_model_string.replace('delta2*delta*','delta3*')
    
    #
    this_code_string += this_model_string + '\n\n\t// Return answer\n\treturn %s;'%k + '\n\n} // END of %s fit implementation\n\n'%k.upper()

    #
    code_string.append( this_code_string )
        

# Write fit equations to file 
codedir = package_dir+'xcp/'
code_path = codedir+'parameter_space_fits.c'
alert('Write fit equations to file at %s'%magenta(code_path))
f = open(code_path,'w+')
f.writelines(code_string)
f.close()

# Write fit equations to file 
headerdir = package_dir+'xcp/'
header_path = headerdir+'parameter_space_fits.h'
alert('Write fit equations HEADER to file at %s'%magenta(header_path))
f = open(header_path,'w+')
f.writelines(header_string)
f.close()

#
alert('All done.')
