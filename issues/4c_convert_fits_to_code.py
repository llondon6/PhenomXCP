#!/usr/bin/env python3

# Setup python environment
from matplotlib.pyplot import *
from numpy import *
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import xcp
from xcp import determine_data_fitting_region,calibration_catalog,metadata_dict,advanced_gmvx_plot,template_amp_phase, gc



# --------------------------------------- #
# Moment independent preliminaries
# --------------------------------------- #

# Load and unpack physical parameter space
package_dir = parent( xcp.__path__[0] )
datadir = package_dir + 'data/version4/'
raw_domain = loadtxt(datadir+'fit_initial_binary_parameters_l2m2.txt')
theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z,Mf,Xf = raw_domain.T

# --------------------------------------- #
# Initialize code strings
# --------------------------------------- #

#
py_code_string = ['\n\n#\ndef generate_model_params(theta,eta,a1):\n\n',
            '\t\'\'\'\n\tHola, soy un codigo escribido por "4b_document_fits.py". \n\t~lionel.london@kcl.ac.uk/pilondon2@gmail.com\n\t\'\'\'  \n\n',
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
header_string_preamble = '''
#ifndef _LALSIM_IMR_PHENOMX_PNR_DEVIATIONS_H
#define _LALSIM_IMR_PHENOMX_PNR_DEVIATIONS_H
/*
* Copyright (C) 2021 The University of Amsterdam
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/


/**
* \\author Lionel London
*/

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


// 
#include <math.h>

'''

code_string_preamble = '''
/*
* Copyright (C) 2021 The University of Amsterdam
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifdef __cplusplus
extern "C" {
#endif

/**
* \\author Lionel London
*/

//
#include "LALSimIMRPhenomX_PNR_deviations.h"

#ifndef _OPENMP
#define omp ignore
#endif

#ifndef PHENOMXHMDEBUG
#define DEBUG 0
#else
#define DEBUG 1
#endif

'''

code_string_ending = '''
#ifdef __cplusplus
}
#endif
'''

header_string_ending = '''
#ifdef __cplusplus
}
#endif

#endif
'''

#
header_string = header_string_preamble+'\n/* Header file for IMRPhenomXCP\'s tuned parameters */\n\n' 

#
code_string = [code_string_preamble+'\n// Import usefuls\n', '#include <math.h>\n\n']

#
full_list_of_python_param_names = []
      
# For all pairs of l and m in the global config file
append_flag = False
for ll,mm in gc.lmlist:

    # --------------------------------------- #
    # Moment dependent preliminaries
    # --------------------------------------- #

    #Load parameter space fit data
    foo_path = datadir+'parameter_space_fits_l%im%i.pickle'%(ll,mm)
    alert('Loading parameter space fit data: %s'%cyan(foo_path))
    foo = pickle.load( open( foo_path, "rb" ) )

    # --------------------------------------- #
    # Load fit objects
    # --------------------------------------- #

    #
    scarecrow = template_amp_phase(0.5, 0.5,zeros(3),zeros(3),lm=(ll,mm),include_nu0=True,floor_dphi=False)
    parameter_names_in_order = scarecrow.__code__.co_varnames[1:scarecrow.__code__.co_argcount]
    fit_object = { k:foo[k] for k in parameter_names_in_order }
    
    # --------------------------------------- #
    # Generate fit python code 
    # --------------------------------------- #

    #
    for k in parameter_names_in_order:
        
        # Store python code for fit
        py_code_string.append( '\t# %s, %i terms\n'%(k,len(foo[k].basis_symbols)) )
        alert('# %s for (%i,%i) has %i terms ... '%(k,ll,mm,len(foo[k].basis_symbols)))

        #
        this_py_code_string = foo[k].__str_python__()
        
        this_py_param_name = k+'_l%im%i'%(ll,mm)
        this_py_code_string = this_py_code_string.replace(k,this_py_param_name)
        full_list_of_python_param_names.append(this_py_param_name
        )
        
        this_py_code_string = this_py_code_string.replace('lambda u,eta,delta,a1: ','')
        this_py_code_string = this_py_code_string.replace('u*u*','u2*')
        this_py_code_string = this_py_code_string.replace('u2*u*','u3*')
        this_py_code_string = this_py_code_string.replace('u3*u*','u4*')
        this_py_code_string = this_py_code_string.replace('u2*u2*','u4*')
        this_py_code_string = this_py_code_string.replace('eta*eta*','eta2*')
        this_py_code_string = this_py_code_string.replace('eta2*eta*','eta3*')
        this_py_code_string = this_py_code_string.replace('delta*delta*','delta2*')
        this_py_code_string = this_py_code_string.replace('delta2*delta*','delta3*')

        #
        py_code_string.append( '\t'+this_py_code_string+'\n\n' )
    
    # --------------------------------------- #
    # Generate fit C code 
    # --------------------------------------- #

    #
    for k in parameter_names_in_order:
        
        #
        this_model_string = foo[k].__str_python__()+';' 
        this_model_string = this_model_string.replace('lambda u,eta,delta,a1: ','')
        
        #
        this_header_string = '// Header formatting for %s of (l,m)=(%i,%i) multipole\ndouble IMRPhenomXCP_%s_l%im%i( double theta, double eta, double a1 )'%(k.upper(),ll,mm,k.upper(),ll,mm)
        header_string += this_header_string + ';\n'
        
        #
        this_code_string = ''.join(['\n\n// %s fit implementation \n%s{ \n\n'%(k.upper(),this_header_string),
                            '\t/*\n\tHola, soy un codigo escribido por "4b_document_fits.py". Dat script is geschreven door een mens die niet kan worden genoemd.\n\t*/  \n\n',
                            '\t// Preliminaries\n',
                            '\tdouble u = cos(theta);\n' if 'u' in this_model_string else '',
                            '\tdouble u2 = u*u;\n' if 'u*u' in this_model_string else '', 
                            '\tdouble u3 = u2*u;\n' if 'u*u*u' in this_model_string else '', 
                            '\tUNUSED double u4 = u3*u;\n' if 'u*u*u*u' in this_model_string else '', 
                            '\tdouble a12 = a1*a1;\n' if 'a1*a1' in this_model_string else '',
                            '\tdouble a13 = a12*a1;\n' if 'a1*a1*a1' in this_model_string else '',
                            '\tdouble eta2 = eta*eta;\n' if 'eta*eta' in this_model_string else '', 
                            '\tdouble eta3 = eta2*eta;\n' if 'eta*eta*eta' in this_model_string else '', 
                            '\tdouble delta = sqrt(1-4*eta);\n' if 'delta' in this_model_string else '',
                            '\tdouble delta2 = delta*delta;\n' if 'delta*delta' in this_model_string else '', 
                            '\tUNUSED double delta3 = delta2*delta;\n' if 'delta*delta*delta' in this_model_string else '' ,
                            '\tdouble %s;\n'%k,
                            '\n\t// Evaluate fit for this parameter\n\t'
                        ])

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
        this_code_string += this_model_string + '\n\n\t// Return answer\n\treturn %s;'%k + '\n\n} // END of %s (%i,%i) fit implementation\n\n'%(k.upper(),ll,mm)

        #
        code_string.append( this_code_string )
                 
#
py_code_string.append( '\t#\n' )
py_code_string.append( '\treturn (%s)\n'%(','.join(full_list_of_python_param_names)) )

# Write fit equations to file 
codedir = package_dir+'xcp/'
code_path = codedir+'parameter_space_fits.py'
alert('Write fit equations to file at %s'%magenta(code_path))
f1 = open(code_path,'w+' if not append_flag else 'a' )
f1.writelines(py_code_string)
f1.close() 

#
code_string += code_string_ending
header_string += header_string_ending

# Write fit equations to file 
codedir = package_dir+'xcp/'
code_path = codedir+'parameter_space_fits.c'
alert('Write fit equations to file at %s'%magenta(code_path))
f2 = open(code_path,'w+' if not append_flag else 'a')
f2.writelines(code_string)
f2.close()

# Write fit equations to file 
headerdir = package_dir+'xcp/'
header_path = headerdir+'parameter_space_fits.h'
alert('Write fit equations HEADER to file at %s'%magenta(header_path))
f3 = open(header_path,'w+' if not append_flag else 'a')
f3.writelines(header_string)
f3.close()

#
if not append_flag:
    append_flag = True

#
alert('All done.')
