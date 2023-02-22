#!/usr/bin/env python3

# Setup python environment
from matplotlib.pyplot import *
from numpy import *
import dill,pickle
from positive import *
from nrutils import scsearch, gwylm
from glob import glob
import xcp
from xcp import determine_data_fitting_region,calibration_catalog,metadata_dict,template_amp_phase,advanced_gmvx_plot,gc

# --------------------------------------- #
# Preliminaries
# --------------------------------------- #

#
alert('Loading parameter space data.')

# Define data location
package_dir = parent( xcp.__path__[0] )
datadir = package_dir + 'data/version4/'

# Load and unpuack physical parameter space
raw_domain = loadtxt(datadir+'fit_initial_binary_parameters_l2m2.txt')
theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2,chi1_x,chi1_y,chi1_z,chi2_x,chi2_y,chi2_z,Mf,Xf = raw_domain.T


# Define desired model domain variables and array 
u = cos(theta)
raw_model_domain = array( [ u, eta, delta, a1 ] ).T

# For all pairs of l and m in the global config file
for ll,mm in gc.lmlist:

    # Load and unpack physical parameter space
    opt_parameter_range = loadtxt(datadir+'fit_opt_parameters_l%im%i.txt'%(ll,mm))
    scarecrow = template_amp_phase(0.5, 0.5,zeros(3),zeros(3),lm=(ll,mm),include_nu0=True,floor_dphi=False)
    parameter_names_in_order = scarecrow.__code__.co_varnames[1:scarecrow.__code__.co_argcount]
    model_range = {  parameter_names_in_order[k]:var for k,var in enumerate(opt_parameter_range.T) }
    # if ll==2:
    #     del model_range['mu4']
    
    #
    if mod(mm,2) == 0:
        model_domain = raw_model_domain
    else:
        mask = around( m1/m2 ) != 1
        model_domain = raw_model_domain[ mask, : ]
        

    #
    foo = {}
    foo['model_domain_header'] = '# columns are: cos(theta), eta, delta, a1 '
    foo['model_domain'] = model_domain

    # --------------------------------------- #
    # Fit dphase parameters 
    # --------------------------------------- #
    '''

    * NOTE that function setting have been determined by manual play

    * NOTE that some parameters are modeled with a rational function and other with a   polynomial

    * NOTE that each model allows a maximum polynomial degree determined by the number of destinct data points in the respective dimension

    '''
    alert('Fitting parameters ...',header=True)


    # # mu1 
    # # ---
    # key = 'mu1'
    # labels={'python':[key,('u', 'eta', 'delta', 'a1'),''],'latex':['\\'+key,(r'\cos(\theta)', r'\eta', r'\delta', r'a_1'),'']}
    # foo[key] = gmvpfit( model_domain, mu1,fitatol=0.0001,verbose=True,maxdeg_list=[4,3,0,3],center=True,labels=labels,estatol=0.01)
    # alert( 'The fit has '+bold(str(len(foo[key].basis_symbols)))+' terms.', header=True )
    # advanced_gmvx_plot(foo[key]); show()

    #
    for key in model_range:
        
        if mod(mm,2) == 0:
            #
            u_order = 3
            a1_order = 2
            mass_ratio_quantity_order = 3
        else:
            #
            u_order = 4 # 3
            a1_order = 3
            mass_ratio_quantity_order = 2
        
        #
        maxdeg_list = [u_order,mass_ratio_quantity_order,0,a1_order]
        temper = not True
        
        # if key in ('nu4'): # if key in ('nu4','mu4'):
        #     maxdeg_list = [u_order,0,mass_ratio_quantity_order,a1_order]
        
        # if key in ('mu3'):
        #     maxdeg_list = [3,mass_ratio_quantity_order,0,a1_order]
        
        alert('Now fitting %s'%red(key),header=True)
        
        #
        labels={'python':[key,('u', 'eta', 'delta', 'a1'),''],'latex':['\\'+key[:-1]+'_'+key[-1],(r'\cos(\theta)', r'\eta', r'\delta', r'a_1'),'']}
        
        # # --- Implement very basic outlier detection and removal --- #
        # if (key == 'mu4') and (ll==3):
        #     mu,sig = mean(model_range[key]),std(model_range[key])
        #     a_ = 2.0
        #     mask =  (model_range[key] >= mu+a_*sig) | (model_range[key] <= mu-a_*sig)
        #     model_range[key][ mask ] = mu
        #     # print(len(mask),sum(mask),mu-a_*sig,mu+a_*sig,lim(model_range[key]),mu,sig)
        #     # error('debugging')
        # # ---------------------------------------------------------- #
        
        foo[key] = gmvpfit( model_domain, model_range[key],fitatol=0.0001,verbose=True,maxdeg_list=maxdeg_list,center=True,estatol=0.0001,labels=labels,temper=temper)
        alert( 'The fit for %s has '%bold(yellow(key))+bold(yellow(str(len(foo[key].basis_symbols))))+' terms.', header=True )

    # # nu4
    # # ---
    # key = 'nu4'
    # foo[key] = gmvpfit( model_domain, nu4,fitatol=0.0001,verbose=True,maxdeg_list=[3,1,1,1],center=True,estatol=0.015)

    # # nu5
    # # ---
    # key = 'nu5'
    # foo[key] = gmvpfit( model_domain, nu5,fitatol=0.0001,verbose=True,maxdeg_list=[2,1,1,1],center=True)

    # # nu6
    # # ---
    # key = 'nu6'
    # foo[key] = gmvpfit( model_domain, nu6,fitatol=0.0001,verbose=True,maxdeg_list=[4,3,0,1],center=True,estatol=0.03)

    # # zeta2
    # # ---
    # key = 'zeta2'
    # foo[key] = gmvpfit( model_domain, zeta2,fitatol=0.0001,verbose=True,maxdeg_list=[4,3,0,1],center=True,estatol=0.05)

    # # --------------------------------------- #
    # # Fit amplitude parameters 
    # # --------------------------------------- #
    # '''

    # * NOTE that function setting have been determined by manual play

    # * NOTE that some parameters are modeled with a rational function and other with a   polynomial

    # * NOTE that each model allows a maximum polynomial degree determined by the number of destinct data points in the respective dimension

    # '''
    # alert('Fitting amplitude parameters ...',header=True)

    # # # mu1 
    # # # ---
    # # key = 'mu1'
    # # foo[key] = gmvpfit( model_domain, mu1,fitatol=0.001,verbose=True,maxdeg_list=[4,0,3,1],center=True,estatol=0.015)
    # # mu2
    # # ---
    # key = 'mu2'
    # foo[key] = gmvpfit( model_domain, mu2,fitatol=0.001,verbose=True,maxdeg_list=[4,3,0,1],center=True)
    # # # mu3
    # # # ---
    # # key = 'mu3'
    # # foo[key] = gmvpfit( model_domain, mu3,fitatol=0.0001,verbose=True,maxdeg_list=[4,3,0,1],center=True,estatol=0.005)
    # # mu4
    # # ---
    # key = 'mu4'
    # foo[key] = gmvpfit( model_domain, mu4,fitatol=0.0001,verbose=True,maxdeg_list=[3,1,1,1],center=True)

    # --------------------------------------- #
    # Saving fit data
    # --------------------------------------- #
    data_path = datadir+'parameter_space_fits_l%im%i.pickle'%(ll,mm)
    alert('Saving '+yellow('parameter_space_fits_l%im%i'%(ll,mm))+' to %s'%magenta(data_path))
    pickle.dump( foo, open( data_path, "wb" ) )

#
alert('All done!',header=True)