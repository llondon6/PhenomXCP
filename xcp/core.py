

#
import pickle
from positive import alert, magenta, parent
from numpy import loadtxt, load
from os.path import exists
import xcp
import json

#
package_dir = parent(xcp.__path__[0])
data_dir = package_dir + 'data/'

# Always load catalog list for calibration runs 
xcp_catalog_path = data_dir+'xcp_catalog.pickle'
xcp_catalog = pickle.load( open( xcp_catalog_path, "rb" ) )
alert('Catalog of calibration runs stored to %s'%magenta('"xcp.xcp_catalog"'),fname='xcp.core')

# Always load curated metadata for calibration runs 
metadata_dict_path = data_dir+'metadata_dict.pickle'
metadata_dict = load(metadata_dict_path,allow_pickle=True)
alert('Metadata dictionary for calibration runs stored to %s'%magenta('"xcp.metadata_dict"'),fname='xcp.core')

#
catalog_paper_md_path = data_dir+'catalog_paper_metadata.json'
with open(catalog_paper_md_path, 'r') as f:
    catalog_paper_metadata = json.load(f)
alert('Metadata dictionary for Ed\'s catalog paper stored to %s'%magenta('"xcp.catalog_paper_metadata"'),fname='xcp.core')


# Function to determine version2 data fitting region
def determine_data_fitting_region(data, fmin=0.03, fmax=0.12):
    '''
    Given version2 data array, determine fitting region dynamically.
    This function assumes clean data within the default or input values of fmin and fmax, and then uses the phase derivate to determine new fmin and fmax values that are ultimately used to define a fitting region.
    '''

    # Import usefuls
    from numpy import argmin, log
    from positive import smooth, find, lim
    from matplotlib.pyplot import figure, plot, show, axhline, xlim, ylim

    # Extract data
    f, amp_td, amp_fd, dphi_td, dphi_fd, phi_td, phi_fd = data

    # Use default domain bounds to determine a mask
    mask = (f >= fmin) & (f <= fmax)

    # Determine the minimum dphi
    # Smooth dphi using postiive's savgol filter
    x = log(f[mask])
    y = smooth(dphi_td[mask]).answer
    knot = argmin(y)
    y_knot = y[knot]
    data[3] = dphi_td - smooth(dphi_td[mask]).answer[knot] + y_knot
    data[4] = dphi_fd - smooth(dphi_fd[mask]).answer[knot] + y_knot

    # Determine new fmin and max using heuristic
    f_knot = f[mask][knot]
    new_fmin = max(f_knot * 0.22, 0.018)  # 0.5 # 0.325
    new_fmax = f_knot + 0.020  # 0.025

    #
    new_mask = (f >= new_fmin) & (f <= new_fmax)
    new_data = data.T[new_mask, :]

    #
    new_knot = find(f >= fmin)[0]+knot

    #
    return new_data, new_knot, new_fmin, new_fmax, f_knot


#
def select_scenty_metadata( sceo ):
    
    '''
    Given nrutils' scentry object, collect metedata useful for generating model waveforms 
    '''
    
    #
    from numpy.linalg import norm
    from numpy import arccos,dot,pi,array
    from positive.physics import calc_chi_eff,calc_chi_p
    
    #
    a = sceo
    
    #
    l = a.L/norm(a.L)
    if (abs(a.m1-a.m2)<1e-3) and (norm(a.X1)<norm(a.X2)):
        a.X1,a.X2 = [ array(k) for k in (a.X2,a.X1) ]
        a.m1,a.m2 = [ float(k) for k in (a.m2,a.m1) ]

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
    return theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2 

   
# Advanced gloss atop mvpolyfit.plot and mvrfit.plot
def advanced_gmvx_plot( fit_object ):
    
    '''
    Advanced gloss atop mvpolyfit.plot and mvrfit.plot
    '''
    
    from matplotlib.pyplot import subplots, plot, xlabel, ylabel, title, sca, gca, figaspect, tight_layout
    from numpy import cos,sin,array,around,ones_like,sort,pi,linspace
    from positive import eta2q,q2eta,eta2delta
    from glob import glob
    from pwca import determine_data_fitting_region,pwca_catalog,metadata_dict
    
    # Load and unpuack physical parameter space
    raw_domain = loadtxt(data_dir+'version2/fit_intial_binary_parameters.txt')
    theta,m1,m2,eta,delta,chi_eff,chi_p,chi1,chi2,a1,a2 = raw_domain.T


    # Define desired model domain variables and array 
    u = cos(theta)
    v = sin(theta)
    q = 1.0/eta2q(eta)

    # --------------------------------------- #
    # Plot ans save fits 
    # --------------------------------------- #

    # Collect set of unique a1 values
    a1_point = around(a1,2)
    a1_set = array(sort(list( set(a1_point) )))

    # Collect set of unique angle values
    degree_point = (theta*180/pi).astype(int)
    theta_point = degree_point*pi/180
    theta_set = array(sort(list( set(theta_point) )))

    # Collect set of unique mass-ratio values
    q_point = around(array([eta2q(n) for n in eta]),2)
    q_set = array(sort(list( set(q_point) )))

    # Collect set of unique eta values
    eta_point = q2eta( q_point )
    eta_set = q2eta(q_set)

    # Summary figure for internal diagnostics 
    summary_fig = fit_object.plot(size_scale=1.5)
    ax = summary_fig.axes
    
    #
    sca(ax[0])
    title(fit_object.labels['python'][0])
    
    #
    num_figs = len(a1_set)*len(theta_set)
    eta_set_figs,set_fig_ax = subplots( len(a1_set), len(theta_set), figsize=5*array([ len(theta_set),len(a1_set) ]) )
    set_fig_ax = set_fig_ax.flatten();
    tight_layout(4,4)
    ax_counter = 0

    #
    for _a1 in a1_set:
        for _theta in theta_set:

            #
            theta_mask = (_theta==theta_point)
            a1_mask = (_a1==a1_point)
            mask = a1_mask & theta_mask

            #
            _eta = eta_point[mask]
            _u = cos(_theta) 

            #
            case_eta   = linspace( min(_eta),max(_eta),1000 ) 
            case_delta = eta2delta( case_eta )
            case_q     = 1.0/eta2q(case_eta)  
            case_theta = _theta * ones_like(case_eta)
            case_u     = cos(case_theta)
            case_a1    = _a1    * ones_like(case_eta)

            #
            case_domain = array([case_u,case_eta,case_delta,case_a1]).T
            case_range = fit_object.eval(case_domain)
            opt_range  = fit_object.eval(fit_object.domain[mask,:])

            #
            sca(ax[0])
            ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color = 'tab:blue' if _a1==a1_set[0] else 'red' )
            
            #
            sca( set_fig_ax[ax_counter] ); ax_counter += 1
            plot( eta[mask], fit_object.range[mask] if hasattr(fit_object,'range') else fit_object.scalar_range[mask], marker='o',ls='none'  )
            plot( eta[mask], opt_range, marker='o',ms=10,mfc='none', color='r',ls='none'  )
            plot( case_eta, case_range, ls='-', color='r' )
            title( r'$a_1=%1.2f$, $\theta=%1.2f$'%(_a1,round(_theta*180.0/pi,0)) )
            xlabel(r'$\eta$')
            ylabel(r'$\%s$'%fit_object.labels['python'][0])
            

    
    #
    num_figs = len(a1_set)*len(eta_set)
    theta_set_figs,set_fig_ax = subplots( len(a1_set), len(eta_set), figsize=5*array([ len(eta_set),len(a1_set) ]) )
    set_fig_ax = set_fig_ax.flatten();
    tight_layout(4,4)
    ax_counter = 0

    #
    for _a1 in a1_set:
        for _eta in eta_set:

            #
            eta_mask = (_eta==eta_point)
            a1_mask = (_a1==a1_point)
            mask = a1_mask & eta_mask

            #
            _theta = theta_point[mask]
            _u = cos(_theta) 

            #
            case_theta   = linspace( min(_theta),max(_theta),1000 ) # 
            case_u     = cos(case_theta)
            case_eta   = _eta * ones_like(case_theta)
            case_delta = eta2delta( case_eta )
            case_a1    = _a1  * ones_like(case_theta)

            #
            case_domain = array([case_u,case_eta,case_delta,case_a1]).T
            case_range = fit_object.eval(case_domain)
            opt_range  = fit_object.eval(fit_object.domain[mask,:])

            #
            sca(ax[0])
            ax[0].plot3D( case_u, case_eta, case_range, lw=1, alpha=1, color = 'tab:blue' if _a1==a1_set[0] else 'red' )
            
            #
            sca( set_fig_ax[ax_counter] ); ax_counter += 1
            plot( cos(theta[mask]), fit_object.range[mask] if hasattr(fit_object,'range') else fit_object.scalar_range[mask], marker='o',ls='none',color='r'  )
            plot( cos(theta[mask]), opt_range, marker='o',ms=10,mfc='none', color='b',ls='none'  )
            plot( cos(case_theta), case_range, ls='-', color='b' )
            title( r'$a_1=%1.2f$, $q=%1.2f$'%(_a1,eta2q(_eta)) )
            xlabel(r'$\cos(\theta)$')
            ylabel(r'$\%s$'%fit_object.labels['python'][0])
            
    #
    return summary_fig, eta_set_figs, theta_set_figs
 