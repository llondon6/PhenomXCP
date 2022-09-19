

#
from positive import *
from matplotlib.pyplot import *

#
def collect_nr_data_plotting_helper( ll, gwylmo, euler_alpha_beta_gamma, unwrap_dphi=True ):
    
    '''
    '''
    
    #
    from numpy import unwrap,pi,mean,array
    # from matplotlib.pyplot import subplots,plot,xlabel,ylabel,xscale,yscale,xlim,ylim,legend,axvline,axhline

    #
    fig,ax = subplots( 3,2, figsize=6*figaspect(0.618*1.5) )
    ax = ax.flatten()

    #
    gwfo_strain = gwylmo[ll,ll]['strain']
    gwfo_psi4 = gwylmo[ll,ll]['psi4']

    #
    domain = gwylmo.f
    est_domain_min = (ll/2)*gwylmo.wstart/(2*pi)
    est_domain_max = gwfo_strain.qnm_prograde_fring
    
    #
    output_data = []
    format_tags = []
    
    # *
    output_data.append(domain)
    format_tags.append('f')

    '''
    Plot FD amplitude
    '''

    # ------------ #
    sca(ax[0])
    # ------------ #

    codomain1 = gwfo_strain.fd_amp
    codomain2 = abs(gwylmo[ll,-ll]['strain'].fd_amp)
    
    # *
    output_data.append(codomain1)
    format_tags.append('fd_amp')

    plot( domain, codomain1, label=r'$|\tilde{h}_{%i%i}|$'%(gwfo_strain.l,gwfo_strain.m) )
    plot(-domain, codomain2, color = 'k', zorder=-20, label='$m=%i$'%(-ll), ls='--', lw=3, alpha=0.4 )
    axvline( est_domain_max )
    axvline( est_domain_min )

    xscale('log')
    yscale('log')

    xlim( est_domain_min*0.75, est_domain_max*1.25 )
    ylim( lim(limy(domain, codomain1),dilate=2.0,dilate_with_multiply=True) )
    legend()
    xlabel('$f\,M$')


    '''
    Plot FD phase derivative. NOTE that here we use psi4.
    '''

    # ------------ #
    sca(ax[1])
    # ------------ #

    #
    x0,x1 = est_domain_min*2, est_domain_max*1.2
    center = lambda D,X: X # - mean(X[(D>=x0) & (D<=x1)])
    
    # #
    # unwrap_function = unwrap
    
    #
    if unwrap_dphi:
        unwrap_function = unwrap 
    else:
        unwrap_function = lambda x,**kwargs: x

    # 
    codomain1 = center( domain, unwrap_function(gwfo_psi4.fd_dphi,discont=4000) )
    codomain2 = center(-domain, unwrap_function(gwylmo[ll,-ll]['psi4'].fd_dphi,discont=4000) )
    
    #
    avg_codomain = 0.5*( codomain1 + spline(domain,codomain2)(-domain) )
    
    # *
    output_data.append( avg_codomain )
    format_tags.append('fd_dphi')

    plot( domain, codomain1, label=r'$d\phi_{%i%i}/df$'%(gwfo_strain.l,gwfo_strain.m) )
    plot(-domain, codomain2, color = 'k', zorder=-20, label='$m=%i$'%(-ll), ls='--', lw=3, alpha=0.4 )
    plot( domain, avg_codomain, label=r'Avg.', color='k' )
    axvline( est_domain_max )
    axvline( est_domain_min )

    xscale('log')

    xlim( x0,x1 )
    ylim( lim(limy(domain, codomain1),dilate=0.1) )
    legend()
    xlabel('$f\,M$')


    '''
    Plot FD alpha
    '''

    # ------------ #
    sca(ax[2])
    # ------------ #

    #
    x0,x1 = est_domain_min*2, est_domain_max*1.1
    center = lambda D,X: X - mean(X[(D>=x0) & (D<=x1)])

    #
    codomain1 = center(domain,-euler_alpha_beta_gamma[0])
    codomain2 = center(-domain,-euler_alpha_beta_gamma[0])
    
    # *
    output_data.append( -euler_alpha_beta_gamma[0] )
    format_tags.append('fd_alpha')

    plot( domain, codomain1, label=r'$\alpha_{%i%i}$'%(gwfo_strain.l,gwfo_strain.m) )
    plot(-domain, codomain2, color = 'k', zorder=-20, label='$f<0$', lw=3, ls='--', alpha=0.4 )
    axvline( est_domain_max )
    axvline( est_domain_min )

    xscale('log')

    xlim( x0,x1 )
    ylim( lim(limy(domain, codomain1),dilate=0.1) )
    legend()
    xlabel('$f\,M$')


    '''
    Plot FD beta
    '''

    # ------------ #
    sca(ax[3])
    # ------------ #

    #
    x0,x1 = est_domain_min*2, est_domain_max*1.1

    #
    codomain1 = euler_alpha_beta_gamma[1]
    codomain2 = euler_alpha_beta_gamma[1]
    
    # *
    output_data.append( codomain1 )
    format_tags.append('fd_beta')

    plot( domain, codomain1, label=r'$\beta_{%i%i}$'%(gwfo_strain.l,gwfo_strain.m) )
    plot(-domain, codomain2, color = 'k', zorder=-20, label='$f<0$', lw=3, ls='--', alpha=0.4 )
    axvline( est_domain_max )
    axvline( est_domain_min )

    xscale('log')

    xlim( x0,x1 )
    ylim( lim(limy(domain, codomain1),dilate=0.1) )
    legend()
    xlabel('$f\,M$')


    '''
    Plot FD gamma
    '''

    # ------------ #
    sca(ax[4])
    # ------------ #

    #
    x0,x1 = est_domain_min*2, est_domain_max*1.1
    center = lambda D,X: X - mean(X[(D>=x0) & (D<=x1)])

    #
    codomain1 = center(domain,-euler_alpha_beta_gamma[2])
    codomain2 = center(-domain,-euler_alpha_beta_gamma[2])
    
    # *
    output_data.append( codomain1 )
    format_tags.append('fd_gamma')

    plot( domain, codomain1, label=r'$\gamma_{%i%i}$'%(gwfo_strain.l,gwfo_strain.m) )
    plot(-domain, codomain2, color = 'k', zorder=-20, label='$f<0$', lw=3, ls='--', alpha=0.4 )
    axvline( est_domain_max )
    axvline( est_domain_min )

    xscale('log')

    xlim( x0,x1 )
    ylim( lim(limy(domain, codomain1),dilate=0.1) )
    legend()
    xlabel('$f\,M$')

    '''
    Axis not used
    '''

    # ------------ #
    sca(ax[5])
    # ------------ #
    plot([0,1],[0,1],color='k')
    plot([0,1],[1,0],color='k')
    xlim(0,1)
    ylim(0,1)
    xticks([])
    yticks([])
    
    #
    output_data = array(output_data).T
    
    #
    return fig, ax, output_data, format_tags

#
def collect_nr_data_helper_v3( ll, gwylmo_fd, gwylmo_td, euler_alpha_beta_gamma_fd, euler_alpha_beta_gamma_td ):
    
    '''
    '''
    
    #
    from numpy import unwrap,pi,mean,array

    #
    domain = gwylmo_fd.f
    
    #
    output_data = []
    format_tags = []
    
    # * 0
    output_data.append(domain)
    format_tags.append('f')
    
    # Amplitudes
    # ---
    
    # * 1
    codomain = gwylmo_fd[ll,ll]['strain'].fd_amp
    output_data.append(codomain)
    format_tags.append('fd_amp')
    
    # * 2
    codomain = gwylmo_td[ll,ll]['strain'].fd_amp
    output_data.append(codomain)
    format_tags.append('td_amp')
    
    # Strain phases
    # ---
    
    # * 3
    codomain = gwylmo_fd[ll,ll]['strain'].fd_phi
    output_data.append(codomain)
    format_tags.append('strain_fd_phase')
    
    # * 4
    codomain = gwylmo_td[ll,ll]['strain'].fd_phi
    output_data.append(codomain)
    format_tags.append('strain_td_phase')
    
    # Psi4 phases
    # ---
    
    # * 5
    codomain = gwylmo_fd[ll,ll]['psi4'].fd_phi
    output_data.append(codomain)
    format_tags.append('psi4_fd_phase')
    
    # * 6
    codomain = gwylmo_td[ll,ll]['psi4'].fd_phi
    output_data.append(codomain)
    format_tags.append('psi4_td_phase')
    
    # FD Angles
    # ---
    
    # * 7
    codomain = -euler_alpha_beta_gamma_fd[0]
    output_data.append(codomain)
    format_tags.append('fd_alpha')
    
    # * 8
    codomain = euler_alpha_beta_gamma_fd[1]
    output_data.append(codomain)
    format_tags.append('fd_beta')
    
    # * 9
    codomain = euler_alpha_beta_gamma_fd[2]
    output_data.append(codomain)
    format_tags.append('fd_gamma')
    
    #
    output_data = array(output_data).T
    
    #
    return output_data, format_tags