'''
Functions for visualising abnormalities and D_RS in continuous iEEG data.

Gabrielle M. Schroeder
CNNP Lab, Newcastle University
May 2023
'''

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sklearn.metrics as sk_metrics
import scipy.stats
import pandas as pd

def heatmap_abnormalities(roi_data,t_days,roi_names,roi_is_resect=None,ax=None,
                          center=None,cmap='viridis',resect_line_clr=(0,0,0),
                          cbar_kws=None,vmin=None,vmax=None):
    '''
    Plot a heatmap of ROI abnormalities (ROIs x time heatmap). Time units are
    hours. Will also sort ROIs by resected/spared and divide categories with 
    a horizontal line if labels are provided.
    
    There must be an integer number of time windows in an hour.

    Parameters
    ----------
    roi_data : 2D numpy array, float
        ROI abnormalities, size number of ROIs x time of time windows.
    t_days : 1D numpy array, float
        Times corresponding to roi_data columns, in days.
    roi_names : list of strings
        List of ROI names corresponding to roi_data rows.
    roi_is_resect : 1D numpy array, boolean, optional
        Whether each ROI was resected. The default is None. If None, ROIs are 
        not sorted by spared/resected in the heatmap.
    ax : matplotlib axes, optional
        Axes in which to plot heatmap. The default is None. If None, new figure
        is generated.
    center : int, float, or None, optional
        Value to put at center of heatmap colorbar. The default is None. If None,
        default colorbar limits and centering are used.
    cmap : str, optional
        Heatmap colormap. The default is 'viridis'.
    resect_line_clr : matplotlib color (e.g., tuple), optional
        Color of line separating resected and spared ROIs. 
        The default is (0,0,0).

    Returns
    -------
    ax : axes object
        Axes in which heatmap is plotted.
    idx : 1D numpy array, int
        indices for sorting ROIs 
    cax: axes object
        handle to colourbar

    '''
    # number of ROIs
    n_roi,n_win = roi_data.shape
    
    # create new figure if no axes passed to function
    if ax is None:
        fig,ax = plt.subplots()
    
    # indices for sorting ROIs by resected and spared, if specified
    # resected ROIs will be at top of plot
    if roi_is_resect is not None:
        idx = np.argsort(np.invert(roi_is_resect))
    else:
        idx = np.arange(n_roi)
        
    # sort ROIs
    sorted_roi = [roi_names[i] for i in idx]
    
    # axes for colorbar
    cax = inset_axes(ax,width='5%',height='100%',
                      loc=3,
                      bbox_to_anchor=(1,0,1,1),
                      bbox_transform=ax.transAxes)
    # heatmap
    sns.heatmap(roi_data[idx,:],center=center,cmap=cmap,
                      yticklabels=sorted_roi,ax=ax,
                      cbar_ax=cax,rasterized=True,
                      cbar_kws=cbar_kws,vmin=vmin,vmax=vmax)

    # labels
    win_per_hr = 1/((t_days[1] - t_days[0])*24)        # number of windows per hr
    if win_per_hr.is_integer():
        win_per_hr = int(win_per_hr)
    else:
        raise ValueError('There must be an integer number of windows per hour')
        
    tick_x = np.arange(0,n_win,win_per_hr*12)   # make tick mark every 12 hrs
    ax.set_xticks(tick_x)                       # set x ticks
    ax.set_xticklabels(t_days[tick_x],rotation=0)    # x tick labels
    ax.set_xlabel('time (days)')                 # x axis label
    ax.set_ylabel('ROI')                        # y axis label
    # frame
    for _,spine in ax.spines.items():
        spine.set_visible(True)
        
    # mark resected vs spared, if specified
    if roi_is_resect is not None:
        n_resect = np.sum(roi_is_resect)
        ax.axhline(n_resect,color=resect_line_clr)
        plt.setp(ax.get_yticklabels()[(n_resect):(n_roi)],color='grey') # spared ROIs labelled in grey

    # return axis 
    return ax, idx, cax
    

def plot_D_RS(t_days,d_rs,ax=None,clr=(92/255,52/255,127/255),alph=0.5):
    '''
    Plot D_RS vs time as lineplot. 

    Parameters
    ----------
    t_days : 1D numpy array, float
        Elapse time (in days) of each time window.
    d_rs : 1D numpy array, float
        D_RS of each time window (same length as t_days).
    ax : axes object, optional
        Axes object for plot; if none provided, plots new figure. The default is None.
    clr : tuple, length 3, optional
        Colour of D_RS line plot. The default is (92/255,52/255,127/255).

    Raises
    ------
    ValueError
        Raises error if there is not an integer number of windows per hour
        (e.g., 120 windows/hr, equivalent to a 30s time window)

    Returns
    -------
    ax : axes object
        Axes object for plot.

    '''
    
    # number of time points (windows)
    n_win = len(t_days)
    
    # make new figure if axes not supplied
    if ax is None:
        fig,ax = plt.subplots()
    
    # plot
    ax.plot(t_days,d_rs,linewidth=0.5,color=clr,alpha=alph)
    
    # limits
    ax.set_xlim(min(t_days),max(t_days))
    ax.set_ylim(0,1)
    
    # labels
    win_per_hr = 1/((t_days[1] - t_days[0])*24)        # number of windows per hr
    if win_per_hr.is_integer():
        win_per_hr = int(win_per_hr)
    else:
        raise ValueError('There must be an integer number of windows per hour')
        
    ax.set_ylabel(r'$D_{RS}$')
    ax.set_title(r'$D_{RS}$')
    tick_x = np.arange(0,n_win,win_per_hr*12)   # make tick mark every 12 hrs
    ax.set_xticks(t_days[tick_x])                # set x ticks
    ax.set_xlabel('time (days)')                 # x axis label
    
    # horizontal line
    ax.axhline(0.5,color=(0,0,0),linewidth=0.25)
    
    return ax


def plot_outcome_auc(good_outcome,pnt_meas, ylim=[0,1],meas_name='$D_{RS}$',
                          outcome_clrs=None,meas_ref=0.5,ax=None,size=7,
                          p_side='less',overlay_violin=True,
                          violin_tint=0.9):
    '''
    Compares measure in patients with good vs poor surgical outcome.
    Plots beeswarm with option violin plot.
    Computes AUC (positive label = bad outcome)

    Parameters
    ----------
    good_outcome : 1D numpy array, boolean
        Which patients had a good surgical outcome.
    pnt_meas : 1D numpy array, float
        Measure for each patient, e.g., D_RS.
    ylim : 1D numpy array, length 1, optional
        y axis limits for plot. The default is [0,1].
    meas_name : string, optional
        Name of measure, used as label for plot y axis. The default is '$D_{RS}$'.
    outcome_clrs : 2x3 numpy array, optional
        Colors used to label good and bad outcome. The default is None.
    meas_ref : float, optional
        Where to draw reference line on y-axis. The default is 0.5.
    ax : axes object, optional
        Axis in which to plot beewswarm; will make new plot if not provided. The default is None.
    size : int, optional
        Size of beeswarm points. The default is 7.
    p_side : string, optional
        H0 alternative type for wilcoxon rank sum test. The default is 'less'.
    overlay_violin : boolean, optional
        Whether to overlay a violin plot. The default is True.
    violin_tint : float, optional
        Number between 0 and 1 to lighten colour for violin plot. The default is 0.9.

    Returns
    -------
    ax : axes object
        Plot axis.
    auc : float
        AUC from using measure to classify outcome.
    p : float
        p-value of AUC (using wilcoxon rank sum test).

    '''
    # first remove any nans
    is_nan = np.isnan(pnt_meas)
    pnt_meas = pnt_meas[np.invert(is_nan)]
    good_outcome = good_outcome[np.invert(is_nan)]
    print(f'removed {np.sum(is_nan)} patients due to NaN measure')
    
    # define bad outcome patients
    bad_outcome = good_outcome == 0
    
    if outcome_clrs is None:
        cmap = cm.get_cmap('Set1',8)
        outcome_clrs = np.zeros((2,3))
        outcome_clrs[0,:] = cmap(1)[:3]
        outcome_clrs[1,:] = cmap(0)[:3]
    
    # compute AUC
    if len(pnt_meas)>0:
        auc = sk_metrics.roc_auc_score(bad_outcome,pnt_meas)
        _,p = scipy.stats.ranksums(pnt_meas[bad_outcome==0],
                                  pnt_meas[bad_outcome==1],
                                  alternative=p_side)
    else:
        auc = np.nan
        p = np.nan
        
    # PLOT
    if ax is None:
        fig,ax = plt.subplots(1,1,figsize=(4,5))
    
    ax.set_ylim(ylim)
    ax.axhline(meas_ref,color = (0.5,0.5,0.5))
    if len(pnt_meas)>0:
        
        # make dataframe
        dct = {'bad_outcome':bad_outcome,'pnt_meas':pnt_meas}
        df = pd.DataFrame(dct)
            
        # violin
        if overlay_violin:
            # lighten clrs
            v_clrs = outcome_clrs*255
            v_clrs = ((255-v_clrs)*violin_tint)+v_clrs
            v_clrs = v_clrs/255
            
            sns.violinplot(data=df,x='bad_outcome',y='pnt_meas',
                            palette=v_clrs,ax=ax,saturation=1,
                            inner='quartiles',cut=0.5)  

        
        # beeswarm
        sns.swarmplot(data=df,x='bad_outcome',y='pnt_meas',hue='bad_outcome',
                      size=size,
                      palette=outcome_clrs,ax=ax)            
            
            
        ax.get_legend().remove()
        ax.set_xticklabels(['ILAE 1','ILAE 2+'])


    ax.set_title(f'{meas_name} \n AUC = {round(auc,2)} \n p = {round(p,4)} ({p_side}) \n n = {len(pnt_meas)}')
    ax.set_ylabel(meas_name)
    ax.set_xlabel('surgical outcome')

    return ax,auc,p