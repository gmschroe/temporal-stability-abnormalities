#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates across patient plots in Figure 4 from Wang et al. (2023) Temporal 
stability of intracranial EEG abnormality maps for localizing epileptogenic 
tissue.

Calculates time-varying D_RS from abnormalities before plotting.

Gabrielle M. Schroeder
CNNP Lab, Newcastle University
May 2023
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import math
import sklearn.metrics as sk_metrics
import scipy.stats
import seaborn as sns


# custom functions
import pyabnorm.common as com
import pyabnorm.compute_abnorm as ab_comp
import pyabnorm.vis_abnorm as ab_vis

#%% directories and list of patients

# directories
data_dir,plot_dir = com.set_paths()

# add subfolder for figure 4 plots
plot_dir = os.path.join(plot_dir,'fig4_across_patients')

# list of patient files
all_fname = os.listdir(data_dir)
n_patients = len(all_fname)

# sort patients by number
pnt_num = [int(''.join(filter(str.isdigit,pnt))) for pnt in all_fname]
pnt_idx = np.argsort(pnt_num)
all_fname = [all_fname[idx] for idx in pnt_idx]

del pnt_num, pnt_idx
#%% compute D_RS for each patient; also save outcomes and interictal/periictal win

# initialise arrays/lists for saving across-patient info
all_pnt_id = []
all_d_rs = [] 
all_good_outcome = np.zeros(n_patients)
all_win_ii = [] # interictal windows
all_win_peri = [] # periictal windows

for i in range(n_patients):

    # load patient data
    pnt_dict = com.mat_load_as_dict(os.path.join(data_dir,all_fname[i]),
                                    'pnt_id','roi_ab','roi_is_resect','good_outcome',
                                    'win_interictal','win_periictal')
    
    print(f'Computing DRS of {pnt_dict["pnt_id"]}')
    
    # resection info as boolean
    pnt_dict['roi_is_resect'] = pnt_dict['roi_is_resect']==1
    
    # compute D_RS 
    d_rs = ab_comp.compute_d_rs(pnt_dict['roi_ab'],pnt_dict['roi_is_resect']) # D_RS vs time

    # save info for plots and outcome calculations
    all_pnt_id.append(pnt_dict['pnt_id'])
    all_good_outcome[i] = pnt_dict['good_outcome']
    all_win_ii.append(pnt_dict['win_interictal'])
    all_win_peri.append(pnt_dict['win_periictal'])
    all_d_rs.append(d_rs)
    
    del pnt_dict, d_rs

# as boolean
all_good_outcome = all_good_outcome == 1

#%% sort patients by outcome
sort_idx = np.argsort(np.invert(all_good_outcome),kind='stable')

all_pnt_id = [all_pnt_id[i] for i in sort_idx]
all_d_rs = [all_d_rs[i] for i in sort_idx]
all_good_outcome = all_good_outcome[sort_idx]
all_win_ii = [all_win_ii[i] for i in sort_idx]
all_win_peri = [all_win_peri[i] for i in sort_idx]

#%% set colours and other vis settings

clr_ii = (178/255,178/255,178/255) # interictal periods
clr_peri = (227/255,124/255,29/255) # periictal periods

outcome_clrs = np.zeros((2,3))
outcome_clrs[1,:] = np.array([0/255,118/255,192/255]) # blue, good outcome 
outcome_clrs[0,:] = np.array([163/255,2/255,52/255]) # red, poor outcome

# fixes text embedding in pdfs
matplotlib.rcParams['pdf.fonttype']=42

# format for saving plots
plot_formats = ['pdf'] 

#%% get D_RS and median D_RS of interictal and periictal periods
# note that one patient does not have periictal data and will therefore be 
# excluded from the across patients analysis

all_d_rs_ii = [] # interictal
all_d_rs_ii_med = np.zeros(n_patients)
all_d_rs_peri = [] # periictal
all_d_rs_peri_med = np.zeros(n_patients)

for i in range(n_patients):
    
    # interictal and periictal D_RS
    d_rs_ii = all_d_rs[i].copy()
    d_rs_ii[all_win_ii[i]==0] = np.nan # remove data that is not interical
    d_rs_peri = all_d_rs[i].copy()
    d_rs_peri[all_win_peri[i]==0] = np.nan # remove data that is not periictal
    
    all_d_rs_ii.append(d_rs_ii)
    all_d_rs_peri.append(d_rs_peri)
    
    # median D_RS of each time period
    all_d_rs_ii_med[i] = np.nanmedian(d_rs_ii)
    all_d_rs_peri_med[i] = np.nanmedian(d_rs_peri)
    

#%% plot interictal and periictal D_RS distributions with median D_RS values marked with vertical line

# number of columns and rows in vis
n_col = 8
n_row = math.ceil(n_patients/n_col)

# start figure
fig = plt.figure(figsize=(12,10))
alph = 0.4 # opacity of distribution
my_bins = np.arange(0,1.1,0.1)

for i in range(n_patients):
    ax = plt.subplot(n_row,n_col,i+1)
    
    ax.hist(all_d_rs_ii[i],bins=my_bins,
        facecolor=clr_ii,alpha=alph)
    ax.hist(all_d_rs_peri[i],bins=my_bins,
        facecolor=clr_peri,alpha=alph)
    ax.set_title(f'{all_pnt_id[i]}\nmedian inter = {round(all_d_rs_ii_med[i],2)}'+
                 f'\nmedian peri = {round(all_d_rs_peri_med[i],2)}',fontsize=8)

    ax.set_xlabel('$D_{RS}$')
    ax.set_xlim(0,1)
    ax.set_xticks((0,0.5,1))
    ax.axvline(all_d_rs_ii_med[i],color=clr_ii,linewidth=2)
    ax.axvline(all_d_rs_peri_med[i],color=clr_peri,linewidth=2)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
fig.tight_layout(rect = [0.03,0.03,0.96,0.96])
fig.suptitle('Interictal and periictal $D_{RS}$ distributions')

com.save_plot(plot_dir,'interical_and_periictal_d_rs_distributions',plot_formats)

#%% scatter plot of median interictal vs periictal D_RS

# first remove interictal data if missing periictal data
# (only compare in patients with both interictal and periictal data)
nan_peri = np.isnan(all_d_rs_peri_med)
all_d_rs_ii_med[nan_peri] = np.nan

# scatter plot
fig,axs = plt.subplots(1,1,figsize=(5,5))
fig.tight_layout(rect = [0.03,0.03,0.96,0.9],h_pad=3,w_pad=3)

# reference (identity line)
sns.lineplot(x=np.array((0,1)),y=np.array((0,1)),color=(0.5,0.5,0.5),
             ax=axs)
# plot median D_RS
sns.scatterplot(x=all_d_rs_ii_med, y=all_d_rs_peri_med,
                hue=all_good_outcome,
                palette=outcome_clrs,ax=axs,legend=False,
                markers=['o'],alpha=0.5)
axs.set_title('Median interictal vs. periictal $D_{RS}$',fontsize = 16)
axs.set_ylabel('median periictal $D_{RS}$',fontsize = 14)
axs.set_xlabel('median interictal $D_{RS}$',fontsize = 14)
axs.set_ylim((0,1))
axs.set_xlim((0,1))
axs.set_aspect('equal')

com.save_plot(plot_dir,'median_interical_and_periictal_d_rs',plot_formats)

#%% compare (signed rank test)
_,p = scipy.stats.wilcoxon(x=all_d_rs_ii_med,
                             y=all_d_rs_peri_med,
                             alternative='two-sided',mode='auto')
_,p_good = scipy.stats.wilcoxon(x=all_d_rs_ii_med[all_good_outcome==1],
                             y=all_d_rs_peri_med[all_good_outcome==1],
                             alternative='two-sided',mode='auto')
_,p_bad = scipy.stats.wilcoxon(x=all_d_rs_ii_med[all_good_outcome==0],
                             y=all_d_rs_peri_med[all_good_outcome==0],
                             alternative='two-sided',mode='auto')


#%% compute AUC using median interictal and periictal D_RS to classify outcomes
# plot beeswarms comparing good/poor outcome

fig,axs = plt.subplots(1,2,figsize=(6,4))
fig.tight_layout(rect = [0.03,0.03,0.96,0.96],h_pad=3,w_pad=3)

# interictal 
ab_vis.plot_outcome_auc(all_good_outcome,all_d_rs_ii_med,
                        meas_name='median interictal $D_{RS}$',ax=axs[0],
                        overlay_violin=True,outcome_clrs=np.flipud(outcome_clrs))

# periictal 
ab_vis.plot_outcome_auc(all_good_outcome,all_d_rs_peri_med,
                        meas_name='median periictal $D_{RS}$',ax=axs[1],
                        overlay_violin=True,outcome_clrs=np.flipud(outcome_clrs))

com.save_plot(plot_dir,'interictal_periictal_d_rs_auc',plot_formats)

