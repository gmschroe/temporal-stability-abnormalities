#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates Figure 3 from Wang et al. (2023) Temporal stability of intracranial EEG
abnormality maps for localizing epileptogenic tissue.

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


# custom functions
import pyabnorm.common as com
import pyabnorm.compute_abnorm as ab_comp
import pyabnorm.vis_abnorm as ab_vis

#%% directories and list of patients

# directories
data_dir,plot_dir = com.set_paths()

# add subfolder for figure 3 plots
plot_dir = os.path.join(plot_dir,'fig3')

# list of patient files
all_fname = os.listdir(data_dir)
n_patients = len(all_fname)

# sort patients by number
pnt_num = [int(''.join(filter(str.isdigit,pnt))) for pnt in all_fname]
pnt_idx = np.argsort(pnt_num)
all_fname = [all_fname[idx] for idx in pnt_idx]

del pnt_num, pnt_idx
#%% compute D_RS for each patient and save outcomes

# initialise arrays/lists for saving across-patient info
all_pnt_id = []
all_d_rs = [] 
all_good_outcome = np.zeros(n_patients)
all_d_rs_med = np.zeros(n_patients)
all_n_day = np.zeros(n_patients)

for i in range(n_patients):

    # load abnormalities, resection labels, outcomes, and patient ID
    pnt_dict = com.mat_load_as_dict(os.path.join(data_dir,all_fname[i]),
                                    'pnt_id','roi_ab','roi_is_resect','good_outcome',
                                    'n_win_per_hr')
    print(f'Computing DRS of {pnt_dict["pnt_id"]}')
    # resection info as boolean
    pnt_dict['roi_is_resect'] = pnt_dict['roi_is_resect']==1
    
    # compute D_RS and median D_RS
    d_rs = ab_comp.compute_d_rs(pnt_dict['roi_ab'],pnt_dict['roi_is_resect']) # D_RS vs time
    all_d_rs_med[i] = np.nanmedian(d_rs)
    
    
    # compute amount of non-missing data (in days)
    n_day = (sum(np.invert(np.isnan(d_rs))))/(pnt_dict['n_win_per_hr']*24)
    all_n_day[i] = n_day

    # save info for plots and outcome calculations
    all_pnt_id.append(pnt_dict['pnt_id'])
    all_good_outcome[i] = pnt_dict['good_outcome']
    all_d_rs.append(d_rs)
    
    del pnt_dict, d_rs, n_day

# as boolean
all_good_outcome = all_good_outcome == 1

#%% sort patients by outcome
sort_idx = np.argsort(np.invert(all_good_outcome),kind='stable')

all_pnt_id = [all_pnt_id[i] for i in sort_idx]
all_d_rs = [all_d_rs[i] for i in sort_idx]
all_good_outcome = all_good_outcome[sort_idx]
all_d_rs_med = all_d_rs_med[sort_idx]
all_n_day = all_n_day[sort_idx]

#%% set colours and other vis settings

outcome_clrs = np.zeros((2,3))
outcome_clrs[1,:] = np.array([0/255,118/255,192/255]) # blue, good outcome 
outcome_clrs[0,:] = np.array([163/255,2/255,52/255]) # red, poor outcome

# fixes text embedding in pdfs
matplotlib.rcParams['pdf.fonttype']=42

# format for saving plots
plot_formats = ['pdf'] 

#%% plot D_RS distributions with median D_RS marked with vertical line

# number of columns and rows in vis
n_col = 8
n_row = math.ceil(n_patients/n_col)

# start figure
fig = plt.figure(figsize=(12,10))
alph = 0.4 # opacity of distribution

for i in range(n_patients):
    ax = plt.subplot(n_row,n_col,i+1)
    ax.hist(all_d_rs[i],bins=np.arange(0,1.1,0.1),
            facecolor=outcome_clrs[int(all_good_outcome[i]),:],alpha=alph)
    ax.set_title(f'{all_pnt_id[i]}'
                 f'\n {round(all_n_day[i],1)} days',fontsize=10)
    ax.set_xlabel('$D_{RS}$')
    ax.set_xlim(0,1)
    ax.set_xticks((0,0.5,1))
    ax.axvline(all_d_rs_med[i],color=outcome_clrs[int(all_good_outcome[i]),:],linewidth=2)
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
fig.tight_layout(rect = [0.03,0.03,0.96,0.96])
fig.suptitle('$D_{RS}$ distributions')

com.save_plot(plot_dir,'d_rs_distributions',plot_formats)

#%% compute percentage of D_RS <= 0.5 ("localising percentage") and visualise

perc_localising = np.zeros(n_patients)

# percentage of non-nan time windows with D_RS <= 0.5
for i in range(n_patients):
    perc_localising[i] = (
        np.sum(all_d_rs[i]<= 0.5)/
        np.sum(np.invert(np.isnan(all_d_rs[i])))
        )*100
    
# plot 
fig,axs = plt.subplots(2,1,figsize = (8,8),sharey = True)
my_bins = np.arange(0,110,10)
alph=0.5
fig.tight_layout(rect = [0.03,0.03,0.96,0.96],h_pad=8,w_pad=3)

# ILAE 1 patients
axs[0].hist(perc_localising[all_good_outcome],bins=my_bins,
                  color=outcome_clrs[1,:],alpha=alph,
                  edgecolor=[0,0,0],linewidth=4)
axs[0].set_title('ILAE = 1',fontsize=16)
axs[0].set_xticks(my_bins)
axs[0].set_ylabel('number of patients',fontsize=14)
axs[0].set_xlabel('percentage $D_{{RS}}$ $\leq$ 0.5',fontsize=14)

# ILAE 2+ patients
axs[1].hist(perc_localising[np.invert(all_good_outcome)],bins=my_bins,
                  color=outcome_clrs[0,:],alpha=alph,
                  edgecolor=[0,0,0],linewidth=4)
axs[1].set_title('ILAE 2+',fontsize=16)
axs[1].set_xticks(my_bins)
axs[1].set_ylabel('number of patients',fontsize=14)
axs[1].set_xlabel('percentage $D_{{RS}}$ $\leq$ 0.5',fontsize=14)

com.save_plot(plot_dir,'localising_percentage',plot_formats)

#%% compute AUC using median D_RS to classify outcomes; plot beeswarms comparing good/poor outcome

fig,axs = plt.subplots(1,1,figsize=(3,4))
fig.tight_layout(rect = [0.03,0.03,0.96,0.96],h_pad=3,w_pad=3)

ab_vis.plot_outcome_auc(all_good_outcome,all_d_rs_med,
                        meas_name='median $D_{RS}$',ax=axs,
                        overlay_violin=True,outcome_clrs=np.flipud(outcome_clrs))

com.save_plot(plot_dir,'d_rs_auc',plot_formats)

#%% ROC curve for median D_RS

# bad outcome = positive label
fpr,tpr,_ = sk_metrics.roc_curve(all_good_outcome == 0, all_d_rs_med)

fig = plt.figure()
plt.plot((0,1),(0,1),color=(0.5,0.5,0.5))
plt.plot(fpr,tpr,color='black',linewidth=2)

ax=plt.gca()
ax.set_aspect('equal')
ax.set_ylabel('true positive rate')
ax.set_xlabel('false positive rate')
ax.set_ylim((0,1))
ax.set_xlim((0,1))
ax.set_title('ROC curve of median $D_{RS}$')

com.save_plot(plot_dir,'d_rs_roc',plot_formats)

#%% wilcoxon one sample signed rank test to test if median D_RS of good/bad outcome 
# patients is </> 0.5, respectively 

good_diff = 0.5 - all_d_rs_med[all_good_outcome==1]
_,p_good = scipy.stats.wilcoxon(good_diff,alternative='greater') # greater = 0.5>patient D_RS values
print(f'p for if good outcome median D_RS < 0.5 (Wilcoxon signed rank): {round(p_good,5)}')

bad_diff = 0.5 - all_d_rs_med[all_good_outcome==0]
_,p_bad = scipy.stats.wilcoxon(bad_diff,alternative='less') # less = 0.5<patient D_RS values
print(f'p for if bad outcome median D_RS > 0.5 (Wilcoxon signed rank): {round(p_bad,5)}')

del good_diff, bad_diff, p_good, p_bad
