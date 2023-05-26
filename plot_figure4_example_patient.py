#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates example patient plot in Figure 4 from Wang et al. (2023) Temporal 
stability of intracranial EEG abnormality maps for localizing epileptogenic 
tissue.

Calculates time-varying D_RS from abnormalities before plotting.

Choose example patient at beginning of script using variable 'pnt'

Gabrielle M. Schroeder
CNNP Lab, Newcastle University
May 2023
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

# custom functions
import pyabnorm.common as com
import pyabnorm.compute_abnorm as ab_comp
import pyabnorm.vis_abnorm as ab_vis

#%% Choose one patient - example in paper is patient 3 

pnt = 'patient 3' 

# no more variable modifications needed after this point

#%% define directories and add patient data

# replace space with underscore 
pnt = pnt.replace(' ','_')

# directories
data_dir,plot_dir = com.set_paths()

# add subfolder for figure 4 plots
plot_dir = os.path.join(plot_dir,'fig4_example')

# patient data
fname = os.path.join(data_dir,f'{pnt}.mat')

# load patient data
pnt_dict = com.mat_load_as_dict(fname)

# read out important variables
good_outcome = pnt_dict['good_outcome'] # patient outcome
pnt_id = pnt_dict['pnt_id'] # patient id
patient_string = pnt_dict['patient_string'] # label for plots
roi_ab = pnt_dict['roi_ab'] # time-varying ROI abnormalities
roi_is_resect = pnt_dict['roi_is_resect'] == 1 # which ROIs were resected
roi_names = np.char.strip(pnt_dict['roi_names']) # ROI names
t_days = pnt_dict['t_days'] # time elapsed (in days) for each time window in roi_ab
win_interictal = pnt_dict['win_interictal']
win_periictal = pnt_dict['win_periictal']
win_ictal = pnt_dict['win_ictal']
#%% set colours and other vis settings

clr_ii = (178/255,178/255,178/255) # interictal periods
clr_peri = (227/255,124/255,29/255) # periictal periods
clr_sz = (206/255,128/255,128/255) # seizures

# fixes text embedding in pdfs
matplotlib.rcParams['pdf.fonttype']=42

# format for saving plots
plot_formats = ['pdf'] 

#%% time windows containing seizures (for plot)

sz_win = np.where(win_ictal)[0]
n_sz_win = len(sz_win)

#%% compute D_RS 

d_rs = ab_comp.compute_d_rs(roi_ab,roi_is_resect) # D_RS vs time

#%% plot time-varying abnormalities and D_RS with seizures and periictal/interictal periods marked

fig,axs = plt.subplots(2,1,figsize=(8,8),gridspec_kw={'height_ratios': [3,1]})
fig.tight_layout(rect = [0.03,0.03,0.96,0.90],h_pad=3,w_pad=3)
plt.suptitle(f'{patient_string} abnormalities and $D_{{RS}}$, periictal vs. interictal')

# plot abnormalities
# first normalise so each window sums to 1
roi_ab_norm = roi_ab.copy() 
roi_ab_norm = roi_ab_norm/np.tile(np.nansum(roi_ab_norm,0),(pnt_dict['n_roi'],1))
axs[0].set_title('abnormalities (normalized)')
ax = ab_vis.heatmap_abnormalities(roi_ab_norm,t_days,roi_names,
                           roi_is_resect=roi_is_resect,ax=axs[0])
# mark windows containing seizures
for i in range(n_sz_win):
    axs[0].axvline(sz_win[i],color=clr_sz,linewidth=2)

# plot D_RS
d_rs_interictal = d_rs.copy()
d_rs_interictal[win_interictal==0] = np.nan # remove data that is not interical
d_rs_periictal = d_rs.copy()
d_rs_periictal[win_periictal==0] = np.nan # remove data that is not periictal

ab_vis.plot_D_RS(t_days,d_rs_interictal,ax=axs[1],clr=clr_ii)
ab_vis.plot_D_RS(t_days,d_rs_periictal,ax=axs[1],clr=clr_peri)
axs[1].set_yticks((0,0.25,0.5,0.75,1))

# mark windows with seizures
for i in range(n_sz_win):
   axs[1].axvline(t_days[sz_win[i]],color=clr_sz,linewidth=0.5)

# save
com.save_plot(plot_dir,f'{pnt}_interictal_vs_periictal',plot_formats)


