#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates Figure 2 from Wang et al. (2023) Temporal stability of intracranial EEG
abnormality maps for localizing epileptogenic tissue.

Change the variable "pnt" at the beginning of the script to plot a different 
example patient.

Calculates time-varying D_RS from abnormalities before plotting.

Brain surface plots not included (created using MATLAB code - see Zenodo file 
https://zenodo.org/record/6090368#.ZC7iGibTWfZ for how to plot this type of 
figure).

Distribution of D_RS values used in Fig. 2 is created in Fig. 3 script.

Gabrielle M. Schroeder
CNNP Lab, Newcastle University
May 2023
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import os

# custom functions
import pyabnorm.common as com
import pyabnorm.compute_abnorm as ab_comp
import pyabnorm.vis_abnorm as ab_vis

#%% Choose one patient - examples in paper are patient 1 and patient 2

pnt = 'patient 1' # good surgical outcome patient
# pnt = 'patient 2' # poor surgical outcome patient 

# no more variable modifications needed after this point

#%% define directories and add patient data

# replace space with underscore 
pnt = pnt.replace(' ','_')

# directories
data_dir,plot_dir = com.set_paths()

# add subfolder for figure 2 plots
plot_dir = os.path.join(plot_dir,'fig2')

# patient data
fname = os.path.join(data_dir,f'{pnt}.mat')

# load patient data
pnt_dict = com.mat_load_as_dict(fname)

# read out important variables
good_outcome = pnt_dict['good_outcome'] # patient outcome
pnt_id = pnt_dict['pnt_id'] # patient id
patient_string = pnt_dict['patient_string'] # label for plos
roi_ab = pnt_dict['roi_ab'] # time-varying ROI abnormalities
roi_is_resect = pnt_dict['roi_is_resect'] == 1 # which ROIs were resected
roi_names = np.char.strip(pnt_dict['roi_names']) # ROI names
t_days = pnt_dict['t_days'] # time elapsed (in days) for each time window in roi_ab

#%% set colours and other vis settings

# blue for good outcome
if good_outcome == 1: 
    clr_d_rs = (0/255,118/255,192/255)
    
# red for poor outcome
else: 
    clr_d_rs = (163/255,2/255,52/255)

# fixes text embedding in pdfs
matplotlib.rcParams['pdf.fonttype']=42

# format for saving plots
plot_formats = ['pdf'] 

#%% compute D_RS and median D_RS

d_rs = ab_comp.compute_d_rs(roi_ab,roi_is_resect) # D_RS vs time
d_rs_med = np.nanmedian(d_rs) # median D_RS

#%% select example time window with median D_RS to visualise

# indices of time windows with median D_RS
idx = np.where(d_rs==d_rs_med)
idx=idx[0]    
n_idx = len(idx) # number of time windows

# select one time window (can change scale to see more examples)
if pnt_id == 'patient 1':
    idx_scale = 0.8 #  0.8 --> index selected is later in recording than 80% of the other indices with median D_RS
elif pnt_id == 'patient 2':
    idx_scale = 0.4
else:
    idx_scale = 0.3
idx = idx[int(np.ceil(n_idx*idx_scale))]
    
del idx_scale

#%% plot time-varying abnormalities and D_RS

fig,axs = plt.subplots(2,1,figsize=(8,6),gridspec_kw={'height_ratios': [3,1]})
fig.tight_layout(rect = [0.03,0.03,0.96,0.90],h_pad=3,w_pad=3)
plt.suptitle(f'{patient_string} abnormalities and $D_{{RS}}$')

# plot abnormalities
# first normalise so each window sums to 1
roi_ab_norm = roi_ab.copy() 
roi_ab_norm = roi_ab_norm/np.tile(np.nansum(roi_ab_norm,0),(pnt_dict['n_roi'],1))
axs[0].set_title('abnormalities (normalized)')
ax = ab_vis.heatmap_abnormalities(roi_ab_norm,t_days,roi_names,
                           roi_is_resect=roi_is_resect,ax=axs[0])

# plot D_RS
ab_vis.plot_D_RS(t_days,d_rs,ax=axs[1],clr=clr_d_rs)

#  mark example window (will be plotted in next figure)
axs[1].scatter(t_days[idx],d_rs[idx],zorder=2,color=clr_d_rs)

# save
com.save_plot(plot_dir,f'{pnt}_abnormalities_and_d_rs',plot_formats)

#%% plot example time window

fig,axs = plt.subplots(1,1,figsize=(3,4))
fig.tight_layout(rect = [0.03,0.15,0.96,0.85])
plt.suptitle(f'{patient_string} example time window with median $D_{{RS}}$')

# time window abnormalities
ab = roi_ab_norm[:,idx]

# colours
cmap_scale = 1.125
cmap = sns.color_palette('viridis',as_cmap=True)
norm = plt.Normalize(vmin=0,vmax=np.max(ab)*cmap_scale)
palette = {h: cmap(norm(h)) for h in ab}

# plot
axs.set_ylim((0,np.max(ab)*cmap_scale))
sns.swarmplot(x=roi_is_resect,y=ab,hue=ab,
                 size=10,palette=palette,ax=axs)
sns.violinplot(x=roi_is_resect,y=ab,
                 color=(0.95,0.95,0.95),ax=axs,saturation=1,
                 inner='quartiles',cut=0.5,width=0.8) 
axs.get_legend().remove()
plt.colorbar(matplotlib.cm.ScalarMappable(cmap=cmap,norm=norm),ax=axs)
axs.set_xticklabels(['spared','resected'])
axs.set_title(f'$D_{{RS}} = ${round(d_rs_med,2)}')
axs.set_ylabel('normalized abnormality')

# save
com.save_plot(plot_dir,f'{pnt}_abnormalities_example',plot_formats)
