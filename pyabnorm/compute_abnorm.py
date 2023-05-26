'''
Functions for computing abnormalities and abnormality-based measures from 
continuous iEEG data.

Gabrielle M. Schroeder
CNNP Lab, Newcastle University
May 2023
'''

import numpy as np
import sklearn.metrics as sk_metrics

def compute_d_rs(roi_data,roi_is_resect):
    '''
    Compute D_RS (AUC) of resected vs. spared ROIs based on ROI measure
    
    D_RS is <0.5 when measure is higher for resected ROIs 

    Parameters
    ----------
    roi_data : 2D numpy array, float
        ROI data, size # ROIs x # time windows
    roi_is_resect : 1D numpy array, boolean
        Whether each ROI was resected (True if resected, False if spared).

    Returns
    -------
    d_rs : 1D numpy array, float
        AUC at each time window, size # time windows 

    '''
    # data dimensions
    n_roi,n_win = roi_data.shape
    
    # initalise d_rs array
    d_rs = np.zeros(n_win)
    
    # check that have both classes
    if np.sum(roi_is_resect)>0 and np.sum(np.invert(roi_is_resect))>0:
        
        # compute d_rs at each time window
        for i in range(n_win):
            
            # meaure from specified time window
            data_i = roi_data[:,i].copy() 
            
            # keep non-nan ROI
            keep_roi = np.invert(np.isnan(data_i))
            data_i = data_i[keep_roi]
            spared_i = ~roi_is_resect[keep_roi] # also invert so boolean indicates spared
            
            # check that still have both classes (will not if window is completely nans)
            if (np.sum(spared_i)>0) and (np.sum(np.invert(spared_i))>0):
                d_rs[i] = sk_metrics.roc_auc_score(spared_i, data_i) # compute D_RS
            else:
                d_rs[i] = np.nan    
    else:
        d_rs = d_rs*np.nan
    
    return d_rs
