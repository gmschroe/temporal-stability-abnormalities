# temporal-stability-abnormalities
Code for Wang et al. (2023) Temporal stability of intracranial EEG abnormality maps for localizing epileptogenic tissue

## Requirements

All code was developed in Python version 3.9.12. 

Dependencies are provided in `requirements.txt`. The following shell code creates a Python environment that uses the same package versions:

    # creates a Python virtual environment
    python -m venv env 
    
    # activate the virtual environment
    source env/bin/activate
    
    # install the dependencies
    python -m pip install -r requirements.txt
    
When you you are finished running the code, you can deactivate the environment by either closing the shell/terminal or using the shell command

    deactivate 

## Data

Patient data is provided at the level of region of interest (ROI) abnormalities vs. time (see paper for methods details). Each patient's data is stored in a MATLAB workspace in the `data` folder. The scripts load each workspace as a dictionary with the following key/value pairs:

* `good_outcome`: (int) whether the patient had a good surgical outcome, defined as ILAE 1 (value = `1`) or poor surgical outcome, defined as ILAE 2+ (value = `0`)
* `n_roi`: (int) the number of ROIs analysed 
* `n_win`: (int) the number of 30s time windows in the patient's iEEG recording
* `n_win_per_hr`: (int) the number of time windows in one hour of data (value = `120`, same for all patients)
* `patient_string`: (string) string containing the patient's ID and outcome info; used for plots
* `pnt_id`: (string) patient ID; matches IDs used in the paper
* `roi_ab`: (2D numpy array, float) time-varying ROI abnormalities, size `n_roi` x `n_win`
* `roi_is_resect`: (1D numpy array, int) vector, length `n_roi`, indicating whether each ROI was resected ( = `1`) or spared ( = `0`)
* `roi_names`: (1D numpy array, string) vector, length `n_roi`, indicating name of each ROI
* `roi_xyz`: (2D numpy array, float) xyz coordinates of each ROI, size `n_roi` x 3
* `t_days`: (1D numpy array, float) vector, length `n_win`, indicating amount of time in days of each time window from the start of the iEEG recording
* `win_ictal`, `win_interictal`, `win_periictal`: (1D numpy array, int) each vector is length `n_win` and indicates which time windows contain seizures, periictal data, and interictal data, respectively (non-overlapping definitions - i.e., each window can only be in one category)

Note that booleans are saved as integers in the mat files.

## Plots

All plots created by the code are already saved in `plots_reference`. When run, the scripts create a new plots folder, `plots`, so that the provided plots are not overwritten.

## Scripts

Run the `plot_figure#.py` scripts to reproduce the plots and analyses for the corresponding figure in the paper. These scripts can be run in any order.

An example patient to plot must be specified using the variable `pnt` at the beginning of `plot_figure2.py` and `plot_figure4_example_patient.py` - change the patient ID to plot a different example patient. Patient IDs are the same as the filenames in `data` (e.g., use either `'patient 1'` or `'patient_1'` to specify patient 1).

Custom functions are provided in the `pyabnorm` package.
