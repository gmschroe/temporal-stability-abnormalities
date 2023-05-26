'''
Function for handling paths, data, and figures.

Gabrielle M. Schroeder
CNNP Lab, Newcastle University
May 2023

'''

import scipy.io as sio
import os
import matplotlib.pyplot as plt
    
def mat_load_as_dict(filename,*args):
    '''
    Loads mat file as python dictionary (each variable is a dictionary entry).

    Parameters
    ----------
    filename : string
        Filename (including path) of mat file to load.
    *args : string(s) (optional)
        Variables to load from mat file, specified as strings.

    Returns
    -------
    mat_dict : dictionary
        Dictionary containing variables from mat file.

    '''
    
    # load all variables, formatted for python
    d=sio.loadmat(filename,squeeze_me=True,chars_as_strings=True,mat_dtype=False)
    
    if not args:
        # load all variables
        mat_dict = d
    else:
        # only keep specified variables in dictionary
        mat_dict = dict()
        for k in args:
            mat_dict[k] = d[k]
    
    return mat_dict


def save_plot(plot_dir,fname,plot_formats=None):
    '''
    Saves plot in desired format(s)

    Parameters
    ----------
    plot_dir : string
        Directory in which to save plot.
    fname : string
        Desired plot filename (excluding extension).
    plot_formats : list of strings (optional)
        List of plot formats, as strings. The default is ['png'].

    '''
    # plot_formats is list of deired plot formats
    # ensure plot_dir directory exists
    
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    
    # default format is png
    if plot_formats is None:
        plot_formats = ['png']
    
    # save plot in each requested format
    for format_i in plot_formats:
        plt.savefig(os.path.join(plot_dir,f'{fname}.{format_i}'),
                    format=format_i,bbox_inches='tight')


def set_paths(data_dir='data',plot_dir='plots'):
    '''
    Specify folders containing the data and plots. Plot folder will be created
    if it does not yet exist. Returns folder names as variables so they can be
    used in the script.

    Parameters
    ----------
    data_dir : string, optional
        Name of folder that contains patient mat files. The default is 'data'.
    plot_dir : string, optional
        Folder in which to save plots. The default is 'plots'.

    Returns
    -------
    data_dir : string
        Name of folder containing data.
    plot_dir : string
        Name of folder in which to save plots.

    '''
    
    print('Setting data and plot directories...')

    # directory containing patient data
    data_dir = 'data'

    print(f'Data directory: {data_dir}')

    # directory for plots (makes directory if doesn't exist)
    os.makedirs(plot_dir,exist_ok=True)

    print(f'Plot directory: {plot_dir}')
    
    return data_dir,plot_dir