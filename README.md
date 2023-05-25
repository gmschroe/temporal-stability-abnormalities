# temporal-stability-abnormalities
Code for Wang et al. (2023) Temporal stability of intracranial EEG abnormality maps for localizing epileptogenic tissue

## Requirements

All code was developed in Python version 3.9.12. 

Dependencies are provided in requirements.txt. The following shell code creates a Python environment that uses the same package versions:

    # creates a Python virtual environment
    $ python -m venv env 
    
    # activate the virtual environment
    $ source env/bin/activate
    
    # install the dependencies
    $ python -m pip install -r requirements.txt
    
When you you are finished running the code, you can deactivate the environment by either closing the shell/terminal or using the shell command

    $ deactivate 
