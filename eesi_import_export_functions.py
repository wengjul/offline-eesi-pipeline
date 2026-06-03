# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 11:42:30 2023

@author: weng_j
"""

import os
import pandas as pd
from tkinter import filedialog as fd # module for interactive file selection
from datetime import datetime # library used for date format converting




def import_file (file_specification = 'Unspecified', header_row = 'infer', spaces = False):
    '''
    Import a datafile via interactive selection by user; 
    a window is promted for this.

    Parameters
    ----------
    file_specification : String, optional
        Specifies which file should be imported (e.g. Autosampler etc). 
        The default is 'Unspecified'.
    header_row : String, optional
        Specifies which row the header starts from. 
        The default is 'infer', i.e. the read_csv function will infer it 
        automatically.
    spaces : Bool, optional
        Specifies whether to take spaces after delimiter into account or not. 
        The default is False.

    Returns
    -------
    df : DataFrame
        The imported data in dataframe format.

    '''
    # Interactive selection of file
    filename = fd.askopenfilename(title = f'Please select the {file_specification} file.')
    print(f'Selected -{file_specification}- file: {filename}')
    file_extension = os.path.splitext(filename)[1]
    if file_extension == '.xlsx': 
        df = pd.read_excel(filename)
    elif file_extension == '.csv':
        df = pd.read_csv(filename, header = header_row, skipinitialspace = spaces)
    elif file_extension == '.txt':
        df = pd.read_csv(filename, sep=' ', header = header_row, skipinitialspace = spaces)
    elif file_extension == '.h5':
        df = pd.read_hdf(filename)
    return df




def export_tofile (df, filename):
    '''
    Exports a dataframe to a csv file.

    Parameters
    ----------
    df : DataFrame
        Dataframe to export.
    filename : String
        Define filename of exported file, including file extension .csv.

    Returns
    -------
    None.

    '''
    # Interactive selection of path
    path = fd.askdirectory(title = 'Please select the path for exporting the .csv file.')
    print(f'Selected path: {path}')
    df.to_csv(path + '/' + filename ,index= False)






###############################################################################
################################ Test Area ####################################
###############################################################################









