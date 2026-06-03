# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 16:24:02 2026

@author: weng_j
"""

# library importing

import numpy as np
import pandas as pd



# Returns only MS part of autosampler df

def get_autosampler_ms_part(auto_df):
    '''
    For autosampler dataframe, cuts of sample information part and returns only MS part,
    i.e. ions as columns.

    Parameters
    ----------
    auto_df : DataFrame
        Autosampler df with sample information and MS part.

    Returns
    -------
    auto_df_ms : DataFrame
        The MS part of the input autosampler df.

    '''
    # Copy data to make sure it's not overwritten
    auto_df_input = auto_df.copy()
    
    # Get column index of last column before MS data
    index = auto_df_input.columns.get_loc('back_rn')
    # Get MS data part of auto_df 
    auto_df_ms = auto_df_input.iloc[:,(index+1):]
    
    return auto_df_ms



# Returns only Sample information part of autosampler df

def get_autosampler_sample_info_part(auto_df):
    '''
    For autosampler dataframe, cuts of MS part returns only sample information part.

    Parameters
    ----------
    auto_df : DataFrame
        Autosampler df with sample information and MS part.

    Returns
    -------
    auto_df_ms : DataFrame
        The sample info part of the input autosampler df.

    '''
    # Copy data to make sure it's not overwritten
    auto_df_input = auto_df.copy()
    
    # Get column index of last column before MS data
    index = auto_df_input.columns.get_loc('back_rn')
    # Get MS data part of auto_df 
    auto_df_ms = auto_df_input.iloc[:,:(index+1)]
    
    return auto_df_ms






# Sorting intermediate dataframe -> Sample type and data type
def ams_intermediate_sorter (auto_df, sample_type='None', data_type='None'):
    '''
    Sorts the intermediate dataframe; cuts it to a selected filter 'type'
    (e.g. Sample, Fieldblank, BB Sample, BB Fieldblank ...)
    and extracts only specified data (default: normalized-WB subtracted).

    Parameters
    ----------
    auto_df : DataFrame
        The EESI dataframe after averaging.
    sample_type : List of Strings; Default: 'None'
        Specification of the 'type' to select samples; if 'None', none are cut.
        Has to be a list of strings in the form [a,b,...], even if only one keyword.
        e.g. ['Sample'] or ['Sample', 'Fieldblank']
    data_type : String ; Default: 'None'
        Specification of which datatype to select; default is 'None', i.e. no
        sorting of datatype.
        e.g. '_norm_wb' for normalized waterblank.

    Returns
    -------
    auto_df_u_sort : DataFrame
        Filtered and sorted dataframe.
    '''
    # Copy data to make sure it's not overwritten
    auto_df_uncut = auto_df.copy()
    
    ## 1. Select only columns that are specified by data_type
    if data_type == 'None':
        auto_df_u_sort = auto_df_uncut # No selection, all columns are passed
    else:
        # create list of organic ion names
        org_spec = [i.replace(data_type, '') for i in auto_df_uncut.columns if i.endswith(data_type)]
        # Get sample info part of autosampler df
        auto_df_u_sort = auto_df_uncut.iloc[:,:(auto_df_uncut.columns.get_loc('back_rn') +1)]
        # Get columns specified by data_type
        auto_df_u_data_type = auto_df_uncut[[str(sp) + data_type for sp in org_spec]]
        # Concat
        auto_df_u_sort = pd.concat([auto_df_u_sort, auto_df_u_data_type], axis = 1)
    
    ## 2. Cutting to selected sample 'type'
    if sample_type == 'None':
        auto_df_u_sort = auto_df_u_sort # No sorting according to sample type
    else:
        # Select samples by 'type'
        auto_df_u_sort = auto_df_u_sort[auto_df_u_sort["type"].isin(sample_type)].copy()
    
    # Reset index of cut df
    auto_df_u_sort = auto_df_u_sort.reset_index(drop=True)
    
    return auto_df_u_sort