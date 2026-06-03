# -*- coding: utf-8 -*-
"""
Created on Tue Mar 18 15:30:25 2025

@author: weng_j
"""

import math
import numpy as np
import pandas as pd



import re # regex parsing
from collections import defaultdict


# Import modules from .py files
import eesi_import_export_functions as impexp # import and export functions

import sys

###############################################################################
## Insert path for offline-ams-pipeline here, to import functions from there ##
###############################################################################
path = ''
sys.path.append(path)  # Set absolute path of that directory

import eesi_visualization_functions as eesi_vis
import eesi_ion_formula_parsing as parsing # EESI ion formula regex parsing


###############################################################################
############################### Functions #####################################
###############################################################################


# Function to parse chemical formulas, including isotopes in square brackets
def regex_parse_formula(formula):
    '''
    Regex (regular expressions) parsing of a chemical formula from the EESI
    amus dataframe. Output is dictionary with counts of each identified chemical 
    element.

    Parameters
    ----------
    formula : String
        EESI amus ion formula.

    Returns
    -------
    element_counts : defaultdict
        Dictionary with counts of each identified chemical 
        element.

    '''
    element_counts = defaultdict(int)
    # Match isotopes in brackets ([15N]) and normal elements (C, H, O, etc.)
    matches = re.findall(r'(\[?\d*[A-Z][a-z]*\]?\d*)', formula)
    
    for match in matches:
        # Extract element and count (handles cases like "[15N]2" and "C6")
        element_match = re.match(r'(\[\d*[A-Z][a-z]*\]|\d*[A-Z][a-z]*)(\d*)', match)
        if element_match:
            element, count = element_match.groups()
            element_counts[element] += int(count) if count else 1
            
    return element_counts







def parse_eesi_amus (amus):
    '''
    On EESI ion formula df (amus), applies regex_parse_formula
    to parse elements & counts for each formula. Adds the counts to element columns
    in amus.
    
    Parameters
    ----------
    amus : Dataframe
        The EESI ion formula df.

    Returns
    -------
    amus : Dataframe
        The EESI ion formula df with element counts.

    '''

    
    # For each amus formula, parse elements included and counts for each element 
    # Extract all element counts
    all_elements = set()
    parsed_data = amus['ion'].apply(lambda x: regex_parse_formula(x) if isinstance(x, str) else {})

    # Collect unique elements
    for elem_counts in parsed_data:
        all_elements.update(elem_counts.keys())

    # Add element columns to amus dataframe
    for element in sorted(all_elements):
        amus[element] = parsed_data.apply(lambda x: x.get(element, 0))
    
    return amus






def filter_amus_ions(amus, allowed_elements):
    '''
    Filters the EESI ion formulas by defined elements, i.e. only keep formulas
    that contain those.

    Parameters
    ----------
    amus : DataFrame
        EESI ion formula amus dataframe, after parsing element counts.
    allowed_elements : set
        Set of allowed element symbols. E.g. {'C', '[15N]'} for carbon and 15N-nitrogen.

    Returns
    -------
    amus : DataFrame
        Filtered DataFrame.

    '''
    # Create a mask that checks if any of the allowed elements are present in each row
    mask = amus[allowed_elements].sum(axis=1) > 0

    # Apply the mask to filter rows
    return amus[mask]



# Filter EESI autosampler by amus ions

def filter_auto_df_eesi_by_amus (auto_df_eesi, amus_filter):
    '''
    Filter ion columns of EESI autosampler df by the amus list.

    Parameters
    ----------
    auto_df_eesi : Dataframe
        The EESI autosampler df.
    amus_filter : Dataframe
        The amus dataframe to be applied for filtering.

    Returns
    -------
    auto_df_filter : Dataframe
        The EESI autosampler df after filtering

    '''
    # Copy input to avouid overwriting
    auto_df_filter = auto_df_eesi.copy()

    # Separate sample info and ms part of the autosampler df
    auto_df_filter_sample_info = general.get_autosampler_sample_info_part(auto_df_filter)
    auto_df_filter_ms = general.get_autosampler_ms_part(auto_df_filter)

    # Get data type in autosampler df
    datatype = '_' + auto_df_filter_ms.columns[0].split("_", 1)[1]

    # Filter autosampler ms part by amus
    auto_df_filter_ms = auto_df_filter_ms.loc[:, auto_df_filter_ms.columns.isin(amus_filter['ion'] + datatype)]

    # Re-combine autosampler parts
    auto_df_filter = pd.concat([auto_df_filter_sample_info, auto_df_filter_ms], axis=1)

    return auto_df_filter






# Retrieve total ion intensity for EESI Ions and save as column

def eesi_total_ion_intensity (auto_df_eesi):
    '''
    Calculate the total ion intensity for each sample (row) in the EESI autosampler
    df and add it as a column in the sample info part.

    Parameters
    ----------
    auto_df_eesi : Dataframe
        The EESI autosampler df, only one datatype.

    Returns
    -------
    df : Dataframe
        The EESI autosampler df with total ion intensity column.

    '''
    # Copy input to avouid overwriting
    df = auto_df_eesi.copy()
    
    # Separate sample info and ms part of the autosampler df
    auto_df_sample_info = general.get_autosampler_sample_info_part(df)
    auto_df_ms = general.get_autosampler_ms_part(df)

    # Get total ion intensities for each sample
    intensities = auto_df_ms.sum(axis=1)

    # Add as a column to sample info part of autosampler df
    auto_df_sample_info.insert(auto_df_sample_info.columns.get_loc('sample_rn'), 'Total ion intensity', intensities)
    
    # Re-combine autosampler parts
    df = pd.concat([auto_df_sample_info, auto_df_ms], axis=1)
    
    return df
    





# Normalize EESI ions by total intensity of each sample -> ion fractions

def eesi_ion_fractions (auto_df_eesi):
    '''
    Calculate EESI ion fraction for each sample in the eesi autosampler df.
    I.e. normalize each ion by total ion intensity of that sample.

    Parameters
    ----------
    auto_df_eesi : Dataframe
        The EESI autosampler df, only one datatype.

    Returns
    -------
    auto_df_norm : Dataframe
        The EESI autosampler df with ion fractions, and a column for total ion intensity.

    '''
    # Copy input to avouid overwriting
    auto_df_norm = auto_df_eesi.copy()
    
    # Calculate total ion intensitie for each sample using eesi_total_ion_intensity()
    auto_df_norm = eesi_total_ion_intensity (auto_df_norm)
    
    # Separate sample info and ms part of the autosampler df
    auto_df_sample_info = general.get_autosampler_sample_info_part(auto_df_norm)
    auto_df_ms = general.get_autosampler_ms_part(auto_df_norm)

    # Normalize ion intensities by sum of ion intensities of each sample
    auto_df_ms = auto_df_ms.div(auto_df_sample_info['Total ion intensity'], axis=0)

    # Re-combine autosampler parts
    auto_df_norm = pd.concat([auto_df_sample_info, auto_df_ms], axis=1)
    
    return auto_df_norm






def transpose_and_element_counts(df_atomic_ratios):
    '''
    Transposes the mass spectrum part of the input df and parses element counts.

    Parameters
    ----------
    df_atomic_ratios : Dataframe
        The EESI autosampler style dataframe to be processed.

    Returns
    -------
    df_in_trans_elements : Dataframe
        The trasposed EESI ion df including element counts and ion names.

    '''
    # Copy input to avoid overwriting
    df_in = df_atomic_ratios.copy()

    # Get mass spectrum part
    df_in_ms = general.get_autosampler_ms_part(df_in)
    # Set row indices to sample names
    df_in_ms.index = df_in['sample']

    # Transpose and keep ion names as a column
    df_in_trans = df_in_ms.T.reset_index()

    # Rename ion name column
    df_in_trans = df_in_trans.rename(columns={'index': 'ion'})
    # Remove datatype denotations from ion names
    df_in_trans['ion'] = df_in_trans['ion'].str.split('_', 1).str[0]

    # Extract element counts from ion names
    df_in_trans_elements = parsing.parse_eesi_amus(df_in_trans)
    
    return df_in_trans_elements




# Calculate EESI atomic ratios
def eesi_atomic_ratios(auto_df_eesi):
    '''
    Calculated atomic ratios for each sample in a EESI autosampler df.
    Following the formula:
    n_element_i = sum(Intensity_j * number_of_element_i_j) for each EESI ion j and element i
    
    atomic ratio i1_to_i2 = n_element_i1 / n_element_i2

    Parameters
    ----------
    auto_df_eesi : DataFrame
        The EESI autosampler df, containing organic ions and after filtering bad ions 
        (especially negative intensity ions after wb or fb subtraction).

    Returns
    -------
    df_atomic_ratios : DataFrame
        The EESI autosampler df with atmoic ratios in its sample info part.

    '''
    # Copy input to avoid overwriting
    df_atomic_ratios = auto_df_eesi.copy()
    
    # Parse element counts and intensities for each EESI ion
    df_trans_elements = transpose_and_element_counts(df_atomic_ratios)

    # Get list of all parsed elements (integer dtype columns)
    int_cols = list(df_trans_elements.select_dtypes(include=['int']).columns)

    # Calculate relative molar intensity for each sample and each element
    # N_element = sum(Intensity_j * number_of_element_j) for each EESI ion j
    df_n_element = pd.DataFrame()
    for i in int_cols:
        # Foe each EESI ion j, calculate Intensity_j * number_of_element_j for each element i
        N_element_i = df_trans_elements.select_dtypes(include=['float']).mul(df_trans_elements[i], axis=0)
        N_element_i = N_element_i.add_suffix(f'_{i}')
        # Add to parsed elements dataframe
        df_n_element = pd.concat([df_n_element,N_element_i], axis=1)
        
    # sum over ions j
    df_n_element = df_n_element.sum()


    # Define the elements to calculate ratios for
    elements = ['O', 'H', 'N', 'P', 'S']

    # Store ratios in a dictionary
    ratios = {}

    # Carbon values for denominator
    c_values = df_n_element[df_n_element.index.str.endswith('_C')].values

    for element in elements:
        key = f'{element.lower()}_to_c'
        element_values = df_n_element[df_n_element.index.str.endswith(f'_{element}')].values
        
        if element_values.size > 0:
            ratios[key] = element_values / c_values
        else:
            ratios[key] = np.full(c_values.size, np.nan) # Return an array of nan if there are no ions with that element at all

    # Insert to the original df
    df_atomic_ratios.insert(df_atomic_ratios.columns.get_loc('sample_rn'), 'O_to_C', ratios['o_to_c'])
    df_atomic_ratios.insert(df_atomic_ratios.columns.get_loc('sample_rn'), 'H_to_C', ratios['h_to_c'])
    df_atomic_ratios.insert(df_atomic_ratios.columns.get_loc('sample_rn'), 'N_to_C', ratios['n_to_c'])
    df_atomic_ratios.insert(df_atomic_ratios.columns.get_loc('sample_rn'), 'P_to_C', ratios['p_to_c'])
    df_atomic_ratios.insert(df_atomic_ratios.columns.get_loc('sample_rn'), 'S_to_C', ratios['s_to_c'])
    
    return df_atomic_ratios






###############################################################################
############################# Execution Panel #################################
###############################################################################
















