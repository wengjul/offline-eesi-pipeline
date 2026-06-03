# -*- coding: utf-8 -*-
"""
Created on Thu Mar 20 13:51:24 2025

@author: weng_j
"""



# library importing
import os
import sys
import shutil
import h5py
import re
import numpy as np
import pandas as pd
import math
import statsmodels.api as sm
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from tkinter import filedialog as fd # module for interactive directory selection

from scipy import stats # used for linear regression


# Import modules from .py files
import eesi_import_export_functions as impexp # import and export functions
import eesi_ion_formula_parsing as parsing # EESI ion formula regex parsing



###############################################################################
## Insert path for offline-ams-pipeline here, to import functions from there ##
###############################################################################
path = ''
sys.path.append(path)  # Set absolute path of that directory
import ams_main_visualization_functions as vis # import visualization functions
import ams_main_general_functions as general # import general functions

###############################################################################
############################### Functions #####################################
###############################################################################





'''
Strategies for ion filtering that can be done row per row individually 

- m/z range filter
- elements filter (e.g. only carbon species)

'''


# elements filter
## Can be applied e.g. to filter out all non-carbonaceous species
def filter_by_elements_eesi(auto_df_eesi, elements):
    '''
    Fiters the EESI ions in the autosampler df by elements. All formulas are filtered
    out that are NOT included in elements. This always includes the unassigned formulas.
    Filtered ions are saved to a csv file. Statistics are displayed.

    Parameters
    ----------
    auto_df_eesi : Dataframe
        The EESI autosampler df 
        (has to hold datatype of normalized, wb subtracted sample intensities).
    elements : set
        Set of element that have to be part of the formulas to not filter out. 
        E.g. {'C', '[15N]'} for carbon and 15N-nitrogen.

    Returns
    -------
    None.

    '''
    # Copy input to avoid overwriting
    auto_df_filter = auto_df_eesi.copy()
    
    # Get all unique ions
    ions = pd.DataFrame({'ion': [col.replace('_norm_wb_sub', '') for col in auto_df_filter.columns if '_norm_wb_sub' == ('_' + col.split('_', 1)[1] if '_' in col else None) ]})
    # Parse element counts
    ions_elements = parsing.parse_eesi_amus(ions)

    # Create a mask that checks if any of the allowed elements are present in each row
    mask = ions[elements].sum(axis=1) > 0
    # Apply the mask and retrieve all ion names that don't contain the specified elements
    filtered_ions = ions[~mask]['ion']
    # Create dataframe for export
    filtered_ions = pd.DataFrame({'EESI Ion':filtered_ions})

    # Export filtered ion list to .csv
    impexp.export_tofile(filtered_ions, f'EESI_ion_filtering_by_elements_{elements}.csv')

    # Filter statistics
    print('\nION FILTER STATISTICS\n')
    print(f'Filtered {len(filtered_ions)} ions out of {len(ions)}')
    print(f'That are {100 * (len(filtered_ions) / len(ions)):.2f} %')










'''
Strategies for ion filtering that depend on a feature matrix of quality control samples
(QC filters, QC repeat samples, Extracted MilliQ) 
For the main sample data, filter list can be applied individually row by row

- filter based on noise/relative deviation of each ion within QC filter samples throughout the campaign 
(no trends/directions of noise considered)


'''






# Filter EESI ions by their relative signal deviation between the QC Filter samples
# Note: it would be an option to filter based on ion fractions (relative signals)
# Instead of absolute signals. However, then it's important to pre-select ions
# As e.g. primary ions, internal standard ions, negative ions (high wb) will strongly
# Impact the ion fractions if they're intense. Therefore, for now don't consider ion fractions
def filter_by_qc_sample_eesi (auto_df_qc, threshold_rel_deviation, declustering='all'):
    '''
    Filters EESI ions, based on their relative signal deviation between the repeat
    QC samples. Ions with a median relative deviation above the defined threshold 
    are filtered out and saved into a .csv file.
    Filter statistics are displayed.

    Parameters
    ----------
    auto_df_qc : dataframe
        The EESI autosampler df, filtered for only the quality control filter samples.
    threshold_rel_deviation : Float
        The relative signal deviation threshold to apply for the filtering.
        E.g. 0.5 means below +- 50% median deviation.
    declustering : string, optional
        Specification which samples should be considered with regard to declustering change. 
        'all' - all samples
        'before' - samples before declustering change
        'after' - samples after declustering change
        The default is 'all'.

    Returns
    -------
    None.

    '''
    # Copy input to avoid overwrite
    auto_df_in = auto_df_qc.copy()
    
    # Select samples based on declustering:
    # all samples 
    if declustering=='all':
        auto_df_in = auto_df_in
    # before declustering change
    elif declustering=='before':
        auto_df_in = auto_df_in[auto_df_in['sample_starting_time'].apply(pd.to_datetime) < pd.to_datetime('2023-02-05 00:00:00')].reset_index(drop=True)
    # after declustering change
    elif declustering=='after':
        auto_df_in = auto_df_in[auto_df_in['sample_starting_time'].apply(pd.to_datetime) > pd.to_datetime('2023-02-05 00:00:00')].reset_index(drop=True)
    
    # Normalize EESI ions by total intensity of each sample -> ion fractions
    #auto_df_fractions = parsing.eesi_ion_fractions(auto_df_in)

    # Get the mass spectra part
    auto_df_ms = general.get_autosampler_ms_part(auto_df_in)

    # Set all negative values in the signal intensities to nan
    auto_df_ms[auto_df_ms < 0] = np.nan
    # Set all 0 values in the signal intensities to nan
    auto_df_ms[auto_df_ms == 0] = np.nan


    # For each unique combination of QC samples, calculate the relative deviation of singal intensity
    # abs(1 - singal_QCFilter1/singal_QCFilter2)

    results = []

    for i in range(len(auto_df_ms)):
        for j in range(len(auto_df_ms)):
            if i < j:  # Ensure we don't divide a row by itself
                row_i = auto_df_ms.iloc[i]
                row_j = auto_df_ms.iloc[j]
                division = np.abs(1- row_i / row_j)
                result_row = pd.Series({'Description': f'QC Sample {i + 1} vs. {j + 1}'}).append(division)
                results.append(result_row)


    ms_reldev = pd.DataFrame(results)


    # Get median of the relative deviations
    ms_reldev_median = ms_reldev[ms_reldev.select_dtypes(include=[int, float]).columns].median().reset_index()
    ms_reldev_median.columns = ['EESI Ion', 'Relative deviation']


    # Filter EESI ions by threshold defined for relative deviation
    filtered_ions = ms_reldev_median[ms_reldev_median['Relative deviation'] > threshold_rel_deviation]

    # Export filtered ion list to .csv
    impexp.export_tofile(filtered_ions, f'EESI_ion_filtering_reldev_qc_filtersamples_{threshold_rel_deviation *100}%_declustering_{declustering}.csv')

    # Filter statistics
    print('\nION FILTER STATISTICS\n')
    print(f'Filtered {len(filtered_ions)} ions out of {len(ms_reldev_median)}')
    print(f'That are {100 * (len(filtered_ions) / len(ms_reldev_median)):.2f} %')










'''
Strategies for ion filtering that directly depend on a feature matrix of samples to be studied

- positive concentration/intensity filter: filter ions that appear in negative concentrations across X% of samples
- sample to waterblank comparison: keep only ions that have for X% of samples Y% higher intensity for sample compared to waterblanks
- sample to fieldblank comparison: keep only ions that have for X% of samples Y% higher intensity for sample compared to fieldblanks


'''










# positive concentration/intensity filter
# Note:
## -> unlike in AMS, when applied to normalized sample data before wb or fb subtraction
## there seem to be no negative ion intensities in EESI

def filter_negative_ions_eesi(auto_df_cut_eesi, threshold, data_type, sample_type=None):
    '''
    Filter ions that appear in negative intensity (or concentration) across X% of samples 
    (X = threshold).
    Filtered ions are saved to a csv file. Statistics are displayed.

    Parameters
    ----------
    auto_df_cut_eesi : Dataframe
        The EESI autosampler dataframe, after applying EESI sample mask 
        (filter double measurements and bad runs etc.).
    threshold : float
        The threshold percentage of samples with positive intensity 
        to determine if an ion is retained or filtered out.
        Orientation can be 90%.
    data_type : String
        The data type to consider for filtering.
        normalized sample, '_norm_s', makes sense 
        if wb & fb effects should not be considered (there are separate filters for those).
        
    sample_type : List of Strings, optional
        The sample types to be considered in filtering. 
        E.g. [Sample, BB Event Sample]
        The default is None, i.e. all types.

    Returns
    -------
    None.

    '''
    # Consider sample types as defined in argument
    if sample_type is not None:
        df_filtered = auto_df_cut_eesi.query('type in @sample_type')
    else:
        df_filtered = auto_df_cut_eesi

    # Get df of all ions with defined data_type
    ions = pd.DataFrame({'EESI Ion': [col.replace(data_type, '') for col in df_filtered.columns if data_type == ('_' + col.split('_', 1)[1] if '_' in col else None) ]})


    # Initialize list of percentage of positive samples for each ion
    percent_pos = []

    for ion in ions['EESI Ion']:
        ion = ion + data_type
        positive_count = (df_filtered[ion] >= 0).sum()
        total_samples = len(df_filtered)
        percent_pos.append((positive_count / total_samples) * 100)

    ions['Percentage of positive samples'] = percent_pos

    # Filter AMS ions by threshold defined for percentage of positive sampels
    filtered_ions = ions[ions['Percentage of positive samples'] < threshold]

    # Export filtered ion list to .csv
    impexp.export_tofile(filtered_ions, f'EESI_ion_filtering_negative_ions_{threshold}%_{data_type}_{sample_type}.csv')



    # Filter statistics
    print('\nION FILTER STATISTICS\n')
    print(f'Filtered {len(filtered_ions)} ions out of {len(ions)}')
    print(f'That are {100 * (len(filtered_ions) / len(ions)):.2f} %')







# sample to waterblank comparison filter


def filter_ions_wb_criteria_eesi(auto_df_cut_eesi, threshold_percentage, factor, sample_type=None):
    '''
    Filter out ions that don't fulfill the requirement: at least threshold_percentage %
    of samples are a factor higher than the according waterblanks (internal standard normalized data).
    
    Filtered ions are saved to a csv file. Statistics are displayed.

    Parameters
    ----------
    auto_df_cut_eesi : Dataframe
        The EESI autosampler dataframe, after applying EESI sample mask 
        (filter double measurements and bad runs etc.).
    threshold_percentage : float
        The percentage threshold to determine if an ion meets the criteria.
    factor : float
        The factor by which the ion concentration should be higher than the waterblank.
    sample_type : List of Strings, optional
        The sample types to be considered in filtering. 
        E.g. [Sample, BB Event Sample]
        The default is None, i.e. all types.

    Returns
    -------
    None.

    '''
    # Consider sample types as defined in argument
    if sample_type is not None:
        df_filtered = auto_df_cut_eesi.query('type in @sample_type')
    else:
        df_filtered = auto_df_cut_eesi

    # Get df of all ions with defined data_type
    ions = pd.DataFrame({'EESI Ion': [col.replace('_norm_wb_sub', '') for col in df_filtered.columns if '_norm_wb_sub' == ('_' + col.split('_', 1)[1] if '_' in col else None) ]})



    # Initialize list of factor above wb for each AMS ion
    factor_above_wb = []


    for ion in ions['EESI Ion']:
        # Considering normalized data
        # Get intensities for samples of this ion
        ion_s = df_filtered[ion + '_norm_s']
        # Get intensities for wb of this ion
        ion_wb = df_filtered[ion + '_norm_wb']
        
        # Set negative values to nan
        ion_s.loc[ion_s < 0] = np.nan
        #ion_s[ion_s == 0] = np.nan
        ion_wb.loc[ion_wb < 0] = np.nan
        #ion_wb[ion_wb == 0] = np.nan
        
        # Calculate factors s/wb
        ion_factors = list(ion_s/ion_wb)
        factor_above_wb.append(ion_factors)



    # Add the factor lists to the ion df
    ions['Factor s/wb'] = factor_above_wb

    # For each ion, calculate percentage of samples with s/wb > factor

    percentage_above_factor = []

    for ion in ions['EESI Ion']:
        # Get number of samples
        number_samples = len(list(ions[ions['EESI Ion'] == ion]['Factor s/wb'])[0])

        # Get number of samples above factor
        number_samples_above_factor = len([x for x in list(ions[ions['EESI Ion'] == ion]['Factor s/wb'])[0] if x > factor])
        # Get the percentage of samples above factor and append to list
        percentage_above_factor.append((number_samples_above_factor/number_samples) * 100)


    # Add the percentages to the ions df
    ions['Percentage of s/wb > factor'] = percentage_above_factor


    # Filter ions based on percentage feature and specified threshold
    filtered_ions = ions[ions['Percentage of s/wb > factor'] < threshold_percentage]


    # Export filtered ion list to .csv
    impexp.export_tofile(filtered_ions, f'EESI_ion_filtering_sample_to_wb_{threshold_percentage}%_{factor}_{sample_type}.csv')

    # Filter statistics
    print('\nION FILTER STATISTICS\n')
    print(f'Filtered {len(filtered_ions)} ions out of {len(ions)}')
    print(f'That are {100 * (len(filtered_ions) / len(ions)):.2f} %')



































'''
Filter sequence and filter performance metrics

'''



# Filtering and dropping ions specified in excel spreadsheet

def drop_ion_by_list(auto_df_eesi, ion_filter):
    '''
    Drops ions from the autosampler dataframe based on an imported ion list 
    (excel spreadsheet).

    Parameters
    ----------
    auto_df_eesi : Dataframe
        The autosampler dataframe with ion names as columns.
    ion_filter : String
        Identifier of spreadsheet ion filter column to select for the filtering.

    Returns
    -------
    auto_ion_drop_filtered : Dataframe
        The autosampler dataframe without the columns corresponding to the filtered ions.

    '''
    
    # Copy input df to make sure it's not overwritten
    auto_ion_drop = auto_df_eesi.copy()

    # Import ion filtering lists, select column specified in ion_filter parameter and get it as list
    ion_filter_list =  impexp.import_file('3_EESI_ion_filtering_lists')[ion_filter].tolist()

    # Remove all nan values from the list
    ion_filter_list = [x for x in ion_filter_list if not isinstance(x, float) or not math.isnan(x)]

    # Add a underscore to each ion name (to make it them unique to other ions that contain these atoms)
    ion_filter_list = [i + '_' for i in ion_filter_list]

    # Check all column names if they start with any of the substrings in ion_filter_list
    # i.e. ionname_. If so, get the column index number and appand it to list ion_filter_columns.
    # Using startswith() to exclude that other ions get selected that contain the 
    # substrings in ion_filter_list only as part of their sum formula
    ion_filter_columns = [index for index, col in enumerate(auto_ion_drop.columns) for sub in ion_filter_list if col.startswith(sub)]

    # Drop the columns with idices from ion_filter_columns list to get filtered autosampler AMS dataframe
    auto_ion_drop_filtered = auto_ion_drop.drop(auto_ion_drop.columns[ion_filter_columns], axis=1)
    
    return auto_ion_drop_filtered























