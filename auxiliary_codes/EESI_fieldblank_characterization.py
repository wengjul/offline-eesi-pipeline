# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 11:24:25 2025

@author: weng_j
"""



import math
import numpy as np
import pandas as pd

import ast


import plotly.express as px # interactive java script based visualization
import plotly.io as pio # setting of default renderer of plotly figures (has to be set to 'browsers' to get interactive plots)
import plotly.graph_objects as go

# Import modules from .py files
import eesi_import_export_functions as impexp # import and export functions
import eesi_ion_filtering as eesi_filter # application of EESI ion filters
import eesi_ion_formula_parsing as parsing # EESI ion formula regex parsing



import sys

###############################################################################
## Insert path for offline-ams-pipeline here, to import functions from there ##
###############################################################################
path = ''

sys.path.append(path)  # Set absolute path of that directory
import ams_main_visualization_functions as vis # import visualization functions
import ams_main_general_functions as general # import general functions
import ams_ion_filtering as ion_filter # filtering AMS ions




'''
Characterize EESI fieldblanks to come up with a strategy for fieldblank subtraction
'''

###############################################################################
############################### Functions #####################################
###############################################################################



# Copy from EESI_quality_control.py

# Scatter matrix of QC filter samples

def qc_filter_scatter_matrix_eesi(auto_df_eesi, axis_range, data_type, axis_scale='linear', diagonal_line=False):
    '''
    Comparison of EESI ions of repeate QC filters or QC repeat samples by 
    scatter matrix.

    Parameters
    ----------
    auto_df_eesi : Dataframe
        The autosampler EESI df.
        Containing only the QC sampes to compare and
        only one datatype of these samples (e.g. normalized, wb subtracted
        and semiquantified data: '_norm_wb_sub_semiquant').
    axis_range : List of floats
        Range of axis in format: [start,end].
    data_type : String
        The datatype that was already selected for auto_df_eesi,
        e.g. '_norm_wb_sub_semiquant'.
    axis_scale : string, optional
        Spcifies linear or log scale. The default is 'linear'.
        log scale: 'log' (range selection doesn't apply for log scale, should be
        scaled by autoscale in web app).
    diagonal_line : Boolean, optional
        Activate a diagonal line through the origin, slope 1, by setting this to True.
        The default is False.
    
    Returns
    -------
    None.

    '''
    
    # Copy input to make sure it's not overwritten
    auto_df_in = auto_df_eesi.copy()
    
    # Get MS data part of auto_df 
    auto_df_in_ms = general.get_autosampler_ms_part(auto_df_in)

    # Transpose dataframe and get column names (ions) as a column
    auto_df_trans = auto_df_in_ms.T.reset_index()
    # Cut off datatype specifier from AMS ion names 
    auto_df_trans['index'] = auto_df_trans['index'].str.replace(data_type, '')
    # Set column names to 'EESI Ion' for EESI ions and the QC sample names for ion fractions
    auto_df_trans.columns = ['EESI Ion'] + list(auto_df_in['sample'])




    # Get a dictionary with all QC filter dimensions and the corresponding data
    dimensions =[]
    for i in list(auto_df_in['sample']):
        dimensions.append(dict(label=i, values=auto_df_trans[i]))



    # Scatter matrix visualization

    fig = go.Figure(data=go.Splom(
                    dimensions = dimensions,
                    #diagonal_visible=False, # remove plots on diagonal
                    showupperhalf=False, # Remove mirrored plots on upper half
                    text=auto_df_trans['EESI Ion'],
                    ))

    if diagonal_line==True:
        n = len(dimensions)
        for i in range(n):
            for j in range(i+1, n):
                max_val = max(dimensions[i]['values'].max(), dimensions[j]['values'].max())
                fig.add_trace(go.Scatter(
                    x=[0, max_val],
                    y=[0, max_val],
                    mode='lines',
                    line=dict(color='red', width=1),
                    showlegend=False,
                    xaxis=f'x{i+1}',
                    yaxis=f'y{j+1}'
                ))


    # Change scale of all x-axes
    for i in range(1, len(auto_df_in['sample']) + 1):
        fig.update_layout(**{f'xaxis{i}': dict(range=axis_range)})
        

    # Change scale of all y-axes
    for i in range(1, len(auto_df_in['sample']) + 1):
        fig.update_layout(**{f'yaxis{i}': dict(range=axis_range)})

    fig.update_xaxes(type=axis_scale)
    fig.update_yaxes(type=axis_scale)

    fig.show()






###############################################################################
############################# Execution Panel #################################
###############################################################################





'''
Import Data
'''

# Autosampler dataframe with EESI ions after averaging and scaling to ambient conditions(semi-quantification)
auto_df_eesi = impexp.import_file('Intermediate autosampler EESI data (after semi-quantification)')


'''
Process Data
'''

# Mark EESI_note according to general samples (i.e. excluding double measuremnts)
auto_df_cut = general.sample_selection(auto_df_eesi, 'EESI_note_general')
# Cut to selected samples in AMS_note
auto_df_cut = general.sample_cut_to_selection(auto_df_cut)
# Filter sampel type and data type
auto_df_type = general.ams_intermediate_sorter(auto_df_cut, 
                                                     sample_type=[
                                                         #'Sample',
                                                         #'SDE Event Sample',
                                                         #'BB Event Sample',
                                                         #'HR Sample',
                                                         #'Fieldblank',
                                                         'BB Fieldblank',
                                                         'SDE Fieldblank',
                                                         'HR Fieldblank'
                                                         ], 
                                                     data_type='_norm_wb_sub_semiquant')


## Extract Sampling dates, months and seasons
auto_df_season = general.sampling_season_month_parser_per_sample(auto_df_type)



# Filtering carbonaceous ion species
auto_df_carbon = eesi_filter.drop_ion_by_list(auto_df_season, ion_filter = 'EESI Ion')


# Filtering based on sample/wb ratios
auto_df_wb_filter = eesi_filter.drop_ion_by_list(auto_df_carbon, ion_filter = 'EESI Ion')

## Calculate total EESI ion intensity
auto_df_total_ion = parsing.eesi_total_ion_intensity(auto_df_wb_filter)







###############################################################################
############################## Visualization ##################################
###############################################################################







'''
Total ion intensity Fieldblanks
'''

# By sampling date
vis.quick_plot_box(auto_df_total_ion, 'Total ion intensity', 'sample_starting_time', 
                   hover_data='sample_starting_time',
                   color='Sampling Season')


# By type
vis.quick_plot_box(auto_df_total_ion, 'Total ion intensity', 'type', 
                   hover_data='sample_starting_time',
                   color='Sampling Season')






'''
Total ion intensity Samples vs Fieldblanks
'''




# Create a categorical column to mark declustering change
# Define the cutoff date
cutoff_date = pd.Timestamp('2023-02-05')

# Convert time column in autosampler df to datetime format
date_cols = ['sample_starting_time']
auto_df_total_ion[date_cols] = auto_df_total_ion[date_cols].apply(pd.to_datetime) # Convert to datetime format

# Check if any of the sampling date columns have a date after the cutoff
auto_df_total_ion['Declustering Category'] = auto_df_total_ion[['sample_starting_time']].gt(cutoff_date).any(axis=1)
# Convert boolean to categorical labels
auto_df_total_ion['Declustering Category'] = auto_df_total_ion['Declustering Category'].map({True: 'After Declustering Change', False: 'Before Declustering Change'})


# By type
vis.quick_plot_box(auto_df_total_ion, 'C6H10NaO5+_norm_wb_sub_semiquant', 'type', 
                   hover_data='sample_starting_time',
                   #color='Declustering Category',
                   category_order= ['Sample',
                       'SDE Event Sample',
                       'BB Event Sample',
                       'HR Sample',
                       'Fieldblank',
                       'BB Fieldblank',
                       'SDE Fieldblank',
                       'HR Fieldblank'
                       ])












'''
Relative ion intensities
'''

# Normalize EESI ions by total intensity of each sample -> ion fractions
auto_df_fractions = parsing.eesi_ion_fractions(auto_df_wb_filter)

# Scatter matrix of Levoglucosan repeats
qc_filter_scatter_matrix_eesi(auto_df_fractions, [0, 0.01], '_norm_wb_sub_semiquant', 
                              axis_scale='linear', 
                              diagonal_line=True
                              )















