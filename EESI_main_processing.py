# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 17:59:53 2026

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
import statsmodels.api as sm
from tkinter import filedialog as fd # module for interactive directory selection

import plotly.express as px # interactive java script based visualization
import plotly.io as pio # setting of default renderer of plotly figures (has to be set to 'browsers' to get interactive plots)


# Import modules from .py files
import eesi_import_export_functions as impexp # import and export functions
import eesi_main_workflow_functions as workflow # main workflow functions
import eesi_main_visualization_functions as vis # import visualization functions




###############################################################################
############################### Functions #####################################
###############################################################################




###############################################################################
################################Execution Panel################################
###############################################################################


#### Import data ####

## Tofware imput
# Sum formulas and corresponding m/z ratios
amus = impexp.import_file('amus .h5 data')

# Time series (after converting to datetime format of python time domain)
#t_series = impexp.import_file('T series in datetime format data')
t_series = impexp.import_file('T series in datetime format example csv data')

# Run number series
#rn_series = impexp.import_file('RN series .h5 data')
rn_series = impexp.import_file('RN series example csv data')


# EESI ion sticks
#mx_data = impexp.import_file('Mx_data .h5 data')
mx_data = impexp.import_file('Mx_data example csv data')


## Autosampler and sample info (incl. extraction volume)
# Processed autosampler file (.xlsx; as .csv file caused problem with date format in excel after editing AMS_note!)
auto_df = impexp.import_file('Autosampler')
# Import sample info file containing the effective spiking volumes
sample_info = impexp.import_file('Sample info')



#### Formating of Input ####
# Initial formating of input dataframes
workflow.raw_df_formater(auto_df, sample_info, t_series, rn_series, amus, mx_data)

# Merging of auto_df and sample_info (left merge on auto_df)
auto_df_sample_info = workflow.sample_info_assigner(auto_df, sample_info)

# Mergin of EESI Tofware output: data, run numbers, times into one df
eesi_df = pd.concat([t_series, rn_series, mx_data], axis=1)






#### Execution of workflow functions ####

## code 0: remove bad runs

# Visualize primary ion Na2I+ to define cutoff threshold values for outliers,
# by visual inspection
# Outliers can be points where the TPS had obvious glitches etc
# atm, this is done based on whole campaign. Should be adjusted to smaller time
# scale, as strong variation of Na2I+ total signal
vis.quick_plot_line(eesi_df, 'Na2I+', 'ms_starting_time', y_title='Na2I+ signal intensity', x_title='Measurement time', font_size = 20, line_width=2)

# Filtering of outliers in primary ion with defined thresholds
eesi_df_c_0 = workflow.bad_run_number_remover_eesi (eesi_df, lower_cutoff = 6500, upper_cutoff = 89000)
# Export to csv, this can be used to start with later
impexp.export_tofile(eesi_df_c_0,'intermediate_eesi_raw.csv')

# Primary ion signal after filtering
vis.quick_plot_line(eesi_df_c_0, 'Na2I+', 'ms_starting_time', y_title='Na2I+ signal intensity', x_title='Measurement time', font_size = 20, line_width=2)





## code 1.1: run number assignment to sample & wb 
auto_df_sample_info_c_1_1 = workflow.auto_df_run_number_assinger_eesi(auto_df_sample_info, eesi_df_c_0)


## Visualization of sample & wb assignment
vis.assigned_rn_visualizer(auto_df_sample_info_c_1_1,eesi_df_c_0,'Na2I+', x_axis='ms_starting_time')
vis.assigned_rn_visualizer(auto_df_sample_info_c_1_1,eesi_df_c_0,'C6H10NaO5+', x_axis='ms_starting_time') # Levoglucosan


## code 1.2: Filtering of outlier datapoints based on median absolute deviation filter
auto_df_sample_info_c_1_2 = workflow.filtering_run_number_remover_eesi(auto_df_sample_info_c_1_1, eesi_df_c_0, ion = 'Na2I+')

## Visualization of sample & wb assignment after filtering
vis.assigned_rn_visualizer(auto_df_sample_info_c_1_2,eesi_df_c_0,'Na2I+', x_axis='ms_starting_time')
vis.assigned_rn_visualizer(auto_df_sample_info_c_1_2,eesi_df_c_0,'C6H10NaO5+', x_axis='ms_starting_time')




## code 2 normalization and averaging


# Code 2.1 Normalization by primary ions
# eesi_df_c_0.iloc[0:250, :]
eesi_df_c_2_1_norm_primion = workflow.norm_by_ions(eesi_df_c_0, ions =  ['Na2I+'])

# Export to csv, this can be used to start with later
impexp.export_tofile(eesi_df_c_2_1_norm_primion,'intermediate_eesi_primion.csv')
### Re-Import intermediate eesi datafiles (after normalization)
eesi_df_c_2_1_norm_primion = impexp.import_file('Intermediate EESI data (after normalization with primary ion)')




# Code 2.2 Normalization by internal standard ions

# Calculate mass concentration of nitrate and sulfate spikes
# 15N-Nitrate concentration (from NH4NO3)
NO3_conc = 3 * (62.9983 / 81.01923) 
# 34-Sulfate concentration (from (NH4)2SO4^2-)
SO4_conc = 3 * (97.94753 / 134.0163)

# Normalize by nitrate ions
eesi_df_c_2_2_norm_no3 = workflow.norm_by_spike_eesi(eesi_df_c_2_1_norm_primion, ions = ['INa3[15N]O3+'], spi_conc=NO3_conc)
# Export to csv, this can be used to start with later
impexp.export_tofile(eesi_df_c_2_2_norm_no3,'intermediate_eesi_no3_INa3[15N]O3+.csv')

# Normalize by sulfate ions
eesi_df_c_2_2_norm_so4 = workflow.norm_by_spike_eesi(eesi_df_c_2_1_norm_primion, ions = ['Na3[34S]O4+', 'INa4[34S]O4+'], spi_conc=SO4_conc)
# Export to csv, this can be used to start with later
impexp.export_tofile(eesi_df_c_2_2_norm_so4,'intermediate_eesi_so4_Na3[34S]O4+_INa4[34S]O4+.csv')


### Re-Import intermediate eesi datafiles (after normalization)
eesi_df_c_2_2_norm_no3 = impexp.import_file('Intermediate EESI data (after normalization with NO3)')
eesi_df_c_2_2_norm_so4 = impexp.import_file('Intermediate EESI data (after normalization with SO4)')





# Code 2.3 Averaging and water blank subtraction

# Code 2.3.1 Averaging 
auto_df_average = workflow.average_per_sample_eesi(auto_df_sample_info_c_1_2, eesi_df_c_2_1_norm_primion, eesi_df_c_2_2_norm_so4)

# Code 2.3.2 Effective extraction volume correction 
auto_df_effectvolume = workflow.sample_extraction_volume_correction_eesi(auto_df_average)

# Code 2.3.3 Water blank subtraction 
auto_df_u = workflow.sample_wb_subtraction_eesi(auto_df_effectvolume)


# Code 3 'Semi-Quantification'
# Correct relative signal intensities for samples with respect to different 
# sampling or measurement conditions (e.g. sampling rate, punch size, 
# number of punches, extraction volume ...)

auto_df_semiquant = workflow.eesi_semi_quantification (auto_df_u)


# Export to csv, this can be used to start with later
impexp.export_tofile(auto_df_semiquant,'auto_df_eesi_semiquant.csv')

### Re-Import intermediate datafiles (after semi-quant)
auto_df_semiquant = impexp.import_file('Intermediate autosampler EESI data (after semi-quantification)')




# Code 4 Fieldblank subtraction

# Fieldblank subtraction with grouping samples and blanks before/after declustering change
# simple median absolute fb intensity approach for subtracted fb mass spectrum

# 'semi-quantified' data
auto_df_fb_sub = workflow.fb_subtraction_eesi_wrap (auto_df_semiquant, 
                                           fieldblank_types = [
                                               'Fieldblank',
                                               'BB Fieldblank',
                                               'SDE Fieldblank',
                                               'HR Fieldblank'
                                               ], 
                                           data_type = '_norm_wb_sub_semiquant')


# Export to csv, this can be used to start with later
impexp.export_tofile(auto_df_fb_sub,'auto_df_eesi_fb_sub.csv')

### Re-Import intermediate datafiles (after semi-quant)
auto_df_fb_sub = impexp.import_file('Intermediate autosampler EESI data (after fb subtraction)')










