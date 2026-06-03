# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 14:43:03 2025

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
import eesi_ion_formula_parsing as parsing # EESI ion formula regex parsing
import eesi_visualization_functions as eesi_vis
import eesi_ion_filtering as eesi_filter



import sys


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







###############################################################################
########################### Execution of Filters ##############################
###############################################################################



'''
Import Data 
'''

# EESI autosampler dataframe after averaging and scaling to ambient conditions (semi-quantification)
auto_df_semiquant = impexp.import_file('Intermediate autosampler EESI data (after semi-quantification)')

# Mark EESI_note according to general samples (i.e. excluding double measuremnts)
auto_df_cut = general.sample_selection(auto_df_semiquant, 'EESI_note_general')
# Cut to selected samples in AMS_note
auto_df_cut = general.sample_cut_to_selection(auto_df_cut)

## Extract Sampling dates, months and seasons
auto_df_season = general.sampling_season_month_parser_per_sample(auto_df_cut)







'''
Filter ions by elements in sum formulas
'''


# EESI elements filter
eesi_filter.filter_by_elements_eesi(auto_df_season, elements= {'C'})


'''
Filter ions by percentage of positive intensity/concentration among selection of samples
'''

eesi_filter.filter_negative_ions_eesi(auto_df_season, threshold = 90, 
                                      data_type = '_norm_s', 
                                      sample_type = ['Sample', 'SDE Event Sample'])

'''
Filter ions by percentage of samples that have intensity a factor above waterblank level
'''

# EESI sample to waterblank comparison ion filter
eesi_filter.filter_ions_wb_criteria_eesi(auto_df_season, threshold_percentage = 40, factor = 3, 
                                         sample_type=[
                                         #'Sample',
                                         'SDE Event Sample', 
                                         'HR Sample', 
                                         #'Fieldblank',
                                         #'BB Fieldblank',
                                         #'SDE Fieldblank',
                                         #'HR Fieldblank',
                                         ])




'''
Filter ions by their relative signal deviation between the QC Filter samples
'''

# Select only QC filter samples. Normalized, WB subtracted, ambient scaled data
auto_df_type_qc = general.ams_intermediate_sorter(auto_df_season, 
                                                     sample_type=[
                                                         'QC Filter'
                                                         #'QC Sample Repeat'
                                                         ], 
                                                     data_type='_norm_wb_sub_semiquant')


# Only consider summer qc samples
auto_df_type_qc_summer = auto_df_type_qc[auto_df_type_qc['Sampling Season'] == 'Summer'].reset_index(drop=True)




# Filter EESI ions by their relative signal deviation between the QC Filter samples 
eesi_filter.filter_by_qc_sample_eesi (auto_df_type_qc_summer, threshold_rel_deviation=0.5, declustering='before')



###############################################################################
############################ Evaluation Panel #################################
###############################################################################









###############################################################################
############################### Test Area #####################################
###############################################################################
































































