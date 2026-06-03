# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 16:12:54 2025

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


import sys

###############################################################################
## Insert path for offline-ams-pipeline here, to import functions from there ##
###############################################################################
path = ''

sys.path.append(path)  # Set absolute path of that directory
import ams_main_visualization_functions as vis # import visualization functions
import ams_main_general_functions as general # import general functions
import ams_ion_filtering as ion_filter # filtering AMS ions





###############################################################################
############################### Functions #####################################
###############################################################################






###############################################################################
############################# Execution Panel #################################
###############################################################################

'''
Import Data
'''

# Autosampler dataframe with EESI ions after averaging and scaling to ambient conditions(semi-quantification)
auto_df_eesi = impexp.import_file('Intermediate autosampler EESI data (after semi-quantification)')
# Some formating
date_cols = ['Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4']
auto_df_eesi[date_cols] = auto_df_eesi[date_cols].apply(pd.to_datetime) # Convert to datetime format


# EMPA IC data
ic_empa = impexp.import_file('EMPA IC and gravimetry data')
# Some formating
date_cols = ['Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4']
ic_empa[date_cols] = ic_empa[date_cols].apply(pd.to_datetime) # Convert to datetime format

# Grenoble IC data
ic_ige = impexp.import_file('IGE IC data')
# Convert entries with conc. below detection and quantification limit to nan
ic_ige.replace(['<QL', '<DL'], np.nan, inplace=True)
# Some formating
date_cols = ['Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4']
ic_ige[date_cols] = ic_ige[date_cols].apply(pd.to_datetime) # Convert to datetime format



# Grenoble LC-MSMS data
organic_markers = impexp.import_file('IGE LC-MSMS data')
# Convert entries with conc. below detection and quantification limit to nan
organic_markers.replace(['<QL', '<LQ', '<DL', '<LD'], np.nan, inplace=True) # Convert entries with conc. below detection and quantification limit to nan



## SDE flags
# Time series SD flags
SD_flags_11_years = impexp.import_file('11 year series SD flags')
SD_flags_HR = impexp.import_file('HR year series SD flags')
# Some formating
date_cols = ['Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4']
SD_flags_11_years[date_cols] = SD_flags_11_years[date_cols].apply(pd.to_datetime) # Convert to datetime format
SD_flags_HR[date_cols] = SD_flags_HR[date_cols].apply(pd.to_datetime) # Convert to datetime format

# Merge the flag data
SD_flags_all = pd.concat([SD_flags_11_years, SD_flags_HR], ignore_index=True)


'''
Formating
'''

## EESI
# Mark EESI_note according to general samples (i.e. excluding double measuremnts)
auto_df_eesi_cut = general.sample_selection(auto_df_eesi, 'EESI_note_general')
# Cut to selected samples in AMS_note
auto_df_eesi_cut = general.sample_cut_to_selection(auto_df_eesi_cut)
# Filter sampel type and data type
auto_df_eesi_cut_type = general.ams_intermediate_sorter(auto_df_eesi_cut, 
                                                     sample_type=[
                                                         'Sample'
                                                         ], 
                                                     data_type='_norm_wb_sub_semiquant')




# EMPA IC
# Remove all samples with sampling data after 01/01/2020, as only monthly composites measured after then
# Define the cutoff date
cutoff_date = pd.Timestamp('2020-01-01')
# Filter out rows where any sampling date column has a date after the cutoff
ic_empa_cut = ic_empa[~(ic_empa[['Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4']] > cutoff_date).any(axis=1)]
# Filter by sample type
sample_type = ['Sample',
               #'BB Event Sample',
               #'SDE Event Sample'
               ]
# Cut to defined sample types
ic_empa_type = ic_empa_cut.query('type in @sample_type').reset_index(drop=True) 


# IGE IC
# Convert from ng/m3 to ug/m3
ic_ige.iloc[:, 6:] = ic_ige.iloc[:, 6:] / 1000
# Filter by sample type
sample_type = ['Sample',
               #'BB Event Sample',
               #'SDE Event Sample'
               ]
# Cut to defined sample types
ic_ige_type = ic_ige.query('type in @sample_type').reset_index(drop=True) 


# IGE LC-MSMS
# Convert from ng/m3 to ug/m3
organic_markers.iloc[:, 6:] = organic_markers.iloc[:, 6:] / 1000
# Filter by sample type
sample_type = ['Sample',
               #'BB Event Sample',
               #'SDE Event Sample'
               ]
# Cut to defined sample types
organic_markers_type = organic_markers.query('type in @sample_type').reset_index(drop=True) 



###############################################################################
################################ Merging data #################################
###############################################################################


data = auto_df_eesi_cut_type.merge(ic_empa_type, how = 'left', on=['sample', 'type', 'Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4'])
data = data.merge(ic_ige_type, how = 'left', on=['sample', 'type', 'Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4'])
data = data.merge(organic_markers_type, how = 'left', on=['sample', 'type', 'Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4'])
data = data.merge(SD_flags_all, how = 'left', on=['sample', 'type', 'Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4'])


# Sort dataframe in chronological time order
data = data.sort_values(by='Sampling Date 1')

## Calculate some combinations of ions that are important for visualization
# all C6H12O6 isomers
data['LC-MSMS Glucose+Fructose+Galactose'] = data['Glucose'].fillna(0) + data['Fructose'].fillna(0) + data['Galactose'].fillna(0)

test = data['LC-MSMS Glucose+Fructose+Galactose']

# Remove all samples with sampling data after 01/01/2019, as declustering change in EESI data
# Define the cutoff date
cutoff_date = pd.Timestamp('2019-01-01')
# Filter out rows where any sampling date column has a date after the cutoff
data_cut = data[~(data[['Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4']] > cutoff_date).any(axis=1)]

# Alternative: create a categorical column to mark declustering change
# Check if any of the sampling date columns have a date after the cutoff
data['Declustering Category'] = data[['Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4']].gt(cutoff_date).any(axis=1)
# Convert boolean to categorical labels
data['Declustering Category'] = data['Declustering Category'].map({True: 'After Declustering Change', False: 'Before Declustering Change'})




# Filter out SD events from 11 year timeseries if necessary
data = data[data['Titanium dust flag'] != 'Dust dominated']
data_cut = data_cut[data_cut['Titanium dust flag'] != 'Dust dominated']


###############################################################################
############################## Visualization ##################################
###############################################################################


'''
Validation of IC data

EMPA IC vs. IGE IC
'''
# Scatter: EMPA IC vs. IGE IC sulfate
vis.quick_plot_scatter(data=data,
                       x='SO4[ug/m3]_mean',
                       y= 'SO42-',
                       x_title = 'EMPA IC sulfate [μg/m<sup>3</sup>]',
                       y_title = 'IGE IC sulfate [μg/m<sup>3</sup>]',
                       #color='type',
                       log=False)

# Scatter: EMPA IC vs. IGE IC nitrate
vis.quick_plot_scatter(data=data,
                       x='NO3[ug/m3]_mean',
                       y= 'NO3-',
                       x_title = 'EMPA IC nitrate [μg/m<sup>3</sup>]',
                       y_title = 'IGE IC nitrate [μg/m<sup>3</sup>]',
                       #color='type',
                       log=False)

# Scatter: EMPA IC vs. IGE IC chloride
vis.quick_plot_scatter(data=data,
                       x='Cl[ug/m3]_mean',
                       y= 'Cl-',
                       x_title = 'EMPA IC chloride [μg/m<sup>3</sup>]',
                       y_title = 'IGE IC chloride [μg/m<sup>3</sup>]',
                       #color='type',
                       log=False)


'''
Sulfate

EESI vs. IC
'''
# Individual EESI ions
eesi_ion = 'Na3SO4+' # Na3SO4+
eesi_ion = 'INa4O4S+' # Na3SO4(NaI)+

# Scatter: EESI Sulfate vs IC Sulfate
vis.quick_plot_scatter(data=data,
                       x='SO42-',
                       y= eesi_ion + '_norm_wb_sub_semiquant',
                       x_title = 'IGE IC sulfate [μg/m<sup>3</sup>]',
                       y_title = eesi_ion + '[a.u.]',
                       #color='type',
                       log=False)


# Combination of relevant ambient quantified sulfate ions
data['EESI combined sulfate'] = data['Na3SO4+_norm_wb_sub_semiquant'] + data['INa4O4S+_norm_wb_sub_semiquant']

# Scatter: EESI combined Sulfate vs IC Sulfate
vis.quick_plot_scatter(data=data,
                       x='EESI combined sulfate',
                       y= 'SO42-',
                       x_title = 'EESI sulfate ambient conc. [μg/m<sup>3</sup>]',
                       y_title = 'IC sulfate ambient conc. [μg/m<sup>3</sup>]',
                       #color='type',
                       log=False)


# Linear regression

vis.quick_linear_regression(data['EESI combined sulfate'], data['SO42-'])




# Time series, one scale



vis.plot_line_multiple_one_scale(x = data['Sampling Date 1'],
                             y_traces = [data['EESI combined sulfate'], data['SO42-']],
                             y_legend_entries= ['Offline EESI', 'Ion Chromatography'],
                             x_axis_title = 'Sampling Date',
                             y_axis_title = 'Sulfate ambient conc. [μg/m<sup>3</sup>]')




'''
Nitrate

EESI vs. IC
'''
          
# Individual EESI ions
eesi_ion = 'Na2NO3+' # Na2NO3+ ## Known to be heavily affected by declustering change in JFJ dataset
eesi_ion = 'INNa3O3+' # Na2NO3(NaI)+
eesi_ion = 'Na2NO3(H2O)+' # Na2NO3(H2O)+
eesi_ion = 'Na2NO3(NaI)(H2O)+' # Na2NO3(NaI)(H2O)+


# Maybe water clusters relevant: to be fit and checked for

# Scatter: EESI Nitrate vs IC Nitrate
vis.quick_plot_scatter(data=data,
                       x='NO3-',
                       y= eesi_ion + '_norm_wb_sub_semiquant',
                       x_title = 'IGE IC nitrate [μg/m<sup>3</sup>]',
                       y_title = eesi_ion + '[a.u.]',
                       #color='type',
                       log=False)



# Combination of relevant ambient quantified nitrate ions ## to be checked if water clusters are relevant once they're fit
data['EESI combined nitrate'] = data['Na2NO3+_norm_wb_sub_semiquant'] + data['INNa3O3+_norm_wb_sub_semiquant']

# Scatter: EESI combined Nitrate vs IC Nitrate
vis.quick_plot_scatter(data=data,
                       x='EESI combined nitrate',
                       y= 'NO3-',
                       x_title = 'EESI nitrate ambient conc. [μg/m<sup>3</sup>]',
                       y_title = 'IC nitrate ambient conc. [μg/m<sup>3</sup>]',
                       #color='type',
                       log=False)


# Linear regression

vis.quick_linear_regression(data.dropna(subset=['NO3-'])['EESI combined nitrate'], data.dropna(subset=['NO3-'])['NO3-'])


# Time series, one scale

vis.plot_line_multiple_one_scale(x = data['Sampling Date 1'],
                             y_traces = [data['EESI combined nitrate'], data['NO3-']],
                             y_legend_entries= ['Offline EESI', 'Ion Chromatography'],
                             x_axis_title = 'Sampling Date',
                             y_axis_title = 'Nitrate ambient conc. [μg/m<sup>3</sup>]')









'''
3-MBTCA (C8H12O6)

EESI vs. IC
'''

# Time series, multiple scales

vis.plot_line_two_scale(x = data_cut['Sampling Date 1'],
                        y_traces = [data_cut['C8H12NaO6+_norm_wb_sub_semiquant'], data_cut['3-MBTCA']],
                        y_legend_entries= ['Offline EESI', 'Ion Chromatography'],
                        x_axis_title = 'Sampling Date',
                        y_axis_titles = ['EESI C8H12O6Na<sup>+</sup> [a.u.]', 'IC 3-MBTCA [μg/m<sup>3</sup>]']
                        )




# Scatter: EESI C8H12O6Na+ vs IC 3-MBTCA
vis.quick_plot_scatter(data=data_cut,
                       x='C8H12NaO6+_norm_wb_sub_semiquant',
                       y= '3-MBTCA',
                       x_title = 'EESI C8H12O6Na<sup>+</sup> [a.u.]',
                       y_title = 'IC 3-MBTCA [μg/m<sup>3</sup>]',
                       #color='Declustering Category',
                       log=False)



# Linear regression
vis.quick_linear_regression(data.dropna(subset=['3-MBTCA'])['C8H12NaO6+_norm_wb_sub_semiquant'], data.dropna(subset=['3-MBTCA'])['3-MBTCA'])







'''
Oxalate/Oxalaic acid (C2H2O4)

EESI vs. IC
'''

# Time series, multiple scales

vis.plot_line_two_scale(x = data_cut['Sampling Date 1'],
                        y_traces = [data_cut['156.985810_norm_wb_sub_semiquant'], data_cut['Oxalate']],
                        y_legend_entries= ['Offline EESI', 'Ion Chromatography'],
                        x_axis_title = 'Sampling Date',
                        y_axis_titles = ['EESI C<sub>2</sub>O<sub>4</sub>Na<sup>3+</sup> [a.u.]', 'IC Oxalate [μg/m<sup>3</sup>]']
                        )




# Scatter: EESI C8H12O6Na+ vs IC 3-MBTCA
vis.quick_plot_scatter(data=data_cut,
                       x='156.985810_norm_wb_sub_semiquant',
                       y= 'Oxalate',
                       x_title = 'EESI C<sub>2</sub>O<sub>4</sub>Na<sup>3+</sup> [a.u.]',
                       y_title = 'Oxalate [μg/m<sup>3</sup>]',
                       #color='Declustering Category',
                       log=False)





# Linear regression
vis.quick_linear_regression(data.dropna(subset=['Oxalate'])['156.985810_norm_wb_sub_semiquant'], data.dropna(subset=['Oxalate'])['Oxalate'])








'''
Mannitol (C6H14O6)

EESI vs. LC-MSMS
'''
#### NOT FIT IN EESI ATM ####

# Time series, multiple scales
vis.plot_line_two_scale(x = data['Sampling Date 1'],
                        y_traces = [data['Mannitol'], data['C6H14O6Na+_norm_wb_sub_semiquant']],
                        y_legend_entries= ['LC-MSMS Mannitol [μg/m<sup>3</sup>]', 'EESI C6H14O6Na+ [a.u.]'],
                        x_axis_title = 'Sampling Date',
                        y_axis_titles = ['LC-MSMS Mannitol [μg/m<sup>3</sup>]', 'EESI C6H14O6Na+ [a.u.]']
                        )




# Scatter: EESI C8H12O6Na+ vs IC 3-MBTCA
vis.quick_plot_scatter(data=data,
                       x='Mannitol',
                       y= 'C6H14O6Na+_norm_wb_sub_semiquant',
                       x_title = 'LC-MSMS Mannitol [μg/m<sup>3</sup>]',
                       y_title = 'EESI C6H14O6Na+ [a.u.]',
                       color='Declustering Category',
                       log=False)




'''
Glucose & isomers (C6H12O6)

EESI vs. LC-MSMS
'''


## just glucose

# Time series, multiple scales
vis.plot_line_two_scale(x = data_cut['Sampling Date 1'],
                        y_traces = [data_cut['Glucose'], data_cut['C6H12NaO6+_norm_wb_sub_semiquant']],
                        y_legend_entries= ['LC-MSMS Glucose [μg/m<sup>3</sup>]', 'EESI C6H12O6Na+ [a.u.]'],
                        x_axis_title = 'Sampling Date',
                        y_axis_titles = ['LC-MSMS Glucose [μg/m<sup>3</sup>]', 'EESI C6H12O6Na+ [a.u.]']
                        )




# Scatter: EESI C8H12O6Na+ vs IC 3-MBTCA
vis.quick_plot_scatter(data=data_cut,
                       x='Glucose',
                       y= 'C6H12NaO6+_norm_wb_sub_semiquant',
                       x_title = 'LC-MSMS Glucose [μg/m<sup>3</sup>]',
                       y_title = 'EESI C6H12O6Na+ [a.u.]',
                       #color='Declustering Category',
                       log=False)



## all C6H12O6 isomers

# Time series, multiple scales
vis.plot_line_two_scale(x = data_cut['Sampling Date 1'],
                        y_traces = [data_cut['LC-MSMS Glucose+Fructose+Galactose'], data_cut['C6H12NaO6+_norm_wb_sub_semiquant']],
                        y_legend_entries= ['LC-MSMS Glucose+Fructose+Galactose [μg/m<sup>3</sup>]', 'EESI C6H12O6Na+ [a.u.]'],
                        x_axis_title = 'Sampling Date',
                        y_axis_titles = ['LC-MSMS Glucose+Fructose+Galactose [μg/m<sup>3</sup>]', 'EESI C6H12O6Na+ [a.u.]']
                        )











'''
Arabitol & isomers (C5H12O5)

EESI vs. LC-MSMS
'''




# Remove all samples with sampling data after 01/01/2019, as declustering change in EESI data
# Define the cutoff date
cutoff_date = pd.Timestamp('2015-01-01')
# Filter out rows where any sampling date column has a date after the cutoff
data_cut = data[~(data[['Sampling Date 1', 'Sampling Date 2', 'Sampling Date 3', 'Sampling Date 4']] > cutoff_date).any(axis=1)]






# Time series, multiple scales
vis.plot_line_two_scale(x = data_cut['Sampling Date 1'],
                        y_traces = [data_cut['Arabitol'], data_cut['C5H12NaO5+_norm_wb_sub_semiquant']],
                        y_legend_entries= ['LC-MSMS Arabitol [μg/m<sup>3</sup>]', 'EESI C5H12O5Na+ [a.u.]'],
                        x_axis_title = 'Sampling Date',
                        y_axis_titles = ['LC-MSMS Arabitol [μg/m<sup>3</sup>]', 'EESI C5H12O5Na+ [a.u.]']
                        )




# Scatter
vis.quick_plot_scatter(data=data_cut,
                       x='Arabitol',
                       y= 'C5H12NaO5+_norm_wb_sub_semiquant',
                       x_title = 'LC-MSMS Arabitol [μg/m<sup>3</sup>]',
                       y_title = 'EESI C5H12NaO5+ [a.u.]',
                       #color='Declustering Category',
                       log=False)




'''
Levoglucosan (C6H10O5)

EESI vs. LC-MSMS
'''





# Time series, multiple scales
vis.plot_line_two_scale(x = data_cut['Sampling Date 1'],
                        y_traces = [data_cut['Levoglucosan'], data_cut['C6H10NaO5+_norm_wb_sub_semiquant']],
                        y_legend_entries= ['LC-MSMS Levoglucosan [μg/m<sup>3</sup>]', 'EESI C6H10O5Na+ [a.u.]'],
                        x_axis_title = 'Sampling Date',
                        y_axis_titles = ['LC-MSMS Levoglucosan [μg/m<sup>3</sup>]', 'EESI C6H10O5Na+ [a.u.]']
                        )




# Scatter
vis.quick_plot_scatter(data=data_cut,
                       x='Levoglucosan',
                       y= 'C6H10NaO5+_norm_wb_sub_semiquant',
                       x_title = 'LC-MSMS Levoglucosan [μg/m<sup>3</sup>]',
                       y_title = 'EESI C6H10NaO5+ [a.u.]',
                       #color='Declustering Category',
                       log=False)









'''
2-Methyl-Tetrol (C5H12O4)

EESI vs. LC-MSMS
'''



# Time series, multiple scales
vis.plot_line_two_scale(x = data_cut['Sampling Date 1'],
                        y_traces = [data_cut['C5H12NaO4+_norm_wb_sub_semiquant'], data_cut['2-Methyltetrols']],
                        y_legend_entries= ['Offline EESI', 'UHPLC-MS/MS'],
                        x_axis_title = 'Sampling Date',
                        y_axis_titles = ['EESI C<sub>5</sub>H<sub>12</sub>O<sub>4</sub>Na<sup>+</sup> [a.u.]', 'UHPLC-MS/MS 2-Methyltetrols [μg/m<sup>3</sup>]']
                        )




# Scatter
vis.quick_plot_scatter(data=data_cut,
                       x='C5H12NaO4+_norm_wb_sub_semiquant',
                       y= '2-Methyltetrols',
                       x_title = 'EESI C<sub>5</sub>H<sub>12</sub>O<sub>4</sub>Na<sup>+</sup> [a.u.]',
                       y_title = 'UHPLC-MS/MS 2-Methyltetrols [μg/m<sup>3</sup>]',
                       #color='Declustering Category',
                       log=False)




 

# Linear regression
vis.quick_linear_regression(data.dropna(subset=['2-Methyltetrols'])['C5H12NaO4+_norm_wb_sub_semiquant'], data.dropna(subset=['2-Methyltetrols'])['2-Methyltetrols'])
























