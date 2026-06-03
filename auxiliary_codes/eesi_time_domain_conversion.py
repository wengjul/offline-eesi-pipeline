# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 13:15:19 2024

@author: weng_j
"""


# Import modules from .py files
import eesi_import_export_functions as impexp # import and export functions
from datetime import datetime
import pandas as pd
import pytz # Conversion of time zones



'''
Conversion of Tofware output t_series.txt time domain.
This domain is privided in seconds from start of Igor time domain start (1.1.1904).
Has to be converted to datetime format.
'''

# Import times .txt file
# The unit there is s from Igor time domain start (1.1.1904) 
times = impexp.import_file('Tofware t_series')
times.rename(columns={'tseries,': 't_series'}, inplace=True)
times['t_series'] = times['t_series'].str.lstrip(',').astype(float)


# IMPORTANT! Python timestamp start from 1970-01-01; Igor timestamp start from  1904-01-01
# Calculate difference in seconds to correct for it later
py_start_date = datetime(1970, 1, 1)
igor_start_date = datetime(1904, 1, 1)
difference = (py_start_date - igor_start_date)
delta_seconds = difference.total_seconds()


# Transfer the times to seconds from start of python domain
# Times - time domain difference (Python vs Igor)
times_python = times - delta_seconds


# Convert seconds since the Unix epoch to datetime and then to the desired string format
times_python_convert = pd.DataFrame(pd.to_datetime(times_python['t_series'], unit='s').dt.strftime('%Y-%m-%d %H:%M:%S'))

# Convert from UTC to MEZ
times_python_convert_mez = pd.to_datetime(times_python_convert['t_series']).dt.tz_localize('UTC').dt.tz_convert('Europe/Berlin').dt.tz_localize(None)

# Export to csv
impexp.export_tofile(times_python_convert_mez,'t_series_datetime.csv')




























