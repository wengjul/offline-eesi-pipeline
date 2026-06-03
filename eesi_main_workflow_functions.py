# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 10:54:25 2024

@author: weng_j
"""

# library importing

import sys
import numpy as np
import pandas as pd
import math



# Import some function from 'JFJ-offline-AMS-analysis' directory

import eesi_main_general_functions as general # import general functions





# Formating of  input dataframes

def raw_df_formater(auto_df, sample_info, t_series, rn_series, amus, mx_data):
    '''
    Formats the raw data input dataframes. Changes unspecified 'Object' data types 
    to strings and datetime objects. 
    For NO3 and SO4 dataframes, calculates sums for labeled and unlabeled species.

    Parameters
    ----------
    auto_df : Dataframe
        The input autosampler df.
    sample_info : Dataframe
        The input HRSO4 df.
    t_series: Dataframe
        The tofware time series output.
    rn_series: Dataframe
        The tofware run number series output.
    amus: DataFrame
        The sum formulas and corresponding m/z ratios df
    mx_data: DataFrame
        The tofware eesi MS signal intensity matrix
    Returns
    -------
    None.

    '''
    # Autosampler
    ## Convert 'Sample' and 'type' to string data type
    auto_df['sample'] = auto_df['sample'].astype("string")
    auto_df['type'] = auto_df['type'].astype("string")
    
    ##  Make sure that autosampler dates are in datatime object format (dtype: datetime64[ns])
    auto_df['sample_starting_time'] = pd.to_datetime(auto_df['sample_starting_time'], format='%Y-%m-%d %H:%M:%S')
    auto_df['sample_ending_time'] = pd.to_datetime(auto_df['sample_ending_time'], format='%Y-%m-%d %H:%M:%S')
    auto_df['LWB_starting_time'] = pd.to_datetime(auto_df['LWB_starting_time'], format='%Y-%m-%d %H:%M:%S')
    auto_df['LWB_ending_time'] = pd.to_datetime(auto_df['LWB_ending_time'], format='%Y-%m-%d %H:%M:%S')
    
    
    # Sample Info 
    ## Convert 'Sample' and 'type' to string data type
    sample_info['sample'] = sample_info['sample'].astype("string")
    sample_info['type'] = sample_info['type'].astype("string")
    
    ##  Make sure that sampling dates are in datatime object format (dtype: datetime64[ns])
    sample_info['Sampling Date 1'] = pd.to_datetime(sample_info['Sampling Date 1'], format='%Y-%m-%d')
    sample_info['Sampling Date 2'] = pd.to_datetime(sample_info['Sampling Date 2'], format='%Y-%m-%d')
    sample_info['Sampling Date 3'] = pd.to_datetime(sample_info['Sampling Date 3'], format='%Y-%m-%d')
    sample_info['Sampling Date 4'] = pd.to_datetime(sample_info['Sampling Date 4'], format='%Y-%m-%d')
    
    # Time series of Tofware EESI output
    ## Renaming
    t_series.rename(columns={'t_series': 'ms_starting_time'}, inplace=True)
    ##  Make sure that t_series dates are in datatime object format (dtype: datetime64[ns])
    t_series['ms_starting_time'] = pd.to_datetime(t_series['ms_starting_time'], format='%Y-%m-%d %H:%M:%S')
    
    # Run number series of Tofware EESI output
    ## Renaming
    rn_series.rename(columns={'rn_series': 'rn'}, inplace=True)
    
    # amus i.e. m/z and sum formulas
    amus.rename(columns={'amus.l': 'ion', 'amus.d': 'mz'}, inplace=True)
    
    # mx_data i.e. eesi MS signal intensity matrix
    mx_data.columns = [col.split('_')[1] for col in mx_data.columns]






# Merging of auto df and sample info df

def sample_info_assigner(auto_df, sample_info):
    '''
    Merges the autosampler dataframe with the sample info dataframe; 
    reorders the resulting columns. 

    Parameters
    ----------
    auto_df : TYPE
        DESCRIPTION.
    sample_info : TYPE
        DESCRIPTION.

    Returns
    -------
    auto_df_sample_info : TYPE
        DESCRIPTION.

    '''
    
    # Make sure not to overwrite input
    auto_df_left = auto_df.copy()  
    sample_info_right = sample_info.copy()
    
    
    # Merge the two dataframes: left merge on autosampler_df, on 'sample' and 'type'
    # both keywords required, otherwise ambiguity
    
    auto_df_sample_info = auto_df_left.merge(sample_info_right, how = 'left', on=['sample', 'type'])
    
    # Reorder columns
    columns = ['sample', 'type', 'sample_starting_time', 'sample_ending_time',
               'LWB_starting_time', 'LWB_ending_time', 'Sampling Date 1', 'Sampling Date 2',
               'Sampling Date 3', 'Sampling Date 4', 'Effective Extraction Volume [mL]',
               'AMS_note', 'LTOF_note', 'Orbi_note', 'Notes']
    
    auto_df_sample_info = auto_df_sample_info.loc[:,columns]
    
    return auto_df_sample_info


# Filtering of outliers in primary ion
def bad_run_number_remover_eesi (eesi_df, lower_cutoff, upper_cutoff):
    '''
    Remove 'bad' run numbers, i.e. spectra, based on the primary ion signal; the idea 
    is to track spikes in primary ion, where e.g. the TPS had glitches etc; 
    exclude those RN.

    Parameters
    ----------
    eesi_df : DataFrame
        The EESI Tofware dataframe after merging of time and run number domain.
    lower_cutoff : Float
        Lower cutoff level for Na2I+ identified by visual inspection of Na2I+ time series.
    upper_cutoff : Float
        Upper cutoff level for Na2I+ identified by visual inspection of Na2I+ time series.

    Returns
    -------
    df_out : TYPE
        DESCRIPTION.

    '''
    
    # Copy data to make sure they're not overwritten
    df = eesi_df.copy()

    # Identify RNs that are outside of tresholds set for primary ion
    drop_list = list(df[(df['Na2I+'] < lower_cutoff) | (df['Na2I+'] > upper_cutoff)]['rn'].values)

    # Remove them from the df and reset index
    df_out = df[~df['rn'].isin(drop_list)].reset_index(drop=True)

    print(f'{len(drop_list)} run numbers were removed based on Na2I+ outliers.')
    
    return df_out



## code 1.1: run number assignment to sample & wb 

def auto_df_run_number_assinger_eesi(auto_df_sample_info, eesi_df_c_0):
    '''
    Assigns run numbers (RN) to the corresponding samples/QC and water blanks (WB) in the autosampler file;
    This is solely based on timings in the AMS and autosampler files; so times 
    between instruments have to match.
    Up to now, the first 6 datapoints and the last datapoint for each sample are removed
    in order to assure neat cutoff of background

    Parameters
    ----------
    auto_df_sample_info : DataFrame
        Autosampler file.
    eesi_df_c_0 : DataFrame
        The EESI data file after filtering outliers of rn in code 0.

    Returns
    -------
    df_autosampler : DataFrame
        The autosampler file with two columns containing RN information for i) the samples/QC
        and the ii) the corresponding WB.

    '''
    # Copy data to make sure they're not overwritten
    df_autosampler = auto_df_sample_info.copy()
    df = eesi_df_c_0.copy()

    # assignment of RN to samples
    sample_rn ,back_rn = [],[]

    for i in range(0,len(df_autosampler),1):
        sample_rn_index = df[(df['ms_starting_time']>df_autosampler['sample_starting_time'].iloc[i])&\
                                 (df['ms_starting_time']<(df_autosampler['sample_ending_time'].iloc[i]))]["rn"].values[6:-1]

        back_rn_index = df[(df['ms_starting_time']>df_autosampler['LWB_starting_time'].iloc[i])&\
                                 (df['ms_starting_time']<(df_autosampler['LWB_ending_time'].iloc[i]))]["rn"].values[6:-1]
        sample_rn.append(sample_rn_index)
        back_rn.append(back_rn_index)

    # assigned RN are saved as list for each sample in autosampler file
    df_autosampler['sample_rn'],df_autosampler['back_rn'] = sample_rn,  back_rn

    return df_autosampler




## code 1.2: Filtering of outlier datapoints based on median absolute deviation filter


def bad_rn_based_on_outlier(data,rn, m = 3):
    '''
    Median absolute deviation filter function for RN applied in filtering_run_number_remover_eesi;
    The RN are filtered based on their signal intensity: the  median for each sample/WB is calculated;
    then the difference of each RN to the sample median is calculated; from this population another median
    is calculated and for each RN the ratio of difference/median (differences): 
    if this exceeds a defined value m, the RN is removed

    Parameters
    ----------
    data : List
        Intensities of RNs in one sample/WB.
    rn : List
        Corresponding RNs of that sample.
    m : Float, optional
        Treshold value for excluding RN. The default is 3.

    Returns
    -------
    rn : List
        List of RN that are kept (i.e. not filtered out).

    '''
    if len(data)>3: # only consider samples with more than 3 data points
        d    = np.abs(data - np.median(data)) # for each RN in sample, get absolute value for difference between RN intensity and sample intensity median 
        mdev = np.median(d)
        s    = d/mdev if mdev else 0.
        ind  = [idx for idx, element in enumerate(s) if element<m]
        if len(ind)>0:
            rn   = [rn[i] for i in ind]
    return rn



def filtering_run_number_remover_eesi(auto_df_sample_info_c_1, eesi_df_c_0, ion):
    '''
    Filters out RN based on two criteria for a difined ion signal intensity:
    i) RN that only contain NaN
    ii) RN that deviate to much from the intensity distribution of the corresponding sample or wb
    The criterion for this is a Median absolute deviation filter defined in bad_rn_based_on_outlier.
    Ion specified for filtering should be present in all samples and wb
    (so internal standard or reagent ion).

    Parameters
    ----------
    auto_df_sample_info_c_1 : DataFrame
        The autosampler df after assignemnt of sample and wb run numbers.
    eesi_df_c_0 : DataFrame
        The EESI dataframe after removal of bad RN based on primary ion signal.
    ion : String
        The ion to use for the Median absolute deviation filter.

    Returns
    -------
    df_autosampler : DataFrame
        The autosampler df with filtered rn removed from 'sample_rn' and 'back_rn'.

    '''
    # Copy data to make sure they're not overwritten
    df_autosampler = auto_df_sample_info_c_1.copy()
    df = eesi_df_c_0.copy()
    
    x_index = []

    # Iterate over all samples in autosampler df
    for i in range(0,len(df_autosampler),1):
        # Extract RN for each sample/WB
        s_rn,lwb_rn  = list(df_autosampler['sample_rn'].iloc[i]),list(df_autosampler['back_rn'].iloc[i]) 
        # If the samples has some rn assigned (judged by its wb)
        if len(lwb_rn)>0:
            # Extract signal intensity values for each sample/WB for defined ion
            lwb_value,s_value  = df[df["rn"].isin(lwb_rn)][ion].values, df[df["rn"].isin(s_rn)][ion].values   
            x_index.append(df[df["rn"].isin(lwb_rn)][ion].index[0]) # Index of first value of each of those lists
            
            ### Method 1.: Filtering RN that have nan
            # If signal intensity is not nan, append index to the list
            lwb_ind,s_ind = [idx for idx, element in enumerate(lwb_value) if np.isnan(element)==0],[idx for idx, element in enumerate(s_value) if np.isnan(element)==0]
            # Get run numbers from the indices (non nan samples)
            lwb_rn,s_rn     = [lwb_rn[c] for c in lwb_ind],[s_rn[c] for c in s_ind]
            ### Re-assign those runs that are not nan (others will be droped) 
            df_autosampler['back_rn'].iloc[i],df_autosampler['sample_rn'].iloc[i]  =lwb_rn,s_rn
            
            ### Method 2.: Removal of outliers based on their deviation from the sample median (median absolute deviation)
            df_autosampler['back_rn'].iloc[i]   = bad_rn_based_on_outlier(lwb_value, lwb_rn, m = 3.)
            df_autosampler['sample_rn'].iloc[i] = bad_rn_based_on_outlier(s_value, s_rn, m = 2.)
    
    return df_autosampler







## code 2: normalization and averaging



'''
How we normalize in EESI data processing

Normalization is done by dividing with the signals of following ions

i) Normalization by reagent ion signal (Na2I+):
   This accounts for performance differences of the EESI spray 

ii) Normalization by internal standard ions:
    This accounts for instabilities in the nebulizer performance, generation of
    aerosol and different concentration levels in the samples (higher mass -> more aerosol formed)
    
    Ions to select (from Tianqu's overview presentation):
    (for internal standard, the respective 15N and 34S species)
        Main nitrate related ions:
            Na2NO3+, Na2NO3(NaI)+
            small contribution: Na2NO3(NaI)(H2O)+

        Main sulfate related ions:
            Na3SO4+, Na3SO4(NaI)+
            small contribution: Na3SO4(H2O)+, Na3SO4(NaI)(H2O)+


'''


## code 2.1 normalization by reagent ion


# Normalization by sum of specified ions

def norm_by_ions (eesi_df_c_0, ions):
    '''
    Normalizes all ion signals in the EESI dataframe by the sum of a specified 
    number of ion species (reagent ions).

    Parameters
    ----------
    eesi_df_c_0 : Dataframe
        The EESI dataframe after removing bad run numbers by threshold.
    ions : List of strings
        The ions to normalize with.
        E.g. ['Na2I+'] for first reagent ion;
        ['Na2I+', 'Na2I(H2O)+'] for first reagent ion and its water cluster.

    Returns
    -------
    df_norm : Dataframe
        The EESI dataframe after normalization with specified ions.

    '''
    # Copy input to avoid overwriting
    df_norm = eesi_df_c_0.copy()
    
    # if there were negative values in specified ion intensity, replace with nan (as not physically meaningful then)
    for ion in ions:
        df_norm[ion][df_norm[ion]<0] = np.nan
    
    # Calculate the sum of ion intensities for the specified ions
    norm_intensity = df_norm[ions].sum(axis=1)
    
    # Normalize by dividing with this sum
    df_norm[df_norm.columns[2:]] = df_norm.iloc[:, 2:].div(norm_intensity, axis=0)

    return df_norm



## code 2.2 normalization by internal standard


# Normalization by sum of specified internal standard ions, considering their spiking concentration

def norm_by_spike_eesi (eesi_df_c_2_1, ions, spi_conc):
    '''
    Normalizes all ion signals in the EESI dataframe by the sum of specified
    internal standard ions.
    
    For unlabeled sulfate and nitrate, the following applies:
    conc(species)[mg/L] = signal(species) * conc(spike)[mg/L]/signal(spike)
    
    Parameters
    ----------
    eesi_df_c_2_1 : Dataframe
        The EESI dataframe after normalization with primary ions.
    ions : List of strings
        The ions to normalize with.
    spi_conc : float
        The concentration of the spike (considering only 34SO4 or 15NO3).

    Returns
    -------
    df_norm : Dataframe
        The EESI dataframe after normalization with specified internal standard ions.

    '''
    
    # Copy input to avoid overwriting
    df_norm = eesi_df_c_2_1.copy()
    
    # if there were negative values in specified ion intensity, replace with nan (as not physically meaningful then)
    for ion in ions:
        df_norm[ion][df_norm[ion]<0] = np.nan
    
    # Calculate the sum of ion intensities for the specified ions
    norm_intensity = df_norm[ions].sum(axis=1)
    
    # Normalize by dividing with this sum
    df_norm[df_norm.columns[2:]] = df_norm.iloc[:, 2:].mul(spi_conc / norm_intensity, axis=0)
    
    return df_norm
    





# Code 2.3 Averaging and water blank subtraction

# Code 2.3.1 Averaging 

# Function to extract and average values
def extract_and_average(row, eesi_norm, columns_to_average, datatype):
    '''
    For a row of the autosampler dataframe, extracts the run numbers assigned to
    it. Then calcualtes the average signal intensities of the associated row 
    in the eesi dataframe (i.e. the ones holding these run numbers).
    As specified in datatype, this is either done for samples or corresponding waterblanks.

    Parameters
    ----------
    row : pd.Series
        A row of the autosampler dataframe.
    eesi_norm : Dataframe
        The eesi dataframe, after normalization steps.
    columns_to_average : List
        List of all column names (i.e. ion names) that are considered in averaging per sample/wb.
    datatype : String
        Specifier for datatype, sample vs. waterblank: 'sample_rn', 'back_rn'.

    Returns
    -------
    pd.Series
        The averaged signal intensities.

    '''
    
    rn_list = row[datatype]
    matching_rows = eesi_norm[eesi_norm['rn'].isin(rn_list)]
    if not matching_rows.empty:
        # Average only the specified columns
        return matching_rows[columns_to_average].mean()
    return pd.Series(index=columns_to_average)



# Function to extract the prefix from a column name
def get_prefix(col):
    '''
    Helper function for ordering different datatypes of each ion next
    to each other, after averaging of signals and assignment to autosampler df.

    Parameters
    ----------
    col : String
        Column name.

    Returns
    -------
    String
        The prefix (before underscore) of column name.

    '''
    parts = col.split('_', 1)
    return parts[0] if len(parts) > 1 else col



# Wrapping function for averaging
def average_per_sample_eesi(auto_df_sample_info_c_1_2, eesi_df_c_2_1, eesi_df_c_2_2_norm):
    '''
    Averaging of sample and waterblank signal intensities of the points (run numbers)
    associated to each sample.

    Parameters
    ----------
    auto_df_sample_info_c_1_2 : DataFrame
        The autosampler df after assignment of run numbers to each sample.
    eesi_df_c_2_1 : DataFrame
        The normalized eesi data, by primary ion.
    eesi_df_c_2_2_norm : DataFrame
        The normalized eesi data, by primary ion and internal standard.

    Returns
    -------
    
    auto_df_u : DataFrame
        The autosamper df with signal intensities assigned to each sample.
        The following datatype specifiers are used:
            primary ion normalized sample: _prim_norm_s
            primary ion normalized wb: _prim_norm_wb
            primary ion & internal standard normalized sample: _norm_s
            primary ion & internal standard normalized wb: _norm_wb

    '''
    # copy data to make sure it is not overwritten
    auto_df = auto_df_sample_info_c_1_2.copy()
    eesi_primion_norm = eesi_df_c_2_1.copy()
    eesi_internalstand_norm = eesi_df_c_2_2_norm.copy()

    # Initialize the final autosampler df to store all the average data in
    auto_df_u = auto_df.copy()

    # Columns from eesi df that we want to average (& keep in the correct order)
    columns_to_average = [col for col in eesi_primion_norm.columns if col not in ['rn', 'ms_starting_time']]

    # Primary ion normalized data
    for datatype in ['sample_rn', 'back_rn']:
        # apply the averaging function to each row in auto_df and collect the results in a DataFrame
        averaged_values_df = auto_df.apply(lambda row: extract_and_average(row, eesi_primion_norm, columns_to_average, datatype=datatype), axis=1)
        
        # rename according to datatype
        if datatype == 'sample_rn':
            averaged_values_df.columns = [i + '_prim_norm_s'  for i in averaged_values_df.columns]
        elif datatype == 'back_rn':
            averaged_values_df.columns = [i + '_prim_norm_wb'  for i in averaged_values_df.columns]

        # Concatenate the original auto_df with the new averaged values
        auto_df_u = pd.concat([auto_df_u, averaged_values_df], axis=1)

    # Internal standard normalized data
    for datatype in ['sample_rn', 'back_rn']:
        # apply the averaging function to each row in auto_df and collect the results in a DataFrame
        averaged_values_df = auto_df.apply(lambda row: extract_and_average(row, eesi_internalstand_norm, columns_to_average, datatype=datatype), axis=1)
        
        # rename according to datatype
        if datatype == 'sample_rn':
            averaged_values_df.columns = [i + '_norm_s'  for i in averaged_values_df.columns]
        elif datatype == 'back_rn':
            averaged_values_df.columns = [i + '_norm_wb'  for i in averaged_values_df.columns]

        # Concatenate the original auto_df with the new averaged values
        auto_df_u = pd.concat([auto_df_u, averaged_values_df], axis=1)

    # Get the columns in the right order, i.e. different datatypes of each ion next to each other
    # Extract unique prefixes from column names
    unique_prefixes = []
    prefix_seen = set()

    for col in general.get_autosampler_ms_part(auto_df_u).columns.tolist():
        prefix = get_prefix(col)
        if prefix not in prefix_seen and '_' in col:
            unique_prefixes.append(prefix)
            prefix_seen.add(prefix)

    # Group columns by their prefix while preserving original order
    prefix_groups = {prefix: [] for prefix in unique_prefixes}
    non_prefix_columns = []

    for col in general.get_autosampler_ms_part(auto_df_u).columns.tolist():
        prefix = get_prefix(col)
        if prefix in prefix_groups:
            prefix_groups[prefix].append(col)
        else:
            non_prefix_columns.append(col)

    # Flatten the groups into a single list maintaining the original order
    ordered_columns = non_prefix_columns
    for prefix in unique_prefixes:
        ordered_columns.extend(prefix_groups[prefix])

    # Add autosampler sample info part
    ordered_columns = general.get_autosampler_sample_info_part(auto_df_u).columns.tolist() + ordered_columns

    # Order autosampler df with that ordered column name list
    auto_df_u = auto_df_u[ordered_columns]
    
    return auto_df_u







# Code 2.3.2 Effective extraction volume correction 

def sample_extraction_volume_correction_eesi(auto_df_average):
    '''
    Applies the correction factor for the effective extraction volume to the 
    average intensities of internal standard normalized samples '*_norm_s' 
    in the autosampler dataframe.

    Parameters
    ----------
    auto_df_average : DataFrame
        The autosampler dataframe after addition of averaged ion intensities .

    Returns
    -------
    auto_df : DataFrame
        The autosampler dataframe with effective extraction volume corrected intensities.

    '''
    # copy autosampler data to make sure it is not overwritten
    auto_df = auto_df_average.copy()

    # get column names for internal standard normalized sample data that should be corrected
    norm_s_names = [i for i in auto_df.columns if '_norm_s' in i if not 'prim' in i]

    # get correction factor for effective volume correction
    correction_factor = 12/auto_df['Effective Extraction Volume [mL]']
    # make sure columns with NaN values (i.e. quality control samples that
    # had appropriate volume) get a correction factor of 1 (i.e. intensities unchanged)
    correction_factor.fillna(1, inplace = True)

    # for each normalized sample column... 
    for name in norm_s_names:
        # ... apply the correction
        auto_df[name] = auto_df[name] * correction_factor

    return auto_df








# Code 2.3.3 Water blank subtraction 



def sample_wb_subtraction_eesi(auto_df_effectvolume):
    '''
    Performs water blank subtraction for each ion in the autosampler dataframe for 
    the internal standard normalized dataset. Columns are reordered to appear next to the other data
    for each ion.

    Parameters
    ----------
    auto_df_effectvolume : DataFrame
        Autosampler dataframe after effective extraction volume correction.

    Returns
    -------
    auto_df : DataFrame
        Autosampler dataframe with columns for water blank subtracted, normalized data.

    '''
    # copy autosampler data to make sure it is not overwritten
    auto_df = auto_df_effectvolume.copy()

    # get a list of all unique ion names listed in the autosampler dataframe
    ions = [i.replace('_norm_s', '') for i in auto_df.columns if '_norm_s' in i and not 'prim' in i]

    # for each ion...
    for sp in ions:
        # ... perform WB subtraction for normalized data and create a new column 
        auto_df[sp + '_norm_wb_sub'] = auto_df[sp + '_norm_s'] - auto_df[sp + '_norm_wb']
        
        # reorder the newly created column next to the existing ones of the same ion
        # get column index of last column for that ion
        index_family = auto_df.columns.get_loc(sp + '_norm_wb')
        # get list of reordered column names
        cols = list(auto_df.iloc[:,:index_family + 1].columns) + [sp + '_norm_wb_sub'] + list(auto_df.iloc[:,index_family +1 :-1].columns)
        # reorder dataframe according to this list
        auto_df = auto_df[cols]
    return auto_df









# Code 3 'Semi-Quantification'
# Correct relative signal intensities for samples with respect to different 
# sampling or measurement conditions (e.g. sampling rate, punch size, 
# number of punches, extraction volume ...)

def eesi_atmo_transform(auto_df, sample_type, ion, punch_diameter, number_punches, 
                           filter_diameter, sampling_time, sampling_rate,
                           extraction_volume):
    '''
    Scales the specified ion in samples with specified sample_type with sampling
    & sample preparation parameters. This results are saved into the 
    autosampler column _norm_wb_sub_semiquant.
    
    (Adapted from the AMS code)

    Parameters
    ----------
    auto_df : DataFrame
        The autosampler dataframe with normalized, waterblank subtracted data, 
        i.e. 'ion_norm_wb_sub' and columns for ambient conc. data '_norm_wb_sub_ambient'.
    ion : String
        The column title of the AMS ion to be processed. This has to be normalized,
        waterblank subtracted data, i.e. 'ion_norm_wb_sub'
    sample_type : String
        The filter type that should be considered for the calculation.
    punch_diameter : Float
        Punch diameter [mm].
    number_punches : Float
        Number of punches for this filter type.
    filter_diameter : Float
        Effective full filter diameter [mm]. I.e. considering only the actually loaded area. 
    sampling_time : Float
        Sampling time per filter [h].
    sampling_rate : Float
        Air sampling rate of the HV sampler [m3/min].
    extraction_volume : Float
        Sample ectraction volume [ml].

    Returns
    -------
    None.

    '''
    # A
    ## Calculate area of filter that was used for punching
    punch_area = number_punches*math.pi*(punch_diameter/2)**2 #mm2
    
    ## Calculate total effective area of filter, i.e. area that is loaded
    filter_area = math.pi*(filter_diameter/2)**2 #mm2
    
    ## Calculate percentage area of filter that was used for punching
    area_ratio = punch_area/filter_area
    #print(f'{area_ratio*100} percent of filter were taken for punching of {sample_type} type.')
    
    # B
    ## Calculate air sampled per filter
    sampling_volume = sampling_time * 60 * sampling_rate # m3 of air sampled per filter
    #print(f'{sampling_volume} m3 of air sampled per filter for {sample_type} type.')
    
    # C
    ## From the water dissolved concentration (_norm_wb_sub), calculate the ambient concentration (_norm_wb_sub_ambient)
    ## ambient mass concentration = ((extract conc [ug/mL] * extraction volume [ml]) / percent of filter area ectracted) / sampled air volume
    auto_df[ion.replace('_norm_wb_sub', '_norm_wb_sub_semiquant')].loc[auto_df['type'] == sample_type] = ((auto_df[ion].loc[auto_df['type'] == sample_type] * extraction_volume) / area_ratio) / sampling_volume



def eesi_semi_quantification (auto_df_u):
    '''
    Scales the EESI ion signals with the according sampling & sample preparation
    conditions, so that all samples are on the same relative scale. Factors are
    for example the number of punches, extraction volume, sampling rate.

    Parameters
    ----------
    auto_df_u : DataFrame
        The EESI dataframe after averaging and wb subtraction.

    Returns
    -------
    
    
    auto_df : DataFrame
        The EESI dataframe after correcting relative signal intensities for sampling
        & sample preparation discrepancies.

    '''
    # Copy input to avoid overwriting
    auto_df = auto_df_u.copy()

    # Get a list of all ion columns to transform
    ions = [col for col in auto_df.columns if '_norm_wb_sub' in col]


    # Create new columns for ambient concentration transformed data
    for ion in ions:
        # get a column name for this ion
        column_name = ion.replace('_norm_wb_sub', '') + '_norm_wb_sub_semiquant'
        auto_df[column_name] = np.nan
        # reorder the newly created column next to the existing ones of the same ion
        # get column index of last column for that ion
        index_family = auto_df.columns.get_loc(ion)
        # get list of reordered column names
        cols = list(auto_df.iloc[:,:index_family + 1].columns) + [column_name] + list(auto_df.iloc[:,index_family +1 :-1].columns)
        # reorder dataframe according to this list
        auto_df = auto_df[cols]  



    # Normal samples and fieldblanks, HR samples and HR fieldblanks, Quality control repeat 2011 samples, 'Blind ab Paket' and 'Stempel'
    sample_types = ['Sample', 'Fieldblank', 'HR Sample', 'HR Fieldblank', 'QC Sample Repeat', 'Blind ab Paket', 'Stempel']

    # Run ams_atmo_transform for each specified sample type
    for i in sample_types:
        # RUn ams_atmo_transform for each ion in the autosampler dataframe 
        for ion in ions:
            eesi_atmo_transform(auto_df, sample_type = i, ion = ion, punch_diameter = 14, number_punches = 8,
                               filter_diameter = 140, sampling_time = 24, sampling_rate = 0.75,
                               extraction_volume = 12)

    # Special event samples and fieldblanks
    sample_types = ['SDE Event Sample', 'SDE Fieldblank', 'BB Event Sample', 'BB Fieldblank', 'BC Event Sample', 'BC Fieldblank']

    # Run ams_atmo_transform for each specified sample type
    for i in sample_types:
        # RUn ams_atmo_transform for each ion in the autosampler dataframe 
        for ion in ions:
            eesi_atmo_transform(auto_df, sample_type = i, ion = ion, punch_diameter = 14, number_punches = 4,
                               filter_diameter = 140, sampling_time = 24, sampling_rate = 0.75,
                               extraction_volume = 12)

    # Daily QC filters
    sample_types = ['QC Filter']

    # Run ams_atmo_transform for each specified sample type
    for i in sample_types:
        # RUn ams_atmo_transform for each ion in the autosampler dataframe 
        for ion in ions:
            eesi_atmo_transform(auto_df, sample_type = i, ion = ion, punch_diameter = 14, number_punches = 2,
                               filter_diameter = 140, sampling_time = 24, sampling_rate = 0.75,
                               extraction_volume = 12)
    
    return auto_df













## code 4: Fieldblank subtraction



'''
Considerations for fieldblank subtraction of JFJ EESI dataset

From characterization of JFJ fieldblanks follows for subsetting fieldblanks

    - Definitly differentiate conditions before and after declustering change, as not
    only ion intensities but also relative ion intensities change

    - effect of sampling season is observed for some sample types, however not large (~20%), for now don't consider

    - effect of sample type also observed, but only significant for BB Fieldblanks vs others;
    we don't know if that transfers to BB samples or not; for now don't consider


How to perform fieldblank subtraction

    Different options

    - calculate median intensity for each ion, directly from absolute intensity fieldblank spectra
    Advantage: no normalization by ion fractions, no need to think about what to do with negative ions
    Disadvantage: difference of total ion intensity between fieldblanks has influence on median spectra
        therefore not an ideal metric; e.g. relative ion intensity of intermediate intensity fieldblank
        could determine median relative ion intensities for all

    - calculate total ion intensities of all fieldblanks in the (sub-)group; calculate
    ion fractions for fb, then take median; scale median spectrum with median intensities
    Advantage: all fb contribute equally to median intensity calculation, provided no strong negative ion outliers
    Disadvantage: if there are strongly negative ions (wb higher than sample), this will
        lead to distorted ion fractions and therefore distorted median fb unit spectra
    Ion filtering could be an option (with respect to highly negative ions), 
    but then one needs to make assumptions and looses ions at this step which might be problematic for certain sample types
    (e.g. some ions can be of importance for only small fraction of samples, like dust or BB samples and then get filtered)
    Setting all negative ion intensities to 0 before calculating fractions could be an option
    but then some fb that have most values below 0 and some outliers above 0 will get large biases


In conclusion:
    either use direct ensamble metric (best median to avoid outlier influences) 
and accept that some median fieldblanks will determine relative composition of fb spectrum
in the end
    or perform some well characterized ion filtering before fb subtraction steps,
    making sure that balance is hold between excluding highly negative ions that 
    largely impact fb ensemble metric calculation and excluding valuable ions (even if only for subset of samples)

For mow, select option 1 for sake of simplicity, option 2 could be more valuable for time series analysis (PMF)
where low intensity samples play a bigger role

'''




# Fieldblank subtraction for a group of samples

def fb_subtraction_eesi (auto_df_cut, auto_df_fb):
    '''
    For the input autosampler df, perform fb subtraction with median fb spectrum
    that is calculated from simple median of all spectra in auto_df_fb (absolute intensities).

    Parameters
    ----------
    auto_df_cut : DataFrame
        The EESI autosampler df with samples to be considered for fb subtraction.
    auto_df_fb : DataFrame
        The EESI autosampler df with only fb to be considered as set to calcuate median fb spectrum.

    Returns
    -------
    auto_df_cut_ms_fb_sub : DataFrame
        The EESI autosampler df mass spectra part with only fb subtracted data.

    '''
    # Copy input to avoid overwriting
    auto_df_in = auto_df_fb.copy()
    
    # Copy input to avoid overwriting
    auto_df_cut_in = auto_df_cut.copy()
    auto_df_fb_in = auto_df_fb.copy()
    
    ## 1.) Calculate median of fieldblank spectra (absolute intensities)
    
    # Separate sample info and ms part of the fb autosampler df
    auto_df_fb_ms = general.get_autosampler_ms_part(auto_df_fb_in)
    
    # Compute median fb spectrum
    fb_median_spectrum = auto_df_fb_ms.median()
    
    ## 2.) Perform fb subtraction for each sample: sample ion intensities - fb_media_spectrum
    
    # Separate sample info and ms part of all samples
    auto_df_cut_sample_info = general.get_autosampler_sample_info_part(auto_df_cut_in)
    auto_df_cut_ms = general.get_autosampler_ms_part(auto_df_cut_in)
    
    # fb subtraction
    auto_df_cut_ms_fb_sub = auto_df_cut_ms - fb_median_spectrum
    
    # Rename columns by adding _fb_sub to specified datatype
    auto_df_cut_ms_fb_sub.columns = [col + "_fb_sub" for col in auto_df_cut_ms_fb_sub.columns]
    
    return auto_df_cut_ms_fb_sub


    

# Wrapping function for fb subtraction

def fb_subtraction_eesi_wrap (auto_df_eesi, fieldblank_types, data_type):
    '''
    Wrapping function for performing fieldblank subtraction for different groups
    of samples and fieldblanks, sub-grouped based on time before or after declustering
    change. It has to specified which data types and which fieldblank types should be
    considered. Adds fieldblank subtracted data columns next to the specified data type
    columns of the same ions.

    Parameters
    ----------
    auto_df_eesi : DataFrame
        The EESI autosampler dataframe.
    fieldblank_types : List of strings
        The fieldblank types to incorporate for calculation of median fb spectrum
        to be used for fb subtraction; e.g. ['Fieldblank', 'BB Fieldblank', 'SDE Fieldblank', 'HR Fieldblank']
    data_type : string
        The data type to be fieldblank subtracted, e.g. '_norm_wb_sub' for normalized
        & waterblank subtracted data.

    Returns
    -------
    auto_df_out : DataFrame
        The EESI autosampler df with fieldblank subtracted data columns added.

    '''
    # Copy input to avoid overwriting
    auto_df_in = auto_df_eesi.copy()
    
    # Create a categorical column to mark declustering change
    # Define the cutoff date, i.e. when declustering change happened
    cutoff_date = pd.Timestamp('2023-02-05')
    
    # Convert time column in autosampler df to datetime format
    date_cols = ['sample_starting_time']
    auto_df_in[date_cols] = auto_df_in[date_cols].apply(pd.to_datetime) # Convert to datetime format
    
    if 'Declustering Category' in auto_df_in.columns:
        pass
    else:
        # Insert categorical declustering column to sample info part of autosampler df
        auto_df_in.insert(auto_df_in.columns.get_loc('sample_rn'), 'Declustering Category', 
                           auto_df_in[['sample_starting_time']].gt(cutoff_date).any(axis=1))
        # Convert boolean values to categorical labels
        auto_df_in['Declustering Category'] = auto_df_in['Declustering Category'].map({True: 'After Declustering Change', False: 'Before Declustering Change'})
    
    
    # Filter specified datatype
    auto_df_cut = general.ams_intermediate_sorter(auto_df_in, data_type=data_type)
    
    # Filter fieldblanks
    auto_df_fb = general.ams_intermediate_sorter(auto_df_cut, sample_type=fieldblank_types)
    
    # Perform fb subtraction for each subset
    # set1 <-> before declustering change
    auto_df_cut_ms_fb_sub_set1 = fb_subtraction_eesi (auto_df_cut[auto_df_cut['Declustering Category'] == 'Before Declustering Change'], auto_df_fb[auto_df_fb['Declustering Category'] == 'Before Declustering Change'])
    # set2 <-> after declustering change
    auto_df_cut_ms_fb_sub_set2 = fb_subtraction_eesi (auto_df_cut[auto_df_cut['Declustering Category'] == 'After Declustering Change'], auto_df_fb[auto_df_fb['Declustering Category'] == 'After Declustering Change'])
    
    # Re-combine subsets
    auto_df_cut_ms_fb_sub = pd.concat([auto_df_cut_ms_fb_sub_set1, auto_df_cut_ms_fb_sub_set2], axis=0)
    
    # Add _fb_sub columns into original dataframe, for each ion next to the other datatype columns of that ion 
    auto_df_out = auto_df_in.copy()
    # Get the columns that were specified for fb subtraction
    columns_to_insert = general.get_autosampler_ms_part(auto_df_cut).columns
    
    # Insert the columns right after the original ones
    # Get the correct order of columns
    original_cols = list(auto_df_out.columns)
    for col in reversed(columns_to_insert):  # Reverse to maintain order after insertion
        idx = original_cols.index(col) + 1
        original_cols[idx:idx] = [col + '_fb_sub']
    
    # Reorder the dataframe using pd.concat()
    auto_df_out = pd.concat([auto_df_out, auto_df_cut_ms_fb_sub], axis=1)[original_cols]
    
    return auto_df_out






































