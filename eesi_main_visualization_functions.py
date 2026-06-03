# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 16:16:13 2026

@author: weng_j
"""


# library importing

import numpy as np
import pandas as pd


from scipy import stats # used for linear regression

import plotly.express as px # interactive java script based visualization
import plotly.io as pio # setting of default renderer of plotly figures (has to be set to 'browsers' to get interactive plots)
import plotly.graph_objects as go
from plotly.subplots import make_subplots # import make_subplots function from plotly.subplots to make grid of plots


config = {
    "toImageButtonOptions": {
        "format": "png",
        "scale": 3
    }}








## code 1.1 visualization of sample & wb assignment
def assigned_rn_visualizer (autosampler_data, df_data, ion, x_axis='rn', by_sample_type = None):
    '''
    Creates an interactive plot of signal intensity vs. run number (RN) for a selected ion;
    the RN assigned to samples/QC and their waterblanks (WB) in the previous step are highlighted;
    therefore this funtion can be used for quality control of the assignment.

    Parameters
    ----------
    autosampler_data : DataFrame
        The autosampler dataframe after assignment of RN to samples/QC/WB.
    df_data : DataFrame
        A MS dataframe after removal of 'bad' RN in previous step.
    ion : String
        The ion of interest (or ion group if assigned prevously in the df).
    x_axis: String, optinal
        The domain to plot on the x-axis. For eesi, can be either rn or time.
    by_sample_type: String, optional
        A string specifiing a certain sample type that should be colored only.
        If this is None, no sample type is specified and all sample and wb assigned rn
        will be colored.
        The default is None.

    Returns
    -------
    None.

    '''
    # Copy data to make sure they're not overwritten
    df_autosampler = autosampler_data.copy()
    df = df_data.copy()
    
    # Flaging of samples in df
    df["sample"] = np.nan
    df["type"] = np.nan
    for i in range(0,len(df_autosampler)):
        rn = df_autosampler["sample_rn"][i]
        df['sample'].loc[df['rn'].isin(rn)] = df_autosampler['sample'][i]
        df['type'].loc[df['rn'].isin(rn)] = df_autosampler['type'][i]

    for i in range(0,len(df_autosampler)):
        rn = df_autosampler["back_rn"][i]
        df['sample'].loc[df['rn'].isin(rn)] = ('WB_' + df_autosampler['sample'][i])
        df['type'].loc[df['rn'].isin(rn)] = ('WB_' + df_autosampler['type'][i])
    
    
    # Interactive visualization wiht plotly
    
    if by_sample_type:
        pio.renderers.default='browser' # change default renderer to web browser
        fig = px.scatter(df, x=x_axis, y=ion, color=df['type'] == by_sample_type)
        fig.show(config=config)
    else:
        pio.renderers.default='browser' # change default renderer to web browser
        fig = px.scatter(df, x=x_axis, y=ion, color=df['sample'].isnull())
        fig.show(config=config)
    







#### Quick plot functions ####




def quick_plot_line(data, y_name, x_name, y_title=None, x_title=None, font_size = 46, line_width=10):
    '''
    Generates a simple line plot with plotly.

    Parameters
    ----------
    data : Dateframe
        The dataframe that contains the data to plot.
    y_name : String
        Column name of the y-axis variable.
    x_name : String
        Column name of the x-axis variable.
    y_titel : String, optional
        A optional titel to specify for the y-axis. Default is None.
    x_title : String, optional
        A optional titel to specify for the x-axis. Default is None.
    font_size : Integer, optional
        The font size of axis titles. The default is 46.
    line_width: Integer, optional
        The line width. The default is 10.

    Returns
    -------
    None.

    '''
    
    
    # change default renderer to web browser
    pio.renderers.default='browser' 
    
    # plot
    fig = px.line(data, y = y_name, x = x_name, 
                 template="simple_white")
    
    # some aesthetics
    if x_title:
        fig.update_xaxes(title=x_title, showline=True, linewidth=2, linecolor='black', mirror=True)
    else:
        fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    
    if y_title:
        fig.update_yaxes(title=y_title, showline=True, linewidth=2, linecolor='black', mirror=True)
    else:
        fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    
    fig.update_traces(line={'width': line_width})
    fig.update_layout(font_size = font_size)
    fig.show(config=config)


