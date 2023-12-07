# -*- coding: utf-8 -*-

"""
Module for averaging AAQ data.

"""

# ======================================================================

# general modules
import numpy as np
import pandas as pd
import csv


# ======================================================================


def read_aaq_data(dir_name, core_name, bottom_id):
    
    """
    Read AAQ csv file

    Arguments:
    ----------
    dir_name: string
        Directory of the folder containing the csv file

    core_name: string
        File name of the csv file without ".csv"

    bottom_id: integer
        Bottom cell ID (row number) in the csv file, below which the data is 
        cutoff from the analysis
    
    Returns:
    --------
    headers: list
        Headers of the AAQ variables
    
    data_array: array
        AAQ raw data cut-off at bottom_id
        
    """
    
    with open(dir_name+'/'+core_name+'.csv', 'r') as fp:
        count = 0
        for line in fp:
            count += 1
            #print(count)
            if count == 66:
                headers = line.split(',')
                ncols = len(headers)
                headers = headers[1:ncols]  # Remove "date" column
            if count > 66:
                data_line = line.split(',')
                data_line = data_line[1:ncols]  # Remove "date" column
                data_layer = np.array(data_line)  # Data of the layer or row
                data_layer = data_layer.astype(np.float64)  # Convert string to float
                if count == 67:
                    data_array = data_layer
                elif count <= bottom_id:
                    data_array = np.vstack((data_array, data_layer))
                elif count == bottom_id:
                    break
    
    return headers, data_array


def depth_bin(u, dep, dz, i_opt):
    

    """
    Binning data along water depth with a thickness dz

    Arguments:
    ----------
    u: array
        data for binning

    dep: array
        water depth corresponding to the data array

    dz: float
        bin thicnkess
    
    i_opt: integer
        Option to treat NaN values in bins
        1: remove NaNs from the array, then the length of array becomes shorter
        2: linear interpolation of NaN values
        otherwise: leave NaN values as is

    Returns:
    --------
    u_bin: array
        binned data
        
    dep_bin: array
        depths of bins
        
    """
    
    # Number of bins
    max_dep = max(dep)
    n_bin = round(max_dep/dz) + 1
    
    # Height of each bin (at the center of bin)
    dep_bin = np.zeros([n_bin])
    for z in range(0,n_bin):
        dep_bin[z] = 0.5*dz + (float(z))*dz
        
    # Binning data
    n_data = len(u)
    u_bin = np.zeros([n_bin])
    u_bin[:] = np.nan
    for bin in range(0,n_bin):
        
        h_i = dz * bin - 1.e-06
        
        count = 0
        u_bar = 0
        for i in range(0,n_data):
            if dep[i] >= h_i and dep[i] < h_i+dz:
                count = count + 1
                u_bar = u_bar + u[i]
        if count > 0:
            u_bar = u_bar / count
        else:
            u_bar = np.nan
        u_bin[bin] = u_bar
        
    if i_opt == 1:
        dep_bin2 = np.zeros([n_bin])
        u_bin2 = np.zeros([n_bin])
        count = 0
        for z in range(0,n_bin):
            if not np.isnan(u_bin[z]):
                dep_bin2[count] = dep_bin[z]
                u_bin2[count] = u_bin[z]
                count = count + 1
        dep_bin = dep_bin2[0:count]
        u_bin = u_bin2[0:count]

    elif i_opt == 2:
        # Linear interporation when there is NaN data in bins
        df = pd.DataFrame({'col1': u_bin})
        df_interp = pd.DataFrame(df.interpolate())
        u_bin = df_interp.to_numpy()

    return dep_bin, u_bin


def aaq_bin(data_array, dz, i_opt):
    
    """
    Create binned AAQ data array

    Arguments:
    ----------
    data_array: array
        AAQ raw data cut-off at bottom_id

    dz: float
        bin thicnkess
    
    i_opt: integer
        Option to treat NaN values in bins
        1: remove NaNs from the array, then the length of array becomes shorter
        2: linear interpolation of NaN values
        otherwise: leave NaN values as is
    
    Returns:
    --------
    bin_array: array
        binned AAQ data array
        
    """
    
    # Depth data
    dep = data_array[:,0]
    
    # Number of AAQ water quality variables
    nvar = data_array.shape[1]-1
    for i in range(0,nvar):
        [dep_bin, u_bin] = depth_bin(data_array[:,i+1], dep, dz, i_opt)
        if i == 0:
            bin_array = np.vstack((dep_bin, u_bin))
        else:
            bin_array = np.vstack((bin_array, u_bin))
    bin_array = np.transpose(bin_array)
    
    return dep_bin, bin_array


def csv_write_aaq(write_flag, core_name, headers, bin_array):
    
    """
    Write binned data to csv file

    Arguments:
    ----------
    write_flag: integer
        Option to write output; 1: yes, otherwise: not    
    
    core_name: string
        File name of the csv file without ".csv"

    headers: list
        Headers of the AAQ variables
    
    bin_array: array
        binned AAQ data array
    
    Returns:
    --------
    none
        
    """
    
    if write_flag == 1:
        file_to_write = core_name+'_average.csv'
        with open(file_to_write, 'w', newline="") as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(bin_array)


def fig_param(var_id, bin_array, xlim_temp, xlim_sal, xlim_cond, xlim_EC25, 
              xlim_density, xlim_sigmaT, xlim_Chl_flu, xlim_Chl_a, 
              xlim_turb, xlim_ph, xlim_DO_sat, xlim_DO, xlim_PAR, 
              dx_temp, dx_sal, dx_cond, dx_EC25, dx_density, dx_sigmaT, 
              dx_Chl_flu, dx_Chl_a, dx_turb, dx_ph, dx_DO_sat, dx_DO, dx_PAR):

    """
    Parameters for plotting for a specified variable

    Arguments:
    ----------
    var_id: integer
        Index of a variable to plot
        
    bin_array: array
        binned AAQ data array
        
    xlim_: list of floats
        X-axis range limit for each variable
        
    dx_: list of floats
        X-axis ticks interval for each variable

    Returns:
    --------
    bin_data: array
        binned data of the selected variable
        
    x_lim: array
        x-axis range
        
    x_ticks: array
        x-axis ticks
        
    x_label: string
        x-axis label
        
    plt_title: string
        plot title
        
    """
    
    # Binned data of the selected variable
    bin_data = bin_array[:,var_id]
    
    if var_id == 1:
        x_lim  = xlim_temp
        dx     = dx_temp
        xlabel = 'Temperature ($^\circ$C)'
    elif var_id == 2:
        x_lim  = xlim_sal
        dx     = dx_sal
        xlabel = 'Salinity (psu)'
    elif var_id == 3:
        x_lim  = xlim_cond
        dx     = dx_cond
        xlabel = 'Conductivity (ms/cm)'
    elif var_id == 4:
        x_lim  = xlim_EC25
        dx     = dx_EC25
        xlabel = 'EC25 (\u03BC'+'s/cm)'
    elif var_id == 5:
        x_lim  = xlim_density
        dx     = dx_density
        xlabel = 'Water density (kg/m\u00b3)'
    elif var_id == 6:
        x_lim  = xlim_sigmaT
        dx     = dx_sigmaT
        xlabel = '\u03C3T'
    elif var_id == 7:
        x_lim  = xlim_Chl_flu
        dx     = dx_Chl_flu
        xlabel = 'Chl-Flu (ppb)'
    elif var_id == 8:
        x_lim  = xlim_Chl_a
        dx     = dx_Chl_a
        xlabel = 'Chl-\u03B1 (\u03BC'+'g/L)'
    elif var_id == 9:
        x_lim  = xlim_turb
        dx     = dx_turb
        xlabel = 'Turbidity (FTU)'
    elif var_id == 10:
        x_lim  = xlim_ph
        dx     = dx_ph
        xlabel = 'pH'
    elif var_id == 11:
        x_lim  = xlim_DO_sat
        dx     = dx_DO_sat
        xlabel = 'DO (%)'
    elif var_id == 12:
        x_lim  = xlim_DO
        dx     = dx_DO
        xlabel = 'DO (mg/L)'
    elif var_id == 13:
        x_lim  = xlim_PAR
        dx     = dx_PAR
        xlabel = 'PAR (\u03BC'+'mol '+"$\mathregular{m^{-2} s^{-1}}$"+')'
        
    x_ticks = np.arange(x_lim[0], x_lim[1]+0.5*dx, dx)

    return bin_data, x_lim, x_ticks, xlabel

