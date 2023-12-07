# -*- coding: utf-8 -*-

"""
Script that bins and averages AAQ data.

"""

# ======================================================================

import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt

# --- Configurations

# --- File name specification --- #
core_name = 'AAQ176-CAD_0434_20221128_075426_0001'
bottom_id = 163   # Bottom cell ID (row number) in the csv file, below which the data is cutoff from the analysis
append_name = 'average'
write_flag = 0    # Option to write output; 1: yes, otherwise: not
# --- Averaging parameters --- #
dz = 0.10         # Bin thickness for averaging (m)
i_opt = 1         # Option to treat NaN values in bins
                  # 1: remove NaNs from the array, then the length of array becomes shorter
                  # 2: linear interpolation of NaN values
                  # otherwise: leave NaN values as is
# --- Option to make a plot --- #
fig_flag = 1      # 1: yes, otherwise: not
# --- Choice of 4 variables to plot --- #
# 1: temp; 2: sal; 3: cond; 4: EC25; 5: density
# 6: sigmaT; 7: Chl-Flu; 8: Chl-a; 9: turb
# 10: pH; 11: DO (%); 12: DO (mg/L); 13: PAR
fig_var_id      = [1, 2, 8, 12]
xaxis_autoscale = [1, 1, 1, 1]  # Apply auto-scale for x-axes of the four variables; 1: yes, otherwise: no
yaxis_autoscale = 1     # Apply auto-scale for y-axis; 1: yes, otherwise: no
# --- Figure adjustment parameters --- #
fig_size = [6.0, 7.0]   # Figure size [width, height]
fig_space = [0.4, 0.3]  # Spaces between the subplots [wspace, hspace]
data_linethck  = 1.2    # Data line thickness
ax_linethck    = 0.7    # Axis line thickness
fontsize_ax    = 9      # Axis ticks font size
fontsize_label = 12     # Axis label font size
marker_option  = 'on'   # Include markers at data points; on: yes, otherwise: no
marker_sybol   = 'o'    # Marker symbol; e.g. o, ^, v, s, x
marker_size    = 5      # Marker size
# --- Y-axis (depth) range limit and tickes interval (when autoscale is off) --- #
dep_lim     = [0, 3]
dep_dy      = 1
# --- X-axis range limit for each variable (when autoscale is off) --- #
temp_lim    = [26.0, 31.0]
sal_lim     = [30, 34]
cond_lim    = [0, 100]
EC25_lim    = [40000, 50000]
density_lim = [1000, 1030]
sigmaT_lim  = [0, 1030]
Chl_flu_lim = [0, 50]
Chl_a_lim   = [0, 50]
turb_lim    = [0, 50]
ph_lim      = [7, 8.5]
DO_sat_lim  = [0, 100]
DO_lim      = [0, 8]
PAR_lim     = [0, 1000]
# --- X-axis ticks interval for each variable (when autoscale is off) --- #
temp_dx    = 2
sal_dx     = 2
cond_dx    = 20
EC25_dx    = 10000
density_dx = 2
sigmaT_dx  = 2
Chl_flu_dx = 10
Chl_a_dx   = 10
turb_dx    = 10
ph_dx      = 0.2
DO_sat_dx  = 10
DO_dx      = 2
PAR_dx     = 100
# ------------------

# ======================================================================

# Read AAQ data
with open(core_name+'.csv', 'r') as fp:
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

# Depth data
dep = data_array[:,0]

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

# Number of AAQ water quality variables
nvar = data_array.shape[1]-1
for i in range(0,nvar):
    [dep_bin, u_bin] = depth_bin(data_array[:,i+1], dep, dz, i_opt)
    if i == 0:
        bin_array = np.vstack((dep_bin, u_bin))
    else:
        bin_array = np.vstack((bin_array, u_bin))
bin_array = np.transpose(bin_array)

# Write binned data to csv file
if write_flag == 1:
    file_to_write = core_name+'_'+append_name+'.csv'
    header = ['Depth [m]','Temp. [degC]','Sal.','Cond. [mS/cm]','EC25 [µS/cm]'
              , 'Density [kg/m3]', 'SigmaT', 'Chl-Flu. [ppb]', 'Chl-a [µg/l]'
              , 'Turb. [FTU]', 'pH', 'DO [%]', 'DO [mg/l]'
              , 'Quant. [µmol/(m2*s)]', 'Mark']
    with open(file_to_write, 'w', newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(bin_array)

def fig_param(var_id):

    """
    Parameters for making a plot for a specified variable

    Arguments:
    ----------
    var_id: integer
        Index of a variable to plot

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
        x_lim = temp_lim
        dx    = temp_dx
        xlabel = 'Temperature ($^\circ$C)'
    elif var_id == 2:
        x_lim = sal_lim
        dx    = sal_dx
        xlabel = 'Salinity (psu)'
    elif var_id == 3:
        x_lim = cond_lim
        dx    = cond_dx
        xlabel = 'Conductivity (ms/cm)'
    elif var_id == 4:
        x_lim = EC25_lim
        dx    = EC25_dx
        xlabel = 'EC25 (\u03BC'+'s/cm)'
    elif var_id == 5:
        x_lim = density_lim
        dx    = density_dx
        xlabel = 'Water density (kg/m\u00b3)'
    elif var_id == 6:
        x_lim = sigmaT_lim
        dx    = sigmaT_dx
        xlabel = '\u03C3T'
    elif var_id == 7:
        x_lim = Chl_flu_lim
        dx    = Chl_flu_dx
        xlabel = 'Chl-Flu (ppb)'
    elif var_id == 8:
        x_lim = Chl_a_lim
        dx    = Chl_a_dx
        xlabel = 'Chl-\u03B1 (\u03BC'+'g/L)'
    elif var_id == 9:
        x_lim = turb_lim
        dx    = turb_dx
        xlabel = 'Turbidity (FTU)'
    elif var_id == 10:
        x_lim = ph_lim
        dx    = ph_dx
        xlabel = 'pH'
    elif var_id == 11:
        x_lim = DO_sat_lim
        dx    = DO_sat_dx
        xlabel = 'DO (%)'
    elif var_id == 12:
        x_lim = DO_lim
        dx    = DO_dx
        xlabel = 'DO (mg/L)'
    elif var_id == 13:
        x_lim = PAR_lim
        dx    = PAR_dx
        xlabel = 'PAR (\u03BC'+'mol '+"$\mathregular{m^{-2} s^{-1}}$"+')'
        
    x_ticks = np.arange(x_lim[0], x_lim[1]+0.5*dx, dx)

    return bin_data, x_lim, x_ticks, xlabel

# Make plots
plt.close('all')
if fig_flag == 1:
    # Make plots in a seperate window
    %matplotlib qt
    
    fig = plt.figure(figsize=(fig_size))
    
    axes= fig.subplots(2, 2)

    plt.subplots_adjust(wspace=fig_space[0], hspace=fig_space[1])
    
    for i in range(0,4):
        [bin_data, x_lim, x_ticks, xlabel] = fig_param(fig_var_id[i])
        if i == 0:
            ax = axes[0][0]
        elif i == 1:
            ax = axes[0][1]
        elif i == 2:
            ax = axes[1][0]
        elif i == 3:
            ax = axes[1][1]

        if marker_option == 'on':
            ax.plot(bin_data, dep_bin, marker=marker_sybol, linestyle='-',
                    linewidth=data_linethck, markersize=marker_size)
        else:
            ax.plot(bin_data, dep_bin, linewidth=data_linethck)
                
        if xaxis_autoscale[i] == 0:
            ax.set_xlim(x_lim[0],x_lim[1])
            ax.set_xticks(x_ticks) 
            ax.set_xticklabels(x_ticks,fontsize=fontsize_ax)
        else:
            xticklabels = ax.get_xticklabels()
            xticks = ax.get_xticks()
            ax.set_xticks(xticks) 
            ax.set_xticklabels(xticklabels,fontsize=fontsize_ax)
            
        if yaxis_autoscale == 0:
            dep_ticks = np.arange(dep_lim[0], dep_lim[1]+0.5*dep_dy, dep_dy)    
            ax.set_ylim(dep_lim[0],dep_lim[1])
            ax.set_yticks(dep_ticks) 
            ax.set_yticklabels(dep_ticks,fontsize=fontsize_ax)
        else:
            yticklabels = ax.get_yticklabels()
            #yticks = ax.get_yticks()
            #ax.set_yticks(yticks) 
            ax.set_yticklabels(yticklabels,fontsize=fontsize_ax)
        
        ax.invert_yaxis()
        ax.set_xlabel(xlabel,fontsize=fontsize_label)
        ax.set_ylabel('Depth (m)',fontsize=fontsize_label)
        
        # change all spines
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(ax_linethck)
        # increase tick width
        ax.tick_params(width=ax_linethck)
    
    plt.show()
