# -*- coding: utf-8 -*-
"""
This program calculates the distances where the relative footprint contributions
(50%, 60%, 70%, 80%) are obtained for the Eddypro full output data. The program
uses Kljun et al. (2015) and Kormann & Meixner (2001) methods to estimate
the footprint contributions.

This program needs to be configured before running it, see the parameters and
more detailed descriptions below.

#References

Kljun2015:
Kljun, N., Calanca, P., Rotach, M. W., & Schmid, H. P. (2015). A simple 
two-dimensional parameterisation for Flux Footprint Prediction (FFP). 
Geoscientific Model Development, 8(11), 3695–3713. 
https://doi.org/10.5194/gmd-8-3695-2015

KM2001:
Kormann, R., & Meixner, F. X. (2001). An analytical footprint model for 
non-neutral stratification. Boundary-Layer Meteorology, 99(2), 207–224. 
https://doi.org/10.1023/A:1018991015119


@author: Mika Korkiakoski (mika.korkiakoski@fmi.fi)
"""

import os
import sys
import warnings
import time
import math
import numpy as np
import pandas as pd

#Location of the eddypro full output datafile
DATA_LOC = '/home/m/Desktop/eddypro_full_output_adv.csv'

#Location of the sector-wise d and z0 file.
#Set to None if no sector-wise values are used.
#If sector-wise d and z0 are used, they with the aerodynamic measurement height
#and recalculated stability parameters will be saved into a separate file.
D_Z0_SECTORS_LOC = None

#Folder to save the footprints
SAVE_LOC = '/home/m/Desktop/'

#Measurement height
MEAS_HEIGHT = 2.5

#Canopy height
#Either canopy height or displacement height is required. Canopy height is used
#to calculate the displacement height, if the displacement height is not given.
#If DISP_HEIGHT is given, CANOPY_HEIGHT is ignored.
#Use the same canopy height as in Eddypro!!!
CANOPY_HEIGHT = 0.3

#Displacement height
#Set to None if not known. Don't use zero.
#If DISP_HEIGHT is None, the displacement height will be calculated from the
#CANOPY_HEIGHT.
#Use the same displacement height as in Eddypro!!!
DISP_HEIGHT = None

#Roughness length
#Set to None if not known. Don't use zero.
#If Z0 is None, it is not used to calculate footprints. Instead, alternative
#methods are used.
#Use the same roughness length as in Eddypro!!!
Z0 = None

#Latitude [deg] of the measurements, needed for Kljun et al. 2015
#Needs to be 0 <= lat <= 90
LAT = 60

#Calculate footprints according to Kljun et al. 2015
#Must be True or False
DO_KLJUN15 = True

#Calculate footprints according to Kormann & Meixner, 2001
#Must be True or False
#KM2001 gives generally larger footprints and is somewhat slower than Kljun15.
DO_KM01 = False



#########################################################################################
VERSION = 'v1.2.1 JAN 2025'
APPNAME = 'Footprints_Kljun15_KM2001'

#Ignore warnings. I know what I'm doing.
warnings.simplefilter('ignore')

#Checks if the directory exists and creates one if it doesn't
def check_create_dir(path):
    """Checks if the directory exists and creates one if it doesn't.
        
    Parameters:
        path (str): Path to the folder that is created. 
    """
    if not os.path.exists(path):
        os.makedirs(path)

def check_params(MEAS_HEIGHT, DISP_HEIGHT, CANOPY_HEIGHT, LAT, SAVE_LOC,
                 DO_KLJUN15, DO_KM01):
    """Checks that the input parameters have valid values. If invalid values are
    encountered, the program is closed with an appropriate error message.
        
    Parameters:
        MEAS_HEIGHT (float): Measurement height.
        DISP_HEIGHT (float): Displacement height.
        CANOPY_HEIGHT (float): Canopy height.
        LAT (float): Latitude.
        SAVE_LOC (str): path to the folder where the results are saved.
        DO_KLJUN15 (bool): Tells if the Kljun15 footprints are calculated.
        DO_KM01 (bool): Tells if the KM01 footprints are calculated.
    """
    #Check the SAVE_LOC parameter, if the final character is not "/" or "\", add it
    if SAVE_LOC[-1] != '/':
        if SAVE_LOC[-1] != '\\':
            #Check the first character, if it is "/", assume linux system and add "/"
            if SAVE_LOC[0] == '/':
                SAVE_LOC = SAVE_LOC+'/'
            #Otherwise, it is a windows system, add "\"
            else:
                SAVE_LOC = SAVE_LOC+'\\'

    #Create the folder for the results
    check_create_dir(SAVE_LOC)

    #If the displacement height is given
    if DISP_HEIGHT is not None:
        #Calculate the measurement height above the displacement height
        zm = MEAS_HEIGHT - DISP_HEIGHT

        if zm <= 0:
            sys.exit('zm (measurement height above the displacement height) must ' +
                     'be larger than zero. Closing the program.')

        if DISP_HEIGHT==0:
            sys.exit('DISP_HEIGHT cannot be zero. Set a proper value or set to ' +
                     'None if not known. If the displacement height was set to ' +
                     'zero in Eddypro, set DISP_HEIGHT to None.')

    #If the displacement height is not given
    if DISP_HEIGHT is None:
        #Check that the canopy height is given
        if isinstance(CANOPY_HEIGHT, (int, float)) is False:
            sys.exit('CANOPY_HEIGHT is required if DISP_HEIGHT is not given. ' +
                     'CANOPY_HEIGHT has an invalid value. Closing the program.')

        #Check that the measurement height is not lower than the canopy height
        if MEAS_HEIGHT < CANOPY_HEIGHT:
            sys.exit('MEAS_HEIGHT cannot be smaller than CANOPY_HEIGHT. Closing the program.')

    #Check that either the displacement height or the canopy height is given
    if DISP_HEIGHT is not None and CANOPY_HEIGHT is not None:
        print('Warning: Both CANOPY_HEIGHT and DISP_HEIGHT are given. ' +
              'CANOPY_HEIGHT will be ignored.')

    #Check that do_kljun is either True or False
    if isinstance(DO_KLJUN15, bool) is False:
        sys.exit('DO_KLJUN15 must be either True or False. Closing the program.')

    #Check that DO_KM01 is either True or False
    if isinstance(DO_KM01, bool) is False:
        sys.exit('DO_KM01 must be either True or False. Closing the program.')

    #Check that latitude is zero or larger and 90 or smaller
    if (LAT < 0 or LAT > 90) and DO_KLJUN15 is True:
        sys.exit('Latitude has to be within 0 <= LAT <= 90. Closing the program.')

def load_d_z0_file(D_Z0_SECTORS_LOC):
    """Loads the file containing sector-wise displacement height (d) and
    roughness length (z0) parameters.
        
    Parameters:
        D_Z0_SECTORS_LOC (str): Path to the sector-wise d and z0 file. 
                
    Returns:
        d_z0_data (dataframe): Dataframe containing sector-wise d and z0 
        parameters.
    """
    try:
        d_z0_data = pd.read_excel(D_Z0_SECTORS_LOC)
    
    #If the file is not an excel file. Try loading it as a csv.
    except (FileNotFoundError, ValueError):
        try:
            d_z0_data = pd.read_csv(D_Z0_SECTORS_LOC)
        except: #If any error after this. Close the program
            if CANOPY_HEIGHT is not None or DISP_HEIGHT is not None:
                print("Sector-wise displacement height and z0 file was not found. " + 
                      "Using one displacement height for the whole site.")
                d_z0_data=np.nan
            
            else:
                sys.exit("Sector-wise displacement height and z0 file was not found. " +
                         "Also, CANOPY_HEIGHT and DISP_HEIGHT is not set! " +
                         "Closing the program.")
    
    return d_z0_data

def disp_z0_height_mngmnt(DISP_HEIGHT, CANOPY_HEIGHT, d_z0_data):
    """Calculate the displacement height if it is not given. Use the same
    method as Eddypro uses.
    For forest canopies, displacement height is estimated to vary between 0.6
    and 0.8 times the height of the canopy (Arya, 1998; Stull, 1988).
        
    Parameters:
        DISP_HEIGHT (float): Displacement height.
        CANOPY_HEIGHT (float): Canopy height.
        
                
    Returns:
        d_z0_data (DataFrane): Dataframe containing sector-wise displacement 
                               and z0 heights.
    """
    #If the file containing sector-wise d and z0 data has been loaded successfully
    if isinstance(d_z0_data, pd.DataFrame):
        return d_z0_data
    
    #If one DISP_HEIGHT value is given, use it for all the wind sectors.
    elif DISP_HEIGHT is not None:
        d_z0_data = pd.DataFrame({'Start_WD [deg]':np.arange(0, 360, 10),
                        'End_WD [deg]':np.arange(10, 370, 10), 'd [m]': DISP_HEIGHT, 
                        'z0 [m]': Z0}, index=[np.arange(0, 36)])
        return d_z0_data
        
    #If sector-wise d and z0 data does not exist and DISP_HEIGHT is not given
    elif DISP_HEIGHT is None:
        try:
            DISP_HEIGHT = 0.67 * CANOPY_HEIGHT
        except TypeError: #If canopy height does not exist, close the program
            raise TypeError('DISP_HEIGHT or CANOPY_HEIGHT is not given! Also, ' +
                            'no sector-wise d or z0 data was not found! One ' +
                            'of the three has to be given to run the program!')
            
        d_z0_data = pd.DataFrame({'Start_WD [deg]':np.arange(0, 360, 10),
                        'End_WD [deg]':np.arange(10, 370, 10), 'd [m]': DISP_HEIGHT, 
                        'z0 [m]': Z0}, index=[np.arange(0, 36)])
        return d_z0_data
    
    
def epro_data_load(DATA_LOC):
    """Loads and formats the eddypro full output data for easier handling by converting
    the time information to datetime and dropping unnecessary columns. Calculates
    30 min averages only to remove possible extra minutes from the timestamp.
        
    Parameters:
        DATA_LOC (str): Path to the Eddypro full output file. 
                
    Returns:
        data (dataframe): Eddypro full output data in pandas dataframe.
        data_units (series): Units of the data in pandas series.
    """
    #Load the Eddypro full output data file
    try:
        data = pd.read_csv(DATA_LOC, index_col=0, skiprows=1, low_memory=False)
    except FileNotFoundError:
        raise FileNotFoundError('Eddypro data file was not found from the ' +
                                'given path. Closing the program.')

    data_units = data.iloc[0,:] #Get units
    data = data.iloc[1:,:] #Remove units

    #Timestamps to datetime
    try:
        data_times = pd.Series(dtype=object)
        for k in range(len(data)):
            if len(data_times)==0:
                data_times = pd.Series(pd.to_datetime(data.date.iloc[k]
                                    + ' ' + data.time.iloc[k], format = '%Y-%m-%d %H:%M'))
            else:
                data_times = pd.concat([data_times, pd.Series(pd.to_datetime(data.date.iloc[k]
                                    + ' ' + data.time.iloc[k], format = '%Y-%m-%d %H:%M'))])
    except AttributeError:
        raise AttributeError('There is something wrong with the datafile. ' +
                             'Check that it is Eddypro full output file! Closing the program.')

    data.index = data_times #Replace index with datetimes

    data = data.drop(['date','time'], axis = 1) #Remove date and time columns
    
    #All data to float format
    cols = data.columns
    for col in cols:
        data[col] = pd.to_numeric(data[col], errors = 'coerce')
    
    #Transform -9999 to nan
    data[data == -9999] = np.nan

    #Make 30 min averaging to remove possible extra minutes in the timestamp
    data = data.resample('30min').mean()

    return data, data_units

def calc_zm_stabparam(data, data_cols, d_z0_data, SAVE_LOC):
    """Calculate height above the measurement height for all the data rows and also 
    recalculate the stability parameter for all the data rows.
        
    Parameters:
        data (dataframe): Eddypro full output data in pandas dataframe.
        data_units (series): Units of the data in pandas series.
        d_z0_data (dataframe): Dataframe containing sector-wise displacement 
                               and z0 heights.
        SAVE_LOC (str): The folder where the result file is saved.
                
    Returns:
        data (dataframe): Eddypro full output data in pandas dataframe with updated
                          height above the measurement height, stability parameter,
                          and z0 values.
    """
    winddir_bins = np.linspace(0, 360, 37) #Bin limits for wind dirs
    labels = np.arange(0,36) #Labels for the bins (index of d_z0_data)
    
    #Bind the data according to wind sectors
    wd_bins = pd.cut(data["wind_dir"], bins = winddir_bins, labels = labels)
    wd_bins = wd_bins.astype(float) #Bin names to float
    
    #If wind dir in the data was nan, set the label to zero. It does not matter
    #what the label actually is as flux is also nan.
    wd_bins[wd_bins.isnull()] = 0
    
    #Calculate the height below the measurement height and get d and z0 for all datarows
    zm = MEAS_HEIGHT - d_z0_data['d [m]'][wd_bins]
    d = d_z0_data['d [m]'][wd_bins]
    z0 = d_z0_data['z0 [m]'][wd_bins]
    
    #Add zm and z0 columns to the data dataframe
    zm.index = z0.index = d.index = data.index
    data['zm'] = zm
    data['d'] = d
    data['z0'] = z0
    
    #Add zm and z0 units to data_cols
    data_cols = pd.concat([data_cols, pd.Series(['[m]','[m]','[m]'], 
                                                index = ['zm', 'd', 'z0'])])
    
    #Recalculate the stability parameter for all datarows and add it to data
    #dataframe, replacing the old values.
    data['(z-d)/L'] = data.zm / data.L
    
    #If using sector-wise values, save the 30-min values of zm, d, z0 and (z-d)/L
    #into a separate file.
    if np.std(d_z0_data['d [m]']) != 0:
        #Get the variables of interest
        data_params = data[['zm','d','z0','(z-d)/L']]
        
        #Reset index
        data_params = data_params.reset_index()
    
        #Change the datetime column (previously the index) name
        data_params = data_params.rename(columns = {'index':'Datetime (Period end)'})
        data_params.to_csv(SAVE_LOC + 'new_zm_d_z0.csv', index = False)

    return data

def bounday_layer_height(Ls, ustars, t_covs, LAT, zLs, air_ts):
    """Calculates the boundary layer height according to Kljun et al. 2015,
    Appendix B. During stable stratification, the height is calculated directly.
    During unstable stratification, the height is calculated iteratively by
    calculating the change in height and adding it into the previous height.
        
    Parameters:
        Ls (series): Series containing the obukhov lengths.
        ustars (series): Series containing the friction velocities.
        t_covs (series): Series containing the temperature covariances.
        LAT (float): Latitude.
        zLs (series): Series containing the stability parameters.
        air_ts (series): Series containing air temperatures.
                
    Returns:
        hs (series): Series containing the calculated boundary layer heights
    """
    omg = 7.2921159e-5 #Angular velocity of Earth's rotation
    f = 2 * omg * np.sin(np.deg2rad(LAT)) #Coriolis parameter

    #If zeta larger than this, startification is considered neutral or stable
    neutral_limit = 0

    #Constants for Eq. B5
    A = 0.2
    B = 2.5
    C = 8
    gamma = 0.01 #Gradient of potential temperature above the convective boundary layer
    g = 9.81 #Gravity

    hs = pd.Series(dtype=float) #Save boundary layer heights here

    #Boundary layer height in neutral or stable stratification
    h_stab = (Ls / 3.8) * (-1 + np.sqrt(1 + 2.28 * (ustars / (f * Ls))))

    #The first index with neutral or stable stratification, because unstable
    #h can only be iterated, not calculated directly.
    ind_ok = h_stab.index[h_stab.index.get_loc(h_stab.first_valid_index())]
    for ind, L, zL, ustar, t_cov, air_t in zip(h_stab.index, Ls, zLs, ustars, t_covs, air_ts):
        if hs.isnull().all() and ind != ind_ok:
            if len(hs)==0:
                hs = pd.Series(np.nan, index = [ind])
            else:
                hs = pd.concat([hs, pd.Series(np.nan, index = [ind])])
            continue

        #If stable or neutral stratification, calculate h directly
        if zL >= neutral_limit:
            #Boundary layer height in stable and neutral stratification
            #Eq. B1
            h_cur = (L / 3.8) * (-1 + np.sqrt(1 + 2.28 * (ustar / (f * L))))
            if len(hs)==0:
                hs = pd.Series(h_cur, index = [ind])
            else:
                hs = pd.concat([hs, pd.Series(h_cur, index = [ind])])

        #In case of unstable stratification, solve h iteratively
        else:
            #Split Eq. B5 to three parts
            D = t_cov / gamma
            E = h_cur**2 / ((1 + 2 * A) * h_cur -2 * B * 0.4 * L)
            F = (C * ustar**2 * air_t) / (gamma * g * ((1 + A) *
                h_cur - B * 0.4 * L))

            #Combine the parts of the equation B5
            dh_dt = D / (E + F) * 1800 #Seconds to 30 min
            h_cur = h_cur + dh_dt #Add the change in height to the previous height
            hs = pd.concat([hs, pd.Series(h_cur, index = [ind])])

    return hs

def get_contribution_dist(r, zm, h, umean, ustar, z0, psi, use_z0):
    """Calculate the distance where the relative contribution of interest is
    obtained in Kljun15 for each 30 min period.
        
    Parameters:
        r (float): The relative contribution [0, 1] that is calculated.
        zm (float): The height above the displacement height.
        h (float): The planetary boundary layer height.
        umean (float): The mean wind speed.
        ustar (float): The friction velocity.
        z0 (float): The roughness length.
        psi (float): The non-dimensional wind shear
        use_z0 (boolean): Which equation is used to calculate the contribution
                
    Returns:
        xr (series): Series containing the distances of the corresponding relative
                     contribution.
    """
    vk = 0.4 #Von Karman

    #Footprint parameterization constants
    #Eq. 17 in Kljun2015
    c = 1.4622
    d = 0.1359
    
    if use_z0 == True:
        xr = (-c / np.log(r) + d) * (zm * (1 - zm / h)**-1) * (np.log(zm / z0) - psi)
   
    else:  #Eq. 25 in Kljun2015
        xr = (-c / np.log(r) + d) * (zm * (1 - zm / h)**-1) * umean / ustar * vk

    return xr

def kljun_2015(zLs, ustars, umeans, hs, zms, Ls, z0s):
    """Calculate the footprints according to Kljun et al. 2015. The
    footprints are calculated by using the mean wind speed. Roughness length
    is not used.
        
    Parameters:
        zLs (float): Series containing the stability parameters.
        ustars (series): Series containing the friction velocities.
        umeans (series): Series containing the mean wind speeds.
        hs (series): Series containing the planetary boundary layer heights.
        zms (series): Series containing the heights above the displacement height.
        Ls (series): Series containing the obukhov lengths
        z0s (series): Series containing the roughness lengths
        
    Returns:
        fps (dataframe): Dataframe containing the footprints calculated with
        different relative contributions for every 30 min period.
    """
    vk = 0.4 #Von Karman

    #If the boundary layer height is smaller than 10 m (not mentioned in the
    #Kuljun2015 paper, but it is used in the attached footprint script)
    #or the measurement height above the displacement height is larger than
    #the height of the entrainment layer (Eq. 27 in Kljun2015), set those
    #boundary layer heights to nan (=makes also resulting footprints nan) as
    #those cases break the assumptions of Kljun2015.
    hs[(hs < 10) | (zms > 0.8 * hs)] = np.nan

    #If ustar < 0.1 m s-1 or the stability parameter is < -15.5,
    #Kljun2015 is not valid, set those cases nan.
    #The ustar limit is not mentioned in Kljun2015 paper, but it is used in
    #the attached footprint script, so it is also applied here.
    ustars[(ustars < 0.1) | (zLs < -15.5)] = np.nan
    
    #If zL is maller than -15.5, set the corresponding boundary layer height
    #to nan (=resulting footprint is nan). Such a low zL value breaks the
    #assumptions of Kljun2015.
    hs[zLs < -15.5] = np.nan
    
    #Print progress
    total = len(zLs)
    p20 = int(total * 0.2)
    p40 = int(total * 0.4)
    p60 = int(total * 0.6)
    p80 = int(total * 0.8)
    line_nums = pd.Series(np.arange(len(zLs)), index=zLs.index)

    fps = pd.DataFrame(dtype=float)
    for ind, ustar, umean, h, L, zm, z0 in zip(ustars.index, ustars, umeans, 
                                               hs, Ls, zms, z0s):
        #Print progress
        print_progress(np.array([p20, p40, p60, p80]), line_nums[ind], 'Kljun2015')
        
        #Maximum footprint contribution
        #If z0 exists, use the format of the equation where it is used 
        #(Eq. 22 in Kljun2015)
        #Otherwise use Eq. 21 in Kljun2015 (based on ratio of umean and ustar)
        if z0 is not None and np.isnan(z0) == False and z0 > 0: #Eq 22
            #Calculate the non-dimensional wind shear, psi
            if L > 0: #Stable stratification
                psi = -5.3 * zm / L
            else: #Unstable stratification
                khi = (1 - 19 * zm / L)**0.25
                psi = (np.log((1 + khi**2) / 2) + 2 * np.log((1 + khi) / 2) - 2 * 
                np.arctan(khi) + np.pi/2)
            
            #Calculate the peak footprint distance
            xr_peak = 0.87 * zm / (1. - (zm / h)) * (np.log(zm / z0) - psi)
            use_z0 = True
        else: #Eq. 21
            xr_peak = 0.87 * zm / (1. - (zm / h)) * (umean / ustar * vk)
            psi = np.nan
            use_z0 = False

        #Footprint offset
        xr_offset = get_contribution_dist(0.01, zm, h, umean, ustar, z0, psi, use_z0)

        #Calculate the distances where the relative contributions are obtained
        xr_50 = get_contribution_dist(0.5, zm, h, umean, ustar, z0, psi, use_z0)
        xr_60 = get_contribution_dist(0.6, zm, h, umean, ustar, z0, psi, use_z0)
        xr_70 = get_contribution_dist(0.7, zm, h, umean, ustar, z0, psi, use_z0)
        xr_80 = get_contribution_dist(0.8, zm, h, umean, ustar, z0, psi, use_z0)

        temp = pd.DataFrame({'x_offset':xr_offset, 'x_peak':xr_peak,
                           'x_50%':xr_50, 'x_60%':xr_60, 'x_70%':xr_70,
                           'x_80%':xr_80}, index = [ind])
        fps = pd.concat([fps, temp])

    return fps

def print_progress(lims, ind, method):
    """Prints a progress message during K&M2001 calculations on every 20% coverage.
    
    Parameters:
        lims (np array): List of percentages of the array length when a progress
        message is printed
        ind (int): Current looping index.
    """
    if ind in lims:
        if ind == lims[0]:
            print ("Calculating " + method + " footprints. 20% done...")
        if ind == lims[1]:
            print ("Calculating " + method + " footprints. 40% done...")
        if ind == lims[2]:
            print ("Calculating " + method + " footprints. 60% done...")
        if ind == lims[3]:
            print ("Calculating " + method + " footprints. 80% done...")

def korm_meix(zLs, ustars, umeans, zms):
    """Calculate the footprints according to Kormann & Meixner (2001). The
    footprints are calculated by using the mean wind speed. Roughness length
    is not used.
        
    Parameters:
        zLs (float): Series containing the stability parameters.
        ustars (series): Series containing the friction velocities.
        umeans (series): Series containing the mean wind speeds.
        zm (float): The height above the displacement height.
        
    Returns:
        fps (dataframe): Dataframe containing the footprints calculated with
        different relative contributions for every 30 min period.
    """
    vk = 0.4 #Von Karman
    
    #If zL is larger than 3 or smaller than -3, set it to nan, making the resulting
    #footprint to nan.
    zLs[zLs.abs() > 3] = np.nan

    fps = pd.DataFrame(dtype = float) #Save footprints here.
    
    #Print progress
    total = len(zLs)
    p20 = int(total * 0.2)
    p40 = int(total * 0.4)
    p60 = int(total * 0.6)
    p80 = int(total * 0.8)
    line_nums = pd.Series(np.arange(len(zLs)), index=zLs.index)
    
    for ind, zL, ustar, umean, zm in zip(zLs.index, zLs, ustars, umeans, zms):
        #Print progress
        print_progress(np.array([p20, p40, p60, p80]), line_nums[ind], 'K&M2001')
        
        #Similarity relations (Paulson, 1970)
        #Eqs. 33, 34 and 35 in KM2001
        if zL > 0:
            phi_m = 1 + 5 * zL
            phi_c = 1 + 5 * zL
        else:
            phi_m = (1 - 16 * zL)**(-1/4)
            phi_c = (1 - 16 * zL)**(-1/2)

        #Intermediate parameters for K&M2001
        #Exponent of the diffusivity power law
        #Eq. 36 in KM2001
        if zL > 0:
            n = 1 / phi_m
        else:
            n = (1 - 24 * zL) / (1 - 16 * zL)

        #Proportionality constant of the diffusivity power law (Eqs. 11 and 32)
        #Eqs. 11 and 32 in KM2001
        kappa = (vk * ustar * zm / phi_c) / zm**n
        
        #Check if zL, ustar or kappa is zero. 
        #If any of them is, set footprints to nan.
        if np.isnan(zL) == True or np.isnan(ustar) == True or kappa == 0:
            temp = pd.DataFrame({'x_offset':np.nan, 'x_peak':np.nan,
                               'x_50%':np.nan, 'x_60%':np.nan, 'x_70%':np.nan,
                               'x_80%':np.nan}, index = [ind])
            fps = pd.concat([fps, temp])
            continue

        #exponent of the wind speed power law
        #Eq. 36 in KM2001
        m = ustar * phi_m / (vk * umean)

        #Proportionality constant of the wind speed power law
        #Eqs. 11 and 31 in KM2001
        UU = umean / zm**m

        #Shape factor
        #Table 1 in KM2001
        r = 2 + m - n

        #Constant
        #Table 1 in KM2001
        mmu = (1 + m) / r

        #The length scale
        #Eq. 19 in KM2001
        zeta = UU * zm**r / (r**2 * kappa)

        #Calculate the distances where the relative contributions are obtained
        do_offset = True
        do50 = True
        do60 = True
        do70 = True
        do80 = True
        int_foot = 0
        di = 1
        for i in range(10000):
            if i == 0:
                i = 1

            #Cross-wind integrated 1D function
            int_foot = int_foot + di * (zeta**mmu * np.exp(-zeta / (i * di)) /
                                        ((i * di)**(1 + mmu) * math.gamma(mmu)))

            if do_offset is True and int_foot > 0.01:
                foot_offset = i * di
                do_offset = False

            if do50 is True and int_foot > 0.5:
                fetch50 = i * di
                do50 = False

            if do60 is True and int_foot > 0.6:
                fetch60 = i * di
                do60 = False

            if do70 is True and int_foot > 0.7:
                fetch70 = i * di
                do70 = False

            if do80 is True and int_foot > 0.8:
                fetch80 = i * di
                do80 = False
                break

        #Peak location (meters from the tower)
        #Eq. 22 in KM2001
        foot_peak = zeta / (1 + mmu)

        temp = pd.DataFrame({'x_offset':foot_offset, 'x_peak':foot_peak,
                           'x_50%':fetch50, 'x_60%':fetch60, 'x_70%':fetch70,
                           'x_80%':fetch80}, index = [ind])
        fps = pd.concat([fps, temp])

    return fps

#Save the dataframes containing the footprint results
def save_files(fp_data, SAVE_LOC, fp_type):
    """Save the datafile as csv to the folder given by the user. The datafile
    contains the footprints calculated either with Kljun2015 or KM2001.
        
    Parameters:
        fp_data (dataframe): Dataframe that includes the footprints for different
                             relative contributions for either Kljun15 or KM2001.
        SAVE_LOC (str): The folder where the result file is saved.
        fp_type (str): String of the footprint method used. Used as a part of the
                       filename.
    """
    #Get start and end times
    st = fp_data.index[0].strftime('%Y-%m-%d')
    et = fp_data.index[-1].strftime('%Y-%m-%d')

    #Reset index
    fp_data = fp_data.reset_index()

    #Change the datetime column (previously the index) name
    fp_data = fp_data.rename(columns = {'index':'Datetime (Period end)'})
    fp_data.to_csv(SAVE_LOC + 'Footprints_' + fp_type + '_' + st + '-' + et +
                   '.csv', index=False)

def main(DISP_HEIGHT):
    """Main function of the EddyproFiltering program.s
    """
    print('Starting ' + APPNAME + ' ' + '(' + VERSION + ').\n')
    
    #Check that the given parameters make sense
    check_params(MEAS_HEIGHT, DISP_HEIGHT, CANOPY_HEIGHT, LAT, SAVE_LOC,
                     DO_KLJUN15, DO_KM01)
    
    #Load sector-wise displacement height and z0 datafile
    d_z0_data = load_d_z0_file(D_Z0_SECTORS_LOC)
    
    #Prepare displacement height and roughness length data for the program, 
    #depending on the input.
    d_z0_data = disp_z0_height_mngmnt(DISP_HEIGHT, CANOPY_HEIGHT, d_z0_data)

    #Load and format the Eddypro data
    data, data_cols = epro_data_load(DATA_LOC)

    print("\nEddypro datafile loaded and formatted.\n")
    
    #Calculate zm and recalculate stability parameters for all wind sectors
    #and save the results into a separate file.
    data = calc_zm_stabparam(data, data_cols, d_z0_data, SAVE_LOC)

    if DO_KLJUN15 is True:
        #Start timer
        start = time.time()

        #Calculate boundary layer heights for Kljun et al. 2015
        print("Calculating planetary boundary layer heights.")
        hs = bounday_layer_height(data.L, data['u*'], data['w/ts_cov'], LAT, 
                                  data['(z-d)/L'], data.air_temperature)

        #Calculate cross-wind integrated footprint according to Kljun et al. 2015
        print("Calculating Kljun et al. (2015) footprints.")
        fps_kljun = kljun_2015(data['(z-d)/L'], data['u*'], data.wind_speed, 
                               hs, data.zm, data.L, data.z0)

        #Save the footprints
        save_files(fps_kljun, SAVE_LOC, 'Kljun2015')

        #End timer
        end = time.time()
        elapsed_time="{:.2f}".format(end-start)

        print("Kljun et al. (2015) footprints saved.")
        print('Time elapsed: ' + elapsed_time + ' seconds.\n')

    else:
        fps_kljun = np.nan

    if DO_KM01 is True:
        #Start timer
        start = time.time()

        #Calculate cross-wind integrated footprint according to Kormann & Meixner 2001
        print("Calculating Kormann & Meixner (2001) footprints.")
        fps_km = korm_meix(data['(z-d)/L'], data['u*'], data.wind_speed, data.zm)

        #Save the footprints
        save_files(fps_km, SAVE_LOC, 'K&M2001')

        #End timer
        end = time.time()
        elapsed_time="{:.2f}".format(end-start)

        print("Kormann & Meixner (2001) footprints saved.")
        print('Time elapsed: ' + elapsed_time + ' seconds.\n')

    else:
        fps_km = np.nan

    print("Program finished!")

    return fps_kljun, fps_km

#If the file is not imported, start the program from function main()
if __name__ == "__main__":
    fps_kljun, fps_km = main(DISP_HEIGHT)
   
