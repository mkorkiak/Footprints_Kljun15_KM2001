# -*- coding: utf-8 -*-
"""
This program calculates the distances where the relative footprint contributions
(50%, 60Q%, 70%, 80%) are obtained for the Eddypro full output data. The program
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

#Folder to save the footprints
SAVE_LOC = '/home/m/Desktop/'

#Measurement height
MEAS_HEIGHT = 26

#Canopy height
#Either canopy height or displacement height is required. Canopy height is used
#to calculate the displacement height, if the displacement height is not given.
#If DISP_HEIGHT is given, CANOPY_HEIGHT is ignored.
#Use the same canopy height as in Eddypro!!!
CANOPY_HEIGHT = 21

#Displacement height
#Set to None if not known. Don't use zero.
#If DISP_HEIGHT is None, the displacement height will be calculated from the
#CANOPY_HEIGHT.
#Use the same displacement height as in Eddypro!!!
DISP_HEIGHT = 9

#Latitude of the measurements, needed for Kljun et al. 2015
#Needs to be 0 <= lat <= 90
LAT = 60

#Calculate footprints according to Kljun et al. 2015
#Must be True or False
DO_KLJUN15 = True

#Calculate footprints according to Kormann & Meixner, 2001
#Must be True or False
#KM2001 gives generally larger footprints and is 6-7 times slower than Kljun15.
DO_KM01 = True



#########################################################################################
VERSION = 'v1.0 MAR 2023'
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

def calc_DISP_HEIGHT(CANOPY_HEIGHT):
    """Calculate the displacement height if it is not given. Use the same
    method as Eddypro uses.
    For forest canopies, displacement height is estimated to vary between 0.6
    and 0.8 times the height of the canopy (Arya, 1998; Stull, 1988).
        
    Parameters:
        CANOPY_HEIGHT (float): Canopy height.
                
    Returns:
        DISP_HEIGHT (float): Displacement height.
    """
    DISP_HEIGHT = 0.67 * CANOPY_HEIGHT #Same as in Eddypro

    return DISP_HEIGHT

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
            data_times = pd.concat([data_times, pd.Series(pd.to_datetime(data.date[k]
                                    + ' ' + data.time[k], format='%Y-%m-%d %H:%M'))])
    except AttributeError:
        raise AttributeError('There is something wrong with the datafile. ' +
                             'Check that it is Eddypro full output file! ' +
                             'Closing the program.')

    data.index = data_times #Replace index with datetimes

    data = data.drop(['date','time'], axis = 1) #Remove date and time columns
    data = data.astype(float) #all data to float format

    #Make 30 min averaging to remove possible extra minutes in the timestamp
    data = data.resample('30min').mean()

    return data, data_units

def bounday_layer_height(Ls, ustars, t_covs, LAT, zm, air_ts):
    """Calculates the boundary layer height according to Kljun et al. 2015,
    Appendix B. During stable stratification, the height is calculated directly.
    During unstable stratification, the height is calculated iteratively by
    calculating the change in height and adding it into the previous height.
        
    Parameters:
        Ls (series): Series containing the obukhov lengths.
        ustars (series): Series containing the friction velocities.
        t_covs (series): Series containing the temperature covariances.
        LAT (float): Latitude.
        zm (float): The height above the displacement height.
        air_ts (series): Series containing air temperatures.
                
    Returns:
        hs (series): Series containing the calculated boundary layer heights
    """
    omg = 7.2921159e-5 #Angular velocity of Earth's rotation
    f = -2 * omg * np.sin(LAT) #Coriolis parameter
    zLs = zm / Ls #Stability parameter

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
        if len(hs)==0 and ind != ind_ok:
            hs = pd.concat([hs, pd.Series(np.nan, index = [ind])])
            continue

        #If stable or neutral stratification, calculate h directly
        if zL > neutral_limit:
            #Boundary layer height in stable and neutral stratification
            #Eq. B1
            h_cur = (L / 3.8) * (-1 + np.sqrt(1 + 2.28 * (ustar / (f * L))))
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

def get_contribution_dist(r, zm, h, umean, ustar):
    """Calculate the distance where the relative contribution of interest is
    obtained in Kljun15 for each 30 min period.
        
    Parameters:
        r (float): The relative contribution [0, 1] that is calculated.
        zm (float): The height above the displacement height.
        h (series): Series containing the planetary boundary layer heights.
        umean (series): Series containing the mean wind speeds.
        ustar (series): Series containing the friction velocities.
                
    Returns:
        xr (series): Series containing the distances of the corresponding relative
                     contribution.
    """
    vk = 0.4 #Von Karman

    #Footprint parameterization constants
    #Eq. 17 in Kljun2015
    c = 1.4622
    d = 0.1359

    #Eq. 25 in Kljun2015
    xr = (-c / np.log(r) + d) * (zm * (1 - zm / h)**-1) * umean / ustar * vk

    return xr

def kljun_2015(zLs, ustars, umeans, hs, zm):
    """Calculate the footprints according to Kljun et al. 2015. The
    footprints are calculated by using the mean wind speed. Roughness length
    is not used.
        
    Parameters:
        zLs (float): Series containing the stability parameters.
        ustars (series): Series containing the friction velocities.
        umeans (series): Series containing the mean wind speeds.
        hs (series): Series containing the planetary boundary layer heights.
        zm (float): The height above the displacement height.
        
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
    hs[(hs < 10) | (zm > 0.8 * hs)] = np.nan

    #If ustar < 0.1 m s-1 or the stability parameter is < -15.5,
    #Kljun2015 is not valid, set those cases nan.
    #The ustar limit is not mentioned in Kljun2015 paper, but it is used in
    #the attached footprint script, so it is also applied here.
    ustars[(ustars < 0.1) | (zLs < -15.5)] = np.nan
    
    #If zL is maller than -15.5, set the corresponding boundary layer height
    #to nan (=resulting footprint is nan). Such a low zL value breaks the
    #assumptions of Kljun2015.
    hs[zLs < -15.5] = np.nan

    fps = pd.DataFrame(dtype=float)
    for ind, ustar, umean, h in zip(ustars.index, ustars, umeans, hs):
        #Maximum footprint contribution
        #Eq. 22 in Kljun2015
        xr_peak = 0.87 * zm / (1. - (zm / h)) * (umean / ustar * vk)

        #Footprint offset
        xr_offset = get_contribution_dist(0.01, zm, h, umean, ustar)

        #Calculate the distances where the relative contributions are obtained
        xr_50 = get_contribution_dist(0.5, zm, h, umean, ustar)
        xr_60 = get_contribution_dist(0.6, zm, h, umean, ustar)
        xr_70 = get_contribution_dist(0.7, zm, h, umean, ustar)
        xr_80 = get_contribution_dist(0.8, zm, h, umean, ustar)

        temp = pd.DataFrame({'x_offset':xr_offset, 'x_peak':xr_peak,
                           'x_50%':xr_50, 'x_60%':xr_60, 'x_70%':xr_70,
                           'x_80%':xr_80}, index = [ind])
        fps = pd.concat([fps, temp])

    return fps

def print_progress_km(lims, ind):
    """Prints a progress message during K&M2001 calculations on every 20% coverage.
    
    Parameters:
        lims (np array): List of percentages of the array length when a progress
        message is printed
        ind (int): Current looping index.
    """
    if ind in lims:
        if ind == lims[0]:
            print ("Calculating KM2001 footprints. 20% done...")
        if ind == lims[1]:
            print ("Calculating KM2001 footprints. 40% done...")
        if ind == lims[2]:
            print ("Calculating KM2001 footprints. 60% done...")
        if ind == lims[3]:
            print ("Calculating KM2001 footprints. 80% done...")
            
def korm_meix(zLs, ustars, umeans, zm):
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
    p20 = np.int(total * 0.2)
    p40 = np.int(total * 0.4)
    p60 = np.int(total * 0.6)
    p80 = np.int(total * 0.8)
    line_nums = pd.Series(np.arange(len(zLs)),index=zLs.index)

    for ind, zL, ustar, umean in zip(zLs.index, zLs, ustars, umeans):
        #Print progress
        print_progress_km(np.array([p20, p40, p60, p80]), ind)
        
        #Similarity relations (Paulson, 1970)
        #Eqs. 33, 34 and 35 in KM2001
        if zL > 0:
            phi_m = 1 + 5 * zL
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
        
        #Check if zL is zero and if kappa is zero (due to ustar being zero). 
        #If either of them is, set footprints to nan.
        if np.isnan(zL) == True or kappa == 0:
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

    #Calculate the displacement height if it is not given
    if DISP_HEIGHT is None:
        DISP_HEIGHT = calc_DISP_HEIGHT(CANOPY_HEIGHT)

    #Load and format the Eddypro data
    data, _ = epro_data_load(DATA_LOC)

    print("\nEddypro datafile loaded and formatted.\n")

    #Measurement height above the displacement height
    zm = MEAS_HEIGHT - DISP_HEIGHT

    if DO_KLJUN15 is True:
        #Start timer
        start = time.time()

        #Calculate boundary layer heights for Kljun et al. 2015
        print("Calculating planetary boundary layer heights.")
        hs = bounday_layer_height(data.L, data['u*'], data['w/ts_cov'], LAT, zm,
                               data.air_temperature)

        #Calculate cross-wind integrated footprint according to Kljun et al. 2015
        print("Calculating Kljun et al. (2015) footprints.")
        fps_kljun = kljun_2015(data['(z-d)/L'], data['u*'], data.wind_speed, hs, zm)

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
        fps_km = korm_meix(data['(z-d)/L'], data['u*'], data.wind_speed, zm)

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
   
