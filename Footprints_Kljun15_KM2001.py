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

import os, sys, warnings, time, math
import numpy as np
import pandas as pd

#Location of the eddypro full output datafile
data_loc = '/home/m/Desktop/eddypro_full_output_adv.csv'

#Folder to save the footprints
save_loc = '/home/m/Desktop/'

#Measurement height
meas_height = 26

#Canopy height
#Either canopy height or displacement height is required. Canopy height is used
#to calculate the displacement height, if the displacement height is not given.
#If disp_height is given, canopy_height is ignored.
#Use the same canopy height as in Eddypro!!!
canopy_height = 21

#Displacement height
#Set to None if not known. Don't use zero.
#If disp_height is None, the displacement height will be calculated from the
#canopy_height.
#Use the same displacement height as in Eddypro!!!
disp_height = 9

#Latitude of the measurements, needed for Kljun et al. 2015
#Needs to be 0 <= lat <= 90
lat = 60

#Calculate footprints according to Kljun et al. 2015
#Must be True or False
do_kljun15 = True

#Calculate footprints according to Kormann & Meixner, 2001
#Must be True or False
#KM2001 gives larger footprints and is 6-7 times slower than Kljun15.
do_km = True



###################################################################################################################################################
VERSION = 'v1.0 MAR 2023'
APPNAME = 'Footprints_Kljun15_KM2001'

#Ignore warnings. I know what I'm doing.
warnings.simplefilter('ignore')

#Checks if the directory exists and creates one if it doesn't
def check_create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
        
    return

def check_params(meas_height, disp_height, canopy_height, lat, save_loc,
                 do_kljun15, do_km):
    #Check the save_loc parameter, if the final character is not "/" or "\", add it
    if save_loc[-1]!='/':
        if save_loc[-1]!='\\':
            #Check the first character, if it is "/", assume linux system and add "/"
            if save_loc[0]=='/':
                save_loc=save_loc+'/'
            #Otherwise, it is a windows system, add "\"
            else:
                save_loc=save_loc+'\\'
    
    #Create the folder for the results
    check_create_dir(save_loc)
    
    #If the displacement height is given
    if disp_height is not None:
        #Calculate the measurement height above the displacement height
        zm = meas_height - disp_height
        
        if zm <= 0:
            sys.exit('zm (measurement height above the displacement height) must ' +
                     'be larger than zero. Closing the program.')
            
        if disp_height==0:
            sys.exit('disp_height cannot be zero. Set a proper value or set to ' +
                     'None if not known. If the displacement height was set to ' +
                     'zero in Eddypro, set disp_height to None.')
    
    #If the displacement height is not given
    if disp_height is None:
        #Check that the canopy height is given
        if isinstance(canopy_height, (int, float))==False:
            sys.exit('canopy_height is required if disp_height is not given. ' + 
                     'canopy_height has an invalid value. Closing the program.')
        
        #Check that the measurement height is not lower than the canopy height
        if meas_height < canopy_height:
            sys.exit('meas_height cannot be smaller than canopy_height. Closing the program.')
    
    #Check that either the displacement height or the canopy height is given
    if disp_height is not None and canopy_height is not None:
        print('Warning: Both canopy_height and disp_height are given. ' +
              'canopy_height will be ignored.')
    
    #Check that do_kljun is either True or False
    if isinstance(do_kljun15, bool)==False:
        sys.exit('do_kljun15 must be either True or False. Closing the program.')
    
    #Check that do_km is either True or False
    if isinstance(do_km, bool)==False:
        sys.exit('do_km must be either True or False. Closing the program.')
    
    #Check that latitude is zero or larger and 90 or smaller
    if (lat < 0 or lat > 90) and do_kljun15 == True:
        sys.exit('Latitude has to be within 0 <= lat <= 90. Closing the program.')
    
    return

#Calculate the displacement height if it is not given
#For forest canopies, it is estimated to vary between 0.6 and 0.8 times
#the height of the canopy (Arya, 1998; Stull, 1988)
def calc_disp_height(canopy_height):
    disp_height = 0.67 * canopy_height #Same as in Eddypro
    
    return disp_height

#Eddypro data formatting
def epro_data_load(data_loc):
    #Load the Eddypro full output data file
    try:
        data = pd.read_csv(data_loc, index_col=0, skiprows=1, low_memory=False)
    except FileNotFoundError:
        raise FileNotFoundError('Eddypro data file was not found from the given path. Closing the program.')
    
    data_units = data.iloc[0,:] #Get units
    data = data.iloc[1:,:] #Remove units
    
    #Timestamps to datetime
    try:
        data_times = pd.Series(dtype=object)
        for k in range(len(data)):
            data_times = pd.concat([data_times, pd.Series(pd.to_datetime(data.date[k]
                                    + ' ' + data.time[k], format='%Y-%m-%d %H:%M'))])
    except AttributeError:
        raise AttributeError('There is something wrong with the datafile. Check that it is Eddypro '+
                 'full output file! Closing the program.')
        
    data.index = data_times #Replace index with datetimes
    
    data = data.drop(['date','time'], axis = 1) #Remove date and time columns
    data = data.astype(float) #all data to float format
    
    #Make 30 min averaging to remove possible extra minutes in the timestamp
    data = data.resample('30min').mean()
    
    return data, data_units

#Calculate the boundary layer height according to Kljun et al. 2015, Appendix B
def bounday_layer_height(L, ustar, t_cov, lat, zm, air_t):
    omg = 7.2921159e-5 #Angular velocity of Earth's rotation
    f = -2 * omg * np.sin(lat) #Coriolis parameter
    zL = zm / L #Stability parameter
    
    #If zeta larger than this, startification is considered neutral or stable
    neutral_limit = 0 
    
    #Constants for Eq. B5
    A = 0.2
    B = 2.5
    C = 8
    gamma = 0.01 #Gradient of potential temperature above the convective boundary layer
    g = 9.81 #Gravity
    
    h = pd.Series(dtype=float) #Save boundary layer heights here
    
    #Boundary layer height in neutral or stable stratification
    h_stab = (L / 3.8) * (-1 + np.sqrt(1 + 2.28 * (ustar / (f * L))))
    
    #The first index with neutral or stable stratification, because unstable
    #h can only be iterated, not calculated directly.
    ind = h_stab.index.get_loc(h_stab.first_valid_index())
    for k in range(len(L)):
        if k < ind:
            h = pd.concat([h, pd.Series(np.nan, index=[L.index[k]])])
            continue
            
        #If stable or neutral stratification, calculate h directly
        if zL[k] > neutral_limit:      
            #Boundary layer height in stable and neutral stratification
            #Eq. B1
            h_cur = (L[k] / 3.8) * (-1 + np.sqrt(1 + 2.28 * (ustar[k] / (f * L[k]))))
            h = pd.concat([h, pd.Series(h_cur, index=[L.index[k]])])
        
        #In case of unstable stratification, solve h iteratively
        else:     
            #Split Eq. B5 to three parts
            D = t_cov[k] / gamma
            E = h_cur**2 / ((1 + 2 * A) * h_cur -2 * B * 0.4 * L[k])
            F = (C * ustar[k]**2 * air_t[k]) / (gamma * g * ((1 + A) * 
                h_cur - B * 0.4 * L[k]))
            
            #Combine the parts of the equation B5
            dh_dt = D / (E + F) * 1800 #Seconds to 30 min
            h_cur = h_cur + dh_dt #Add the change in height to the previous height
            h = pd.concat([h, pd.Series(h_cur, index=[L.index[k]])])
            
    return h

#Calculate the distance where the contribution of interest is obtained in Kljun2015
def get_contribution_dist(r, zm, h, umean, ustar):
    vk = 0.4 #Von Karman
    
    #Footprint parameterization constants
    #Eq. 17 in Kljun2015
    c = 1.4622
    d = 0.1359
    
    #Eq. 25 in Kljun2015
    xr = (-c / np.log(r) + d) * (zm * (1 - zm / h)**-1) * umean / ustar * vk
    
    return xr

#Footprints according to Kljun et al. 2015
def kljun_2015(zL, ustar, wind_speed, h, zm):
    vk=0.4 #Von Karman
    
    #If the boundary layer height is smaller than 10 m (not mentioned in the 
    #Kuljun2015 paper, but it is used in the attached footprint script) 
    #or the measurement height above the displacement height is larger than 
    #the height of the entrainment layer (Eq. 27 in Kljun2015), set those 
    #boundary layer heights to nan (=makes also resulting footprints nan) as 
    #those cases break the assumptions of Kljun2015.
    h[(h < 10) | (zm > 0.8*h)] = np.nan
    
    #If ustar < 0.1 m s-1 or the stability parameter is < -15.5, 
    #Kljun2015 is not valid, set those cases nan.
    #The ustar limit is not mentioned in Kljun2015 paper, but it is used in
    #the attached footprint script, so it is also applied here.
    ustar[(ustar < 0.1) | (zL < -15.5)]=np.nan
    
    fps = pd.DataFrame(dtype=float)
    for k in range(len(ustar)):
        #Maximum footprint contribution 
        #Eq. 22 in Kljun2015
        xr_peak = 0.87 * zm / (1. - (zm / h)) * (wind_speed[k] / ustar[k] * vk)
        
        #Footprint offset
        xr_offset = get_contribution_dist(0.01, zm, h[k], wind_speed[k], ustar[k])
        
        #Calculate the distances where the relative contributions are obtained
        xr_50 = get_contribution_dist(0.5, zm, h[k], wind_speed[k], ustar[k])
        xr_60 = get_contribution_dist(0.6, zm, h[k], wind_speed[k], ustar[k])
        xr_70 = get_contribution_dist(0.7, zm, h[k], wind_speed[k], ustar[k])
        xr_80 = get_contribution_dist(0.8, zm, h[k], wind_speed[k], ustar[k])
        
        temp = pd.DataFrame({'x_offset':xr_offset, 'x_peak':xr_peak, 
                           'x_50%':xr_50, 'x_60%':xr_60, 'x_70%':xr_70,
                           'x_80%':xr_80}, index=[ustar.index[k]])
        fps = pd.concat([fps, temp])
        
    return fps

#Footprints according to Kormann & Meixner, 2001
def korm_meix(zL, ustar, wind_speed, zm, disp_height):
    vk = 0.4 #Von Karman
    
    fps = pd.DataFrame(dtype=float)
    for k in range(len(zL)):
        #Similarity relations (Paulson, 1970)
        #Eqs. 33, 34 and 35 in KM2001
        if zL[k] > 0:
            phi_m = 1 + 5 * zL[k]
            phi_c = phi_m
            #psi_m = -5 * data['(z-d)/L'][k]
        else:
            phi_m = (1 - 16 * zL[k])**(-1/4)
            phi_c = (1 - 16 * zL[k])**(-1/2)
            #eta = (1 - 16 * data['(z-d)/L'][k])**(1/4)
            #psi_m = -2*np.log((1 + eta)/2)-np.log((1 + eta**2)/2)
            #+2*np.arctan(eta)-np.pi/2
    
        #Intermediate parameters for K&M2001
        #Exponent of the diffusivity power law
        #Eq. 36 in KM2001
        if zL[k] > 0:
            n = 1 / phi_m
        else:
            n = (1 - 24 * zL[k]) / (1 - 16 * zL[k])
 
        #Proportionality constant of the diffusivity power law (Eqs. 11 and 32)
        #Eqs. 11 and 32 in KM2001
        kappa = (vk * ustar[k] * zm / phi_c) / zm**n

        #exponent of the wind speed power law
        #Eq. 36 in KM2001
        m = ustar[k] * phi_m / (vk * wind_speed[k])

        #Proportionality constant of the wind speed power law 
        #Eqs. 11 and 31 in KM2001
        UU = wind_speed[k] / zm**m

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
            
            if do_offset == True and int_foot > 0.01:
               foot_offset = i * di
               do_offset = False
    
            if do50 == True and int_foot > 0.5:
               fetch50 = i * di
               do50 = False
               
            if do60 == True and int_foot > 0.6:
               fetch60 = i * di
               do60 = False

            if do70 == True and int_foot > 0.7:
               fetch70 = i * di
               do70 = False

            if do80 == True and int_foot > 0.8:
               fetch80 = i * di
               do80 = False
               break

        #Peak location (meters from the tower)
        #Eq. 22 in KM2001
        foot_peak = zeta / (1 + mmu)
        
        temp = pd.DataFrame({'x_offset':foot_offset, 'x_peak':foot_peak, 
                           'x_50%':fetch50, 'x_60%':fetch60, 'x_70%':fetch70,
                           'x_80%':fetch80}, index=[zL.index[k]])
        fps = pd.concat([fps, temp])
        
    return fps

#Save the dataframes containing the footprint results
def save_files(fp_data, save_loc, fp_type):
    #Get start and end times
    st = fp_data.index[0].strftime('%Y-%m-%d')
    et = fp_data.index[-1].strftime('%Y-%m-%d')
    
    #Reset index
    fp_data = fp_data.reset_index()
    
    #Change the datetime column (previously the index) name
    fp_data = fp_data.rename(columns = {'index':'Datetime (Period end)'})
    fp_data.to_csv(save_loc + 'Footprints_' + fp_type + '_' + st + '-' + et + 
                   '.csv', index=False)
    
    return

def main(disp_height):
    print('Starting '+APPNAME+' '+'('+VERSION+').\n')
    
    #Calculate the displacement height if it is not given
    if disp_height<=0:
        disp_height = calc_disp_height(canopy_height)
    
    #Check that the given parameters make sense
    check_params(meas_height, disp_height, canopy_height, lat, save_loc,
                     do_kljun15, do_km)
    
    #Load and format the Eddypro data
    data, data_units = epro_data_load(data_loc)
    
    print("\nEddypro datafile loaded and formatted.\n")
    
    #Measurement height above the displacement height
    zm = meas_height - disp_height
    
    if do_kljun15 == True:
        #Start timer
        start = time.time()
        
        #Calculate boundary layer heights for Kljun et al. 2015
        print("Calculating planetary boundary layer heights.")
        h = bounday_layer_height(data.L, data['u*'], data['w/ts_cov'], lat, zm,
                               data.air_temperature)
        
        #Calculate cross-wind integrated footprint according to Kljun et al. 2015
        print("Calculating Kljun et al. (2015) footprints.")
        fps_kljun = kljun_2015(data['(z-d)/L'], data['u*'], data.wind_speed, h, zm)
        
        #Save the footprints
        save_files(fps_kljun, save_loc, 'Kljun2015')
        
        #End timer
        end = time.time()
        elapsed_time="{:.2f}".format(end-start)
        
        print("Kljun et al. (2015) footprints saved.")
        print('Time elapsed: '+elapsed_time+' seconds.\n')
        
    else:
        fps_kljun = np.nan
    
    if do_km == True:
        #Start timer
        start = time.time()
        
        #Calculate cross-wind integrated footprint according to Kormann & Meixner 2001
        print("Calculating Kormann & Meixner (2001) footprints.")
        fps_km = korm_meix(data['(z-d)/L'], data['u*'], data.wind_speed, zm, disp_height)
        
        #Save the footprints
        save_files(fps_km, save_loc, 'K&M2001')
        
        #End timer
        end = time.time()
        elapsed_time="{:.2f}".format(end-start)
        
        print("Kormann & Meixner (2001) footprints saved.")
        print('Time elapsed: '+elapsed_time+' seconds.\n')
        
    else:
        fps_km = np.nan
        
    print("Program finished!")
    
    return fps_kljun, fps_km

#If the file is not imported, start the program from function main()
if __name__ == "__main__":
   fps_kljun, fps_km = main(disp_height)

"""
epro_data=pd.read_csv('C:/Users/korkiak/OneDrive - Ilmatieteen laitos/Desktop/Eddypro demo/Lettosuo/eddypro_letto_fp_korr_meix_full_output_2023-03-16T095027_adv.csv',index_col=0,skiprows=1)
epro_data, _ = epro_data_mod(epro_data)

plt.figure()
plt.plot(fps_kljun['x_70%'])
plt.plot(fps_km['x_70%'])
plt.plot(epro_data['x_70%'])

plt.figure()
plt.plot(fps_kljun.x_peak)
plt.plot(fps_km.x_peak)
plt.plot(epro_data.x_peak)
"""