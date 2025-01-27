**Footprints_Kljun15_KM2001**
=
This program calculates the distances where the relative footprint contributions
(50%, 60%, 70%, 80%) are obtained for the Eddypro full output data. The program
uses Kljun et al. (2015) and Kormann & Meixner (2001) methods to estimate
the footprint contributions.

This program needs to be configured before running it, see the parameters and
more detailed descriptions below.

Kormann & Meixner (2001) gives larger footprints and is somewhat slower than 
Kljun et al. (2015) in this program.

No fancy libraries needed. Check that Numpy and Pandas are installed.

Author: Mika Korkiakoski

Bug reports: mika.korkiakoski@fmi.fi

Parameters
=

**DATA_LOC**  
Location of the eddypro full output datafile.  
The location needs to be in hyphens.  
The folders can be separated by "/" (Linux and windows) or by backslash (Windows), but if the folders are separated by backslash, one needs to add "r" before the first hyphen like this: data_loc = r'C:\m\Desktop\results'.

**D_Z0_SECTORS_LOC**  
Location of the sector-wise d and z0 file.  
Set to None if no sector-wise values are used.  
If sector-wise d and z0 are used, they with the aerodynamic measurement height and recalculated stability parameters will be saved into a separate file.

**SAVE_LOC**  
Folder to save the footprints. Same rules for path formatting as for data_loc.

**MEAS_HEIGHT**  
Measurement height in meters. Typically means the sonic anemometer height.   
This variable is mandatory.

**CANOPY_HEIGHT**  
Canopy height.  
Either canopy height or displacement height is required. Canopy height is used to calculate the displacement height, if the displacement height is not given. If disp_height is given, canopy_height is ignored and a warning is raised.  
Use the same canopy height as in Eddypro!!!

**DISP_HEIGHT**  
Displacement height.  
Set to None if not known. Don't use zero.  
If disp_height is None, the displacement height will be calculated from the canopy_height similarly as in Eddypro.  
Use the same displacement height as in Eddypro!!!  

**Z0**  
Roughness length  
Set to None if not known. Don't use zero.  
If Z0 is None, it is not used to calculate footprints. Instead, alternative methods are used. See the references for more details.
Use the same roughness length as in Eddypro!!!  

**LAT**  
Latitude of the measurements  .
Needed for Kljun et al. 2015 method for boundary layer height calculation.  
Needs to be 0 <= lat <= 90.

**DO_KLJUN15**  
Calculate footprints according to Kljun et al. (2015).  
Must be True or False.

**DO_KM01**  
Calculate footprints according to Kormann & Meixner (2001).  
Must be True or False.


**References**
=
Kljun2015:
Kljun, N., Calanca, P., Rotach, M. W., & Schmid, H. P. (2015). A simple 
two-dimensional parameterisation for Flux Footprint Prediction (FFP). 
Geoscientific Model Development, 8(11), 3695–3713. 
https://doi.org/10.5194/gmd-8-3695-2015

KM2001:
Kormann, R., & Meixner, F. X. (2001). An analytical footprint model for 
non-neutral stratification. Boundary-Layer Meteorology, 99(2), 207–224. 
https://doi.org/10.1023/A:1018991015119
