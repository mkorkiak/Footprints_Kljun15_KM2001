**Footprints_Kljun15_KM2001**
=
This program calculates the distances where the relative footprint contributions
(50%, 60%, 70%, 80%) are obtained for the Eddypro full output data. The program
uses Kljun et al. (2015) and Kormann & Meixner (2001) methods to estimate
the footprint contributions.

This program needs to be configured before running it, see the parameters and
more detailed descriptions below.

Kormann & Meixner (2001) gives larger footprints and is 6-7 times slower than 
Kljun et al. (2015) in this program.

No fancy libraries needed. Check that Numpy and Pandas are installed.
=
**data_loc**
Location of the eddypro full output datafile.
The location needs to be in hyphens
The folders can be separated by "/" (Linux and windows) or by "\" (Windows)
But if the folders are separated by "\", one needs to add "r" before the first hyphen like this:
data_loc=r'C:\m\Desktop\results\'

**save_loc**
Folder to save the footprints. Same rules as for data_loc.

**meas_height**
Measurement height in meters. Typically means the sonic anemometer height. 
This variable is mandatory.

**canopy_height**
Canopy height
Either canopy height or displacement height is required. Canopy height is used
to calculate the displacement height, if the displacement height is not given.
If disp_height is given, canopy_height is ignored and a warning is raised.
Use the same canopy height as in Eddypro!!!

**disp_height**
Displacement height
Set to None if not known. Don't use zero.
If disp_height is None, the displacement height will be calculated from the
canopy_height similarly as in Eddypro.
Use the same displacement height as in Eddypro!!!

**lat**
Latitude of the measurements
Needed for Kljun et al. 2015 method for boundary layer height calculation.
Needs to be 0 <= lat <= 90

**do_kljun15**
Calculate footprints according to Kljun et al. 2015
Must be True or False

**do_km**
Calculate footprints according to Kormann & Meixner, 2001
Must be True or False


**References**

Kljun2015:
Kljun, N., Calanca, P., Rotach, M. W., & Schmid, H. P. (2015). A simple 
two-dimensional parameterisation for Flux Footprint Prediction (FFP). 
Geoscientific Model Development, 8(11), 3695–3713. 
https://doi.org/10.5194/gmd-8-3695-2015

KM2001:
Kormann, R., & Meixner, F. X. (2001). An analytical footprint model for 
non-neutral stratification. Boundary-Layer Meteorology, 99(2), 207–224. 
https://doi.org/10.1023/A:1018991015119
