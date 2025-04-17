# Identifying Activity-Sensitive Lines in the Sun Using K-Means Clustering

This repository is a guide to fitting spectral lines to identify potential activity-sensitive lines that could interfere with radial
velocity detections. NEID files are made into templates, which are fed into the pipeline to pull out line parameters. A .csv file is saved
from the pipeline so the resutls can be read and plotted in the appropriate notebook. This is where correlation coefficients are calculated,
and k-means clustering is performed to identify groups of activity indicators.

# Scientific Motivation

Radial Velocity detections are one of the most robust ways of detecting exoplanets. When a planet orbits around its host star, the pull of the planet causes the star 
to "wiggle" about the center of gravity between the two objects. This creates a Doppler-shifting effect--the star appears blue-shifted as it moves toward us,
and red-shifted when it moves away from us. We can observe this effect through spectroscopy, which is a method of separating light into its component parts
by wavelength using a diffraction grating. We observe the Doppler shift in the star's spectrum, as seen below.

![Radial_velocity_doppler_spectroscopy](https://github.com/user-attachments/assets/defe75d8-135b-4126-a731-2cbc2988de45)

Unfortunately, the star has a very active atmosphere, which can distort spectral lines, and interfere with RV signals. In fact, in order to detect an Earth-like planet around 
a star similar in size to our Sun, we would need to be able to detect a velocity of ~9 cm/s, but stellar activity creates a noise floor around 1 m/s, nearly an order of
magnitude above this signal. Thus, it is important to figure out which spectral lines are magnetically-sensitive, so our RV detections are less-prone to
contamination from the stellar surface. We begin with a line list created for the Sun, which uses physics and radiative-transfer code to calculate
theoretical wavelength, depth, and full-width half-maximum (FWHM) values for spectral lines we should observe in the Sun. This is the motivation for this project. 

## Here is an outline of the project:

1. spectra_pipeline.ipynb : Begin with this notebook. This notebook starts with reading in NEID spectra, which are in the format of .fits files with can be found
at the Google Drive link below:

Link to neid_solar_data on Google Drive: https://drive.google.com/file/d/1sRNEKcAEOc93sTisQP-flMgLUUpcklcV/view?usp=drive_link
Folder containing 3 NEID .fits solar files over the course of 20 days near solar noon. These are used to create
the templates.

There are 3 .fits files containing spectra per day taken near solar noon over the course of 20 days. Using these files, we then
create a template for each day. A template is a daily-averaged spectrum, and increases the SNR of our data while also filling in data gaps or nan values.
We average these 3 files per day to create 20 total templates. An example template spectrum is shown below:

<img width="891" alt="Screenshot 2025-04-17 at 3 12 22 PM" src="https://github.com/user-attachments/assets/9a921aec-c9b2-4be8-a4b7-689caec3ca81" />

Once these templates are created, they are fed into the line-by-line pipeline. This goes through each individual line list wavelength, and finds where it should
be located in the observed NEID spectrum (the template). Once a line is identified, a modified Gaussian is fit to it, and FWHM and Depth are recorded.
Below shows an example of a single absorption line from the NEID spectra, along with the parameters that are measured.

<img width="443" alt="Screenshot 2025-04-17 at 2 42 59 PM" src="https://github.com/user-attachments/assets/54a4eaf4-d6a3-4431-8e37-cdfd77397e96" />

These parameters are written to a .csv file in the user's directory, so they can be read in in the second notebook, plot_results.

3. plot_results.ipynb : This notebook reads in the dataproducts from the pipeline for analysis. Pearson correlation coefficients are
calculated between line parameters and solar activity indicators. Then, k-means clustering is performed to identify groups of activity
indicators.

4. tools.py : Contains functions used in spectra_pipeline.ipynb and plot_results.ipynb. Included in the import statements.

5. final_calcsold.csv : These are dataproducts from the Solaster pipeline, which calculates solar activity indicators using Solar Dynamics 
Observatory data.

6. new_gurtovenko_2015.csv : Theoretical line list pulled from Gurtovenko et al. 2015. Observed lines are matched to theoretical ones from 
this list based upon theoretical wavelength and theoretical line depth. If two observed lines are approximately equidistant from the
theoretical line, the observed line with the depth closest to the theoretical depth is selected.

Link to neid_solar_data on Google Drive: https://drive.google.com/file/d/1sRNEKcAEOc93sTisQP-flMgLUUpcklcV/view?usp=drive_link
Folder containing 3 NEID .fits solar files over the course of 20 days near solar noon. These are used to create
the templates.
