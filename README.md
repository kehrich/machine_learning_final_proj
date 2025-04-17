# machine_learning_final_proj

This repository is a guide to fitting spectral lines to identify potential activity-sensitive lines that could interfere with radial
velocity detections. NEID files are made into templates, which are fed into the pipeline to pull out line parameters. A .csv file is saved
from the pipeline so the resutls can be read and plotted in the appropriate notebook. This is where correlation coefficients are calculated,
and k-means clustering is performed to identify groups of activity indicators.

It is divided into the following files, and notebooks should be run in this order:

1. spectra_pipeline.ipynb : Begin with this notebook. This notebook reads in NEID spectra, creates a template (daily-averaged spectrum) for each day, and uses
a linelist and these templates in the line-by-line analysis to match observed lines to theoretical ones and measure their parameter variations
for each day. These measured dataproducts are saved to a .csv file, so they can be read-in in plot_results.ipynb.

2. plot_results.ipynb : This notebook reads in the dataproducts from the pipeline for analysis. Pearson correlation coefficients are
calculated between line parameters and solar activity indicators. Then, k-means clustering is performed to identify groups of activity
indicators.

3. tools.py : Contains functions used in spectra_pipeline.ipynb and plot_results.ipynb. 

4. final_calcs.csv : These are dataproducts from the Solaster pipeline, which calculates solar activity indicators using Solar Dynamics 
Observatory data.

5. new_gurtovenko_2015.csv : Theoretical line list pulled from Gurtovenko et al. 2015. Observed lines are matched to theoretical ones from 
this list based upon theoretical wavelength and theoretical line depth. If two observed lines are approximately equidistant from the
theoretical line, the observed line with the depth closest to the theoretical depth is selected.

6. neid_solar_data : Folder containing 3 NEID .fits solar files over the course of 20 days near solar noon. These are used to create
the templates.
