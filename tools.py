import numpy as np
from astropy.io import fits
import glob
import matplotlib.pyplot as plt
# %matplotlibinline
import os
from pathlib import Path
import sys
plt.style.use('seaborn-v0_8-darkgrid')
import pandas
from astropy.modeling import models, fitting
from astropy import modeling
from scipy.signal import savgol_filter
from math import floor
from random import choice
import csv
import pandas as pd
import matplotlib.pyplot as plt
from scipy import constants
from scipy import interpolate
from scipy import optimize
from scipy import stats
from scipy import ndimage
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

SPEEDOFLIGHT = 3e5 #km/s

#----------------------------------------
def findel(num, arr):
    '''
    Finds index of nearest element in array to a given number
    Parameters
    ----------
    num : int/float
        number to search for nearest element in array
    arr : :obj:`ndarray` of :obj:`float`
        array to search for ‘num’
    Returns
    -------
    idx : int
        index of array element with value closes to ‘num’
    S Halverson - JPL - 29-Sep-2019
    '''
    # arr = np.array(arr)
    # idx = (np.abs(arr - num)).argmin()
    arr = np.asarray(arr)
    diff = np.abs(arr - num)
    idx = (np.where(diff == np.nanmin(diff))[0][0])
    return idx
#----------------------------------------

# convert from air wavelengths to vacuum
def airtovac(wvl_air):
    s = 104 / wvl_air
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.0001599740894897 / (38.92568793293 - s**2)
    wvl_vac = wvl_air * n  
    return wvl_vac


def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def gauss_line(x, H, A, x0, sigma, B):
    return H + B*x + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
                                
def parabola(x, A, B, C):
    return A * x**2 + B*x + C
    
def compare_lines1(run, line_list, wavelengths, flux, obs_jd, maxima_dict = {}, init_window_wid=0.03, blaze_correct=False, vac=True):
    '''
    Compares line lists to empirical data and returns parameters and characteristics for each line.
    
    Inputs
    ------
    run : 'int'
        run #, used to create unique filenames
    line_list : .csv file
        list of lines to compare to empirical data
    wavelengths : list
        list of empirical stellar spectral wavelength data, parsed to find 
        corresponding lines in line_list
    flux : list
        list of empirical stellar spectral flux data, corresponds to wavelength list
    maxima_dict : dictionary
        dictionary in the format as follows: {empirical_line : left_bound, right bound};
        if empty, boundaries will be computed by peak finding, else provided bounds will be used
    init_window_wid : 'float'
        sets initial window width to match theoretical line with empirical line in terms of 
        wavelength in Angstrom, default = 0.03
    blaze_correct : 'False' or list
        if list provided, flux will be flattened using blaze function; default = 'False'
    vac : Boolean--'True' or 'False'
        if 'False', line_list will be converted using air_to_vac; default = 'True'

    Returns
    -------
    fwhms : list
        list of FWHM vals from Gaussian fits of data
    gauss_meanlambda: list
        list of data wavelengths that have a corresponding line list wavelength
    plotted_lines : list
        list of lines from line list which were successfully matched to an empirical data lineand plotted
    species : list
        list of strings of species that correspond to plotted_lines
    depth : list
        list of amplitudes from Gaussian fit of each data wavelength corresponding to list wavelength
    summary_dict : dict
        dictionary summarizing returns, key is line list wavelength from plotted_lines, and values are fwhms,
        Gaussian mean wavelengths (corresponding data wavelengths), species, and depths corresponding to
        each line list wavelength
    
    K Ehrich - JPL - 15-Jun-2021
    '''

    #find median wavelengths of each order, in Angstroms
    med_wls = np.array([])
    n_ord = len(wavelengths)
    for order_ind in range(n_ord):
        median = np.nanmedian(wavelengths[order_ind, :])
        med_wls = np.append(med_wls, median)

    # convert to appropriate units, current template gives wavelengths in nm 
#     line_wls_air = line_wls_nm * 10 #convert to Angstroms

    if vac == False:
        line_wls = airtovac(line_wls_air)
    else:
        line_wls = line_wls_air
        
    if blaze_correct != False:
        print('Blaze inputted, blaze correct ON')
    elif blaze_correct == False:
        print('Blaze correct OFF')
    
    #initialize arrays
    flux_orders = []
    winning_flux_mins = []
    winning_wl_mins = []
    fwhms = []
    gauss_meanlambda = []
    gauss_meanlambda2 = []
    plotted_lines = []
    species = []
    depths = []
    depths2 = []
    summary_dict = {}
    std_fits = []
    std_fits2 = []
    left_maximas = []
    right_maximas = []
    ords_for_lines = []
    dep_count = 0

    if maxima_dict == {}:
        print('***** LEFT/RIGHT MAXIMAS NOT PROVIDED *****')
        for ind, line in enumerate(line_wls):
            print(f'*** Index is {ind}, line is {line} ***')
            #make absolute depth cut to eliminate noise for fitting purposes
            if cent_line_depth[ind]>0.2:
                #find order with wavelengths closest to each line
                ord_for_line = findel(line, med_wls)
                flux_orders.append(ord_for_line)

                #create boundaries for narrowed window to find local minima
                order_wl_bounds_ind_left = [] #corresponds with flux_orders 
                order_wl_bounds_ind_right = []
                left_wave_ind = findel(line-init_window_wid, wavelengths[ord_for_line])
                right_wave_ind = findel(line+init_window_wid, wavelengths[ord_for_line])
                order_wl_bounds_ind_left.append(left_wave_ind)
                order_wl_bounds_ind_right.append(right_wave_ind)
                wlbound_left = wavelengths[ord_for_line][left_wave_ind]
                wlbound_right = wavelengths[ord_for_line][right_wave_ind]

                #create arrays specific to each order
                flux_ord = flux[ord_for_line,:] / np.nanmax(flux[ord_for_line,:])
                wave_ord = wavelengths[ord_for_line,:]
                if blaze_correct != False:
                    blaze = blaze_correct
                    blaze_ord = blaze[ord_for_line,:] / np.nanmax(blaze[ord_for_line,:])
                    flat_flux = flux_ord / blaze_ord
                elif blaze_correct == False:
                    flat_flux = flux_ord

                #narrow wavelengths to within window
                window_wls = []
                window_fluxes = [] #flattened flux, zoomed in
#                 window_tell = []
                for wave in wave_ord:
                    if wave <= wlbound_right and wave >= wlbound_left:
                        window_wls.append(wave)
                        wave_ind = findel(wave, wave_ord)
                        window_fluxes.append(flat_flux[wave_ind])
#                         window_tell.append(tell_ord[wave_ind])

                nan_inds = []
                for k, i in enumerate(window_fluxes):
                    if np.isnan(i) == True:
                        nan_inds.append(k)

                if len(nan_inds) > 0:
                    max_nanind = max(nan_inds)
                    print(max_nanind)
                    window_fluxes = window_fluxes[max_nanind+1:]
                    window_wls = window_wls[max_nanind+1:]
#                     window_tell = window_tell[max_nanind+1:]



                print('window flux len is', len(window_fluxes))
                

                if len(window_fluxes) > 5 and len(nan_inds) == 0:
                    #find all local minima within general window
                    mins = (np.diff(np.sign(savgol_filter(np.diff(window_fluxes), 11, 5))) > 0).nonzero()[0] + 1 
                    maxes = (np.diff(np.sign(savgol_filter(np.diff(window_fluxes), 11, 5))) < 0).nonzero()[0] + 1 
                    if len(maxes) == 0:
                        print(f'***** No maxima found, line {line} skipped *****')
                    elif len(mins) == 0:
                        print(f'***** No minima found, line {line} skipped *****')
                    else:
                        print('maxes', maxes)
                        print('mins', mins)

                        #find best actual minima, a.k.a. "winning" minima, by comparing to expected line depth
                        initmins_depthcomp = []
                        initmins_wlcomp = []
                        for i in mins:
                            initmins_depthcomp.append(abs(window_fluxes[i] - (1-cent_line_depth[ind])))
                            initmins_wlcomp.append(abs(window_wls[i] - line))

                        mins_depthinds = np.argsort(np.asarray(initmins_depthcomp))
                        mins_depthcomp = np.asarray(initmins_depthcomp)[mins_depthinds]
                        mins_wlcomp = np.asarray(initmins_wlcomp)[mins_depthinds]
                        mins_weightedscoreinds = np.argsort(mins_wlcomp * mins_depthcomp)
                        mins_weightedscore = (mins_wlcomp * mins_depthcomp)[mins_weightedscoreinds]

                        min_flux_subs = []
                        min_wl_subs = []   
                        for m in range(len(mins_depthcomp)):
                            if m <= 5:
                                min_flux_subs.append(mins_depthcomp[m])
                                min_wl_subs.append(mins_wlcomp[m])

                        win_flux_mins = []
                        win_wl_mins = []
                        for minn in mins:
                            if abs(np.asarray(window_fluxes)[minn] - (1-cent_line_depth[ind])) * abs(np.asarray(window_wls)[minn] - line) == mins_weightedscore[0] and np.asarray(window_fluxes)[minn]<0.9:
                                win_flux_mins.append(np.asarray(window_fluxes)[minn])
                                win_wl_mins.append(np.asarray(window_wls)[minn])
            

                        #refined "winning" minima arrays, shortened to 4 max data points
                        ref_winflux_ind = np.argsort(win_flux_mins)
                        ref_win_fluxmins = np.asarray(win_flux_mins)[ref_winflux_ind]
                        ref_win_wlmins = np.asarray(win_wl_mins)[ref_winflux_ind]
                        ref_wl_inds = np.argsort(ref_win_wlmins)
                            
                        
                        if len(ref_win_fluxmins) > 0:
                            flux_maxes = np.asarray(window_fluxes)[maxes]
                            wl_maxes = np.asarray(window_wls)[maxes]
                            closest = findel(ref_win_wlmins[0], wl_maxes)
                            if wl_maxes[closest] - ref_win_wlmins[0] > 0: ## if +
                                if wl_maxes[closest] - ref_win_wlmins[0] > 0.005: ## if legit max
                                    print('winner is', wl_maxes[closest])
                                    right_wl_bound = wl_maxes[closest]
                                    left_wl_bound = ref_win_wlmins[0] - abs(wl_maxes[closest] - ref_win_wlmins[0])
                            elif wl_maxes[closest] - ref_win_wlmins[0] < 0: ## if -
                                if abs(wl_maxes[closest] - ref_win_wlmins[0]) > 0.005: ## if legit max
                                    left_wl_bound = wl_maxes[closest]
                                    right_wl_bound = ref_win_wlmins[0] + abs(ref_win_wlmins[0] - wl_maxes[closest])

                            left_waveord_ind = findel(left_wl_bound, wave_ord)
                            right_waveord_ind = findel(right_wl_bound, wave_ord)
                            winning_flux_mins.append(ref_win_fluxmins[0])
                            winning_wl_mins.append(ref_win_wlmins[0])
                            winwlsub = line - ref_win_wlmins[0]
                            
                            if abs(ref_win_fluxmins[0] - (1-cent_line_depth[ind]))<0.5:
                                #----------create Gaussian fit------------#
                                #create a new window centered around minimum 
                                new_center_ind = findel(ref_win_wlmins[0], wave_ord)
                                flux_for_fit = np.asarray(window_fluxes)
                                wl_for_fit = np.asarray(window_wls)
#                                 tell_for_fit = tell_ord[left_waveord_ind:right_waveord_ind]

                                #create arrays for gaussian fit
                                schunk = flux_for_fit
                                upside = (schunk * -1)+1
                                xchunk = wl_for_fit

                                #initialize and create Gaussian model
                                init = [0, np.nanmax(upside),np.nanmean(xchunk), (wl_for_fit[-1] - wl_for_fit[0])/7, 0]
                                fit_params, covar = curve_fit(gauss_line, wl_for_fit, upside, p0=init, maxfev=50000, bounds=([0,0.01,-np.inf,0,-np.inf],[1,1,np.inf,1,np.inf]))
                                model = gauss_line(xchunk, *fit_params)

                                #find FWHM, and append arrays:
                                fwhm = 2*np.sqrt(2*np.log(2)) * fit_params[3] #in Angstrom
                                fwhms.append(fwhm)
                                print('FWHM =', fwhm)
                                gauss_meanlambda.append(fit_params[2])
                                plotted_lines.append(line)
                                species.append(full_species[ind])
                                depths.append(fit_params[1])
                                summary_dict[line] = fit_params[2], fwhm, fit_params[1], full_species[ind]
                                maxima_dict[line] = wave_ord[left_waveord_ind], wave_ord[right_waveord_ind]
                                left_maximas.append(wave_ord[left_waveord_ind])
                                right_maximas.append(wave_ord[right_waveord_ind])
                                ords_for_lines.append(ord_for_line)
                                
                                #Make narrower fit
                                fit2_boundind_left = findel(fit_params[2]- fit_params[2]*(2/3e5), wl_for_fit)
                                fit2_boundind_right = findel(fit_params[2]+ fit_params[2]*(2/3e5), wl_for_fit)
                                wl_for_fit2 = wl_for_fit[fit2_boundind_left:fit2_boundind_right+1]
                                flux_for_fit2 = flux_for_fit[fit2_boundind_left:fit2_boundind_right+1]
                                if len(wl_for_fit2) >= 5:
                                    print(len(wl_for_fit2))
                                    schunk2 = flux_for_fit2
                                    upside2 = (schunk2 * -1) + 1
                                    xchunk2 = wl_for_fit2
                                    init2 = [1, 0.5, 0.2]
                                    fit_params2, covar2 = curve_fit(parabola, wl_for_fit2, upside2, p0=init2, maxfev=50000, bounds=([-np.inf,-np.inf,0],[np.inf, np.inf,1.5])) 
                                    model2 = parabola(xchunk2, *fit_params2)
                
                                    gauss_meanlambda2.append(-fit_params2[1]/(2*fit_params2[0])) #-b/2a

                                    std_fit2 = np.std(((model2*-1)+1) - flux_for_fit2)
                                    gof2 = std_fit2 / fit_params[1]
                                    std_fits2.append(std_fit2)
                                else:
                                    gauss_meanlambda2.append('NA')
                                    std_fits2.append('NA')


                        
                                #calculate 'goodness of fit' 
                                std_fit = np.std(((model*-1)+1) - flux_for_fit)
                                gof = std_fit / fit_params[1]
                                std_fits.append(std_fit)

                                if gof > 0.05:
                                    print(f'Goodness of Fit is {gof}')
                                    plt.clf()
                                    #Create plots to compare line list to empirical wavelengths and
                                    #plot Gaussian fit, only plot BAD fits
                                    plt.plot(wl_for_fit,flux_for_fit, label='True Observed Line')
                                    plt.plot(window_wls, window_fluxes)
                                    plt.scatter(wl_maxes, flux_maxes)
                                    plt.scatter(wave_ord[left_waveord_ind], flat_flux[left_waveord_ind], c='magenta', label='Maxima for Fit')
                                    plt.scatter(wave_ord[right_waveord_ind], flat_flux[right_waveord_ind], c='magenta')
                                    plt.scatter(ref_win_wlmins[0], ref_win_fluxmins[0], c='r', label='Potential Line Matches')
        #                             plt.plot(window_wls, window_tell,alpha=0.5, label='Telluric Model')
                                    plt.plot(xchunk, (model*-1)+1, linestyle='dashed', label='Gaussian Fit')
                                    plt.ylim(ref_win_fluxmins[0]*0.6, np.nanmax(flat_flux)*1.1)
                                    plt.vlines(fit_params[2], 0, np.nanmax(flux_ord),'r',lw=5,alpha=0.3, label='Line Minimum $\lambda$ Value')
                                    plt.vlines(line, 0, np.nanmax(flux_ord),color='violet',lw=5,alpha=0.3, label='Line List $\lambda$ Value')
                                    plt.hlines(1-cent_line_depth[ind], wl_for_fit[0], wl_for_fit[-1],'teal',lw=5,alpha=0.3, label='Theoretical Depth Value')
                                    plt.ticklabel_format(useOffset=False, style='plain')
                                    plt.title(f'Order {ord_for_line} Wavelength vs. Normalized Flux')
                                    plt.legend(loc=9, fontsize='small')
                                    plt.savefig(f'{spec_path}/parab_fits/{obs_jd}_symm_line_{line}plot_{line_list}_run{run}.png')
                                    if len(wl_for_fit2) >=5:
                                        plt.plot(xchunk2, (model2*-1)+1, linestyle='dashed', label='Narrow Gaussian Fit')
                                    plt.show()
                            
            elif cent_line_depth[ind] < 0.2:
                print('Did not make depth cut')
                dep_count += 1
                

    #if left/right bounds are supplied
    else:
        print('***** LEFT/RIGHT MAXIMAS PROVIDED *****')
        for ind, line in enumerate(line_wls):
            print(f'*** Index is {ind} ***')
            #make absolute depth cut to eliminate noise for fitting purposes
            if cent_line_depth[ind]>0.2 and line in maxima_dict:
                #find order with wavelengths closest to each line
                ord_for_line = findel(line, med_wls)
                flux_orders.append(ord_for_line)

                #create arrays specific to each order
                flux_ord = flux[ord_for_line,:] / np.nanmax(flux[ord_for_line,:])
                wave_ord = wavelengths[ord_for_line,:]

                #narrow wavelengths to within window
                window_wls = []
                window_fluxes = [] #flattened flux, zoomed in

                left_wl_bound = maxima_dict[line][0]
                right_wl_bound = maxima_dict[line][1]
                for wave in wave_ord:
                    if wave <= right_wl_bound and wave >= left_wl_bound:
                        window_wls.append(wave)
                        wave_ind = findel(wave, wave_ord)
                        window_fluxes.append(flux_ord[wave_ind])
                
                nan_inds = []
                for k, i in enumerate(window_fluxes):
                    if np.isnan(i) == True:
                        nan_inds.append(k)

                if len(nan_inds) > 0:
                    max_nanind = max(nan_inds)
                    print(max_nanind)
                    window_fluxes = window_fluxes[max_nanind+1:]
                    window_wls = window_wls[max_nanind+1:]
                

                if len(window_fluxes) > 5 and len(nan_inds) == 0:
                    win_wl_maxes = [left_wl_bound, right_wl_bound]
                    win_flux_maxes = [window_fluxes[0], window_fluxes[-1]]


                    left_waveord_ind = findel(left_wl_bound, wave_ord)
                    right_waveord_ind = findel(right_wl_bound, wave_ord)


                    #----------create Gaussian fit------------#
                    #create a new window centered around minimum 
                    flux_for_fit = np.asarray(window_fluxes)
                    wl_for_fit = np.asarray(window_wls)
                    print('flux for fit', type(flux_for_fit))
                    schunk = flux_for_fit
                    upside = (schunk * -1)+1
                    xchunk = wl_for_fit

                    #initialize and create Gaussian model
                    init = [0, np.nanmax(upside),np.nanmean(xchunk), (wl_for_fit[-1] - wl_for_fit[0])/7, 0]
                    fit_params, covar = curve_fit(gauss_line, wl_for_fit, upside, p0=init, maxfev=50000,bounds=([0,0.01,-np.inf,0,-np.inf],[1,1,np.inf,1,np.inf]))
                    model = gauss_line(xchunk, *fit_params)
                    
                    #find minimum to plot
                    wl_min = fit_params[2]
                    mean_ind = findel(wl_min, window_fluxes)
                    flux_min = window_fluxes[mean_ind]

                    #find FWHM, and append arrays:
                    fwhm = 2*np.sqrt(2*np.log(2)) * fit_params[3] #in Angstrom
                    fwhms.append(fwhm)
                    print('FWHM =', fwhm)
                    gauss_meanlambda.append(fit_params[2])
                    plotted_lines.append(line)
                    species.append(full_species[ind])
                    depths.append(fit_params[1])
                    summary_dict[line] = fit_params[2], fwhm, fit_params[1], full_species[ind]
                    left_maximas.append(left_wl_bound)
                    right_maximas.append(right_wl_bound)
                    maxima_dict[line] = left_wl_bound, right_wl_bound


        #             #calculate 'goodness of fit' 
                    std_fit = np.std(((model*-1)+1) - flux_for_fit)
                    gof = std_fit / fit_params[1]
                    std_fits.append(std_fit)
                    
                    #create second, narrower Gaussian fit 
                    fit2_boundind_left = findel(fit_params[2]- fit_params[2]*(2/3e5), wl_for_fit)
                    fit2_boundind_right = findel(fit_params[2]+ fit_params[2]*(2/3e5), wl_for_fit)
                    wl_for_fit2 = wl_for_fit[fit2_boundind_left:fit2_boundind_right+1]
                    if len(wl_for_fit2) >= 5:
                        flux_for_fit2 = flux_for_fit[fit2_boundind_left:fit2_boundind_right+1]
                        schunk2 = flux_for_fit2
                        upside2 = (schunk2 * -1) + 1
                        xchunk2 = wl_for_fit2
#                         init2 = [0, np.nanmax(upside2),np.nanmean(xchunk2), (wl_for_fit2[-1] - wl_for_fit2[0])/7, 0]
#                         fit_params2, covar2 = curve_fit(gauss_line, wl_for_fit2, upside2, p0=init2, maxfev=50000, bounds=([0,0.01,-np.inf,0,-np.inf],[1,1,np.inf,1,np.inf]))
#                         model2 = gauss_line(xchunk2, *fit_params2)

                        init2 = [1, 0.5, 0.2]
                        fit_params2, covar2 = curve_fit(parabola, wl_for_fit2, upside2, p0=init2, maxfev=50000, bounds=([-np.inf,-np.inf,0],[np.inf, np.inf,1.5])) 
                        model2 = parabola(xchunk2, *fit_params2)
    
                        gauss_meanlambda2.append(-fit_params2[1]/(2*fit_params2[0])) #-b/2a

                        std_fit2 = np.std(((model2*-1)+1) - flux_for_fit2)
                        gof2 = std_fit2 / fit_params[1]
                        std_fits2.append(std_fit2)
                    else:
                        gauss_meanlambda2.append('NA')
                        std_fits2.append('NA')

                    print(f'Line is {line}, left bound is {left_wl_bound} right bound is {right_wl_bound}')

                    if gof > 0.05:
                        print(f'Goodness of Fit is {gof}')
                        plt.clf()
                        #Create plots to compare line list to empirical wavelengths and plot Gaussian fit
                        plt.plot(wl_for_fit,flux_for_fit, label='True Observed Line')
                        plt.scatter(win_wl_maxes, win_flux_maxes, c='magenta', label='Maxima for Fit')
                        plt.scatter(wl_min, flux_min, c='r', label='Potential Line Matches')
                        plt.plot(xchunk, (model*-1)+1, linestyle='dashed', label='Gaussian Fit')
                        if len(wl_for_fit2) >=5:
                            plt.plot(xchunk2, (model2*-1)+1, linestyle='dashed', label='Narrow Gaussian Fit')
                        plt.ylim(flux_min*0.6, np.nanmax(flux_ord)*1.1)
                        plt.vlines(fit_params[2], 0, np.nanmax(flux_ord),'r',lw=5,alpha=0.3, label='Line Minimum $\lambda$ Value')
                        plt.vlines(line, 0, np.nanmax(flux_ord),color='violet',lw=5,alpha=0.3, label='Line List $\lambda$ Value')
                        plt.hlines(1-cent_line_depth[ind], wl_for_fit[0], wl_for_fit[-1],'teal',lw=5,alpha=0.3, label='Theoretical Depth Value')
                        plt.ticklabel_format(useOffset=False, style='plain')
                        plt.title(f'Order {ord_for_line} Wavelength vs. Normalized Flux')
                        plt.legend(loc=9, fontsize='small')
                        plt.savefig(f'{spec_path}/parab_fits/{obs_jd}_symm_line_{line}plot_{line_list}_run{run}.png')
                        plt.show()
        
    fwhms = np.asarray(fwhms)[np.asarray(fwhms) < 1]
    plt.clf()
    n, foo, patches = plt.hist(x=fwhms, bins=50)
    plt.xlabel('FWHM Values')
    plt.ylabel('N')
    plt.title('Histogram of FWHM Values')
    plt.savefig(f'{spec_path}/parab_fits/{obs_jd}_symm_fwhm_histogram_{line_list}_run{run}.png')
    plt.show()
    
    plt.clf()
    n1, foo1, patches1 = plt.hist(x=(np.asarray(std_fits) / np.asarray(depths)), bins=50)
    plt.xlabel('St Dev. / Amplitude')
    plt.ylabel('N')
    plt.title('Histogram of Standard Deviation / Amplitude Values')
    plt.savefig(f'{spec_path}/parab_fits/{obs_jd}_symm_std_histogram_{line_list}_run{run}.png')
    plt.show()
    
    plotted_lines = np.array(plotted_lines)
    depths = np.array(depths)
    species = np.array(species)
    fwhms = np.array(fwhms)
    gauss_meanlambda = np.array(gauss_meanlambda)
    left_maximas = np.array(left_maximas)
    right_maximas = np.array(right_maximas)
    gofs = np.asarray(std_fits) / np.asarray(depths)
    gauss_meanlambda2 = np.array(gauss_meanlambda2)
    
    np.savetxt(f'{spec_path}/parab_fits/{obs_jd}_symm_summarytable_{line_list}_run{run}.csv', np.column_stack([plotted_lines, depths, species, fwhms, gauss_meanlambda, left_maximas, right_maximas, gofs, gauss_meanlambda2]), header='line (Angstroms), depth, species, FWHMs, observed line (mean), left bound, right bound, goodness of fit, narrow observed line', delimiter=',', fmt='%s')
    print('FINAL LENGTH OF DEPTH COUNTER IS', dep_count)


    return plotted_lines, depths, species, fwhms, gauss_meanlambda, summary_dict, maxima_dict, ords_for_lines


'''Template creation code below written by Abigail Burrows'''

def blaze_correct(lamp_spec,test_ord,wave,spec):
    percentile_flux = np.nanpercentile(spec[test_ord,:], 95)
    # 1. Fit with the lamp
    if test_ord == 15: # order 17 if we don't remove the NaN orders
        # for lower
        ratio = spec[test_ord-1,3000:8000]/lamp_spec[test_ord-1,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fitlow = lamp_spec[test_ord-1,:]*scale_fac

        # for upper
        ratio = spec[test_ord+2,3000:8000]/lamp_spec[test_ord+2,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fithigh = lamp_spec[test_ord+2,:]*scale_fac

        # tie together and average
        stacked = np.vstack((fitlow, fithigh))
        fitspec = np.nanmean(stacked, axis=0)
        spec_noblaze = spec[test_ord,:]/fitspec
        wvl_test = wave[test_ord,:]
    elif test_ord == 16: # order 18 if we don't remove the NaN order
        # for lower
        ratio = spec[test_ord-2,3000:8000]/lamp_spec[test_ord-2,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fitlow = lamp_spec[test_ord-2,:]*scale_fac

        # for upper
        ratio = spec[test_ord+1,3000:8000]/lamp_spec[test_ord+1,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fithigh = lamp_spec[test_ord+1,:]*scale_fac

        # tie together and average
        stacked = np.vstack((fitlow, fithigh))
        fitspec = np.nanmean(stacked, axis=0)
        spec_noblaze = spec[test_ord,:]/fitspec
        wvl_test = wave[test_ord,:]
    # elif test_ord > 70:
    #     lamp_norm = lamp_spec[test_ord,1000:8000]/np.nanmax(lamp_spec[test_ord,1000:8000])
    #     lamp_norm = lamp_norm * np.nanpercentile(spec[test_ord,1000:8000],95)
    #     spec_norm = spec[test_ord,1000:8000]
    #     wvl_test = wave[test_ord,1000:8000]
    #     spec_noblaze = spec_norm/lamp_norm
    else:
        lamp_norm = lamp_spec[test_ord,:]/np.nanmax(lamp_spec[test_ord,:])
        lamp_norm = lamp_norm * np.nanpercentile(spec[test_ord,:],95)
        spec_norm = spec[test_ord,:]
        wvl_test = wave[test_ord,:]
        spec_noblaze = spec_norm/lamp_norm

    # ok try to do a various polynomial fit
    if test_ord > 10 and test_ord < 115:
        wave_portions = np.array_split(wvl_test, 9)
        flux_portions = np.array_split(spec_noblaze, 9)
        maxs_flux = []
        maxs_wvl = []
        for i in range(len(flux_portions)):
            portion = flux_portions[i]
            # maxes
            maxs_flux.append(np.nanmax(portion))
            if len(wave_portions[i][np.where(portion == np.nanmax(portion))]) > 1:
                maxs_wvl.append(float(wave_portions[i][np.where(portion == np.nanmax(portion))][0]))
            else:
                maxs_wvl.append(float(wave_portions[i][np.where(portion == np.nanmax(portion))]))
        # add first -- max in the first 
        maxs_flux.append(np.nanmax(spec_noblaze[0:25]))
        maxs_wvl.append(wave_portions[0][0])

        # add last
        maxs_flux.append(np.nanmax(spec_noblaze[-25:-1]))
        maxs_wvl.append(wave_portions[-1][-1])

        # spline
        spline = interpolate.interp1d(maxs_wvl, maxs_flux)
        spline_lines = spline(wvl_test)

        spec_noblaze = spec_noblaze/spline_lines

    return spec_noblaze*percentile_flux

def reverse_blaze(lamp_spec,test_ord,wave,spec):
    percentile_flux = np.nanpercentile(spec[test_ord,:], 95)
    # 1. Fit with the lamp
    if test_ord == 15: # order 17 if we don't remove the NaN orders
        # for lower
        ratio = spec[test_ord-1,3000:8000]/lamp_spec[test_ord-1,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fitlow = lamp_spec[test_ord-1,:]*scale_fac

        # for upper
        ratio = spec[test_ord+2,3000:8000]/lamp_spec[test_ord+2,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fithigh = lamp_spec[test_ord+2,:]*scale_fac

        # tie together and average
        stacked = np.vstack((fitlow, fithigh))
        lamp_norm = np.nanmean(stacked, axis=0)
        spec_noblaze = spec[test_ord,:]/lamp_norm
        wvl_test = wave[test_ord,:]
    elif test_ord == 16: # order 18 if we don't remove the NaN order
        # for lower
        ratio = spec[test_ord-2,3000:8000]/lamp_spec[test_ord-2,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fitlow = lamp_spec[test_ord-2,:]*scale_fac

        # for upper
        ratio = spec[test_ord+1,3000:8000]/lamp_spec[test_ord+1,3000:8000]
        scale_fac= np.nanpercentile(ratio, 95)

        fithigh = lamp_spec[test_ord+1,:]*scale_fac

        # tie together and average
        stacked = np.vstack((fitlow, fithigh))
        lamp_norm = np.nanmean(stacked, axis=0)
        spec_noblaze = spec[test_ord,:]/lamp_norm
        wvl_test = wave[test_ord,:]
    else:
        lamp_norm = lamp_spec[test_ord,:]/np.nanmax(lamp_spec[test_ord,:])
        lamp_norm = lamp_norm * np.nanpercentile(spec[test_ord,:],95)
        spec_norm = spec[test_ord,:]
        wvl_test = wave[test_ord,:]
        spec_noblaze = spec_norm/lamp_norm

    # ok try to do a various polynomial fit
    if test_ord > 10 and test_ord < 115:
        wave_portions = np.array_split(wvl_test, 9)
        flux_portions = np.array_split(spec_noblaze, 9)
        maxs_flux = []
        maxs_wvl = []
        for i in range(len(flux_portions)):
            portion = flux_portions[i]
            # maxes
            maxs_flux.append(np.nanmax(portion))
            if len(wave_portions[i][np.where(portion == np.nanmax(portion))]) > 1:
                maxs_wvl.append(float(wave_portions[i][np.where(portion == np.nanmax(portion))][0]))
            else:
                #print(wave_portions[i][np.where(portion == np.nanmax(portion))])
                maxs_wvl.append(float(wave_portions[i][np.where(portion == np.nanmax(portion))]))
        # add first -- max in the first 
        maxs_flux.append(np.nanmax(spec_noblaze))
        maxs_wvl.append(wave_portions[0][0])

        # add last
        maxs_flux.append(np.nanmax(spec_noblaze))
        maxs_wvl.append(wave_portions[-1][-1])

        # spline
        spline = interpolate.interp1d(maxs_wvl, maxs_flux)
        spline_lines = spline(wvl_test)
        #spline_lines = np.pad(spline_lines, (1000,9216-8000), 'empty')
        #print(spline_lines)

        spec_noblaze = spec_noblaze/spline_lines
        return lamp_norm*spline_lines
    return lamp_norm
    #return wvl_test, spec_noblaze*percentile_flux

def wave_vel_shift(wave, vel):
    '''
    Appropriately shifts wavelength(s) 'wave' for velocity 'vel'

    Parameters
    ----------
    wave : float or array of floats
        Wavelength(s) to be stretched

    vel : float 
        Velocity to shift wave [km/s]
    Returns
    -------
    wave_shifted : float or array of floats
        Wavelength(s) shifted by velocity 'vel'
    '''

    # shift wavelengths to provided velocity
    wave_shifted = wave * (1. + vel / constants.c)

    return wave_shifted

def airtovac(wvl_air):
    s = 104 / wvl_air
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2) + 0.0001599740894897 / (38.92568793293 - s**2)
    wvl_vac = wvl_air * n

    return wvl_vac

def make_temps(sort_files):
    '''Creates a template for each day
    
    Inputs
    --------------
    sort_files: list
        list of files from same day within larger list which will be used to create each template
    
    Outputs
    --------------
    templates: list
        list of templates with one template per day
    masterwaves: list
        list of wavelengths to be used in conjunction with template fluxes
        
    K Ehrich - JPL - 18-Jul-2022
    '''

    for jd_ind, spec_fits_files in enumerate(sort_files):
        #---------------------------Master Wavelength----------------------------------
        # find average measured_rv
        print(jd_ind,'/',len(sort_files))
        avg_ccfrv = 0
        num_files = 0
        for file in spec_fits_files:
            ccf_header = fits.getheader(file, 'CCFS')
            measured_rv = ccf_header['CCFRVMOD']
            avg_ccfrv += measured_rv
            num_files += 1
        avg_ccfrv = avg_ccfrv/num_files
        # first, make the master_wavelengths file
        wave_master = fits.getdata(spec_fits_files[0], 'SCIWAVE') / 10.
        # red-shift each spectrum to account for the barrycentric velocity
        # shift: lambda_rest = lambda_obs / (1+v/c)
        # pull velocity information from header -- this is from neid_measure_act slash this is how we were shifting for velocity offset there
        hdr_spec = fits.getheader(spec_fits_files[0])
        ccf_header = fits.getheader(spec_fits_files[0],'CCFS')

        bc_vel = hdr_spec['SSBRV160'] # km/s
        q_vel = hdr_spec['QRV'] # km/s
        measured_rv = ccf_header['CCFRVMOD']

        # shift wavelengths to stellar rest frame
        vel_factor = -1. * (float(bc_vel) + float(q_vel) + float(measured_rv)) - avg_ccfrv
        # add a 7.55 for the otherstar -- add a variable for adding extra velocity correction?
        # # add .675 for the sun -- trying to get closer from the start (may need to add this to the template but try this first)

        # shift factor
        z = (1-((vel_factor*10**3)/constants.c))

        # shift wave
        wave_master = wave_master * z # or is it wave/z ?
        # print(wave_master[50][3000:3016])
        #----------------------------Initialize Template Array---------------------------------
        template = np.zeros(np.shape(wave_master))
        #----------------------------Loop---------------------------------
        # For each file: blaze correct, cut out the first 1000 and last 1000 pixels
        hasblaze = 0
        noblaze = 0
        for fits_file in spec_fits_files:
            main_header = fits.getheader(fits_file)
            #if main_header['DRIFTFUN'] == 'simultaneous' and main_header['WAVECAL'] == 'LFCplusThAr': # trying to take out bad days
            if True:
                #----------------------------Read In---------------------------------
                # read in wave, spec, and var
                # the actual spectrum (in flux units)
                flux = fits.getdata(fits_file, 'SCIFLUX')

                # the variance or noise of the spectrum
                var = fits.getdata(fits_file, 'SCIVAR')

                # the wavelength solution of the spectrum
                wave = fits.getdata(fits_file, 'SCIWAVE') / 10. # nm

                if np.shape(flux) == (122, 9216):
                    while True:
                        try:
                            blaze = fits.getdata(fits_file, 'SCIBLAZE')
                            hasblaze += 1
                            break
                        except KeyError:
                            #print(fits_file)
                            noblaze += 1
                            blaze = np.array((10,10))
                            break
                    #---------------------------Clean NaN orders----------------------------------
                    # get rid of NaN orders, if any
                    # get order mean wavelengths
                    nord = flux.shape[0]
                    wvls_ord = []

                    # get mean order wavelengths
                    for ind, order in enumerate(range(nord)):
                        wvls_ord.append(np.nanmean(wave[ind,:]))
                    wvls_ord = np.asarray(wvls_ord)

                    # now cut the nan orders
                    nan_ords = np.isnan(wvls_ord)
                    #----------------------------Blaze correct and cut---------------------------------
                    # blaze correct each spectrum, cut out the first and last, and save
                    if np.shape(blaze) != np.shape(flux):
                        for i in range(len(flux)):
                            if not nan_ords[i]:
                                wvl_ord = wave[i]
                                try:
                                    flux_ord = blaze_correct(i,wave,flux)
                                except TypeError:
                                    print('Skipped because missing argument')
                                    print(fits_file)
                                #flux_ord = flux[i] / blaze[i]
                                if i == 50:
                                    plt.scatter(wvl_ord, flux_ord)
                                flux_ord = flux_ord / ndimage.maximum_filter1d(flux_ord, 3, axis=- 1, output=None, mode='reflect', cval=0.0, origin=0)
                            else:
                                wvl_ord = wave[i]
                                flux_ord = flux[i]
                            wave[i] = wvl_ord
                            flux[i] = flux_ord
                    else:
                        for i in range(len(flux)):
                            if not nan_ords[i]:
                                wvl_ord = wave[i]
                                #flux_ord = blaze_correct(i,wave,flux)
                                percentile_flux = np.nanpercentile(flux[i,:], 95)
                                flux_ord = flux[i] / blaze[i]
                                #flux_ord = flux_ord / ndimage.maximum_filter1d(flux_ord, 3, axis=- 1, output=None, mode='reflect', cval=0.0, origin=0)
                                flux_ord = flux_ord*percentile_flux
                            else:
                                wvl_ord = wave[i]
                                flux_ord = flux[i]
                            wave[i] = wvl_ord
                            flux[i] = flux_ord
                    #-------------------------------Shift into Star's RF------------------------------
                    # red-shift each spectrum to account for the barrycentric velocity
                    # shift: lambda_rest = lambda_obs / (1+v/c)
                    # pull velocity information from header -- this is from neid_measure_act slash this is how we were shifting for velocity offset there
                    hdr_spec = fits.getheader(fits_file)
                    ccf_header = fits.getheader(fits_file,'CCFS')

                    bc_vel = hdr_spec['SSBRV160'] # km/s # 'Barycentric corr. (km/s) for order 160' why do you take the 160th order? # barrycentric correction for each order WE CAN USE ONE FOR ALL ORDERS BC SUN
                    q_vel = hdr_spec['QRV'] # km/s # 'Queue RV for Star'
                    measured_rv = ccf_header['CCFRVMOD'] # 'Bary. corrected RV for weighted orders' -- this is the actual RV

                    # shift wavelengths to stellar rest frame
                    vel_factor = -1. * (float(bc_vel)) + float(q_vel) + float(measured_rv) - avg_ccfrv
                    # add a 7.55 for the otherstar -- add a variable for adding extra velocity correction?
                    # add .675 for the sun -- trying to get closer from the start (may need to add this to the template but try this first)


                    # shift factor
                    z = (1-((vel_factor*10**3)/constants.c))
                    z_file = (1+hdr_spec['SSBZ060'])

                    # shift wave
                    wave = wave * z # or is it wave/z ?
                    #wave = wave*z_file
                    #------------------------------Interpolate with Master Wavelengths-------------------------------
                    for ind in range(len(wave)):
                        f = interpolate.interp1d(wave[ind], flux[ind], fill_value="extrapolate")# can only use the fill value if your new values are BARELY out of range
                        flux[ind] = f(wave_master[ind])
                        wave[ind] = wave_master[ind]

                        # cut flux
                        flux[ind][0:1000] = np.nan
                        flux[ind][8000:] = np.nan

                        # cut variance
                        var[ind][0:1000] = np.nan
                        var[ind][8000:] = np.nan
                        #------------------------------Add to template-------------------------------
                        template[ind] += flux[ind]
                        
        # save each template to .csv file
        templatePath = os.path.join(local_path, f'test_temp{uniques[jd_ind]}.csv')
        with open(templatePath,"w+") as my_csv:
            csvWriter = csv.writer(my_csv,delimiter=',')
            csvWriter.writerows(template)
            
        # save each masterwave to a .csv file
        masterwavePath = os.path.join(local_path, f'test_wave{uniques[jd_ind]}.csv')
        with open(masterwavePath,"w+") as my_csv:
            csvWriter = csv.writer(my_csv,delimiter=',')
            csvWriter.writerows(wave_master)    
                        
        # read in template
        template = np.loadtxt(open(templatePath, "rb"), delimiter=",")

        # read in masterwave
        masterwave = np.loadtxt(open(masterwavePath, "rb"), delimiter=",")
        

    return template, masterwave












