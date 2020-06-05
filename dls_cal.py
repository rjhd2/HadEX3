#!/usr/local/sci/bin/python
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 492                                     $:  Revision of last commit
#$Author::                                      $:  Author of last commit
#$Date:: 2020-04-29 16:51:05 +0100 (Wed, 29 Apr#$:  Date of last commit
#------------------------------------------------------------
#  
#  Calculate the Decorrelation Length Scale (DLS) for the data   
#
#------------------------------------------------------------
'''
Calculates the Decorrelation Length Scale from all the index data

Outline of the process:
* Read in the Station ID ist
* Read in all the station information
* Assign station to latitude bands
* Cross Correlate everything within each band
* Bin up the results
* Fit with an exponential decay
* Find the Decorrelation Length Scale (DLS) and plot
* Interpolate between the bands
* Write output files

dls_cal.py invoked by typing::

  python dls_cal.py --index "TX90p" --qc_flags --diagnostics

Input arguments:

--index         Which ETCCDI index to use

--qc_flags      Which QC flags to use

--diagnostics   Output extra info (default False)

'''
import os
import datetime as dt
import calendar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.stats import chi2
import numpy as np

import utils


# set up calendar month names
month_names = calendar.month_abbr[:]
month_names[0] = "Ann"


mdi = -99.9

# constants for DLS calculation
SUFFICIENT_YEARS = 30
MIN_PER_BIN = 10
MAX_SEPARATION = utils.MAX_DLS # perhaps 6000?
BIN_WIDTH = 100

# for fitting the exponential (offset of asymptote from y=zero) - take into an argparse later
C = 0.0

nyears = len(utils.REFERENCEYEARS)

# #*******************************************
# def map_2_points(point0, point1):
#     '''
#     Copied from FORTRAN code in utility_2.5.2.f90
#     Adapted to work on all array at once.
#     Calculates the great circle shortest distance and angle between two points

#     :param array point0: input point (lat, lon) - gridbox centre, repeated to match length of station
#     :param array point1: input point (lat, lon) - station
#     :returns: separation, angle - separation in kilometres and angle in radians (-pi to pi)
#     '''

#     if len(point0) != len(point1):
#         print("lengths do not match")
#         raise ValueError

#     lat0, lon0 = point0[:, 0], point0[:, 1]
#     lat1, lon1 = point1[:, 0], point1[:, 1]

#     R_earth = 6378206.4  # metres

#     coslt1 = np.cos(np.radians(lat1))
#     sinlt1 = np.sin(np.radians(lat1))
#     coslt0 = np.cos(np.radians(lat0))
#     sinlt0 = np.sin(np.radians(lat0))
    
#     cosl0l1 = np.cos(np.radians(lon1 - lon0))
#     sinl0l1 = np.sin(np.radians(lon1 - lon0))

#     cosc = sinlt0 * sinlt1 + coslt0 * coslt1 * cosl0l1

#     gt1 = np.where(cosc > 1.)
#     if len(gt1) > 0: cosc[gt1] = 1.
#     lt1 = np.where(cosc < -1.)
#     if len(lt1) > 0: cosc[lt1] = -1.

#     #if cosc > 1. : cosc=1.
#     #if cosc < -1.: cosc=-1.

#     sinc = np.sqrt(1.0 - cosc**2)
    

#     cosaz = np.ones(len(point0))
#     sinaz = np.zeros(len(point0))

#     nonzero = np.where(np.abs(sinc) > 1.e-7)

#     cosaz[nonzero] = ((coslt0[nonzero] * sinlt1[nonzero]) - \
#                         (sinlt0[nonzero] * coslt1[nonzero] * cosl0l1[nonzero])) / sinc[nonzero]
#     sinaz[nonzero] = sinl0l1[nonzero] * coslt1[nonzero] / sinc[nonzero]

#     #if np.abs(sinc) > 1.e-7:
#     #    cosaz=(coslt0 * sinlt1 - sinlt0*coslt1*cosl0l1) / sinc #Azimuth
#     #    sinaz = sinl0l1*coslt1/sinc
#     #else:		#Its antipodal
#     #    cosaz = 1.d0
#     #    sinaz = 0.d0

#     separation = np.arccos(cosc) * R_earth / 1000. # convert from m to km
#     angle = np.arctan2(sinaz, cosaz)

#     return separation.astype(int), angle.astype(int) #  map_2_points

#*******************************************
def corr2_coeff_matrix(A, B): # DEPRECATED
    """Correlate each n with each m.

    DEPRECATED
    Correlate each n with each m.

    :param np.array A: Shape N x T
    :param np.array B: Shape M x T

    :returns np.array: Shape N x M in which each element is a correlation coefficient.

    RJHD - 1/2/2016
    http://stackoverflow.com/questions/30143417/computing-the-correlation-coefficient-between-two-multi-dimensional-arrays

    This cannot handle the common array step needed for full arraywise operation
    """
     # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1);
    ssB = (B_mB**2).sum(1);

    # Finally get corr coeff
    return np.ma.dot(A_mA, B_mB.T)/np.ma.sqrt(np.ma.dot(ssA[:, None], ssB[None])) # corr2_coeff_matrix

#*******************************************
def corr2_coeff_loop(A): # DEPRECATED
    """Correlate each n with each m.

    DEPRECATED
    Correlate each n with each m.

    :param np.array A: Shape N x T

    :returns np.array: Shape N x M in which each element is a correlation coefficient.
 
    Uses terminology from above subroutine but long-hand, loopy way.
    http://stackoverflow.com/questions/30143417/computing-the-correlation-coefficient-between-two-multi-dimensional-arrays

    scipy.stats.mstats.pearsonr works in this instance as just two, prepared timeseries

    But this is very slow.  Left in code for completeness
    """

    # A is full matrix.
    import scipy.stats

    corrs = np.ma.zeros([A.shape[0], A.shape[0]])
    corrs.mask = np.ones(corrs.shape)

    # copy from UNSW Fortran code.  

    start_time = time.time()    
    for i, raw_source in enumerate(A):
        for j, raw_target in enumerate(A):
            if j <= i: continue
            
            # align masks
            common_mask = np.ma.mask_or(raw_source.mask, raw_target.mask)

            # need to ensure that A.mask isn't amended
            source = np.ma.array(raw_source.data[:], mask=common_mask)
            target = np.ma.array(raw_target.data[:], mask=common_mask)
            
            # ensure sufficient overlap
            if len(source.compressed()) < SUFFICIENT_YEARS: continue

            corrs[i, j] = corrs[j, i] = scipy.stats.mstats.pearsonr(source, target)[0]
            corrs.mask[i, j] = corrs.mask[j, i] = False

            # below is transcribed from Fortran
            # this produces the same as the command scipy.stats.mstats.pearsonr(source, target)[0]
            #  retained for traceability
            # 29 August 2017 RJHD

            # S_mS = source - np.ma.mean(source)
            # T_mT = target - np.ma.mean(target)

            # ssS = np.ma.sum(S_mS**2)
            # ssT = np.ma.sum(T_mT**2)
        
            # corrs[i, j] = corrs[j, i] = np.ma.sum(S_mS * T_mT) / (np.ma.sqrt(ssS) * np.ma.sqrt(ssT))
            # corrs.mask[i, j] = corrs.mask[j, i] = False

    return corrs # corr2_coeff_loop

#*******************************************
def corr2_coeff(A):
    """Correlate each n with each m.

    :param np.array A: Shape N x T

    :returns np.array: Shape N x M in which each element is a correlation coefficient.
 
    Uses terminology from:
    http://stackoverflow.com/questions/30143417/computing-the-correlation-coefficient-between-two-multi-dimensional-arrays
    but applying common mask row-wise.  Testing to ensure match with mstats.pearsonr occurred in loop version

    Trying to generalise corr2_coeff for array-wise processing, and hence faster
    """
    corrs = np.ma.zeros([A.shape[0], A.shape[0]])
    corrs.mask = np.ones(corrs.shape)

    for i, raw_source in enumerate(A):
        raw_targets = np.ma.copy(A)

        raw_sources = np.tile(np.ma.copy(raw_source), (raw_targets.shape[0], 1))

        common_mask = np.ma.mask_or(raw_sources.mask, raw_targets.mask)

        # apply common mask
        sources = np.ma.array(raw_sources.data[:], mask=common_mask)
        targets = np.ma.array(raw_targets.data[:], mask=common_mask)

        good_overlap = np.ma.count(sources, axis=1)

        # get mean
        mS = np.ma.mean(sources, axis=1)
        mT = np.ma.mean(targets, axis=1)

        # make anomalies
        S_mS = sources - np.tile(mS[:, np.newaxis], (1, sources.shape[1]))
        T_mT = targets - np.tile(mT[:, np.newaxis], (1, targets.shape[1]))
        
        ssS = np.ma.sum(S_mS**2, axis=1)
        ssT = np.ma.sum(T_mT**2, axis=1)
        
        corrs[i, :] = np.ma.sum(S_mS * T_mT, axis=1) / (np.ma.sqrt(ssS) * np.ma.sqrt(ssT))

        corrs.mask[i, good_overlap < SUFFICIENT_YEARS] = True
            
    return corrs # corr2_coeff


#*******************************************
def separations_and_correlations(data, separations, names, flatten=True, diagnostics=False):
    """
    Return the separations and Pearson correlations ready to bin and plot

    :param array data: data to correlate (stations x years)
    :param array separations: the separation array
    :param array names: string array of station names being processed
    :param bool flatten: flatten the triangular arrays to 1-D
    :param bool diagnostics: output diagnostic information

    :returns: two 1-D arrays
    seps - separation
    cors - correlations

    """

    # loopwise code - also accounts for sufficient overlap period
    correlations = corr2_coeff(data)

    # have two matrices (seps & corrs)
    # need to obtain the upper triangle

    # and bin up.
    tril = np.tril_indices(correlations.shape[0])

    correlations.mask[tril] = True
    separations = np.ma.array(separations)
    separations.mask = correlations.mask # this also accounts for the overlap of common data

    if diagnostics:
        # check for correlations which are too good.
        spurious_corrs_i, spurious_corrs_j = np.ma.where(np.logical_and(correlations > 0.99, separations > utils.DEFAULT_DLS))

        for i, j in zip(spurious_corrs_i, spurious_corrs_j):
            print("{}-{} r = {}, sep = {}km".format(names[i], names[j], correlations[i, j], separations[i, j]))

    # flatten the triangular matrices
    if flatten:
        seps = np.reshape(separations, [-1]).compressed()
        cors = np.reshape(correlations, [-1]).compressed()

        return seps, cors # separations_and_correlations
    else:
        return separations, correlations # separations_and_correlations

#*******************************************
def exponential_curve(x, A, B, C=0): 
    """ 
    Quick calculation of exponential curve
    
    Input params are: x, A, B, C = 0

    :Returns: y = A * exp(-Bx) + C
    """

    if C == 0:
        return A * np.exp(-B * x)
    else:
        return A * np.exp(-B * x) + C # exponential_curve
        
#*******************************************
def exponential_inverted(y, A, B, C=0): 
    """ 
    Quick calculation of inverse of exponential curve to get x from y
    
    Input params are: y, A, B, C=0

    :Returns: ln((y - C)/A) / -B
    """
    return np.log((y - C)/A) / -B # x at value of y

#*******************************************
def exponential_fit(x, y, sigma, C=0):
    """
    Fit the exponential curve to extract the DLS.

    :param array x: bin centers
    :param array y: bin values
    :param array sigma: error in y values
    :param float C: non-zero asymptote

    :returns: dls & y-values for fitted curve
    """

    # need to have calculated a st.dev for the data to be used
    non_zero_y = y[sigma != 0]
    non_zero_x = x[sigma != 0]
    non_zero_s = sigma[sigma != 0]


    ctr = 0
    while x[ctr] != non_zero_x[0]:
        ctr += 1    

    MAXFEV = 4000

    # C set in file preamble
    
    min_loc = np.argmin(non_zero_y)
    
    # in case there are wiggles, use the minimum point to fit the initial drop
    initial_A = np.mean(non_zero_y[:2])

    # if initial few bins are empty (unlikely), then inflate the intercept 
    #   from the mean of the first few points given likely shape of curve
    if ctr != 0:
        initial_A = initial_A + (ctr * 0.1)

    initial_B = 1./500.
    if 1 < min_loc < int(0.9*len(non_zero_y)) and C == 0: # if not within 10% of end (need at least 2 points to fit with 2 parameters)

        popt, pcov = curve_fit(exponential_curve, non_zero_x[:int(min_loc*1.)], non_zero_y[:int(min_loc*1.)], \
                               p0=[initial_A, initial_B], maxfev=MAXFEV, absolute_sigma=True, sigma=non_zero_s[:int(min_loc*1.)])

        A, B = popt

    elif 2 < min_loc < int(0.9*len(non_zero_y)) and C != 0: # if not within 10% of end (need at least 3 points to fit with 3 parameters)
        initial_C = non_zero_y[int(min_loc*1.)]
        popt, pcov = curve_fit(exponential_curve, non_zero_x[:int(min_loc*1.)], non_zero_y[:int(min_loc*1.)], \
                               p0=[initial_A, initial_B, initial_C], maxfev=MAXFEV, absolute_sigma=True, sigma=non_zero_s[:int(min_loc*1.)])

        A, B, C = popt

    else:
        if C == 0:
            popt, pcov = curve_fit(exponential_curve, non_zero_x, non_zero_y, \
                                   p0=[initial_A, initial_B], maxfev=MAXFEV, absolute_sigma=True, sigma=non_zero_s)

            A, B = popt
        else:
            initial_C = non_zero_y[-1]
            popt, pcov = curve_fit(exponential_curve, non_zero_x, non_zero_y, \
                                   p0=[initial_A, initial_B, initial_C], maxfev=MAXFEV, absolute_sigma=True, sigma=non_zero_s)

            A, B, C = popt

    
    plot_curve = exponential_curve(x, *popt)

    # chi sq
    residuals = y - plot_curve
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r_squared = 1 - (ss_res / ss_tot)
    chisq = np.sum((residuals[sigma != 0])**2/sigma[sigma != 0]**2)
    print("R2 = {:6.4f}, chi-sq = {:12.6f}, dof = {:4.1f}".format(r_squared, chisq, len(non_zero_x)-len(popt)))

    # change zero for C if wanting to 
    dls = exponential_inverted(((A/np.exp(1.))+0), A, B, C) 

    return dls, plot_curve, chisq, r_squared # exponential_fit

#*******************************************
def extrap(interpolator):
    """
    http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate-give-an-extrapolated-result-beyond-the-input-range

    Linearly extrapolate outside of defined range

    :param: interpolator from scipy.interpolate import interp1d
    """

    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(list(map(pointwise, np.array(xs))))

    return ufunclike # extrap

#*******************************************
def assign_to_latitude_bands(stations):
    """
    Go through station list and see which latitude band they sit in

    :param array stations: array of station objects

    :returns: station_bands - array of which band
    """


    station_bands = np.zeros(len(stations))
    station_bands[:] = -1
    
    # spin through latitude bands
    for lb, band in enumerate(utils.LAT_BANDS):
        
        # and then all stations
        for s, stn in enumerate(stations):
            
            if station_bands[s] == -1: # only check unassigned stations
                if band[0] == -90.0:
                    if band[0] <= stn.latitude <= band[1]:
                        # account for station at south pole
                        station_bands[s] = lb
                else:
                    if band[0] < stn.latitude <= band[1]:
                        station_bands[s] = lb  
                
 
    return station_bands # assign_to_latitude_bands

#*******************************************
def get_separations(stations, all_locations):
    """
    Get the sepatation of stations with a list of locations (mainly, all other stations)

    :param array stations: array of station objects - for array size
    :param array all_locations: array of lat/lon locations

    :returns: separation & angle arrays
    """

    separation = np.zeros([len(stations), len(stations)]).astype(int)
    angle = np.zeros([len(stations), len(stations)]).astype(int)
    
    for l, loc in enumerate(all_locations):
        
        loc = np.repeat([loc], len(all_locations), axis=0)
        
        separation[l, :], angle[l, :] = utils.map_2_points(loc, all_locations, angle_int=True)

    return separation, angle # get_separations

#*******************************************
def read_data(stations, index, timescale, nyears, nmonths):
    """
    Read in all the station data into the list - and only keep years of interest

    :param array stations: array of station objects
    :param str index: index to process
    :param str timescale: timescale to process
    :param int nyears: number of years of data use
    :param int nmonths: number of months of data [legacy]

    :returns: Nstations x Nyears x Nmonths array    
    """


    # need to read in all the data to be able to process
    data = np.ma.zeros([len(stations), nyears, nmonths])
    data[:] = utils.HADEX_MDI
    data._FillValue = utils.HADEX_MDI

    for s, stn in enumerate(stations):
        
        # allow for stations which might not exist - usually annual if monthly present
        if os.path.exists(os.path.join(stn.location, stn.id, "{}_{}_{}.csv".format(stn.id, index.lower(), timescale))):

            times, indata = utils.read_station_index(stn, index, timescale)

            match = np.in1d(utils.REFERENCEYEARS, times)
            match_back = np.in1d(times, utils.REFERENCEYEARS)

            if len(indata.shape) == 1:
                indata = indata.reshape([indata.shape[0], 1])
            data[s, match, :] = indata[match_back] # store all the info
        else:
            print("File missing - skipped:")
            print("        {}".format(os.path.join(stn.location, stn.id, "{}_{}_{}.csv".format(stn.id, index.lower(), timescale))))

    return data # read_data

#*******************************************
def get_all_data(stations, index, timescale, nyears, nmonths):
    """
    Reads in all the data across timescales if appropriate

    :param array stations: array of station objects
    :param str index: index to process
    :param str timescale: timescale to process
    :param int nyears: number of years of data use
    :param int nmonths: number of months of data [legacy]

    :returns: Nstations x Nyears x Nmonths array
    """

    print("reading in all data")

    # make sure annual and monthly read if appropriate
    if index in utils.MONTHLY_INDICES:

        all_monthly_data = read_data(stations, index, timescale, nyears, nmonths-1)
        all_annual_data = read_data(stations, index, "ANN", nyears, 1)

        # need to read in all the data to be able to process
        all_data = np.ma.zeros([len(stations), nyears, nmonths])
        all_data[:] = utils.HADEX_MDI
        all_data._FillValue = utils.HADEX_MDI

        all_data[:, :, 1:] = all_monthly_data
        all_data[:, :, 0] = np.squeeze(all_annual_data) # annual at index 0

    else:
        all_data = read_data(stations, index, timescale, nyears, nmonths)

    # and apply the mask
    all_data = np.ma.masked_where(all_data <= -99.9, all_data)

    return all_data # get_all_data

#*******************************************
def interpolate_dls_to_grid(dls, nmonths):
    """
    Interpolate DLS calculated for latitude bands onto grid box centres

    :param array dls: array of DLS values [bands x months]
    :param int nmonths: number of months [legacy]

    :returns: array of DLS on grid [box_centre_lats x months]
    """

    band_centre_latitudes = np.mean(utils.LAT_BANDS, axis=1)
    grid_dls = np.zeros([len(utils.box_centre_lats), nmonths])

    for month in range(nmonths):           
        # scipy.interpolate function interp1d
        f_interp = interp1d(band_centre_latitudes, dls[:, month])

        # extrap is local subroutine, above
        f_extrap = extrap(f_interp)

        grid_dls[:, month] = f_extrap(utils.box_centre_lats)

    # replace those dls < default with default and > max with max
        
    grid_dls[grid_dls < utils.DEFAULT_DLS] = utils.DEFAULT_DLS
    grid_dls[grid_dls > utils.MAX_DLS] = utils.MAX_DLS

    return grid_dls # interpolate_dls_to_grid

#*******************************************
def write_dls_file(filename, grid_dls, nmonths, month_names):
    """
    Write the DLS output file

    :param str filename: output filename
    :param array grid_dls: array of [lat x months] DLS values
    :param int nmonths: number of months
    :param str month_names: names of months to use
    """

    # write header and then list of latitude and DLS values
    with open(filename, "w") as outfile:

        outfile.write("  {} monthly/annual values,\n".format(nmonths))
        outfile.write("  {} latitude values,\n".format(len(utils.box_centre_lats)))
        outfile.write("  {} is default/minimum dls value.\n".format(utils.DEFAULT_DLS))
        if nmonths == 13:
            outfile.write("{:8s} {}\n".format("", " ".join("{:>5s}".format(x) for x in month_names)))
        else:
            outfile.write("{:8s} {:>5s}\n".format("", "Ann"))

        for l, lat in enumerate(utils.box_centre_lats):
            outfile.write("{:8.4f} {}\n".format(lat, " ".join("{:5.0f}".format(x) for x in grid_dls[l, :])))

        outfile.close()

    return # write_dls_file

#*********************************************
def main(index="TX90p", diagnostics=False, qc_flags=""):
    """
    The main DLS function

    :param str index: which index to run
    :param bool diagnostics: extra verbose output
    :param str qc_flags: which QC flags to process W, B, A, N, C, R, F, E, V, M
    """

    if index in utils.MONTHLY_INDICES:
        nmonths = 13
        timescale = "MON"
    else:
        nmonths = 1
        timescale = "ANN"

    # move this up one level eventually?
    all_datasets = utils.get_input_datasets()

    # spin through all datasets
    stations = np.array([])
    for dataset in all_datasets:

        try:
            ds_stations = utils.read_inventory(dataset, subdir="formatted/indices", final=True, timescale=timescale, index=index, qc_flags=qc_flags)
            good_stations = utils.select_qc_passes(ds_stations, qc_flags=qc_flags)
            
            stations = np.append(stations, good_stations)

            print("Adding {} ({} stations), nstations = {}".format(dataset.name, len(good_stations), len(stations)))

        except IOError:
            # file missing
            print("No stations with data for {}".format(dataset.name))
            
    nstations = len(stations)

    # array of lats and lons for calculation of separations
    all_locations = np.array([[stn.latitude, stn.longitude] for stn in stations])

    # get the separations (km, radians)
    stn_separation, stn_angle = get_separations(stations, all_locations)

    # assign stations to bands
    StationBands = assign_to_latitude_bands(stations)

    # read in all the station data
    all_data = get_all_data(stations, index, timescale, nyears, nmonths)

    # set up the DLS defaults
    bins = np.arange(0, MAX_SEPARATION + BIN_WIDTH, BIN_WIDTH)

    all_dls = np.zeros([len(utils.LAT_BANDS), nmonths])
    all_dls[:] = utils.DEFAULT_DLS


    # now spin through all latitude bands and months.
    for lb, band in enumerate(utils.LAT_BANDS):

        stations_in_bands, = np.where(StationBands == lb)

        if len(stations_in_bands) <= 30:
            # insufficient stations within this latitude band, next band
            if diagnostics:
                print("Index {}, Band {} to {}".format(index, band[0], band[1]))
                print("Number of stations {}".format(len(stations_in_bands)))
            print("Ann, Jan -- Dec, DLS = {} km".format(utils.DEFAULT_DLS))

            # spin through months to remove old plots if they exist
            for month in range(nmonths):
                if os.path.exists(os.path.join(utils.PLOTLOCS, "DLS", "DLS_{}_{}_{}to{}.png".format(index, month_names[month], band[0], band[1]))):
                    os.remove(os.path.join(utils.PLOTLOCS, "DLS", "DLS_{}_{}_{}to{}.png".format(utils.PLOTLOCS, index, month_names[month], band[0], band[1])))
            continue

        print("{}, # stations {}".format(band, len(stations_in_bands)))

        # process each month
        for month in range(nmonths):
            print(month_names[month])

            month_data = all_data[stations_in_bands, :, month]

            names = [s.id for s in stations[stations_in_bands]]

            # get the separation and correlation for each cross pair
            # correlations only from 1951 (match HadEX2)
            cor_yr = 1951-utils.STARTYEAR.year 
            seps, cors = separations_and_correlations(month_data[:, cor_yr:], stn_separation[stations_in_bands, :][:, stations_in_bands], names, diagnostics=diagnostics)

            if len(seps) == 0 and len(cors) == 0:
                # then none of the available stations either had sufficient overlapping data
                #  or values at that particular point (correlations of lots of zeros doesn't mean anything)
                #  so escape and go on to next month
                if diagnostics:
                    print("Index {}, Band {} to {}, month {}".format(index, band[0], band[1], month_names[month]))
                    print("Number of stations {}".format(len(stations_in_bands)))
                    print("Likely that all values for this index, month and band are zero\n hence correlations don't mean anything")
                    print("Using default DLS = {}km".format(utils.DEFAULT_DLS))
                else:
                    print("No stations, {} - {} DLS = {} km".format(band, month_names[month], utils.DEFAULT_DLS))
                continue

            # get the bins
            bin_assignment = np.digitize(seps, bins, right=True) # "right" means left bin edge included
            bin_centers = bins - BIN_WIDTH/2.

            # average value for each bin if sufficient correlations to do so.
            means = np.zeros(len(bins))
            sigmas = np.zeros(len(bins))
            for b, bin in enumerate(bins):
                locs, = np.where(bin_assignment == b)

                if len(locs) > MIN_PER_BIN:
    #                means[b] = np.ma.mean(cors[locs])
                    means[b] = np.ma.median(cors[locs])
                    sigmas[b] = np.ma.std(cors[locs])
    #                print(bin, means[b], len(locs), cors[locs])
    #                raw_input("stop")


            filled_bins, = np.where(means != 0)

            # if sufficient bins are filled then fit the curve
            if len(filled_bins)/float(len(bins)) >= 0.5:

                if utils.FIX_ZERO:
                    # fix zero bin to be 1.0, and use bin edges, not centres (HadEX2)
                    means[0] = 1.
                    sigmas[0] = sigmas[1]
                    dls, plot_curve, chisq, R2 = exponential_fit(bins, means, sigmas, C=C)
                else:
                    dls, plot_curve, chisq, R2 = exponential_fit(bin_centers[1:], means[1:], sigmas[1:], C=C)

                # only take fit if greater than minimum set overall
                all_dls[lb, month] = np.max([dls, utils.DEFAULT_DLS])
                
                # test at 5% level and 2 or 3 dofs, as per HadEX2
                if utils.FIX_ZERO and chisq >= chi2.isf(0.05, len(bins[sigmas != 0])-2):
                    print("inadequately good fit")
                    all_dls[lb, month] = utils.DEFAULT_DLS
                elif chisq >= chi2.isf(0.05, len(bins[sigmas != 0])-3):
                    print("inadequately good fit")
                    all_dls[lb, month] = utils.DEFAULT_DLS

                # plot the fit if required
                plt.clf()
                plt.scatter(seps, cors, c='b', marker='.', alpha=0.1, edgecolor=None)

                # calculate the 2D density of the data given
                counts, xbins, ybins = np.histogram2d(seps, cors, bins=50)

                # make the contour plot (5 levels)
                plt.contour(counts.transpose(), 5, extent=[xbins.min(), xbins.max(),
                                                           ybins.min(), ybins.max()], linewidths=1, colors='black',
                            linestyles='solid')

                if utils.FIX_ZERO:
                    plt.plot(bins[sigmas != 0], means[sigmas != 0], 'ro')
                    plt.errorbar(bins[sigmas != 0], means[sigmas != 0], yerr=sigmas[sigmas != 0], fmt="none", ecolor="r")
                    plt.plot(bins, plot_curve, c='cyan', ls='-', lw=2)
                else:                      
                    plt.plot(bin_centers[1:][sigmas[1:] != 0], means[1:][sigmas[1:] != 0], 'ro')
                    plt.errorbar(bin_centers[1:][sigmas[1:] != 0], means[1:][sigmas[1:] != 0], yerr=sigmas[1:][sigmas[1:] != 0], fmt="none", ecolor="r")
                    plt.plot(bin_centers[1:], plot_curve, c='cyan', ls='-', lw=2) # plot curve will have been truncated

                plt.axvline(dls, c='magenta', ls="--", lw=2)
                plt.axvline(utils.DEFAULT_DLS, c='k', ls=":", lw=1)
                plt.axvline(utils.MAX_DLS, c='k', ls=":", lw=1)
                plt.text(dls+10, 0.95, "dls = {:4.0f}km".format(dls))
                plt.text(3010, -0.95, "r2 = {:6.4f}".format(R2))
                plt.text(3010, -0.85, "chi2 = {:6.4f}".format(chisq))
                plt.text(3010, -0.75, "Nstat = {}".format(len(stations_in_bands)))
                plt.xlim([-100, 5000])
                plt.ylim([-1, None])

                plt.xlabel("Separation (km)")
                plt.ylabel("Correlation")
                plt.title("{} - {}; {} to {}".format(index, month_names[month], band[0], band[1]))

                # add text to show what code created this and when
                if utils.WATERMARK:
                    watermarkstring = "/".join(os.getcwd().split('/')[4:])+'/' + os.path.basename(__file__)+"   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
                    plt.figtext(0.01, 0.01, watermarkstring, size=6)

                if utils.FIX_ZERO:
                    plt.savefig(os.path.join(utils.PLOTLOCS, "DLS", \
                                             "DLS_{}_{}_{}_{}to{}_fixzero.png".format(index, "{}-{}".format(str(utils.REF_START)[-2:], \
                                                                                                    str(utils.REF_END)[-2:]), month_names[month], band[0], band[1])), dpi=300)
                else:
                    plt.savefig(os.path.join(utils.PLOTLOCS, "DLS", \
                                             "DLS_{}_{}_{}_{}to{}.png".format(index, "{}-{}".format(str(utils.REF_START)[-2:], \
                                                                                                    str(utils.REF_END)[-2:]), month_names[month], band[0], band[1])), dpi=300)
                    

                print("DLS = {:7.2f} km".format(dls))

            else:
                print("insufficient bins for fit ({}/{})".format(len(filled_bins), float(len(bins))))
                if os.path.exists(os.path.join(utils.PLOTLOCS, "DLS", \
                                             "DLS_{}_{}_{}_{}to{}.png".format(index, "{}-{}".format(str(utils.REF_START)[-2:], \
                                                                                                    str(utils.REF_END)[-2:]), month_names[month], band[0], band[1]))):
                    # remove old plots!
                    os.remove(os.path.join(utils.PLOTLOCS, "DLS", \
                                             "DLS_{}_{}_{}_{}to{}.png".format(index, "{}-{}".format(str(utils.REF_START)[-2:], \
                                                                                                    str(utils.REF_END)[-2:]), month_names[month], band[0], band[1])))

    # replace those dls < default with default and > max with max
    all_dls[all_dls < utils.DEFAULT_DLS] = utils.DEFAULT_DLS

    # interpolate
    grid_dls = interpolate_dls_to_grid(all_dls, nmonths)

    # write output file
    write_dls_file(os.path.join(utils.DLSLOCS, "dls_{}.txt".format(index)), grid_dls, nmonths, month_names)

    return # main
    

#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', dest='index', action='store', default="TX90p",
                        help='Which index to run')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default=False')
    parser.add_argument('--qc_flags', dest='qc_flags', action='store', default="",
                        help='Which QC flags to use when filtering stations, default=""')

    args = parser.parse_args()          

    main(index=args.index, diagnostics=args.diagnostics, qc_flags=args.qc_flags)


#*******************************************
# END
#*******************************************
