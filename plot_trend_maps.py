#!/usr/bin/ python
#***************************************
#
#   Plot the trend maps for the globe.
#
#   21 May 2018 RJHD
#***************************************
#************************************************************************
#                    SVN Info
# $Rev:: 473                                           $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2020-02-27 09:54:10 +0000 (Thu, 27 Feb 2020) $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
"""
Plots trend maps

plot_trend_maps.py invoked by typing::

  python plot_trend_maps.py --index "TX90p" --anomalies --grid --diagnostics

Input arguments:

--index         Which ETCCDI index to use

--anomalies     Run on anomaly outputs

--grid          Set the gridding method (ADW)

--diagnostics   Output extra info (default False)
"""
import os
import calendar
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import iris
import iris.plot
import iris.analysis.stats

# local utilities scripts
import utils
import plt_utils as putils

FLTMDI = -1.e30

# set up calendar month names
month_names = calendar.month_abbr[:]
month_names[0] = "Ann"


#************************************************************************
def GetYears(cube):

    print("this might not give correct results in Py3 - test")
    times = cube.coord('time').points
    years = np.round(np.array([(t - 101)//10000 for t in times]))

    return years

#************************************************************************
def ApplyClimatology(cube):

    years = GetYears(cube)

    clim_years = np.where((years >= utils.CLIM_START) & (years <= utils.CLIM_END))
    
    climatology = np.ma.mean(cube.data[clim_years], axis=0)
    
    cube.data = cube.data - climatology                                
        
    return cube

#***************************************
def TrendingCalculation(cube, verbose=False, plots=False):
    '''
    Calculate trend for index for each gridbox
    return the range, significance and the detrended data

    :param array data: masked array
    :param bool verbose: extra print to screen
    :param bool plots: extra plots of the trend calculation for each gridbox
    '''

    nyears, nlat, nlon = cube.data.shape

    # convert YYYYMMDD to YYYY
    years = np.arange(nyears) + utils.TREND_START

    # set up arrays
    trend = np.zeros([nlat, nlon])
    trend.fill(FLTMDI)
    sigma = np.zeros([nlat, nlon])
    sigma.fill(FLTMDI)
    significance = np.zeros([nlat, nlon])
    significance.fill(FLTMDI)

    detrended_data = np.zeros(cube.data.shape)
    detrended_data.fill(FLTMDI)

    # run through each grid box
    for t, lat in enumerate(range(nlat)):
        if verbose: print("%i/%i" % (t, nlat))
        for n, lon in enumerate(range(nlon)):
            # get the timeseries for this gridbox
            gridbox = cube.data[:, lat, lon]
         
            # need to do completeness check with spread through decades
            if len(gridbox.compressed()) >= np.floor(utils.TREND_COMPLETENESS * nyears):
                # if less than 66% complete, don't calculate trend

                if np.ma.array(years, mask=gridbox.mask).compressed()[-1] < utils.TREND_FINAL_YEAR:
                    # if final year is too early, then don't calculate trend
                    continue

                # median pairwise slopes
                mpw, mpw_l, mpw_u = utils.MedianPairwiseSlopes(years, gridbox.filled(FLTMDI), FLTMDI)

                # get trend, sigma and significance
                trend[lat, lon] = mpw

                sigma[lat, lon] = (mpw_u - mpw_l) / 2. / 2.  # 95% ~ 2 sigma full width

                if ((mpw > 0) and (mpw_l > 0)) or ((mpw < 0 and mpw_u < 0)):
                    # "significant" as both slope is clearly different from zero
                    significance[lat, lon] = 1
                else:
                    significance[lat, lon] = 0

                # detrend

                # find point on line
                ypivot = np.mean(gridbox.compressed())
                xpivot = np.mean(years[np.where(gridbox.mask != True)])

                # get linear trend values at each year
                lineary = [mpw * (yr - xpivot) + ypivot for yr in years[np.where(gridbox.mask != True)]]


                detrended = gridbox.compressed() - lineary

                if plots:
                    
                    plt.clf()
                    plt.plot(years[np.where(gridbox.mask != True)], gridbox.compressed(), 'bo')
                    plt.plot(years[np.where(gridbox.mask != True)], lineary, 'r-')
                    plt.show()

                detrended_data[np.where(gridbox.mask != True), lat, lon] = detrended

    trend = np.ma.masked_where(trend == FLTMDI, trend)*10 # to get per decade
    sigma = np.ma.masked_where(sigma == FLTMDI, sigma)
    significance = np.ma.masked_where(significance == FLTMDI, significance)

    trend = putils.make_iris_cube_2d(trend, cube.coord("grid_latitude").points, cube.coord("grid_longitude").points, "trend", "")
    sigma = putils.make_iris_cube_2d(sigma, cube.coord("grid_latitude").points, cube.coord("grid_longitude").points, "sigma", "")
    significance = putils.make_iris_cube_2d(significance, cube.coord("grid_latitude").points, cube.coord("grid_longitude").points, "trend_significance", "")
#    detrended_data = putils.make_iris_cube_2d(detrended_data, cube.coord("latitude").points, cube.coord("longitude").points, "detrended", "")

    return trend, sigma, significance#, detrended_data # Trending Calculation


#*********************************************************
def periodConstraint(cube, t1, lower=True):
    """
    Return constraint that cube times are >= t1 (or >= t1 if lower=False)
    """
    if lower:
        return iris.Constraint(time=lambda cell: cell.point >= t1) # periodConstraint
    else:
        return iris.Constraint(time=lambda cell: cell.point <= t1) # periodConstraint

#************************************************************************

#*********************************************************
def main(index, diagnostics=False, anomalies="None", grid="ADW"):
    """
    Plot maps of linear trends

    :param str index: which index to run
    :param bool diagnostics: output diagnostic information
    :param str anomalies: run code on anomalies or climatology rather than raw data
    :param str grid: gridding type ADW/CAM
    """


    # get details of index
    index = utils.INDICES[index] 

    # allow for option of running through each month
    if index.name in utils.MONTHLY_INDICES:
        nmonths = 13
    else:
        nmonths = 1

    # sort the colour maps
    RdYlBu, RdYlBu_r = putils.adjust_RdYlBu()
    BrBG, BrBG_r = putils.make_BrBG()

    cube_list = iris.load(os.path.join(utils.FINALROOT, utils.make_filenames(index=index.name, grid=grid, anomalies=anomalies, extra="", month_index="")))    

    names = np.array([cube.var_name for cube in cube_list])

    #*************
    # plot trend map

    for month in range(nmonths):

        if month == 0:
            # annual
            if anomalies != "climatology":
                if index.name in ["TX90p", "TN90p", "SU", "TR", "GSL"]:
                    bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
                    cmap = RdYlBu_r
                elif index.name in ["DTR", "ETR"]:
                    bounds = [-100, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 100]
                    cmap = RdYlBu_r
                elif index.name in ["WSDI"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu_r
                elif index.name in ["TX10p", "TN10p", "FD", "ID"]:
                    bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
                    cmap = RdYlBu
                elif index.name in ["CSDI"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu
                elif index.name in ["TXn", "TNn"]:
                    bounds = [-100, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 100]
                    cmap = RdYlBu_r
                elif index.name in ["TXx", "TNx"]:
                    bounds = [-100, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 100]
                    cmap = RdYlBu_r
                elif index.name in ["Rx1day"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = BrBG
                elif index.name in ["Rx5day"]:
                    bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
                    cmap = BrBG
                elif index.name in ["PRCPTOT"]:
                    bounds = [-100, -20, -10, -5, -2, 0, 2, 5, 10, 20, 100]
                    cmap = BrBG
                elif index.name in ["Rnnmm", "R95p", "R99p"]:
                    bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
                    cmap = BrBG
                elif index.name in ["R95pTOT"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = BrBG
                elif index.name in ["R99pTOT"]:
                    bounds = [-100, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 100]
                    cmap = BrBG
                elif index.name in ["R10mm"]:
                    bounds = [-100, -3, -1.5, -0.75, -0.25, 0, 0.25, 0.75, 1.5, 3, 100]
                    cmap = BrBG
                elif index.name in ["R20mm"]:
                    bounds = [-100, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 100]
                    cmap = BrBG
                elif index.name in ["CWD"]:
                    bounds = [-100, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 100]
                    cmap = BrBG
                elif index.name in ["SDII"]:
                    bounds = [-100, -0.75, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 0.75, 100]
                    cmap = BrBG
                elif index.name in ["CDD"]:
                    bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
                    cmap = BrBG_r
                elif index.name in ["CDDcold18"]:
                    bounds = [-10000, -100, -50, -20, -10, 0, 10, 20, 50, 100, 10000]
                    cmap = RdYlBu_r
                elif index.name in ["HDDheat18"]:
                    bounds = [-10000, -800, -400, -200, -100, 0, 100, 200, 400, 800, 10000]
                    cmap = RdYlBu
                elif index.name in ["GDDgrow10"]:
                    bounds = [-10000, -400, -200, -100, -50, 0, 50, 100, 200, 400, 10000]
                    cmap = RdYlBu
                elif index.name in ["WSDI3"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu_r
                elif index.name in ["CSDI3"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu
                elif index.name in ["TNlt2", "TNltm2", "TNltm20", "TMlt10", "TMlt5"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu
                elif index.name in ["TXge30", "TXge35", "TMge5", "TMge10", "TXge50p"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu_r
                elif index.name in ["TNm", "TXm", "TMm", "TXTN"]:
                    bounds = [-100, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 100]
                    cmap = RdYlBu_r
                elif index.name in ["TXbTNb"]:
                    bounds = [-100, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 100]
                    cmap = RdYlBu
                elif index.name in ["RXday"]:
                    bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
                    cmap = BrBG
                else:
                    bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
                    cmap = RdYlBu_r
                   

            else:
                if index.name in ["TX90p", "TN90p", "WSDI", "SU", "TR", "GSL", "DTR", "ETR"]:
                    bounds = np.arange(0, 60, 5)
                    cmap = plt.cm.YlOrRd
                elif index.name in ["TX10p", "TN10p", "CSDI", "FD", "ID"]:
                    bounds = np.arange(0, 60, 5)
                    cmap = plt.cm.YlOrRd_r
                elif index.name in ["TXx", "TNx"]:
                    bounds = np.arange(-10, 40, 5)
                    cmap = plt.cm.YlOrRd
                elif index.name in ["TXn", "TNn"]:
                    bounds = np.arange(-30, 10, 5)
                    cmap = plt.cm.YlOrRd
                elif index.name in ["Rx1day", "Rx5day"]:
                    bounds = np.arange(0, 100, 10)
                    cmap = plt.cm.YlGnBu
                elif index.name in ["PRCPTOT"]:
                    bounds = np.arange(0, 1000, 100)
                    cmap = plt.cm.YlGnBu
                elif index.name in ["CWD", "R20mm", "R10mm", "Rnnmm", "R95pTOT", "R99pTOT", "R95p", "R99p", "SDII"]:
                    bounds = np.arange(0, 200, 20)
                    cmap = plt.cm.YlGnBu
                elif index.name in ["CDD"]:
                    bounds = np.arange(0, 200, 20)
                    cmap = plt.cm.YlGnBu_r
                else:
                    bounds = np.arange(0, 60, 5)
                    cmap = plt.cm.YlOrRd

        else:
            # monthly
            if anomalies != "climatology":
                if index.name in ["TX90p", "TN90p", "SU", "TR", "GSL", "DTR", "ETR"]:
                    bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
                    cmap = RdYlBu_r
                elif index.name in ["WSDI"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu_r
                elif index.name in ["TX10p", "TN10p", "FD", "ID"]:
                    bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
                    cmap = RdYlBu
                elif index.name in ["CSDI"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu
                elif index.name in ["TXx", "TNx", "TXn", "TNn"]:
                    bounds = [-100, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 100]
                    cmap = RdYlBu_r
                elif index.name in ["Rx1day", "Rx5day"]:
                    bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
                    cmap = BrBG
                elif index.name in ["PRCPTOT"]:
                    bounds = [-100, -10, -5, -2, -1, 0, 1, 2, 5, 10, 100]
                    cmap = BrBG
                elif index.name in ["Rnnmm", "R95pTOT", "R99pTOT", "R95p", "R99p"]:
                    bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
                    cmap = BrBG
                elif index.name in ["R20mm", "R10mm"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = BrBG
                elif index.name in ["CWD"]:
                    bounds = [-100, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 100]
                    cmap = BrBG
                elif index.name in ["SDII"]:
                    bounds = [-100, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 100]
                    cmap = BrBG
                elif index.name in ["CDD"]:
                    bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
                    cmap = BrBG_r
                elif index.name in ["CDDcold18"]:
                    bounds = [-100, -20, -15, -10, -5, 0, 5, 10, 15, 20, 100]
                    cmap = RdYlBu_r
                elif index.name in ["HDDheat18"]:
                    bounds = [-10000, -800, -400, -200, -100, 0, 100, 200, 400, 800, 10000]
                    cmap = RdYlBu
                elif index.name in ["GDDgrow10"]:
                    bounds = [-10000, -400, -200, -100, -50, 0, 50, 100, 200, 400, 10000]
                    cmap = RdYlBu
                elif index.name in ["WSDI3"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu_r
                elif index.name in ["CSDI3"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu
                elif index.name in ["TNlt2", "TNltm2", "TNltm20", "TMlt10", "TMlt5"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu
                elif index.name in ["TXge30", "TXge35", "TMge5", "TMge10", "TXge50p"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = RdYlBu_r
                elif index.name in ["TNm", "TXm", "TMm", "TXTN"]:
                    bounds = [-100, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 100]
                    cmap = RdYlBu_r
                elif index.name in ["TXbTNb"]:
                    bounds = [-100, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 100]
                    cmap = RdYlBu
                elif index.name in ["RXday"]:
                    bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
                    cmap = BrBG
                elif index.name in ["24month_SPI", "12month_SPI", "6month_SPI", "3month_SPI", "24month_SPEI", "12month_SPEI", "6month_SPEI", "3month_SPEI"]:
                    bounds = [-100, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 100]
                    cmap = BrBG
                else:
                    bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
                    cmap = RdYlBu_r

            else:
                if index.name in ["TX90p", "TN90p", "WSDI", "SU", "TR", "GSL", "DTR", "ETR"]:
                    bounds = np.arange(0, 60, 5)
                    cmap = plt.cm.YlOrRd
                elif index.name in ["TX10p", "TN10p", "CSDI", "FD", "ID"]:
                    bounds = np.arange(0, 60, 5)
                    cmap = plt.cm.YlOrRd_r
                elif index.name in ["TXx", "TNx"]:
                    bounds = np.arange(-10, 40, 5)
                    cmap = plt.cm.YlOrRd
                elif index.name in ["TXn", "TNn"]:
                    bounds = np.arange(-30, 10, 5)
                    cmap = plt.cm.YlOrRd
                elif index.name in ["Rx1day", "Rx5day"]:
                    bounds = np.arange(0, 100, 10)
                    cmap = plt.cm.YlGnBu
                elif index.name in ["PRCPTOT"]:
                    bounds = np.arange(0, 1000, 100)
                    cmap = plt.cm.YlGnBu
                elif index.name in ["CWD", "R20mm", "R10mm", "Rnnmm", "R95pTOT", "R99pTOT", "R95p", "R99p", "SDII"]:
                    bounds = np.arange(0, 200, 20)
                    cmap = plt.cm.YlGnBu
                elif index.name in ["CDD"]:
                    bounds = np.arange(0, 200, 20)
                    cmap = plt.cm.YlGnBu_r
                else:
                    bounds = np.arange(0, 60, 5)
                    cmap = plt.cm.YlOrRd

        selected_cube, = np.where(names == month_names[month])

        cube = cube_list[selected_cube[0]]
        try:
            cube.coord('grid_latitude').guess_bounds()
            cube.coord('grid_longitude').guess_bounds()  
        except ValueError:
            pass

        # fix percent -> days issue for these four
        if index.name in ["TX90p", "TN90p", "TX10p", "TN10p"]:
            cube.data = cube.data * 3.65
            index.units = "days"

        # get recent period and trend
        if anomalies != "climatology":
            postYYYY = periodConstraint(cube, utils.TREND_START)
            cube = cube.extract(postYYYY)
            preYYYY = periodConstraint(cube, utils.TREND_END, lower=False)
            cube = cube.extract(preYYYY)
            trend_cube, sigma, significance = TrendingCalculation(cube, verbose=diagnostics)


        if index.units == "degrees_C":
            units = '$^{\circ}$'+"C"
        else:
            units = index.units


        if anomalies != "climatology":
            outname = putils.make_filenames("trend", index=index.name, grid=grid, anomalies=anomalies, month=month_names[month])

            putils.plot_smooth_map_iris("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname), trend_cube, cmap, bounds, "Trend ({}/10 year)".format(units), title="{} - {}, Linear Trend {}-{}".format(index.name, month_names[month], utils.TREND_START, utils.TREND_END), figtext="(a)", significance=significance)

        else:
            outname = putils.make_filenames("climatology", index=index.name, grid=grid, anomalies=anomalies, month=month_names[month])

            putils.plot_smooth_map_iris("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname), cube[0], cmap, bounds, "{}".format(units), title="{} - {}, Climatology {}-{}".format(index.name, month_names[month], utils.CLIM_START.year, utils.CLIM_END.year), figtext="(a)")

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
    parser.add_argument('--anomalies', dest='anomalies', action='store', default="None",
                        help='To use the anomalies or climatology rather than actuals, default="None"') 
    parser.add_argument('--grid', dest='grid', action='store', default="ADW",
                        help='gridding routine to use (ADW or CAM) default="ADW"') 

    args = parser.parse_args()
    
    main(index=args.index, diagnostics=args.diagnostics, anomalies=args.anomalies, grid=args.grid)
    

#*******************************************
# END
#*******************************************
