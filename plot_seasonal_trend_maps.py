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
Plots seasonal trend maps

plot_seasonal_trend_maps.py invoked by typing::

  python plot_seasonal_trend_maps.py --index "TX90p" --anomalies --grid --diagnostics

Input arguments:

--index         Which ETCCDI index to use

--anomalies     Run on anomaly outputs

--grid          Set the gridding method (ADW)

--diagnostics   Output extra info (default False)
"""
import os
import calendar
import copy
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
month_names = np.array(calendar.month_abbr[:])
month_names[0] = "Ann"

SEASON_DICT = {"MAM":["Mar", "Apr", "May"], "JJA":["Jun", "Jul", "Aug"], \
                   "SON":["Sep", "Oct", "Nov"], "DJF":["Dec", "Jan", "Feb"]}
SEASONS = ["DJF", "MAM", "JJA", "SON"]

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
def TrendingCalculation(cube, verbose=False):
    '''
    Calculate trend for index for each gridbox
    return the range, significance and the detrended data

    :param array data: masked array
    :param bool verbose: extra print to screen
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
            if len(gridbox.compressed()) > (utils.TREND_COMPLETENESS * nyears):
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
                if verbose:
                    
                    plt.clf()
                    plt.plot(years[np.where(gridbox.mask != True)], gridbox.compressed(), 'bo')
                    plt.plot(years[np.where(gridbox.mask != True)], lineary, 'r-')
                    plt.show()

                detrended_data[np.where(gridbox.mask != True), lat, lon] = detrended

    trend = np.ma.masked_where(trend == FLTMDI, trend)*10 # to get per decade
    sigma = np.ma.masked_where(sigma == FLTMDI, trend)
    significance = np.ma.masked_where(significance == FLTMDI, trend)

    trend = putils.make_iris_cube_2d(trend, cube.coord("grid_latitude").points, cube.coord("grid_longitude").points, "trend", "")
    sigma = putils.make_iris_cube_2d(sigma, cube.coord("grid_latitude").points, cube.coord("grid_longitude").points, "sigma", "")
    significance = putils.make_iris_cube_2d(significance, cube.coord("grid_latitude").points, cube.coord("grid_longitude").points, "trend_significance", "")
#    detrended_data = putils.make_iris_cube_2d(detrended_data, cube.coord("latitude").points, cube.coord("longitude").points, "detrended", "")

    return trend, sigma, significance#, detrended_data # Trending Calculation


#*********************************************************
def periodConstraint(cube, t1):
    """
    Return constraint that cube times are >= t1
    """

    return iris.Constraint(time=lambda cell: cell.point >= t1) # periodConstraint

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

    # sort the colour maps
    RdYlBu, RdYlBu_r = putils.adjust_RdYlBu()
    BrBG, BrBG_r = putils.make_BrBG()

    cube_list = iris.load(os.path.join(utils.FINALROOT, utils.make_filenames(index=index.name, grid=grid, anomalies=anomalies, extra="", month_index="")))    

    names = np.array([cube.var_name for cube in cube_list])

    #*************
    # plot trend map

    for season in SEASONS:

        three_month_data = []
        months = SEASON_DICT[season]
        for month in months:

            if anomalies != "climatology":
                if index.name in ["TX90p", "TN90p", "SU", "TR"]:
                    bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
                    cmap = RdYlBu_r
                elif index.name in ["TX10p", "TN10p", "FD", "ID"]:
                    bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
                    cmap = RdYlBu
                elif index.name in ["DTR", "ETR"]:
                    bounds = [-100, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 100]
                    cmap = RdYlBu_r
                elif index.name in ["TXx", "TNx", "TXn", "TNn"]:
                    bounds = [-100, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 100]
                    cmap = RdYlBu_r
                elif index.name in ["Rx1day", "Rx5day"]:
                    bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
                    cmap = BrBG
                elif index.name in ["CWD"]:
                    bounds = [-100, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 100]
                    cmap = BrBG
                elif index.name in ["CDD", "PRCPTOT"]:
                    bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
                    cmap = BrBG_r
                elif index.name in ["R10mm", "R20mm"]:
                    bounds = [-100, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 100]
                    cmap = BrBG
                else:
                    bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
                    cmap = RdYlBu_r

            else:
                if index.name in ["TX90p", "TN90p", "SU", "TR"]:
                    bounds = np.arange(0, 60, 5)
                    cmap = plt.cm.YlOrRd
                elif index.name in ["TX10p", "TN10p", "FD", "ID"]:
                    bounds = np.arange(0, 60, 5)
                    cmap = plt.cm.YlOrRd_r
                elif index.name in ["DTR", "ETR"]:
                    bounds = np.arange(0, 30, 5)
                    cmap = plt.cm.YlOrRd
                elif index.name in ["TXx", "TNx"]:
                    bounds = np.arange(-10, 40, 5)
                    cmap = plt.cm.YlOrRd
                elif index.name in ["TXn", "TNn"]:
                    bounds = np.arange(-30, 10, 5)
                    cmap = plt.cm.YlOrRd
                elif index.name in ["Rx1day", "Rx5day"]:
                    bounds = np.arange(0, 100, 10)
                    cmap = plt.cm.YlGnBu
                elif index.name in ["CWD"]:
                    bounds = np.arange(0, 100, 10)
                    cmap = BrBG
                elif index.name in ["CDD", "PRCPTOT", "R10mm", "R20mm"]:
                    bounds = np.arange(0, 100, 10)
                    cmap = BrBG_r
                else:
                    bounds = np.arange(0, 60, 5)
                    cmap = plt.cm.YlOrRd

            selected_cube, = np.where(names == month)

            cube = cube_list[selected_cube[0]]
            try:
                cube.coord('grid_latitude').guess_bounds()
                cube.coord('grid_longitude').guess_bounds()  
            except ValueError:
                pass

            # fix percent -> days issue for these four
            if index.name in ["TX90p", "TN90p", "TX10p", "TN10p"]:
                this_month, = np.where(month_names == month)
                ndays = np.array([calendar.monthrange(y, this_month[0])[1] for y in utils.REFERENCEYEARS])

                cube.data = cube.data * ndays[:, None, None] / 100.
                index.units = "days"

            three_month_data += [cube.data]

        # extracted the three months of the season
        season_cube = copy.deepcopy(cube)
        
        three_month_data = np.ma.array(three_month_data)
        # take appropriate seasonal value
        if index.name in ["TX90p", "TN90p", "TX10p", "TN10p"]:
            season_cube.data = np.ma.sum(three_month_data, axis=0)
        elif index.name in ["FD", "ID", "SU", "TR"]:
            season_cube.data = np.ma.sum(three_month_data, axis=0)
        elif index.name in ["TXx", "TNx", "ETR"]:
            season_cube.data = np.ma.max(three_month_data, axis=0)
        elif index.name in ["TXn", "TNn"]:
            season_cube.data = np.ma.min(three_month_data, axis=0)
        elif index.name in ["Rx1day", "Rx5day"]:
            season_cube.data = np.ma.max(three_month_data, axis=0)
        elif index.name in ["CDD", "CWD"]:
            season_cube.data = np.ma.max(three_month_data, axis=0)
        elif index.name in ["R10mm", "R20mm", "PRCPTOT"]:
            season_cube.data = np.ma.sum(three_month_data, axis=0)
        elif index.name in ["R95pTOT", "R99pTOT", "DTR"]:
            season_cube.data = np.ma.mean(three_month_data, axis=0)
        elif index.name in ["TNlt2", "TNltm2", "TNltm20", "TXge35", "TXge30", "TMlt10", "TMge10", "TMlt5", "TMge5"]:
            season_cube.data = np.ma.sum(three_month_data, axis=0)
        elif index.name in ["TMm", "TXm", "TNm", "TXgt50p"]:
            season_cube.data = np.ma.mean(three_month_data, axis=0)


        # mask if fewer that 2 months present
        nmonths_locs = np.ma.count(three_month_data, axis=0)
        season_cube.data = np.ma.masked_where(nmonths_locs < 2, season_cube.data)

        # get recent period and trend
        if anomalies != "climatology":
            postYYYY = periodConstraint(season_cube, utils.TREND_START)
            season_cube = season_cube.extract(postYYYY)
            trend_cube, sigma, significance = TrendingCalculation(season_cube)

        if anomalies != "climatology":
            figtext=""
            if index.name == "TX90p":
                if season == "DJF":
                    figtext = "(a)"
                elif season == "MAM":
                    figtext = "(b)"
                elif season == "JJA":
                    figtext = "(c)"
                elif season == "SON":
                    figtext = "(d)"
            elif index.name == "TN10p":
                if season == "DJF":
                    figtext = "(e)"
                elif season == "MAM":
                    figtext = "(f)"
                elif season == "JJA":
                    figtext = "(g)"
                elif season == "SON":
                    figtext = "(h)"

            outname = putils.make_filenames("trend", index=index.name, grid=grid, anomalies=anomalies, month=season)

            putils.plot_smooth_map_iris("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname), trend_cube, cmap, bounds, "Trend ({}/10 year)".format(index.units), title="{} - {}, {}-2018".format(index.name, season, utils.TREND_START), figtext=figtext, significance=significance)

        else:
            outname = putils.make_filenames("climatology", index=index.name, grid=grid, anomalies=anomalies, month=season)

            putils.plot_smooth_map_iris("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname), cube[0], cmap, bounds, "{}".format(index.units), title="{} - {}, {}-{}".format(index.name, season, utils.CLIM_START.year, utils.CLIM_END.year))

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
    
    if args.index in ["TX90p", "TN90p", "TX10p", "TN10p", "TXx", "TNx", "TXn", "TNn", "FD", "ID", "SU", "TR", "Rx1day", "Rx5day", "R10mm", "R20mm", "DTR", "ETR", "PRCPTOT", "CDD", "CWD"]:
        main(index=args.index, diagnostics=args.diagnostics, anomalies=args.anomalies, grid=args.grid)
    if args.index in ["TMm", "TXm", "TNm", "TXge35", "TXge30", "TMlt10", "TMge10", "TMlt5", "TMge5", "TXgt50p", "TNlt2", "TNltm2", "TNltm20"]:
        main(index=args.index, diagnostics=args.diagnostics, anomalies=args.anomalies, grid=args.grid)
        

#*******************************************
# END
#*******************************************
