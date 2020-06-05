#!/usr/bin/ python
#***************************************
#
#   Plot the trend maps for the globe.
#
#   21 May 2018 RJHD
#***************************************
#************************************************************************
#                    SVN Info
# $Rev:: 332                                           $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2019-04-16 16:06:06 +0100 (Tue, 16 Apr 2019) $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
"""
Plots maps of differences between average over two 30-year periods

plot_temporal_difference.py invoked by typing::

  python plot_temporal_difference.py --index "TX90p" --anomalies --grid --diagnostics

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


#*********************************************************
def periodConstraint(cube, t1, t2):
    """
    Return constraint that cube times are t1 <= cube < t2
    """

    return iris.Constraint(time=lambda cell: t1 <= cell.point < t2) # periodConstraint

#*********************************************************
def get_climatology(cube, start, length):

    period = periodConstraint(cube, start, start+length)
    cube = cube.extract(period)

    # extract number of years present
    nyears = np.ma.count(cube.data, axis = 0)

    clim_cube = cube.collapsed(['time'], iris.analysis.MEAN)
    sigma_cube = cube.collapsed(['time'], iris.analysis.STD_DEV)

    # mask regions with less than 50% completeness
    clim_cube.data = np.ma.masked_where(nyears <= length/2., clim_cube.data)
    sigma_cube.data = np.ma.masked_where(nyears <= length/2., sigma_cube.data)
    
    return clim_cube, sigma_cube # get_climatology

#*********************************************************
def main(index, first, second, length, diagnostics=False, anomalies="None", grid="ADW"):
    """
    Plot maps of linear trends

    :param str index: which index to run
    :param int first: start of first period
    :param int second: start of second period
    :param int length: length of periods
    :param bool diagnostics: output diagnostic information
    :param str anomalies: run code on anomalies or climatology rather than raw data
    :param str grid: gridding type ADW/CAM
    """

    if first + length -1 > second:
        print("Periods overlap, please re-specify")
        return

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

    print(index.name)

    if index.name in ["TX90p", "TN90p", "SU", "TR", "GSL"]:
        bounds = [-100, -20, -15, -10, -5, 0, 5, 10, 15, 20, 100]
        cmap = RdYlBu_r
    elif index.name in ["DTR", "ETR"]:
        bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
        cmap = RdYlBu_r
    elif index.name in ["WSDI"]:
        bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
        cmap = RdYlBu_r
    elif index.name in ["TX10p", "TN10p", "ID"]:
        bounds = [-100, -20, -15, -10, -5, 0, 5, 10, 15, 20, 100]
        cmap = RdYlBu
    elif index.name in ["FD"]:
        bounds = [-100, -12, -9, -6, -3, 0, 3, 6, 9, 12, 100]
        cmap = RdYlBu
    elif index.name in ["CSDI"]:
        bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
        cmap = RdYlBu
    elif index.name in ["TXn", "TNn"]:
        bounds = [-100, -5, -3, -2, -1, 0, 1, 2, 3, 5, 100]
        cmap = RdYlBu_r
    elif index.name in ["TXx", "TNx"]:
        bounds = [-100, -4, -2, -1, 0.5, 0, 0.5, 1, 2, 4, 100]
        cmap = RdYlBu_r
    elif index.name in ["Rx1day"]:
        bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
        cmap = BrBG
    elif index.name in ["Rx5day"]:
        bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
        cmap = BrBG
    elif index.name in ["PRCPTOT"]:
        bounds = [-100, -40, -20, -10, -5, 0, 5, 10, 20, 40, 100]
        cmap = BrBG
    elif index.name in ["Rnnmm", "R95p", "R99p"]:
        bounds = [-100, -20, -15, -10, -5, 0, 5, 10, 15, 20, 100]
        cmap = BrBG
    elif index.name in ["R95pTOT", "R99pTOT"]:
        bounds = [-100, -10, -5, -2.5, -1, 0, 1, 2.5, 5, 10, 100]
        cmap = BrBG
    elif index.name in ["R10mm"]:
        bounds = [-100, -10, -5, -2.5, -1, 0, 1, 2.5, 5, 10, 100]
        cmap = BrBG
    elif index.name in ["R20mm"]:
        bounds = [-100, -5, -2.5, -1, -0.5, 0, 0.5, 1, 2.5, 5, 100]
        cmap = BrBG
    elif index.name in ["CWD"]:
        bounds = [-100, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 100]
        cmap = BrBG
    elif index.name in ["SDII"]:
        bounds = [-100, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 100]
        cmap = BrBG
    elif index.name in ["CDD"]:
        bounds = [-100, -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10, 100]
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
        bounds = [-100, -20, -15, -10, -5, 0, 5, 10, 15, 20, 100]
        cmap = RdYlBu_r

    cube_list = iris.load(os.path.join(utils.FINALROOT, utils.make_filenames(index=index.name, grid=grid, anomalies=anomalies, extra="", month_index="")))    

    names = np.array([cube.var_name for cube in cube_list])

    #*************
    # plot difference map

    for month in range(nmonths):
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

        # get two cubes to difference 
        first_cube, first_sigma = get_climatology(cube, first, length)

        second_cube, second_sigma = get_climatology(cube, second, length)

        differences = second_cube - first_cube

        # get "significance" by looking at non-overlapping sigmas
        total_sigma = first_sigma + second_sigma

        significance = np.ma.zeros(first_cube.shape)
        significance[differences.data > total_sigma.data] = 1
        significance.mask = differences.data.mask

        significance = putils.make_iris_cube_2d(significance, cube.coord("grid_latitude").points, cube.coord("grid_longitude").points, "difference_significance", "")


        first_string = "{}{}".format(str(first)[-2:], str(first+length-1)[-2:])
        second_string = "{}{}".format(str(second)[-2:], str(second+length-1)[-2:])

        if index.units == "degrees_C":
            units = '$^{\circ}$'+"C"
        else:
            units = index.units

        if anomalies != "climatology":
            outname = putils.make_filenames("diff_{}-{}".format(second_string, first_string), index=index.name, grid=grid, anomalies=anomalies, month=month_names[month])
            
            putils.plot_smooth_map_iris("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname), differences, cmap, bounds, "Difference {}-{} ({})".format(second_string, first_string, units), title="{} - {}, Difference ({}-{}) - ({}-{})".format(index.name, month_names[month], second, second+length-1, first, first+length-1), figtext="(b)") #, significance=significance)
 
    return # main

#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', dest='index', action='store', default="TX90p",
                        help='Which index to run')
    parser.add_argument('--first', dest='first', action='store', default=1951,
                        help='Start of first period')
    parser.add_argument('--second', dest='second', action='store', default=1981,
                        help='Start of second period')
    parser.add_argument('--length', dest='length', action='store', default=30,
                        help='Length of periods')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default=False')
    parser.add_argument('--anomalies', dest='anomalies', action='store', default="None",
                        help='To use the anomalies rather than actuals, default="None"') 
    parser.add_argument('--grid', dest='grid', action='store', default="ADW",
                        help='gridding routine to use (ADW or CAM) default="ADW"') 

    args = parser.parse_args()
    
    if args.anomalies != "climatology":
        main(index=args.index, first=int(args.first), second=int(args.second), length=int(args.length), \
             diagnostics=args.diagnostics, anomalies=args.anomalies, grid=args.grid)
    

#*******************************************
# END
#*******************************************
