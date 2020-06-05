#!/usr/bin/ python
#***************************************
#
#   Plot the trend maps for the globe.
#
#   21 May 2018 RJHD
#***************************************
#************************************************************************
#                    SVN Info
# $Rev:: 495                                           $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2020-06-05 11:19:31 +0100 (Fri, 05 Jun 2020) $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
"""
Plots maps of differences between grids using the two reference periods

plot_maps_interRefPeriod..py invoked by typing::

  python plot_maps_interRefPeriod.py --index "TX90p" --anomalies --grid --diagnostics

Input arguments:

--index         Which ETCCDI index to use

--anomalies     Run on anomaly outputs

--grid          Set the gridding method (ADW)

--diagnostics   Output extra info (default False)
"""
import os
import copy
import calendar
import datetime as dt
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import iris
import iris.plot
import iris.analysis.stats
import cartopy

# local utilities scripts
import utils
import plt_utils as putils

FLTMDI = -1.e30

# set up calendar month names
month_names = calendar.month_abbr[:]
month_names[0] = "Ann"

COLOURS = {"61-90" : plt.cm.Reds, "81-10" : plt.cm.Blues}
coverage = {}


#************************************************************************
def make_iris_cube_2d(data, lats, lons, name, units):
    """
    Make an Iris cube of a single year of data from arrays of lat, lon and data

    :param array data: data (lat x lon)
    :param array lons: longitudes
    :param array lats: latitudes
    :param str name: name for the cube
    :param str units: units of the data

    :returns: cube - iris cube
    """
    

    # create the iris cube
    cube = iris.cube.Cube(data)
    cube.rename(name)
    cube.units = units
    
    # single field, so no time dimension needed
  
    latcoord = iris.coords.DimCoord(lats, standard_name='latitude', units='degrees')
    loncoord = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    cube.add_dim_coord(latcoord, 0)
    cube.add_dim_coord(loncoord, 1)

    try:
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
    except ValueError:
        pass

    # set coordinate system - 07-08-2014
    cs = iris.coord_systems.GeogCS(6371229)

    cube.coord('latitude').coord_system = cs
    cube.coord('longitude').coord_system = cs

    return cube # make_iris_cube_2d

#***************************************
def TrendCoverageCalculation(cube, verbose=False, plots=False):
    '''
    Calculate coverage of the trend map for index for each gridbox
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

                # get trend, sigma and significance
                trend[lat, lon] = 1


    trend = np.ma.masked_where(trend == FLTMDI, trend)

    return trend # TrendCoverageCalculation


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

    cube_list = iris.load(os.path.join(utils.FINALROOT, utils.make_filenames(index=index.name, grid=grid, anomalies=anomalies, extra="", month_index="")))    

    names = np.array([cube.var_name for cube in cube_list])

    #*************
    # plot trend map

    for month, name in enumerate(month_names[:nmonths]):

        for BASEP in ["61-90", "81-10"]:

            # have to do by hand - and amend when adding anomalies
            infilename = "{}/{}_{}_{}-{}_{}_{}_{}.nc".format(utils.FINALROOT, "HadEX3", index.name, utils.STARTYEAR.year, utils.ENDYEAR.year-1, grid, BASEP, "{}x{}deg".format(utils.DELTALAT, utils.DELTALON))

            try:
                cubelist = iris.load(infilename)
                names = np.array([c.var_name for c in cubelist])
                cube = cubelist[np.where(names == name)[0][0]]

            except RuntimeError:
                print("File not found: {}".format(infilename))
                continue

            except IOError:
                print("File not found: {}".format(infilename))
                continue

            except KeyError:
                print("Month not available in {}".format(filename))
                continue

            try:
                cube.coord('grid_latitude').guess_bounds()
                cube.coord('grid_longitude').guess_bounds()  
            except ValueError:
                pass
  
            # get recent period and trend
            postYYYY = periodConstraint(cube, utils.TREND_START)
            cube = cube.extract(postYYYY)
            preYYYY = periodConstraint(cube, utils.TREND_END, lower=False)
            cube = cube.extract(preYYYY)
            coverage[BASEP] = TrendCoverageCalculation(cube, verbose=diagnostics)

        # now have two cubes.  To calculate individual and overlap
        fig = plt.figure(figsize=(8, 5))
        
        plt.clf()
        ax = plt.axes([0.01, 0.05, 0.98, 0.95], projection=cartopy.crs.Robinson())
        ax.gridlines() #draw_labels=True)
        ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
        ax.coastlines()
        
        ext = ax.get_extent() # save the original extent

        parent_cube = make_iris_cube_2d(np.ma.ones(coverage["61-90"].shape), cube.coord("grid_latitude").points, cube.coord("grid_longitude").points, "trend", "")

        # both
        cmap = plt.cm.Set2
        norm = mpl.cm.colors.BoundaryNorm(np.arange(0, 1, 0.125), cmap.N)
        raw_mask = np.ones(coverage["61-90"].shape)
        locs = np.where(np.logical_and(coverage["61-90"].mask == False,  coverage["81-10"].mask == False))
        raw_mask[locs] = 0
        dummy_cube = copy.deepcopy(parent_cube)
        dummy_cube.data.data[:] = 0.1
        dummy_cube.data.mask = raw_mask       
        mesh = iris.plot.pcolormesh(dummy_cube, cmap=cmap, norm=norm)

        # 81-10
        raw_mask = np.ones(coverage["61-90"].shape)
        locs = np.where(np.logical_and(coverage["61-90"].mask == True,  coverage["81-10"].mask == False))
        raw_mask[locs] = 0
        dummy_cube = copy.deepcopy(parent_cube)
        dummy_cube.data.data[:] = 0.7
        dummy_cube.data.mask = raw_mask       
        mesh = iris.plot.pcolormesh(dummy_cube, cmap=cmap, norm=norm)

        # 61-90
        raw_mask = np.ones(coverage["61-90"].shape)
        locs = np.where(np.logical_and(coverage["61-90"].mask == False,  coverage["81-10"].mask == True))
        raw_mask[locs] = 0
        dummy_cube = copy.deepcopy(parent_cube)
        dummy_cube.data.data[:] = 0.4
        dummy_cube.data.mask = raw_mask       
        mesh = iris.plot.pcolormesh(dummy_cube, cmap=cmap, norm=norm)

        ax.set_extent(ext, ax.projection) # fix the extent change from colormesh

        import matplotlib.patches as mpatches
        colors = ["#66c2a5", "#e78ac3", "#ffd92f"]
        texts = ["Both", "61-90", "81-10"]
        patches = [ mpatches.Patch(color=colors[i], label="{:s}".format(texts[i]) ) for i in range(len(texts)) ]
        plt.legend(handles=patches, bbox_to_anchor=(0.5, -0.12), loc='lower center', ncol=3, frameon=False, prop={'size':utils.FONTSIZE}, labelspacing=0.15, columnspacing=0.5)

        plt.title("{} - Trend coverage comparison".format(index.name), fontsize=utils.FONTSIZE)

        if utils.WATERMARK:
            watermarkstring = "{} {}".format(os.path.join("/".join(os.getcwd().split('/')[4:]), os.path.basename(__file__)), dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M"))
            plt.figtext(0.01, 0.01, watermarkstring, size=6)

        if index.name in ["TX90p"]:
            fig.text(0.03, 0.97, "(c)", fontsize=utils.FONTSIZE)
        elif index.name in ["TN10p"]:
            fig.text(0.03, 0.97, "(d)", fontsize=utils.FONTSIZE)
        elif index.name in ["R95p"]:
            fig.text(0.03, 0.97, "(g)", fontsize=utils.FONTSIZE)
        elif index.name in ["R99p"]:
            fig.text(0.03, 0.97, "(h)", fontsize=utils.FONTSIZE)
            
        outname = putils.make_filenames("interRefP_map", index=index.name, grid=grid, anomalies="None", month=name)
        
        plt.savefig("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname), dpi=300)
        plt.close()


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
