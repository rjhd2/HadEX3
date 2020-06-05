#!/usr/local/sci/bin/python
#***************************************
#
#   Plot the timeseries for the globe.
#     compares between existing ClimdEx family
#
#   21 May 2018 RJHD
#***************************************
#************************************************************************
#                    SVN Info
# $Rev:: 371                                           $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2019-05-30 11:46:57 +0100 (Thu, 30 May 2019) $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
"""
Plots timeseries between grids using the two reference periods

plot_timeseries_interRefPeriod..py invoked by typing::

  python plot_timeseries_interRefPeriod.py --index "TX90p" --diagnostics

Input arguments:

--index         Which ETCCDI index to use

--diagnostics   Output extra info (default False)
"""
import os
import calendar
import datetime as dt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import numpy as np

import iris
import iris.coord_categorisation
import cf_units

# local utilities scripts
import utils
import plt_utils as putils
import plot_timeseries as pts

# plot defaults
plt.rcParams["xtick.major.pad"] = "8"
plt.rcParams["xtick.direction"] = "out"
plt.rcParams["ytick.direction"] = "out"
minorLocator = MultipleLocator(1)


TSOUTLOCATION = utils.INFILELOCS

DAYSPERYEAR = 365. # to convert from days in the year to percent (0->365 to 0->100)

# set up calendar month names
month_names = calendar.month_abbr[:]
month_names[0] = "Ann"

#***************************************
def PlotCoverage(plt, ax, valid_boxes, colours, starttime, endtime, index, month, leg_title, ncol=2):
    """
    Plots the coverage plot

    :param object plt: plot instance
    :param object ax: axes instance
    :param array valid_boxes: array of number of valid boxes with time
    :param array colours: array of valid colours
    :param int starttime: start of x-axis
    :param int endtime: end of x-axis
    :param str index: index
    :param str month: month name
    :param str leg_title: legend title
    :param array labels: array of line labels
    :param int ncol: number of columns in legend
    :returns: None
    """
 
    for version, (times, nboxes) in valid_boxes.items():
        years = np.array([int(str(y)[0:4]) for y in times])
        good = np.where(nboxes > 0)
        ax.plot(years[good], nboxes[good], colours[version], label=version, lw=2, ls="-")

    ax.set_ylabel("Valid grid-boxes (%)", fontsize=utils.FONTSIZE)
    ax.set_xlim([starttime-1, endtime])
    ax.set_xticks([i for i in range(starttime-1, endtime+10, 10)])
    ax.tick_params(labelsize=utils.FONTSIZE)
    ax.set_title("{} - {}".format(index, month_names[month]), fontsize=utils.FONTSIZE)
    leg = ax.legend(loc='lower center', ncol=ncol, bbox_to_anchor=(0.46, -0.05), frameon=False, title=leg_title, prop={'size':utils.FONTSIZE}, labelspacing=0.15, columnspacing=0.5)
    if leg_title != "":
        plt.setp(leg.get_title(), fontsize=utils.FONTSIZE)

    if utils.WATERMARK:
        watermarkstring = "{} {}".format(os.path.join("/".join(os.getcwd().split('/')[4:]), os.path.basename(__file__)), dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M"))
        plt.figtext(0.01, 0.01, watermarkstring, size=6)

    return # PlotCoverage

#***************************************
def fix_time_coord(incube):

    try:
        newtimes = np.array([dt.datetime.strptime("{}".format(int(d)), "%Y%m%d") for d in incube.coord("time").points])
    except ValueError:
        newtimes = np.array([dt.datetime.strptime("{}".format(int(d)), "%Y") for d in incube.coord("time").points])
        
    newdiff = np.array([(t - newtimes[0]).days for t in newtimes])
    
    # replace cube time coordinate with new one
    time_unit = cf_units.Unit('days since ' + dt.datetime.strftime(newtimes[0], "%Y-%m-%d %H:%M"), calendar=cf_units.CALENDAR_GREGORIAN)   
    timecoord = iris.coords.DimCoord(newdiff, standard_name='time', units=time_unit, var_name="time") # add bounds?
    incube.remove_coord('time')
    incube.add_dim_coord(timecoord, 0)

    return incube # fix_time_coord

#***************************************
def MaskData(data, mdi, mask):
    """
    Apply the masks to the data

    :param array data: data array
    :param flt mdi: missing data indicator
    :param array mask: np.ma of mask

    :returns: nboxes - number of valid boxes, masked_data - np.ma of data with masks applied

    """

    nyears, nlats, nlons = data.shape
    # the masked regional dataset (per year)
    masked_data = np.ma.zeros([nyears, nlats, nlons])
    masked_data.fill(mdi)
    
    nboxes = np.zeros(nyears)

    for year in np.arange(nyears):

        # apply  mask
        masked_data[year] = np.ma.array(data[year], mask=mask)

        # count number of boxes
        nboxes[year] = len(np.ma.where(masked_data[year] != mdi)[0])

    return nboxes, masked_data # MaskData
#************************************************************************
#************************************************************************
def main(index, comparison=False, diagnostics=False, anomalies="None", grid="ADW", normalise=False, matched=False):
    """
    :param str index: which index to run
    :param bool comparison: compare against other datasets
    :param bool diagnostics: output diagnostic information
    :param str anomalies: run code on anomalies or climatology rather than raw data
    :param str grid: gridding type ADW/CAM
    :param bool normalise: plot as anomalies from 1961-90
    :param bool matched: match HadEX3 to HadEX2 coverage and plot timeseries
    """

    # get details of index
    index = utils.INDICES[index] 

    # allow for option of running through each month
    if index.name in utils.MONTHLY_INDICES:
        nmonths = 13
        timescale = "MON"
    else:
        nmonths = 1
        timescale = "ANN"

    # currently not stored - but ready just in case
    outfilename = "{}/timeseries_{}.dat".format(TSOUTLOCATION, index.name)
    if os.path.exists(outfilename): os.remove(outfilename)

    coveragefilename = "{}/timeseries_{}_boxes.dat".format(TSOUTLOCATION, index.name)
    if os.path.exists(coveragefilename): os.remove(coveragefilename)

    # plot all month versions at once
    for month, name in enumerate(month_names[:nmonths]):

        if diagnostics:
            print(name)

        # set up the figure
        fig = plt.figure(figsize=(8, 5.5))
        plt.clf()
        ax = fig.add_axes([0.17, 0.35, 0.8, 0.6])
        ax2 = fig.add_axes([0.17, 0.1, 0.8, 0.2], sharex=ax)

        # number of valid boxes
        timeseries = {}
        valid_boxes = {}
        land_boxes = {}
        colours = {}

        masks = []
        # get common mask first
        for BASEP in ["61-90", "81-10"]:

            # have to do by hand - and amend when adding anomalies
            infilename = "{}/{}_{}_{}-{}_{}_{}_{}.nc".format(utils.FINALROOT, "HadEX3", index.name, utils.STARTYEAR.year, utils.ENDYEAR.year-1, grid, BASEP, "{}x{}deg".format(utils.DELTALAT, utils.DELTALON))

            try:
                cubelist = iris.load(infilename)
                names = np.array([c.var_name for c in cubelist])
                incube = cubelist[np.where(names == name)[0][0]]
                incube.coord("grid_latitude").standard_name = "latitude"
                incube.coord("grid_longitude").standard_name = "longitude"
        
                incube = fix_time_coord(incube)
            except RuntimeError:
                print("File not found: {}".format(infilename))
                continue

            except IOError:
                print("File not found: {}".format(infilename))
                continue

            except KeyError:
                print("Month not available in {}".format(filename))
                continue

            masks += [incube.data.mask]

            # have to extract box counts before using a common mask
            LABEL = "{}".format(BASEP)
            coord = incube.coord("time")
            years = np.array([c.year for c in coord.units.num2date(coord.points)])

            # find which boxes have x% of years with data - default is 90% for timeseries
            completeness_mask = utils.CompletenessCheckGrid(incube.data, utils.ENDYEAR.year, utils.STARTYEAR.year)

            # apply completeness mask, and obtain box counts
            nboxes_completeness, completeness_masked_data = MaskData(incube.data, incube.data.fill_value, completeness_mask)    
            incube.data = completeness_masked_data

            nboxes = np.zeros(nboxes_completeness.shape[0])
            for year in range(incube.data.shape[0]):
                nboxes[year] = np.ma.count(incube.data[year])

            # for coverage, number of boxes will vary with lat/lon resolution
            max_boxes = np.product(incube.data.shape[1:])
            if diagnostics:
                print("{}: total grid boxes = {}, max filled = {}".format(LABEL, max_boxes, int(np.max(nboxes))))

            # collect outputs
            valid_boxes[LABEL] = [years, 100.*nboxes/max_boxes] # scale by total number as lat/lon resolution will be different
            if month == 0:
                lsm = utils.get_land_sea_mask(incube.coord("latitude").points, incube.coord("longitude").points, floor=False)
                n_land_boxes = len(np.where(lsm == False)[0])
            land_boxes[LABEL] = [years, 100.*nboxes/n_land_boxes] # scale by number of land boxes

        # now merge masks

        final_mask = np.ma.mask_or(masks[0], masks[1])

        for BASEP in ["61-90", "81-10"]:

            # have to do by hand - and amend when adding anomalies
            infilename = "{}/{}_{}_{}-{}_{}_{}_{}.nc".format(utils.FINALROOT, "HadEX3", index.name, utils.STARTYEAR.year, utils.ENDYEAR.year-1, grid, BASEP, "{}x{}deg".format(utils.DELTALAT, utils.DELTALON))

            try:
                cubelist = iris.load(infilename)
                names = np.array([c.var_name for c in cubelist])
                incube = cubelist[np.where(names == name)[0][0]]
                incube.coord("grid_latitude").standard_name = "latitude"
                incube.coord("grid_longitude").standard_name = "longitude"

                incube = fix_time_coord(incube)

                incube.data.mask = final_mask

            except RuntimeError:
                print("File not found: {}".format(infilename))
                continue

            except IOError:
                print("File not found: {}".format(infilename))
                continue

            except KeyError:
                print("Month not available in {}".format(filename))
                continue


            # fix percent -> days issue for these four
            if index.name in ["TX90p", "TN90p", "TX10p", "TN10p"]:
                if name == "Ann":
                    incube.data = incube.data * (DAYSPERYEAR/100.)
                else:
                    incube.data = incube.data * (calendar.monthrange(2019, month)[1]/100.)
                    
                index.units = "days"

            # restrict to times of interest
            time_constraint = iris.Constraint(time=lambda cell: utils.STARTYEAR <= cell <= utils.ENDYEAR)

            incube = incube.extract(time_constraint)

            # find which boxes have x% of years with data - default is 90%
            completeness_mask = utils.CompletenessCheckGrid(incube.data, utils.ENDYEAR.year, utils.STARTYEAR.year)

            # apply completeness mask, and obtain box counts
            nboxes_completeness, completeness_masked_data = MaskData(incube.data, incube.data.fill_value, completeness_mask)    
            incube.data = completeness_masked_data

            if normalise:
                # apply normalisation!
#                clim_constraint = iris.Constraint(time=lambda cell: dt.datetime(utils.REF_START, 1, 1) <= cell <= dt.datetime(utils.REF_END, 1, 1))

                if BASEP == "61-90":
                    clim_constraint = iris.Constraint(time=lambda cell: dt.datetime(1961, 1, 1) <= cell <= dt.datetime(1990, 1, 1))
                elif BASEP == "81-10":
                    clim_constraint = iris.Constraint(time=lambda cell: dt.datetime(1981, 1, 1) <= cell <= dt.datetime(2010, 1, 1))

                norm_cube = incube.extract(clim_constraint)
                norm = norm_cube.collapsed(['time'], iris.analysis.MEAN)

                incube = incube - norm

            # weights for the region
            weights = iris.analysis.cartography.cosine_latitude_weights(incube)
            ts = incube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=weights)

            # only plot where there are non-missing values
            coord = ts.coord("time")
            years = np.array([c.year for c in coord.units.num2date(coord.points)])


            # do the plot, highlighting vanilla HadEX3
            LABEL = "{}".format(BASEP)
            if BASEP == "61-90":
                line = ax.plot(years, ts.data, c='k', ls="-", lw=2, label=LABEL, zorder=10)
            else:
                line = ax.plot(years, ts.data, ls="-", lw=2, label=LABEL, zorder=5)
            timeseries[LABEL] = [years, ts]
            colours[LABEL] = line[0].get_color()

        # once all lines plotted, then tidy up
        putils.SortAxesLabels(plt, ax, index, utils.STARTYEAR.year, utils.ENDYEAR.year, month)
        putils.SortAxesLabels(plt, ax2, index, utils.STARTYEAR.year, utils.ENDYEAR.year, month)
        ax.set_xlim([1900, 2020])

        if normalise:
            # only plot zero line if done as anomalies
            ax.axhline(0, color='0.5', ls='--')
        elif index.name in ["TX90p", "TX10p", "TN90p", "TN10p"]:
            # or plot expected line otherwise
            if name == "Ann":
                ax.axhline(36.5, color='0.5', ls='--')
            else:
                ax.axhline(calendar.monthrange(2019, month)[1]*0.1, color='0.5', ls='--')

        # plot legend below figure
        leg = ax.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.46, -0.05), frameon=False, title='', prop={'size':utils.FONTSIZE}, labelspacing=0.15, columnspacing=0.5)

        # and plot the difference
        ax2.plot(years, timeseries["61-90"][1].data-timeseries["81-10"][1].data, ls="-", lw=2, c="#e41a1c")
        ax2.set_title("")
        ax2.set_ylabel("Difference\n({})".format(index.units), fontsize=utils.FONTSIZE)
#        ax2.axhline(0, color='k', ls='--')


        # extra information
        if utils.WATERMARK:
            watermarkstring = "{} {}".format(os.path.join("/".join(os.getcwd().split('/')[4:]), os.path.basename(__file__)), dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M"))
            plt.figtext(0.01, 0.01, watermarkstring, size=6)

        ax = putils.Rstylee(ax)
        ax2 = putils.Rstylee(ax2)
        plt.setp(ax.get_xticklabels(), visible=False)

        if index.name in ["TX90p"]:
            fig.text(0.03, 0.97, "(a)", fontsize=utils.FONTSIZE)
        elif index.name in ["TN10p"]:
            fig.text(0.03, 0.97, "(b)", fontsize=utils.FONTSIZE)
        elif index.name in ["R95p"]:
            fig.text(0.03, 0.97, "(e)", fontsize=utils.FONTSIZE)
        elif index.name in ["R99p"]:
            fig.text(0.03, 0.97, "(f)", fontsize=utils.FONTSIZE)
            
        # and save
        outname = putils.make_filenames("interRefP_ts", index=index.name, grid=grid, anomalies="None", month=name)

        plt.savefig("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname), dpi=300)

        plt.close()

        #*****************
        # plot coverage - how many grid boxes have values (scaled by total number)

        fig = plt.figure(figsize=(8, 5.5))
        plt.clf()
        ax = fig.add_axes([0.17, 0.35, 0.8, 0.6])
        ax2 = fig.add_axes([0.17, 0.1, 0.8, 0.2], sharex=ax)

#        PlotCoverage(plt, ax, valid_boxes, colours, utils.STARTYEAR.year, utils.ENDYEAR.year, index.name, month, '', ncol=2)
        PlotCoverage(plt, ax, land_boxes, colours, utils.STARTYEAR.year, utils.ENDYEAR.year, index.name, month, '', ncol=2)
        ax.set_xlim([1900, 2020])

        # and plot the difference
#        ax2.plot(years, valid_boxes["81-10"][1].data-valid_boxes["61-90"][1].data, ls="-", lw=2, c="#e41a1c")
        ax2.plot(years, land_boxes["81-10"][1]-land_boxes["61-90"][1], ls="-", lw=2, c="#e41a1c")
        ax2.set_title("")
        ax2.set_ylabel("Difference\n(%)", fontsize=utils.FONTSIZE)
        ax2.tick_params(labelsize=utils.FONTSIZE)
       
        ax = putils.Rstylee(ax)
        ax2 = putils.Rstylee(ax2)
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.xaxis.set_minor_locator(minorLocator)

        outname = putils.make_filenames("interRefP_coverage_ts", index=index.name, grid=grid, anomalies="None", month=name)
        print("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname))
        plt.savefig("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname), dpi=300)

        plt.close("all")

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

    args = parser.parse_args()

    main(index=args.index, diagnostics=args.diagnostics, normalise=False)

#************************************************************************
#                                 END
#************************************************************************
