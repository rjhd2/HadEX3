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
# $Rev:: 435                                           $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2019-11-14 16:42:09 +0000 (Thu, 14 Nov 2019) $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
"""
Plots timeseries with optional normalisation and comparison to other datasets using Iris Cubes

Adds in coverage uncertainty from ERA5

plot_timeseries.py invoked by typing::

  python plot_timeseries.py --index "TX90p" --comparison --anomalies --grid --normalise --matched --diagnostics

Input arguments:

--index         Which ETCCDI index to use

--comparison    Compare to HadEX, HadEX2 and GHCNDEX

--anomalies     Run on anomaly outputs

--grid          Set the gridding method (ADW)

--normalise     Normalise the output over 1961-90

--matched       Match coverage between HadEX3 and HadEX2

--diagnostics   Output extra info (default False)
"""
import os
import copy
import calendar
import datetime as dt
import numpy as np

import iris
import iris.coord_categorisation
import cf_units

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# local utilities scripts
import utils
import plt_utils as putils


# plot defaults
plt.rcParams["xtick.major.pad"] = "8"
plt.rcParams["xtick.direction"] = "out"
plt.rcParams["ytick.direction"] = "out"
minorLocator = MultipleLocator(1)

PLOTERA = False
DATASETS = ["HadEX3", "HadEX2", "HadEX", "GHCNDEX", "ERA5"] 
GHCNDEX_VERSION = "v20190308"

LABELS = {"HadEX3" : "HadEX3", "HadEX2" : "HadEX2", "HadEX" : "HadEX", "GHCNDEX" : "GHCNDEX", "HadEX3-matched" : "HadEX3-matched", "ERA5" : "ERA5"} # , "_H2ADW" : "HadEX2 ADW"}
COLOURS = {"HadEX3" : "k", "HadEX2" : "#e41a1c", "HadEX" : "#4daf4a", "GHCNDEX" : "#3773b8", "HadEX3-matched" : "k", "ERA5" : "#984ea3"}
LS = {"HadEX3" : "-", "HadEX2" : "-", "HadEX" : "-", "GHCNDEX" : "-", "HadEX3-matched" : "--", "ERA5" : "-"}
ZORDER = {"HadEX3" : 10, "HadEX2" : 5, "HadEX" : 4, "GHCNDEX" : 3, "HadEX3-matched" : 9, "ERA5" : 6}

TSOUTLOCATION = utils.INFILELOCS
DAYSPERYEAR = 365.

# set up calendar month names
month_names = calendar.month_abbr[:]
month_names[0] = "Ann"

#***************************************
def CalculateGlobalAverage(masked_data, times, index, nyears, mdi, newlats, normalise=False):
    """
    Calculate the global average of the data

    :param array masked_data: np.ma of the data
    :param array times: time array
    :param str index: index name
    :param int nyears: number of years
    :param flt mdi: missing data indicator
    :param array newlats: latitude array
    :param bool normalise: adjust to climatology period

    :returns: years - time array with years, average - global, cosine weighted average timeseries

    """
    
    # get the average
    average = np.zeros(nyears)
    average.fill(mdi)
 
    for year in np.arange(nyears):
        
        # mask latitudes by the data mask (remove the 360deg longitude bit)
        theselats = np.ma.array(newlats, mask=masked_data[year].mask)
        theselats = np.ma.compressed(theselats)
        # apply the data mask
        theseboxes = np.ma.compressed(masked_data[year])
       # combine the arrays to get the datamask

        # only calculate if there is data
        if len(theselats) != 0:
            if np.sum(np.cos(np.deg2rad(theselats))) > 0.:
                average[year] = np.sum(theseboxes * np.cos(np.deg2rad(theselats))) / np.sum(np.cos(np.deg2rad(theselats)))

    years = np.array([int(str(y)[0:4]) for y in times])
    
    good_avg = np.where(average != mdi)

    # adjust to "climatology period"
    if normalise:
        clim_loc = np.where((years >= utils.CLIM_START) & (years <= utils.CLIM_END) & (average != mdi))[0]

        if len(clim_loc) > 0:
 
            climatology = np.mean(average[clim_loc[0]:clim_loc[-1]+1])
        
            average[good_avg] = average[good_avg]-climatology

    return years, average # CalculateGlobalAverage

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

#************************************************************************
def make_iris_cube_3d(data, times, time_units, lats, lons, name, units):
    """
    Make an Iris cube of data from arrays of years, lat, lon and data

    :param array data: data (years x lat x lon)
    :param array times: times
    :param str time_units: units for time dimension
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
    
    time_unit = cf_units.Unit(time_units, calendar=cf_units.CALENDAR_GREGORIAN)
    timecoord = iris.coords.DimCoord(times, standard_name='time', units=time_unit, var_name="time") # add bounds?
    latcoord = iris.coords.DimCoord(lats, standard_name='latitude', units='degrees')
    loncoord = iris.coords.DimCoord(lons, standard_name='longitude', units='degrees')
    
    cube.add_dim_coord(timecoord, 0)
    cube.add_dim_coord(latcoord, 1)
    cube.add_dim_coord(loncoord, 2)

    cube.coord('time').guess_bounds()
    try:
        cube.coord('latitude').guess_bounds()
    except ValueError:
        # bounds already exist
        pass
    try:
        cube.coord('longitude').guess_bounds()
    except ValueError:
        # bounds already exist
        pass

    # set coordinate system - 07-08-2014
    cs = iris.coord_systems.GeogCS(6371229)

    cube.coord('latitude').coord_system = cs
    cube.coord('longitude').coord_system = cs

    return cube # make_iris_cube_2d

#***************************************
def PlotCoverage(plt, ax, valid_boxes, colours, starttime, endtime, index, month, leg_title, ncol=2, land=False):
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
    :param bool land: true if only land grid boxes (for axis label)

    :returns: None
    """
 
    for dataset, (times, nboxes) in valid_boxes.items():
        if dataset == "ERA5":
            continue
        years = np.array([int(str(y)[0:4]) for y in times])
        good = np.where(nboxes > 0)
        plt.plot(years[good], nboxes[good], COLOURS[dataset], label=LABELS[dataset], lw=2, zorder=ZORDER[dataset], ls=LS[dataset])

    if land:
        plt.ylabel("Valid land grid-boxes (%)", fontsize=utils.FONTSIZE)
    else:
        plt.ylabel("Valid grid-boxes (%)", fontsize=utils.FONTSIZE)
        
    ax.set_xlim([starttime-1, endtime])
    ax.set_xticks([i for i in range(starttime-1, endtime+10, 10)])
    ax.tick_params(labelsize=utils.FONTSIZE)
    ax.set_title("{} - {}".format(index, month_names[month]), fontsize=utils.FONTSIZE)
    leg = plt.legend(loc='lower center', ncol=ncol, bbox_to_anchor=(0.46, -0.31), frameon=False, title=leg_title, prop={'size':utils.FONTSIZE}, labelspacing=0.15, columnspacing=0.5)
    if leg_title != "":
        plt.setp(leg.get_title(), fontsize=utils.FONTSIZE)

    if utils.WATERMARK:
        watermarkstring = "{} {}".format(os.path.join("/".join(os.getcwd().split('/')[4:]), os.path.basename(__file__)), dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M"))
        plt.figtext(0.01, 0.01, watermarkstring, size=6)

    return # PlotCoverage

#************************************************************************
#************************************************************************

def main(index, comparison=False, diagnostics=False, anomalies="None", grid="ADW", normalise=False, matched=False, uncertainties=False):
    """
    :param str index: which index to run
    :param bool comparison: compare against other datasets
    :param bool diagnostics: output diagnostic information
    :param str anomalies: run code on anomalies or climatology rather than raw data
    :param str grid: gridding type ADW/CAM
    :param bool normalise: plot as anomalies from e.g. 1961-90
    :param bool matched: match HadEX3 to HadEX2 coverage and plot timeseries
    :param bool uncertainties: plot ERA5 derived coverage uncertainties
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
        ax = fig.add_axes([0.15, 0.2, 0.82, 0.75])

        # number of valid boxes
        timeseries = {}
        valid_boxes = {}
        land_boxes = {}

        e5cube = 0
        # spin through all comparison datasets
        for ds, dataset in enumerate(DATASETS):
            incube = 0

            print(dataset)
            if not comparison:
                # but if not doing comparisons, skip (just run HadEX3)
                if dataset != "HadEX3":
                    continue

            if dataset == "HadEX":
                if name != "Ann":
                    continue
                else:
                    try:
                        if index.name == "R95pTOT":
                            filename = "{}/HadEX_{}_1951-2003.txt".format(utils.HADEX_LOC, "R95pT")
                        else:
                            filename = "{}/HadEX_{}_1951-2003.txt".format(utils.HADEX_LOC, index.name)


                        all_data, years = [], []
                        data = np.zeros((72, 96))
                        latc = -1
                        with open(filename, "r") as infile:

                            for lc, line in enumerate(infile):

                                if lc == 0:
                                    # skip header
                                    continue

                                # read each line
                                line = line.split()
                                if len(line) < 10:
                                    years += [int(line[0])]
                                    if line[0] != "1951":
                                        all_data += [data]
                                    # reset storage
                                    data = np.zeros((72, 96))
                                    latc = -1 
                                else:
                                    latc += 1
                                    data[latc, :] = line

                        # add final year
                        all_data += [data]
                        all_data = np.array(all_data).astype(float)
                        all_data = np.ma.masked_where(all_data == -999.99, all_data)
                        if index.name == "R95pTOT":
                            all_data *= 100
                        latitudes = np.arange(90, -90, -2.5)
                        longitudes = np.arange(-180, 180, 3.75)

                        incube = make_iris_cube_3d(all_data, years, "unknown", latitudes, longitudes, name, index.units)
                        incube = fix_time_coord(incube)

                    except RuntimeError:
                        print("File not found: {}".format(filename))
                        continue

                    except IOError:
                        print("File not found: {}".format(filename))
                        continue

#                    except IndexError:
#                        print("Month not available in {}".format(filename))
#                        continue
                    
            elif dataset == "HadEX2":
                filename = "{}/HadEX2_{}_1901-2010_h2_mask_m4.nc".format(utils.HADEX2_LOC, index.name)

                try:
                    cubelist = iris.load(filename)
                    names = np.array([c.var_name for c in cubelist])
                    incube = cubelist[np.where(names == name)[0][0]]
                    incube.coord("lat").standard_name = "latitude"
                    incube.coord("lon").standard_name = "longitude"

                    incube = fix_time_coord(incube)

                except RuntimeError:
                    print("File not found: {}".format(filename))
                    continue
                        
                except IOError:
                    print("File not found: {}".format(filename))
                    continue

                except IndexError:
                    print("Month not available in {}".format(filename))
                    continue
                    
            elif dataset == "HadEX3":
                
                filename = os.path.join(utils.FINALROOT, utils.make_filenames(index=index.name, grid=grid, anomalies=anomalies, extra="", month_index=""))


                try:
                    cubelist = iris.load(filename)
                    names = np.array([c.var_name for c in cubelist])
                    incube = cubelist[np.where(names == name)[0][0]]
                    incube.coord("grid_latitude").standard_name = "latitude"
                    incube.coord("grid_longitude").standard_name = "longitude"

                    incube = fix_time_coord(incube)

                    h3_cube = copy.deepcopy(incube)

                    if anomalies == "anomalies":

                        # to make into actuals, add climatology to the anomalies
                        clim_filename = os.path.join(utils.FINALROOT, utils.make_filenames(index=index.name, grid=grid, anomalies="climatology", extra="", month_index=""))

                        clim_cubelist = iris.load(clim_filename)
                        names = np.array([c.var_name for c in cubelist])
                        clim_cube = clim_cubelist[np.where(names == name)[0][0]]
                        try:
                            clim_cube.coord("grid_latitude").standard_name = "latitude"
                            clim_cube.coord("grid_longitude").standard_name = "longitude"
                        except iris.exceptions.CoordinateNotFoundError:
                            pass

                        if clim_cube.coord("time").units.origin != "days since 1901-01-01 00:00":
                            clim_cube = fix_time_coord(clim_cube)

                except RuntimeError:
                    print("File not found: {}".format(filename))
                    continue
                        
                except IOError:
                    print("File not found: {}".format(filename))
                    continue

                except IndexError:
                    print("Month not available in {}".format(filename))
                    continue
                    
            elif dataset == "GHCNDEX":

                filename = "{}/{}/GHCND_{}_1951-2019_RegularGrid_global_2.5x2.5deg_LSmask.nc".format(utils.GHCNDEX_LOC, GHCNDEX_VERSION, index.name)


                try:
                    cubelist = iris.load(filename)
                    names = np.array([c.var_name for c in cubelist])
                    incube = cubelist[np.where(names == name)[0][0]]

                    incube = fix_time_coord(incube)

                except RuntimeError:
                    print("File not found: {}".format(filename))
                    continue
                        
                except IOError:
                    print("File not found: {}".format(filename))
                    continue
                    
                except IndexError:
                    print("Month not available in {}".format(filename))
                    continue
                    
            elif dataset == "ERA5":

                filename = "{}/ERA5_{}_1979-2019.nc".format(utils.ERA5_LOC, index.name)

                try:
                    cubelist = iris.load(filename)
                    names = np.array([c.var_name for c in cubelist])
                    incube = cubelist[np.where(names == name)[0][0]]
                    
                    # match latitude order
                    incube.coord('latitude').points = incube.coord('latitude').points[::-1]
                    incube.data = incube.data[:, ::-1, :]

                    # match to H3 grid
                    try:
                        h3_cube.coord("longitude").guess_bounds()
                        h3_cube.coord("latitude").guess_bounds()
                    except ValueError:
                        # already has bounds
                        pass
                    try:
                        incube.coord("longitude").guess_bounds()
                        incube.coord("latitude").guess_bounds()
                    except ValueError:
                        # already has bounds
                        pass
                    incube = incube.regrid(h3_cube, iris.analysis.Linear(extrapolation_mode="mask"))
                    
                    e5_cube = copy.deepcopy(incube)

                except RuntimeError:
                    print("File not found: {}".format(filename))
                    e5_cube = False
                    continue
                        
                except IOError:
                    print("File not found: {}".format(filename))
                    e5_cube = False
                    continue
                    
                except IndexError:
                    print("Month not available in {}".format(filename))
                    e5_cube = False
                    continue
                    
                print("need to match grids")


            # process data ready for plotting

            # if anomalies from HadEX3, then need to add onto climatology
            if anomalies == "anomalies" and dataset == "HadEX3":
                incube.data = incube.data + clim_cube.data

            # fix percent -> days issue for these four
            if index.name in ["TX90p", "TN90p", "TX10p", "TN10p"]:
                incube.data = incube.data * (DAYSPERYEAR/100.)
                index.units = "days"

            # restrict to times of interest
            time_constraint = iris.Constraint(time=lambda cell: utils.STARTYEAR <= cell <= utils.ENDYEAR)

            incube = incube.extract(time_constraint)

            if matched and (utils.DELTALON == 3.75 and utils.DELTALAT == 2.5):
                # if matching coverage, retain hadex2
                if dataset == "HadEX2":
                    hadex2_mask = incube.data.mask[:]
                if dataset == "HadEX3":
                    # new cube to hold
                    matched_hadex3 = incube.data[:hadex2_mask.shape[0]]


            # find which boxes have x% of years with data - default is 90%
            if dataset == "GHCNDEX":
                completeness_mask = utils.CompletenessCheckGrid(incube.data, utils.ENDYEAR.year, 1951)
            elif dataset == "HadEX":
                completeness_mask = utils.CompletenessCheckGrid(incube.data, 2003, 1951)
            elif dataset == "HadEX2":
                completeness_mask = utils.CompletenessCheckGrid(incube.data, 2010, utils.STARTYEAR.year)
            elif dataset == "HadEX3":
                completeness_mask = utils.CompletenessCheckGrid(incube.data, utils.ENDYEAR.year, utils.STARTYEAR.year)
            elif dataset == "ERA5":
                completeness_mask = utils.CompletenessCheckGrid(incube.data, utils.ENDYEAR.year, 1979)

            # extract number of boxes before applying temporal completeness
            nboxes = np.zeros(incube.data.shape[0])
            for year in range(incube.data.shape[0]):
                nboxes[year] = np.ma.count(incube.data[year])

            # apply completeness mask, and obtain box counts
            nboxes_completeness, completeness_masked_data = MaskData(incube.data, incube.data.fill_value, completeness_mask)    
            incube.data = completeness_masked_data
            
            if normalise:
                # apply normalisation!
                clim_constraint = iris.Constraint(time=lambda cell: dt.datetime(utils.REF_START, 1, 1) <= cell <= dt.datetime(utils.REF_END, 1, 1))
                norm_cube = incube.extract(clim_constraint)
                norm = norm_cube.collapsed(['time'], iris.analysis.MEAN)

                incube = incube - norm

            # weights for the region
            weights = iris.analysis.cartography.cosine_latitude_weights(incube)
            ts = incube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=weights)

            # only plot where there are non-missing values
            coord = ts.coord("time")
            years = np.array([c.year for c in coord.units.num2date(coord.points)])
            if dataset == "ERA5":
                if PLOTERA:
                    # only plot ERA5 if selected
                    plt.plot(years, ts.data, c=COLOURS[dataset], ls=LS[dataset], lw=2, label=LABELS[dataset], zorder=ZORDER[dataset])
            else:
                plt.plot(years, ts.data, c=COLOURS[dataset], ls=LS[dataset], lw=2, label=LABELS[dataset], zorder=ZORDER[dataset])

            if dataset == "HadEX3":
                h3_ts = ts
                h3_years = years

            # save
            max_boxes = np.product(incube.data.shape[1:])
            if diagnostics:
                print("{}: total grid boxes = {}, max filled = {}".format(dataset, max_boxes, int(np.max(nboxes))))

            valid_boxes[dataset] = [years, 100.*nboxes/max_boxes] # scale by total number as lat/lon resolution will be different
            
            store_ts = np.ones(h3_ts.shape) * utils.HADEX_MDI
            match = np.in1d(h3_years, years)
            match_back = np.in1d(years, h3_years)
            store_ts[match] = ts.data[match_back]
            timeseries[dataset] = [h3_years, store_ts]

            # get land sea mask
            if month == 0:
                lsm = utils.get_land_sea_mask(incube.coord("latitude").points, incube.coord("longitude").points, floor=False)
                n_land_boxes = len(np.where(lsm == False)[0])
            land_boxes[dataset] = [years, 100.*nboxes/n_land_boxes] # scale by number of land boxes


        # once all lines plotted
        # get coverage error for HadEX3
        if uncertainties:
            try:
                # test to see if there was an actual cube from this loop (else using stale cube from previous)
                if e5_cube != 0:
                    coverage_offset, coverage_stdev = putils.compute_coverage_error(h3_cube, e5_cube)
                    coverage_stdev *= 2. # 90%, 2s.d.
                    plt.fill_between(h3_years, h3_ts.data-coverage_stdev, h3_ts.data+coverage_stdev, color='0.5', label="ERA5 coverage uncertainty")
            except UnboundLocalError:
                # e5_cube not referenced - i.e. no ERA5 data for this index?
                pass

        # then tidy up
        putils.SortAxesLabels(plt, ax, index, utils.STARTYEAR.year, utils.ENDYEAR.year, month)

        plt.xlim([1900, 2020])
        if normalise:
            # only plot zero line if done as anomalies
            plt.axhline(0, color='k', ls='--')

        # plot legend below figure
        leg = plt.legend(loc='lower center', ncol=3, bbox_to_anchor=(0.46, -0.31), frameon=False, title='', prop={'size':utils.FONTSIZE}, labelspacing=0.15, columnspacing=0.5)
        # extra information
        if utils.WATERMARK:
            watermarkstring = "{} {}".format(os.path.join("/".join(os.getcwd().split('/')[4:]), os.path.basename(__file__)), dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M"))
            plt.figtext(0.01, 0.01, watermarkstring, size=6)
        plt.figtext(0.03, 0.95, "(c)", size=utils.FONTSIZE)

        ax = putils.Rstylee(ax)

        # and save
        if uncertainties:
            outname = putils.make_filenames("timeseries_uncertainties", index=index.name, grid=grid, anomalies=anomalies, month=name)
        else:
            outname = putils.make_filenames("timeseries", index=index.name, grid=grid, anomalies=anomalies, month=name)
            

        plt.savefig("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname), dpi=300)

        plt.close()
        
        # output data file
        if comparison and name == "Ann":
            with open(os.path.join(utils.INFILELOCS, "{}_timeseries.dat".format(index.name)), "w") as outfile:

                outfile.write("{:4s} {:7s} {:7s} {:7s} {:7s}\n".format("Year", "HadEX3", "HadEX2", "HadEX", "GHCNDEX"))
                
                years = timeseries["HadEX3"][0]

                for y, year in enumerate(years):
                    items = [year]

                    for key in ["HadEX3", "HadEX2", "HadEX", "GHCNDEX"]:
                        if key in timeseries.keys():
                            items += [timeseries[key][1][y]]
                        else:
                            items += [utils.HADEX_MDI]

                    outfile.write("{:4d} {:7.2f} {:7.2f} {:7.2f} {:7.2f}\n".format(items[0], items[1], items[2], items[3], items[4]))

        #*****************
        # Plot coverage - how many grid boxes have values (scaled by total number)

        fig = plt.figure(figsize=(8, 5.5))
        plt.clf()
        ax = fig.add_axes([0.15, 0.2, 0.82, 0.75])

        PlotCoverage(plt, ax, valid_boxes, COLOURS, utils.STARTYEAR.year, utils.ENDYEAR.year, index.name, month, '', ncol=2)
        plt.figtext(0.03, 0.95, "(d)", size=utils.FONTSIZE)

        plt.xlim([1900, 2020])

        ax = putils.Rstylee(ax)

        outname = putils.make_filenames("timeseries_coverage", index=index.name, grid=grid, anomalies=anomalies, month=name)

        plt.savefig("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname), dpi=300)

        plt.close("all")

        #*****************
        # plot coverage - how many grid boxes have values (scaled by land fraction)

        fig = plt.figure(figsize=(8, 5.5))
        plt.clf()
        ax = fig.add_axes([0.15, 0.2, 0.82, 0.75])

        PlotCoverage(plt, ax, land_boxes, COLOURS, utils.STARTYEAR.year, utils.ENDYEAR.year, index.name, month, '', ncol=2, land=True)
        plt.figtext(0.03, 0.95, "(d)", size=utils.FONTSIZE)

        plt.xlim([1900, 2020])

        ax = putils.Rstylee(ax)

        outname = putils.make_filenames("timeseries_land_coverage", index=index.name, grid=grid, anomalies=anomalies, month=name)

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
    parser.add_argument('--comparison', dest='comparison', action='store_true', default=False,
                        help='Compare against other datasets if possible')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default=False')
    parser.add_argument('--anomalies', dest='anomalies', action='store', default="None",
                        help='To use the anomalies or climatology rather than actuals, default="None"') 
    parser.add_argument('--grid', dest='grid', action='store', default="ADW",
                        help='gridding routine to use (ADW or CAM), default="ADW"') 
    parser.add_argument('--normalise', dest='normalise', action='store_true', default=False,
                        help='normalise the timeseries to 1961-90, default=False') 
    parser.add_argument('--matched', dest='matched', action='store_true', default=False,
                        help='plot the coverage matched timeseries, default=False') 
    parser.add_argument('--uncertainties', dest='uncertainties', action='store_true', default=False,
                        help='plot the ERA5 derived uncertainties, default=False') 

    args = parser.parse_args()

    if args.anomalies != "climatology":
        main(index=args.index, comparison=args.comparison, diagnostics=args.diagnostics, anomalies=args.anomalies, grid=args.grid, matched=args.matched, normalise=args.normalise, uncertainties=args.uncertainties)
    else:
        print("timeseries not calculable on climatology")

#************************************************************************
#                                 END
#************************************************************************
