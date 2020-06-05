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
# $Rev:: 252                                           $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2018-11-14 12:13:22 +0000 (Wed, 14 Nov 2018) $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import os
import calendar
import netCDF4 as ncdf
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

# local utilities scripts
import utils
import plt_utils as putils

# set up calendar month names
month_names = calendar.month_abbr[:]
month_names[0] = "Ann"


#************************************************************************

def main(index, diagnostics=False, anomalies="None", grid="ADW"):
    """
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
        timescale = "MON"
    else:
        nmonths = 1
        timescale = "ANN"

    # setting up colours
    cmap = plt.cm.viridis
    bounds = np.arange(0, 110, 10)
    norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)


    # plot all month versions at once
    for month, name in enumerate(month_names):

        if diagnostics:
            print(name)

        # set up the figure
        fig = plt.figure(figsize=(8, 6))
        plt.clf()
        ax = fig.add_axes([0.1, 0.1, 0.85, 0.85])

        filename = os.path.join(utils.FINALROOT, utils.make_filenames(index=index.name, grid=grid, anomalies=anomalies, extra="", month_index=""))

        try:
            ncfile = ncdf.Dataset(filename, 'r')
            
            timevar = ncfile.variables['time'] # array of YYYMMDD
            latvar = ncfile.variables['latitude'] # array of lats
            lonvar = ncfile.variables['longitude'] # array of lons
            
            annualvar = ncfile.variables[month_names[month]] # array of arrays

            if anomalies == "anomalies":
                # to make into actuals, add climatology to the anomalies
                clim_filename = os.path.join(utils.FINALROOT, utils.make_filenames(index=index.name, grid=grid, anomalies="climatology", extra="", month_index=""))

                clim_file = ncdf.Dataset(clim_filename, 'r')
                climvar = clim_file.variables[month_names[month]]
        except RuntimeError:
            print("File not found: {}".format(filename))
            
        except IOError:
            print("File not found: {}".format(filename))
            
        except KeyError:
            continue

        # extract the information
        times = timevar[:]
        lats = latvar[:]
        lons = lonvar[:]

        # get land sea mask
        if month == 0:
            lsm = utils.get_land_sea_mask(lats, lons, floor=False)
            n_land_boxes = np.sum(lsm.astype(int), axis=1).astype(float)
            n_land_boxes = np.ma.expand_dims(n_land_boxes, axis=1)

        # if anomalies from HadEX3, then need to add onto climatology
        if anomalies == "anomalies":
            annual_data = annualvar[:] + climvar[:]
        else:
            annual_data = annualvar[:]

        # go through each year and count up
        zonal_boxes = np.zeros(annual_data.shape[:2][::-1])
        for y, year_data in enumerate(annual_data):
            zonal_boxes[:, y] = np.ma.count(year_data, axis=1).astype(float)

        zonal_boxes = np.ma.masked_where(zonal_boxes == 0, zonal_boxes)
            
        normalised_boxes = 100. * zonal_boxes / np.tile(n_land_boxes, [1, annual_data.shape[0]])

        newtimes, newlats = np.meshgrid(np.append(times, times[-1]+1), utils.box_edge_lats)
            
        mesh = plt.pcolormesh(newtimes, newlats, normalised_boxes, cmap=cmap, norm=norm)
        
        cb = plt.colorbar(mesh, orientation='horizontal', pad=0.07, fraction=0.05, \
                            aspect=30, ticks=bounds[1:-1], label="%", drawedges=True)
        
        cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
        cb.outline.set_linewidth(2)
        cb.dividers.set_color('k')
        cb.dividers.set_linewidth(2)

        plt.ylim([-90, 90])
        plt.yticks(np.arange(-90, 120, 30))
        ax.yaxis.set_ticklabels(["{}N".format(i) if i > 0 else "{}S".format(abs(i)) if i < 0 else "{}".format(i) for i in np.arange(-90, 120, 30)])
        plt.ylabel("Latitude")
        plt.title("{} - {}".format(index.name, name))

        outname = putils.make_filenames("latitude_coverage", index=index.name, grid=grid, anomalies=anomalies, month=name)

        plt.savefig("{}/{}/{}".format(utils.PLOTLOCS, index.name, outname))

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
                        help='gridding routine to use (ADW or CAM), default="ADW"') 
    args = parser.parse_args()

    if args.anomalies != "climatology":
        main(index=args.index, diagnostics=args.diagnostics, anomalies=args.anomalies, grid=args.grid)
    else:
        print("timeseries not calculable on climatology")

#************************************************************************
#                                 END
#************************************************************************
