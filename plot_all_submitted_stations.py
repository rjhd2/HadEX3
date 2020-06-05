#!/usr/local/sci/bin/python
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 106                                     $:  Revision of last commit
#$Author::                                      $:  Author of last commit
#$Date:: 2018-05-22 10:06:47 +0100 (Tue, 22 May#$:  Date of last commit
#------------------------------------------------------------
#  
#  Plot station networks  
#
#------------------------------------------------------------
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy

import utils
import plt_utils as putils


COLOURS = {"ghcnd" : "purple", \
               "ecad" : "lime", \
               "lacad" : "lime", \
               "sacad" : "lime", \
               "eobs" : "lime", \
               "laobs" : "lime", \
               "saobs" : "lime", \
               "spain" : "red", \
               "russia" : "blue", \
               "hadex2" : "cyan", \
               "acre" : "yellow", \
               "honduras" : "red", \
               "pacific" : "purple", \
               "nz" : "red", \
               "south_america" : "orange", \
           "west_africa_pptn" : "blue", \
           "west_africa_indices" : "red", \
           "south_africa" : "red", \
           "australia" : "magenta", \
           "arabia" : "orange",\
           "chile" : "purple", \
           "colombia" : "blue", \
           "canada" : "yellow", \
           "decade" : "magenta", \
           "india" : "magenta", \
           "malaysia" : "orange", \
           "brunei" : "orange", \
           "myanmar" : "orange", \
           "vietnam" : "orange", \
           "philippines" : "orange", \
           "singapore" : "orange", \
           "thailand" : "orange", \
           "indonesia" : "orange", \
           "china" : "red", \
           "japan" : "magenta", \
           "mexico" : "blue", \
           "iran" : "yellow", \
           "brazil" : "blue", \
           "brazil_sp" : "red", \
           "ghcndex" : "0.5"\
}

#*********************************************
def get_label(name):
    """
    Fix dataset name appropriately
    """
    
    # specific fixes
    if name == "nz":
        new_name = "New Zealand"
    elif name in ["ghcnd", "ecad", "lacad", "sacad", "decade", "acre", "ghcndex"]:
        new_name = name.upper()
    elif name == "hadex2":
        new_name = "HadEX2"
    elif name == "west_africa_pptn":
        new_name = "West Africa 1"
    elif name == "west_africa_indices":
        new_name = "West Africa 2"
   
    # generic fixes
    elif "_" in name:
        new_name = " ".join([n.capitalize() for n in name.split("_")])

    else:
        new_name = name.capitalize()
        


    return new_name

#*********************************************
def main(diagnostics=False):
    """
    Read inventories and make scatter plot

    :param bool diagnostics: extra verbose output

    """

    # move this up one level eventually?
    all_datasets = utils.get_input_datasets()

    # set up the figure
    fig = plt.figure(figsize=(10, 6.7))
    plt.clf()
    ax = plt.axes([0.025, 0.14, 0.95, 0.90], projection=cartopy.crs.Robinson())
    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
    ax.coastlines()

    # dummy scatters for full extent
    plt.scatter([-180, 180, 0, 0], [0, 0, -90, 90], c="w", s=1, transform=cartopy.crs.Geodetic(), \
                    edgecolor='w', linewidth='0.01')

    # run all datasets
    total = 0
    for dataset in all_datasets:

        try:
            # choose appropriate subdirectory.
            subdir = "formatted/indices"

            ds_stations = utils.read_inventory(dataset, subdir=subdir, final=False, \
                                               timescale="", index="", anomalies="None", qc_flags="")

        except IOError:
            # file missing
            print("No stations with data for {}".format(dataset.name))
            ds_stations = []

        if len(ds_stations) > 0:
            lats = np.array([stn.latitude for stn in ds_stations])
            lons = np.array([stn.longitude for stn in ds_stations])

            # and plot
            scatter = plt.scatter(lons, lats, c=COLOURS[dataset.name], s=15, \
                                      label="{} ({})".format(get_label(dataset.name), len(ds_stations)), \
                                      transform=cartopy.crs.Geodetic(), edgecolor='0.5', linewidth='0.5')

            total += len(ds_stations)

    # make a legend
    leg = plt.legend(loc='lower center', ncol=5, bbox_to_anchor=(0.50, -0.34), \
                         frameon=False, title="", prop={'size':12}, labelspacing=0.15, columnspacing=0.5, numpoints=3)
    plt.setp(leg.get_title(), fontsize=12)

    plt.figtext(0.05, 0.92, "{} Stations".format(total))

    plt.title("HadEX3 stations")

    # and save
    outname = putils.make_filenames("station_locations", index="All", grid="ADW", anomalies="None", month="All")

    plt.savefig("{}/{}".format(utils.PLOTLOCS, outname), dpi=300)

    plt.close()

    return # main


#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default=False')

    args = parser.parse_args()
          
    main(diagnostics=args.diagnostics)


#*******************************************
# END
#*******************************************
