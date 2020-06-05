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
"""
Plots timeseries with optional normalisation and comparison to other datasets using Iris Cubes

Adds in coverage uncertainty from ERA5

plot_stations.py invoked by typing::

  python plot_stations.py --index "TX90p" --qc_flags --anomalies --diagnostics

Input arguments:

--index         Which ETCCDI index to use

--qc_flags      Which QC flags to apply

--anomalies     Run on anomaly outputs

--diagnostics   Output extra info (default False)


"""
import os
import datetime as dt
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
def main(index="TX90p", diagnostics=False, qc_flags="", anomalies="None"):
    """
    Read inventories and make scatter plot

    :param str index: which index to run
    :param bool diagnostics: extra verbose output
    :param str qc_flags: which QC flags to process W, B, A, N, C, R, F, E, V, M
    :param str anomalies: run code on anomalies or climatology rather than raw data

    """

    if index in utils.MONTHLY_INDICES:
        timescale = ["ANN", "MON"]
    else:
        timescale = ["ANN"]

    # move this up one level eventually?
    all_datasets = utils.get_input_datasets()

    for ts in timescale:
        # set up the figure
        fig = plt.figure(figsize=(10, 6.5))
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
                if anomalies == "None":
                    subdir = "formatted/indices"
                elif anomalies == "anomalies":
                    subdir = "formatted/anomalies"
                elif anomalies == "climatology":
                    subdir = "formatted/climatology"

                ds_stations = utils.read_inventory(dataset, subdir=subdir, final=True, \
                                                   timescale=ts, index=index, anomalies=anomalies, qc_flags=qc_flags)
                ds_stations = utils.select_qc_passes(ds_stations, qc_flags=qc_flags)

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
        leg = plt.legend(loc='lower center', ncol=5, bbox_to_anchor=(0.50, -0.3), \
                             frameon=False, title="", prop={'size':12}, labelspacing=0.15, columnspacing=0.5, numpoints=3)
        plt.setp(leg.get_title(), fontsize=12)

        plt.figtext(0.06, 0.91, "{} Stations".format(total))
        plt.title("{} - {}".format(index, ts))

        # extra information
        if utils.WATERMARK:
            watermarkstring = "{} {}".format(os.path.join("/".join(os.getcwd().split('/')[4:]), os.path.basename(__file__)), dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M"))
            plt.figtext(0.01, 0.01, watermarkstring, size=6)
#        plt.figtext(0.03, 0.95, "(c)", size=14)

        # and save
        outname = putils.make_filenames("station_locations", index=index, grid="ADW", anomalies=anomalies, month=ts.capitalize())

        plt.savefig("{}/{}/{}".format(utils.PLOTLOCS, index, outname))
            
        plt.close()

        # write out total station number
        if ts == "ANN":
            with open(os.path.join(utils.INFILELOCS, "{}_stations.txt".format(index)), "w") as outfile:
                outfile.write("{}\n".format(index))
                outfile.write("{}".format(total))
        
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
    parser.add_argument('--anomalies', dest='anomalies', action='store', default="None",
                        help='To use the anomalies or climatology rather than actuals, default="None"') 

    args = parser.parse_args()
          
    main(index=args.index, diagnostics=args.diagnostics, qc_flags=args.qc_flags, anomalies=args.anomalies)


#*******************************************
# END
#*******************************************
