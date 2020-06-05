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
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy
import numpy as np
import time

import utils
import plot_stations as ps
import plt_utils as putils

#*********************************************
def time_presence(station, index, timescale):
    
    # get data
    times, indata = utils.read_station_index(station, index, timescale)

    # match times
    match = np.in1d(utils.REFERENCEYEARS, times)       
    match_back = np.in1d(times, utils.REFERENCEYEARS)       
    
    # sort array
    if timescale == "ANN":
        data = np.zeros([len(utils.REFERENCEYEARS), 1])
    elif timescale == "MON":
        data = np.zeros([len(utils.REFERENCEYEARS), 12])
        
    data[:] = utils.HADEX_MDI

    if len(indata.shape) == 1:
        indata = indata.reshape([indata.shape[0], 1])
        data[match, :] = indata[match_back] # store all the info

    data[data != utils.HADEX_MDI] = 1
    data[data == utils.HADEX_MDI] = 0

    return data.astype(int) # time_presence

#*********************************************
def main(index="TX90p", diagnostics=False, qc_flags="", anomalies="None"):
    """
    Read inventories and make scatter plot

    :param str index: which index to run
    :param bool diagnostics: extra verbose output
    :param str qc_flags: which QC flags to process W, B, A, N, C, R, F, E, V, M
    :param str anomalies: run code on anomalies or climatology rather than raw data

    """
    with open(os.path.join(utils.INFILELOCS, "{}_yearly_stations.txt".format(index)), "w") as outfile:
        outfile.write("{}\n".format(index))

    if index in utils.MONTHLY_INDICES:
        timescale = ["ANN", "MON"] # allow for future!
    else:
        timescale = ["ANN"]

    # move this up one level eventually?
    all_datasets = utils.get_input_datasets()

    for ts in timescale:

        # run all datasets
        for d, dataset in enumerate(all_datasets):
            
            print(dataset)

            try:
                # choose appropriate subdirectory.
                subdir = "formatted/indices"

                ds_stations = utils.read_inventory(dataset, subdir=subdir, final=True, \
                                                       timescale=ts, index=index, anomalies=anomalies, qc_flags=qc_flags)
                ds_stations = utils.select_qc_passes(ds_stations, qc_flags=qc_flags)
                
            except IOError:
                # file missing
                print("No stations with data for {}".format(dataset.name))
                ds_stations = []

            # extract relevant info for this dataset
            if len(ds_stations) > 0:
                
                # extract values for this dataset
                for s, stn in enumerate(ds_stations):
                    presence = time_presence(stn, index, ts) # year/month
                    if s == 0:
                        ds_presence = np.expand_dims(presence, axis=0)[:]
                    else:
                        ds_presence = np.append(ds_presence, np.expand_dims(presence, axis=0), axis=0) # station/year/month

                ds_lats = np.array([stn.latitude for stn in ds_stations])
                ds_lons = np.array([stn.longitude for stn in ds_stations])

                # store in overall arrays
                try:
                    all_lats = np.append(all_lats, ds_lats[:], axis=0)
                    all_lons = np.append(all_lons, ds_lons[:], axis=0)
                    all_presence = np.append(all_presence, ds_presence[:], axis=0) # dataset*station/year/month
                    all_dataset_names = np.append(all_dataset_names, np.array([dataset.name for i in ds_lats]))
                except NameError:
                    # if not yet defined, then set up
                    all_lats = ds_lats[:]
                    all_lons = ds_lons[:]
                    all_presence = ds_presence[:]
                    all_dataset_names = np.array([dataset.name for i in ds_lats])


        for y, year in enumerate(utils.REFERENCEYEARS):

            # set up the figure
            fig = plt.figure(figsize=(10, 6.5))
            plt.clf()
            ax = plt.axes([0.025, 0.10, 0.95, 0.90], projection=cartopy.crs.Robinson())
            ax.gridlines() #draw_labels=True)
            ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
            ax.coastlines()

            # dummy scatters for full extent
            plt.scatter([-180, 180, 0, 0], [0, 0, -90, 90], c="w", s=1, transform=cartopy.crs.Geodetic(), \
                            edgecolor='w', linewidth='0.01')

            total = 0
            for dataset in all_datasets:
                
                ds, = np.where(all_dataset_names == dataset.name)
                locs, = np.where(all_presence[ds, y, 0] == 1)

                if len(locs) > 0:
                    plt.scatter(all_lons[ds][locs], all_lats[ds][locs], c=ps.COLOURS[dataset.name], \
                                    s=15, label="{} ({})".format(ps.get_label(dataset.name), len(locs)), \
                                    transform=cartopy.crs.Geodetic(), edgecolor='0.5', linewidth='0.5')
                    total += len(locs)
                else:
                    # aiming to show all, even if zero
                    plt.scatter([-180], [-90], c=ps.COLOURS[dataset.name], s=15, \
                                    label="{} ({})".format(ps.get_label(dataset.name), len(locs)), \
                                    transform=cartopy.crs.Geodetic(), edgecolor='0.5', linewidth='0.5')
                time.sleep(1)

            # make a legend
            leg = plt.legend(loc='lower center', ncol=6, bbox_to_anchor=(0.50, -0.25), frameon=False, \
                                 title="", prop={'size':10}, labelspacing=0.15, columnspacing=0.5, numpoints=3)
            plt.setp(leg.get_title(), fontsize=12)

            plt.figtext(0.05, 0.92, "{} Stations".format(total))

            plt.title("{} - {} - {}".format(index, ts, year))

            # and save
            outname = putils.make_filenames("station_locations_{}_{}".format(ts.capitalize(), year), index=index, grid="ADW", anomalies=anomalies)

            plt.savefig("{}/{}/{}".format(utils.PLOTLOCS, index, outname))
            
            plt.close()
            plt.clf()
            print("{} done".format(year))

            # write out total station number
            with open(os.path.join(utils.INFILELOCS, "{}_yearly_stations.txt".format(index)), "a") as outfile:
                outfile.write("{} {}\n".format(year, total))

            time.sleep(1)

        # reset namespace
        del all_lats
        del all_lons
        del all_presence
        del all_dataset_names
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
