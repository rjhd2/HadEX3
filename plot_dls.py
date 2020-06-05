#!/usr/local/sci/bin/python
#***************************************
#
#   Code to plot summary of DLS values as a gridded plot
#
#   13 December 2012 RJHD
#************************************************************************
#                    SVN Info
# $Rev:: 446                                           $:  Revision of last commit
# $Author:: rdunn                                      $:  Author of last commit
# $Date:: 2020-01-13 16:18:46 +0000 (Mon, 13 Jan 2020) $:  Date of last commit
#************************************************************************
#                                 START
#************************************************************************
import os
import calendar
import datetime as dt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# RJHD routines
import utils


# Definitions/Subroutines
#***************************************
def read_dls(filename):
    '''
    Read in DLS data from input datafile.

    :param string filename: full path to file
    :returns: np array of latitudes, np array of dls's
    '''

    indata = np.genfromtxt(filename, skip_header=4, dtype=(float))

    latitudes = indata[:, 0]

    dls = indata[:, 1:]

    return latitudes, dls


#***************************************
def plot_dls(lats, dls, index):
    '''
    Plot DLS values as grid using pcolor

    :param lats array: array of latitudes
    :param dls array: array of DLS values - monthly and annual as appropriate
    :param index string: name of index
    :returns: outputs png file
    '''
    
    cmap = plt.cm.YlGnBu


    plt.clf()
    fig = plt.figure()
 
    months = len(dls[0, :])

    # if annual or monthly values
    if months == 1:
        
        ax = plt.axes([0.45, 0.10, 0.08, 0.80])
        xaxis = np.array([0, 1])

    else:
        ax = plt.axes([0.10, 0.10, 0.85, 0.80])
        xaxis = np.arange(0, months+1)
    
    if np.min(dls) == np.max(dls):
        # if all values == 200
        grid = plt.pcolor(xaxis, lats, dls, cmap=cmap, vmin=200.1, vmax=1000.)
    else:
        grid = plt.pcolor(xaxis, lats, dls, cmap=cmap, vmin=200.1)
    grid.cmap.set_under('0.5')
    plt.ylabel("Latitudes")

    plt.ylim([-90, 90])
    plt.yticks(np.arange(-90, 90 + 30, 30))

    if months == 1:
        plt.xticks([1.5], ['Ann'], rotation=45)
        plt.xlim([0, 1])
        cb = plt.colorbar(orientation='horizontal', pad=0.1, fraction=0.05, aspect=25, \
                              extend='max', shrink=10)

    else:
        month_names = calendar.month_abbr[:]
        month_names[0] = 'Ann'
        plt.xticks([x + 0.5 for x in np.arange(0, months + 1)], month_names, rotation=45)
        plt.axvline(x=1, color='k', ls='-', lw=2)

        plt.xlim([0, 13])
        cb = plt.colorbar(orientation='horizontal', pad=0.1, fraction=0.05, aspect=25, extend='max')

    cb.set_label("DLS (km)")

    # add text to show what code created this and when
    if utils.WATERMARK:
        watermarkstring = "/".join(os.getcwd().split('/')[4:])+'/'+os.path.basename(__file__)+ \
                          "   "+dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M")
        plt.figtext(0.01, 0.01, watermarkstring, size=6)

    plt.title("DLS - "+index)
    plt.savefig(os.path.join(utils.PLOTLOCS, "DLS", 'dls_mesh_'+index+'.png'))
    plt.close()

    return  # plot_dls

#***************************************
def main():

    # MAIN PROGRAMME STARTS HERE
    # list files and work through one by one

    for dls_file in os.listdir(utils.DLSLOCS):

        index = dls_file.split("_")[1].split(".")[0]

        latitudes, dls_data = read_dls(os.path.join(utils.DLSLOCS, dls_file))

        plot_dls(latitudes, dls_data, index)

        print(index+' done')

    print('Finished')

    return # main

#************************************************************************
if __name__ == "__main__":

    main()

#************************************************************************
#                                 END
#************************************************************************
