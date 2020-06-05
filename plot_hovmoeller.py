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
import numpy as np
import datetime as dt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

import iris
import iris.quickplot as qplt

# local utilities scripts
import utils
import plt_utils as putils

# set up calendar month names
month_names = calendar.month_abbr[:]
month_names[0] = "Ann"


#************************************************************************

def main(index, diagnostics=False, normalise=True, anomalies="None", grid="ADW"):
    """
    :param str index: which index to run
    :param bool diagnostics: output diagnostic information
    :param str anomalies: run code on anomalies or climatology rather than raw data
    :param str grid: gridding type ADW/CAM
    """
    
    cosine = False

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

    # assign bounds and colormaps
    if index.name in ["TX90p", "TN90p", "SU", "TR", "GSL"]:
        bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
        cmap = RdYlBu_r
    elif index.name in ["DTR", "ETR"]:
        bounds = [-100, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 100]
        cmap = RdYlBu_r
    elif index.name in ["WSDI"]:
        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
        cmap = RdYlBu_r
    elif index.name in ["TX10p", "TN10p", "FD", "ID"]:
        bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
        cmap = RdYlBu
    elif index.name in ["CSDI"]:
        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
        cmap = RdYlBu
    elif index.name in ["TXn", "TNn"]:
        bounds = [-100, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 100]
        cmap = RdYlBu_r
    elif index.name in ["TXx", "TNx"]:
        bounds = [-100, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 100]
        cmap = RdYlBu_r
    elif index.name in ["Rx1day"]:
        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
        cmap = BrBG
    elif index.name in ["Rx5day"]:
        bounds = [-100, -4, -3, -2, -1, 0, 1, 2, 3, 4, 100]
        cmap = BrBG
    elif index.name in ["PRCPTOT"]:
        bounds = [-100, -20, -10, -5, -2, 0, 2, 5, 10, 20, 100]
        cmap = BrBG
    elif index.name in ["Rnnmm", "R95p", "R99p"]:
        bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
        cmap = BrBG
    elif index.name in ["R95pTOT"]:
        bounds = [-100, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 100]
        cmap = BrBG
    elif index.name in ["R99pTOT"]:
        bounds = [-100, -1, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 1, 100]
        cmap = BrBG
    elif index.name in ["R10mm"]:
        bounds = [-100, -3, -1.5, -0.75, -0.25, 0, 0.25, 0.75, 1.5, 3, 100]
        cmap = BrBG
    elif index.name in ["R20mm"]:
        bounds = [-100, -2, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 2, 100]
        cmap = BrBG
    elif index.name in ["CWD"]:
        bounds = [-100, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 100]
        cmap = BrBG
    elif index.name in ["SDII"]:
        bounds = [-100, -0.75, -0.5, -0.25, -0.1, 0, 0.1, 0.25, 0.5, 0.75, 100]
        cmap = BrBG
    elif index.name in ["CDD"]:
        bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
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
        bounds = [-100, -8, -4, -2, -1, 0, 1, 2, 4, 8, 100]
        cmap = RdYlBu_r

    norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

    cube_list = iris.load(os.path.join(utils.FINALROOT, utils.make_filenames(index=index.name, grid=grid, anomalies=anomalies, extra="", month_index="")))    

    names = np.array([cube.var_name for cube in cube_list])

    # plot all month versions at once
    for month, mname in enumerate(month_names):

        if diagnostics:
            print(mname)

        selected_cube, = np.where(names == mname)

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


        # Take the mean over latitude
        cube = cube.collapsed('grid_longitude', iris.analysis.MEAN)

        # if show relative to climatology
        if normalise:            
            clim_constraint = iris.Constraint(time=lambda cell: utils.REF_START <= cell <= utils.REF_END)
            norm_cube = cube.extract(clim_constraint)
            norm_cube = norm_cube.collapsed(['time'], iris.analysis.MEAN)

            cube = cube - norm_cube

        # plot
        # set up the figure
        fig = plt.figure(figsize=(8, 6))
        plt.clf()
        ax = fig.add_axes([0.1, 0.1, 0.85, 0.85])
        ax.patch.set_facecolor("0.8")

        contour = iris.plot.pcolor(cube, cmap=cmap, norm=norm)#, vmax=bounds[-2], vmin=bounds[1])

        cb = plt.colorbar(contour, orientation='horizontal', pad=0.07, fraction=0.05, \
                            aspect=30, ticks=bounds[1:-1], drawedges=True)

        cb.set_label(index.units, size=utils.FONTSIZE)
        # thicken border of colorbar and the dividers
        # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
        cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
        cb.ax.tick_params(labelsize=utils.FONTSIZE, size=0)

#        cb.outline.set_color('k')
        cb.outline.set_linewidth(2)
        cb.dividers.set_color('k')
        cb.dividers.set_linewidth(2)

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(utils.FONTSIZE)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(utils.FONTSIZE)

        ax.set_xlim([1900,2020])

        if cosine:
            ax.set_ylim(np.sin(np.deg2rad(np.array([-90, 90]))))
            ax.set_yticks(np.sin(np.deg2rad(np.array([-90, -60, -30, 0, 30, 60, 90]))))
            ax.set_yticklabels(["-90"+r'$^{\circ}$'+"S", "-60"+r'$^{\circ}$'+"S", \
                                    "-30"+r'$^{\circ}$'+"S", "0"+r'$^{\circ}$'+"", \
                                    "30"+r'$^{\circ}$'+"N", "60"+r'$^{\circ}$'+"N", "90"+r'$^{\circ}$'+"N"], fontsize=utils.FONTSIZE)
        else:
            ax.set_ylim([-90, 90])
            ax.set_yticks([-60, -30, 0, 30, 60])
            ax.set_yticklabels(["-60"+r'$^{\circ}$'+"S", "-30"+r'$^{\circ}$'+"S", \
                                    "0"+r'$^{\circ}$'+"", "30"+r'$^{\circ}$'+"N", "60"+r'$^{\circ}$'+"N"], fontsize=utils.FONTSIZE)


        plt.title("{} - {}, Hovmöller".format(index.name, month_names[month]), fontsize=utils.FONTSIZE)
        fig.text(0.03, 0.95, "(e)", fontsize=utils.FONTSIZE)

        if utils.WATERMARK:
            watermarkstring = "{} {}".format(os.path.join("/".join(os.getcwd().split('/')[4:]), os.path.basename(__file__)), dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M"))
            plt.figtext(0.01, 0.01, watermarkstring, size=6)

        outname = putils.make_filenames("hovmoeller", index=index.name, grid=grid, anomalies=anomalies, month=mname)
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
                        help='gridding routine to use (ADW or CAM), default="ADW"') 
    args = parser.parse_args()

    if args.anomalies != "climatology":
        main(index=args.index, diagnostics=args.diagnostics, anomalies=args.anomalies, grid=args.grid)
    else:
        print("hovmoeller not calculable on climatology")

#************************************************************************
#                                 END
#************************************************************************
