#!/usr/bin/ python
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 332                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-04-16 16:06:06 +0100 (Tue, 16 Apr 2019) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""
Plotting utility routines and definitions
"""
import os
import copy
import datetime as dt
import calendar
import numpy as np

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import iris
import iris.coord_categorisation
import iris.plot
import cartopy

import netcdf_procs as ncdfp
import utils

from matplotlib.ticker import MultipleLocator
minorLocator = MultipleLocator(1)

# set up calendar month names
month_names = calendar.month_abbr[:]
month_names[0] = "Ann"

#************************************************************************
def adjust_RdYlBu():
    '''
    Original RdYlBu has a greenish colour on the blue end, which doesn't contrast
    well with the yellow.  As this is a diverging colourmap, want to have a larger
    distinction for the positives and negatives.

    Take the colourmap, move the blues along a bit, and repeat at the end.
    '''

    cmap = plt.cm.RdYlBu

    cmaplist = [cmap(i) for i in range(cmap.N)]

    # ignore first 20 of the blues
    N = 20    
    retained_colours = np.array(cmaplist[255//2 + N :])

    stretched = np.ones((len(cmaplist[255//2:]), 4))

    # interpolate out
    for i in range(stretched.shape[1]):
        stretched[:, i] = np.interp(np.linspace(0, retained_colours.shape[0], stretched.shape[0]), np.arange(retained_colours.shape[0]), retained_colours[:, i])

    # convert to tuple
    new_colours = []
    for c in stretched:
        new_colours += [tuple(c)]

    # copy back
    cmaplist[255//2:] = new_colours

    # enforce white for the centre of the colour range
    for i in [126, 127, 128, 129]:
        cmaplist[i] = (1.0, 1.0, 1.0, 1.0)

    # and make the colour map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap_r = cmap.from_list('Custom cmap', cmaplist[::-1], cmap.N)

    return cmap, cmap_r # adjust_RdYlBu

#************************************************************************
def make_BrBG():
    '''
    Enforce white at centre of BrBG (not the case using the default)
    '''


    cmap = plt.cm.BrBG
    cmaplist = [cmap(i) for i in range(cmap.N)]

    # enforce white for the centre of the colour range
    for i in [126, 127, 128, 129]:
        cmaplist[i] = (1.0, 1.0, 1.0, 1.0)

    # and make the colour map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    cmap_r = cmap.from_list('Custom cmap', cmaplist[::-1], cmap.N)

    return cmap, cmap_r # make_BrBu

#***************************************
def extract_boxes_scatter(cube):
    '''
    e.g. extract which boxes have "significant" trends

    :param array cube: nlat x nlon cube of significant trends (1/0)

    :returns: slons, slats - arrays of significant lat-lon pairs
    '''

    nlat, nlon = cube.data.shape
    slats = []
    slons = []

    for t, lat in enumerate(cube.coord("latitude").points):
        for n, lon in enumerate(cube.coord("longitude").points):
            # cannot be less than -90 lat or 0 lon in this case
            if cube.data[t, n] == 1 and cube.data.mask[t, n] == False:
                
                slats += [lat]
                slons += [lon]

    return slons, slats # ExtractSignificantBoxes

#************************************************************************
def plot_smooth_map_iris(outname, cube, cmap, bounds, cb_label, scatter=[], figtext="", title="", contour=False, cb_extra="", save_netcdf_filename="", significance=None):
    '''
    Standard scatter map

    :param str outname: output filename root
    :param array cube: cube to plot
    :param obj cmap: colourmap to use
    :param array bounds: bounds for discrete colormap
    :param str cb_label: colorbar label
    :param str save_netcdf_filename: filename to save the output plot to a cube.
    '''

    norm = mpl.cm.colors.BoundaryNorm(bounds, cmap.N)

    fig = plt.figure(figsize=(8, 5.5)) # same width as timeseries, but more landscape

    plt.clf()
    ax = plt.axes([0.01, 0.10, 0.98, 0.95], projection=cartopy.crs.Robinson())
    ax.gridlines() #draw_labels=True)
    ax.add_feature(cartopy.feature.LAND, zorder=0, facecolor="0.9", edgecolor="k")
    ax.coastlines("50m")

    ext = ax.get_extent() # save the original extent

    plot_cube = copy.deepcopy(cube)

        
    mesh = iris.plot.pcolormesh(plot_cube, cmap=cmap, norm=norm)
           
    cb = plt.colorbar(mesh, orientation='horizontal', pad=0.03, fraction=0.05, \
                        aspect=30, ticks=bounds[1:-1], drawedges=True)

    cb.set_label(cb_label, size=utils.FONTSIZE)
    # thicken border of colorbar and the dividers
    # http://stackoverflow.com/questions/14477696/customizing-colorbar-border-color-on-matplotlib
    cb.set_ticklabels(["{:g}".format(b) for b in bounds[1:-1]])
    cb.ax.tick_params(labelsize=utils.FONTSIZE, size=0)

#    cb.outline.set_color('k')
    cb.outline.set_linewidth(2)
    cb.dividers.set_color('k')
    cb.dividers.set_linewidth(2)

    if cb_extra != "":
        fig.text(0.04, 0.05, cb_extra[0], fontsize=12, ha="left")
        fig.text(0.96, 0.05, cb_extra[1], fontsize=12, ha="right")
        

    ax.set_extent(ext, ax.projection) # fix the extent change from colormesh

    # plot "significant" boxes with a scatter dot
    if significance is not None:
        sig_x, sig_y = extract_boxes_scatter(significance)
        ax.scatter(sig_x, sig_y, s=0.5, c="k", marker='.', zorder=20, transform=cartopy.crs.PlateCarree())

    plt.title(title, fontsize=utils.FONTSIZE)
    fig.text(0.03, 0.95, figtext, fontsize=utils.FONTSIZE)

    if utils.WATERMARK:
        watermarkstring = "{} {}".format(os.path.join("/".join(os.getcwd().split('/')[4:]), os.path.basename(__file__)), dt.datetime.strftime(dt.datetime.now(), "%d-%b-%Y %H:%M"))
        plt.figtext(0.01, 0.01, watermarkstring, size=6)

    plt.savefig(outname, dpi=300)
    plt.close()

    return # plot_smooth_map_iris

#*********************************************
def make_filenames(name, index="TX90p", grid="ADW", anomalies="None", month="Ann"):
    """
    Build appropriate filenames for all the information

    :param str name: plot name
    :param str index: the index
    :param str grid: the grid
    :param str anomalies: the anomalies
    :param str month: name of month (or timescale)

    :returns: string of filename
    """
    attribs = ncdfp.read_global_attributes(os.path.join(utils.INFILELOCS, "attributes.dat"))
    dataset_name = attribs["title"]

    baseP_string = "{}-{}".format(str(utils.REF_START)[-2:], str(utils.REF_END)[-2:])
    resol_string = "{}x{}deg".format(utils.DELTALAT, utils.DELTALON)

    outfilename = "{}_{}_{}-{}_{}_{}_{}_{}".format(dataset_name, index, utils.STARTYEAR.year, utils.ENDYEAR.year-1, grid, baseP_string, resol_string, month)

    if anomalies != "None":
        outfilename = "{}_{}".format(outfilename, anomalies)

    outfilename = "{}_{}".format(outfilename, name)
    
    outfmt = "png"

    return "{}.{}".format(outfilename, outfmt) # make_filenames


#*******************************
def compute_coverage_error(observations, reanalysis):
    '''
    Calculate the coverage error on a monthly basis

    Takes each month in observations, find all the corresponding calendar months in reanalysis
    (i.e. for Jan 1973, finds all Jans in ERA).  Then mask by data coverage, and get
    range of residuals in global average.  Use these residuals to estimate error

    '''
    
    offset = np.zeros(len(observations.coord("time").points))
    st_dev = np.zeros(len(observations.coord("time").points))

    # add month names into data
    iris.coord_categorisation.add_month(reanalysis, 'time', name='month')
    iris.coord_categorisation.add_month(observations, 'time', name='month')

    for m, month  in enumerate(observations.coord("time").points):

        # get weightings for cubes
        grid_areas = iris.analysis.cartography.cosine_latitude_weights(observations[m]) 
        area_average = observations[m].collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)

        if area_average == observations.data.fill_value:
            # no data for this month
            offset[m] = observations.data.fill_value
            st_dev[m] = observations.data.fill_value
            
        else:
            # make a copy so can update the mask without worrying
            reanal = copy.deepcopy(reanalysis)

            # get a clean-mean
            grid_areas = iris.analysis.cartography.cosine_latitude_weights(reanal)              

            total_mean = reanal.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)

            # apply the observation coverage mask to all months
            # combine the lsm with the data mask
            combined_mask = np.ma.mask_or(np.array([observations[m].data.mask for i in range(reanal.data.shape[0])]), reanal.data.mask)

            reanal.data.mask = combined_mask

            # get a masked mean
            grid_areas = iris.analysis.cartography.cosine_latitude_weights(reanal)              
            masked_mean = reanal.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)
            
            # calculate residuals and find mean and st_dev
            residuals = masked_mean.data - total_mean.data

            offset[m] = np.mean(residuals)
            st_dev[m] = np.std(residuals, ddof=1) # to match IDL

    return offset, st_dev # compute_coverage_error

#***************************************
def SortAxesLabels(plt, ax, index, starttime, endtime, month):
    """Writes all the axis labels, tickmarks and other bits to the plot

    :param object plt: plot instance
    :param object ax: axes instance
    :param str index: name of index being plotted
    :param int starttime: start of x-axis
    :param int endtime: end of x-axis
    :returns: None

    """

    if index.units == "degrees_C":
        units = '$^{\circ}$'+"C"
    else:
        units = index.units

    ax.set_ylabel("{} ({})".format(index.name, units), horizontalalignment='center', size=utils.FONTSIZE)
    ax.set_title("{} - {}".format(index.name, month_names[month]), size=utils.FONTSIZE)
    ax.yaxis.set_minor_locator(minorLocator)

    ax.set_xlim([starttime-1, endtime])
    ax.set_xticks([i for i in range(starttime-1, endtime+10, 10)])
#    ax.xaxis.set_minor_locator(minorLocator)

    ax.tick_params(labelsize=utils.FONTSIZE)

    return # SortAxesLabels

#***************************************
def Rstylee(ax):
    """
    Other settings to make plots "R" style
    """

    # other settings
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.yaxis.set_ticks_position("left")
    ax.xaxis.set_ticks_position("bottom")
    ax.spines["bottom"].set_position(("axes", -0.05))
    ax.spines["left"].set_position(("axes", -0.05))

    return ax # Rstylee


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
