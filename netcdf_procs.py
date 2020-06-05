#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 450                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2020-01-15 15:58:00 +0000 (Wed, 15 Jan 2020) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""
Routines used for netCDF reading/writing
"""
import os
import calendar
import datetime as dt
import netCDF4 as ncdf
import numpy as np

# RJHD utilities
import utils

#************************************************************************
def read_global_attributes(attr_file):
    '''
    Reads attributes file and returns a dictionary of key/values
    '''
        
    try:
        with open(attr_file, 'r') as infile:        
            lines = infile.readlines()
        
    except IOError:
        print("Attributes file not found at " + attr_file)
    
    
    attributes = {}
    
    for line in lines:
        split_line = line.split()
        
        attributes[split_line[0]] = " ".join(split_line[1:])    
        
    return attributes # read_global_attributes

#*********************************************
def netcdf_read(filename, month=0):
    """
    Read the individual monthly netCDF files
    """

    ncfile = ncdf.Dataset(filename, 'r')
    
    nclon = ncfile.variables['longitude'][:]
    nclat = ncfile.variables['latitude'][:]
    nctime = ncfile.variables['time'][:]

    month_names = calendar.month_abbr[0:]
    month_names[0] = "Ann"

#    indata = ncfile.variables[month_names[month]]
    indata = ncfile.variables[month_names[month]][:]

    ncfile.close()
 
    return indata[:], nclon, nclat, nctime # netcdf_read

#************************************************************************
def write_coordinates(outfile, nc_var, short_name, standard_name, long_name, units, axis, data, coordinate_length=1, do_zip=True, least_significant_digit=None):
    """
    Write coordinates as variables

    :param str outfile: output netcdf file
    :param str short_name: netcdf short_name
    :param str standard_name: netcdf standard_name
    :param str long_name: netcdf long_name
    :param str units: netcdf units
    :param str axis: netcdf axis
    :param flt data: coordinate 
    :param int coordinate_length: length of dimension
    :param bool do_zip: allow for zipping
    :param int least_significant_digit: smallest reliable decimal place
    """

    try:
        nc_var.standard_name = standard_name
    except AttributeError:
        pass
    nc_var.long_name = long_name
    nc_var.units = units
    nc_var.axis = axis

    nc_var[:] = data

    return # write_coordinates

#*********************************************
def netcdf_write(filename, index, data, times, lats, lons, single_month=-1, station_count=False):
    """
    Write the full array to file.
    """

    outfile = ncdf.Dataset(filename, 'w', format='NETCDF4', clobber=True)

    lon_dim = outfile.createDimension('longitude', len(lons))
    lat_dim = outfile.createDimension('latitude', len(lats))
    if data.shape[0] == 1: # climatology array
        time_dim = outfile.createDimension('time', 1)
    else:
        time_dim = outfile.createDimension('time', len(times))
    
    ncbnd = outfile.createDimension('ncbnd', 2) # number of bounds

    lonvar = outfile.createVariable('longitude', np.dtype('float'), ('longitude',))
    latvar = outfile.createVariable('latitude', np.dtype('float'), ('latitude',))
    timevar = outfile.createVariable('time', np.dtype('int'), ('time',))

    # set up the bounds
    lonBvar = outfile.createVariable('longitude_bounds', 'f4', ('longitude', 'ncbnd'))
    latBvar = outfile.createVariable('latitude_bounds', 'f4', ('latitude', 'ncbnd'))
    # and associate with the coordinates
    lonvar.bounds = "longitude_bounds"
    latvar.bounds = "latitude_bounds"

    # do these need to be rolled?

    # write the coordinate values
    write_coordinates(outfile, lonvar, "longitude", "grid_longitude", "longitudes of grid box centres", "degrees_north", "X", lons)
    write_coordinates(outfile, latvar, "latitude", "grid_latitude", "latitude of grid box centres", "degrees_east", "Y", lats)
    if data.shape[0] == 1: # climatology array
        write_coordinates(outfile, timevar, "time", "time", "time", "years", "T", utils.REF_START)
    else:
        write_coordinates(outfile, timevar, "time", "time", "time", "years", "T", times)

    # and the bounds
    lonBvar[:] = np.array(list(zip(utils.box_edge_lons[:-1], utils.box_edge_lons[1:])))
    latBvar[:] = np.array(list(zip(utils.box_edge_lats[:-1], utils.box_edge_lats[1:])))

    # no need to give full information for bounds, but kept for reference
    # write_coordinates(outfile, lonBvar, "longitude_bounds", "grid_longitude_bounds", "longitudes of grid box edges", \
    #                   "degree", "X", np.array(list(zip(utils.box_edge_lons[:-1], utils.box_edge_lons[1:]))))
    # write_coordinates(outfile, latBvar, "latitude_bounds", "grid_latitude_bounds", "latitudes of grid box edges", \
    #                   "degree", "Y", np.array(list(zip(utils.box_edge_lats[:-1], utils.box_edge_lats[1:]))))

    month_names = calendar.month_abbr[0:]
    month_names[0] = "Ann"
    
    # if not a monthly index, just retain the first element
    if index not in utils.MONTHLY_INDICES:
        month_names = month_names[:1]

    # if writing just a single month
    if single_month != -1:
        month_names = [month_names[single_month]] # need to keep as list!

    for m, name in enumerate(month_names):
        datavar = outfile.createVariable(name, 'f4', ('time', 'latitude', 'longitude'), fill_value=utils.HADEX_MDI)
        datavar.missing_value = np.float64(utils.HADEX_MDI)
        datavar.FillValue = np.float64(utils.HADEX_MDI)
        datavar[:] = data[:, m, :, :]  
        datavar.long_name = utils.INDICES[index].long_name
        datavar.standard_name = index

        if station_count:
            datavar.units = "1"
            datavar.definition = "Number of stations within each grid box"
        else:
            datavar.units = utils.INDICES[index].units
            datavar.definition = utils.INDICES[index].definition

        datavar.team = utils.INDICES[index].team
        datavar.coordinates = "latitude longitude time"
    
    # Global Attributes
    # from file
    attribs = read_global_attributes(os.path.join(utils.INFILELOCS, "attributes.dat"))
    for attr in attribs:
        outfile.__setattr__(attr, attribs[attr])
    outfile.history = "Created by " + os.path.basename(__file__)
    outfile.title = index
    outfile.date_created = dt.datetime.strftime(dt.datetime.now(), "%a %b %d, %H:%M %Y")
    outfile.Conventions = 'CF-1.6' 
    outfile.Metadata_Conventions = 'Unidata Dataset Discovery v1.0,CF Discrete Sampling Geometries Conventions'
#    outfile.processing_date = processing_date

    outfile.geospatial_lat_resolution = utils.DELTALAT
    outfile.geospatial_lon_resolution = utils.DELTALON
    if utils.INDICES[index].monthly == True:
        outfile.time_coverage_resolution = "Monthly"
    else:
        outfile.time_coverage_resolution = "Annual"

    outfile.close()

    print("{} written".format(filename))

    return # netcdf_write_monthly

#------------------------------------------------------------
# END
#------------------------------------------------------------
