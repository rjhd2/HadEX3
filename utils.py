#!/usr/bin python
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 489                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2020-04-23 16:38:49 +0100 (Thu, 23 Apr 2020) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""
Utility routines and definitions
"""
import os
import sys
import datetime as dt
import math
import glob
import shutil
import configparser
import json
import numpy as np

import netcdf_procs as ncdfp


#******************************************************************************************
#******************************************************************************************
class cd:
    """
    Context manager for changing the current working directory

    https://stackoverflow.com/questions/431684/how-do-i-cd-in-python/13197763#13197763
    """
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

#*********************************************
class Dataset(object):
    '''
    Class to hold information for the input Dataset

    '''
    def __init__(self, name, version, base_period):
        self.name = name
        self.version = version
        self.base_period = base_period
        self.location = os.path.join(INROOT, name, version)
        
        if not os.path.exists(os.path.join(INROOT, name)):
            os.mkdir(os.path.join(INROOT, name))

        if not os.path.exists(self.location):
            os.mkdir(self.location)
            os.mkdir(os.path.join(self.location, "raw"))
            os.mkdir(os.path.join(self.location, "formatted"))
        

    def __str__(self):     
        return "name: {}, version: {}".format(self.name, self.version)

    __repr__ = __str__
   
#*********************************************
class Station(object):
    '''
    Class to hold information for the station

    '''
    def __init__(self, ID, latitude, longitude, location, dataset):
        self.id = ID
        self.latitude = latitude
        self.longitude = longitude
        self.location = location
        self.source = dataset
        self.qc = ""
        

    def __str__(self):     
        return "station: {}".format(self.id)

    __repr__ = __str__
   
#*********************************************
class Timeseries(object):
    '''
    Class for timeseries
    '''
    
    def __init__(self, name, times, data):
        self.name = name
        self.times = times
        self.data = data
        

    def __str__(self):     
        return "timeseries of {}".format(self.name)

    __repr__ = __str__

#*********************************************
class Index(object):
    '''
    Class for timeseries
    '''
    
    def __init__(self, name, long_name, units, definition, team, monthly, baseperiod):
        self.name = name
        self.long_name = long_name
        self.units = units
        self.definition = definition
        self.team = team
        self.monthly = monthly
        self.baseperiod = baseperiod
        

    def __str__(self):     
        return "{} index".format(self.name)

    __repr__ = __str__

#******************************************************************************************
#******************************************************************************************
def copy(src, dest):
    print(r'{}*'.format(os.path.expanduser(src)))
    for filename in glob.glob(r'{}*'.format(os.path.expanduser(src))):
        if not os.path.exists("{}/{}".format(dest, filename.split("/")[-1])):
            shutil.copy(filename, dest)
        print(filename, dest)

    return # copy
#*********************************************
def get_input_datasets():
    '''
    Automatically read through all the input datasets and store info

    :return: array of datasets
    '''
    baseP_string = "{}-{}".format(str(REF_START)[-2:], str(REF_END)[-2:])
    initial_datasets = []

    # spin through input datasets
    with open(os.path.join(INFILELOCS, "input_datasets.dat"), "r", encoding="latin-1") as infile:
        
        for line in infile:

            split = line.split()

            if split[0] == "singapore":
                # both reference periods exist, so can choose to match
                initial_datasets += [Dataset(split[0], split[1], baseP_string)]

            else:
                initial_datasets += [Dataset(split[0], split[1], split[2])]

    initial_datasets = np.array(initial_datasets)

    # remove those where have two versions, one for each base period
    names = [ds.name for ds in initial_datasets]
    duplicated = set([x for x in names if names.count(x) > 1])

    datasets = []
    for ds in initial_datasets:

        if ds.name in duplicated:
            if ds.base_period == baseP_string:
                # only select if matching base period
                datasets += [ds]
        else:
            datasets += [ds]

    return np.array(datasets)
 
#*********************************************
def write_climpact_inventory_header(filename):
    '''
    Write the inventory header in HadEX format

    :param str filename: name of output file (full path)
    '''

    # makes new file (overwrites)
    with open(filename, "w") as outfile:
        outfile.write("{:20s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s} {:>10s}\n".format("station_file", "latitude", "longitude", "wsdin", "csdin", "Tb_HDD", "Tb_CDD", "Tb_GDD", "rx_ud", "rnnmm_ud", "txtn_ud", "SPEI"))

    outfile.close()

    return # write_climpact_inventory_header

#*********************************************
def write_climpact_inventory(filename, station):
    '''
    Write the inventory in HadEX format

    :param str filename: name of output file (full path)
    :param station station: station to write
    '''

    with open(filename, "a") as outfile:

        outfile.write("{:20s} {:10.2f} {:10.3f} {:10d} {:10d} {:10d} {:10d} {:10d} {:10d} {:10d} {:10d} {:10d}\n".format("{}.txt".format(station.id), station.latitude, station.longitude, WSDI_N, CSDI_N, TB_HDD, TB_CDD, TB_GDD, RX_UD, RNNMM_UD, TXTN_UD, SPEI))

    outfile.close()

    return # write_climpact_inventory


#*********************************************
def standardise_timeseries(station):
    '''
    Merge the parameter individual timeseries together

    :param station station: station object to process
    '''

    #******************
    def overwrite_arrays(all_times, ts):
        '''
        Overwrite the data with the standardised and padded times and data

        :param array all_times: final list of unique times
        :param timeseries ts: timeseries object to amend
        '''
 
        # matching is quicker with arrays of integers, so extract timedeltas
        start = np.min([all_times[0], ts.times[0]])
        all_times_diffs = np.array([(at - start).days for at in all_times])
        ts_diffs = np.array([(t - start).days for t in ts.times])
        
        match = np.in1d(all_times_diffs, ts_diffs)

        newdata = np.zeros(all_times.shape)
        newdata.fill(HADEX_MDI)
        newdata[match] = ts.data

        ts.data = newdata
        ts.times = all_times

        return ts # overwrite_arrays

    
    all_times = station.tmax.times
    all_times = np.append(all_times, station.tmin.times)
    all_times = np.append(all_times, station.prcp.times)

    # returns a unique set of dates
    all_times = np.unique(all_times)

    station.tmax = overwrite_arrays(all_times, station.tmax)
    station.tmin = overwrite_arrays(all_times, station.tmin)
    station.prcp = overwrite_arrays(all_times, station.prcp)

    return station # standardise_timeseries

#*********************************************
def write_climpact_format(location, station):
    '''
    Write the output files in HadEX/Climpact format

    :param str location: path location to station data
    :param station station: station object
    '''

    station = standardise_timeseries(station)

    with open(os.path.join(location, "{}.txt".format(station.id)), "w") as outfile:

        for t in range(len(station.tmax.times)):

            outfile.write("{:4d}{:6d}{:8d}{:11.1f}{:8.1f}{:8.1f}\n".format(station.tmax.times[t].year, station.tmax.times[t].month, station.tmax.times[t].day, station.prcp.data[t], station.tmax.data[t], station.tmin.data[t]))
                    
    outfile.close()

    return # write_climpact_format

#*********************************************
def write_station_index(filename, station, index, doMonthly=False):
    """
    Write data held in station object into climpact style file

    :param str filename: output file to write to
    :param stationObj station: station object to write to tfile
    :param str index: ETCCDI index being written
    :param bool doMonthly: if a monthly index
    """
    
    # open and write header
    with open(filename, "w") as outfile:

        outfile.write('"Description: ","{}"\n'.format(INDICES[index].definition))
        outfile.write('"Station: ","{}"\n'.format(station.id))
        outfile.write('"Latitude: ","{}"\n'.format('{:10.2f}"\n'.format(station.latitude).strip()))
        outfile.write('"Longitude: ","{}"\n'.format('{:10.3f}"\n'.format(station.longitude).strip()))
        outfile.write('"ClimPACT2_version: ","0.0.0\n')
        outfile.write('"Date_of_calculation: ","{}"\n'.format(dt.date.strftime(dt.date.today(), "%Y-%m-%d")))
        outfile.write("time, {}, normalised (all years)\n".format(index))

        # if monthly index
        if doMonthly:
            for m, mnth in enumerate(station.months):
                outfile.write("{:04.0f}-{:02.0f}, {}, {}\n".format(station.myears[m], station.months[m], station.monthly[m], 0))
        else:
            for y, yr in enumerate(station.years):
                outfile.write("{:04.0f}, {}, {}\n".format(station.years[y], station.annual[y], 0))

    return # write_station_index

#*********************************************
def check_directories(dataset, station, anomalies=False):

    # prepare to write the data - check all directories present
    if not os.path.exists(os.path.join(dataset.location, "formatted")):
        try:
            os.mkdir(os.path.join(dataset.location, "formatted"))                               
        except OSError:
            # likely that another instance has made the directory, so just skip
            pass

    # if standard run, place in formatted/indices/stationID
    if not anomalies:
        if not os.path.exists(os.path.join(dataset.location, "formatted", "indices")):
            try:
                os.mkdir(os.path.join(dataset.location, "formatted", "indices"))                               
            except OSError:
                # likely that another instance has made the directory, so just skip
                pass
        if not os.path.exists(os.path.join(dataset.location, "formatted", "indices", station.id)):
            try:
                os.mkdir(os.path.join(dataset.location, "formatted", "indices", station.id))
            except OSError:
                # likely that another instance has made the directory, so just skip
                pass

    # if anomaly/climatology run, place in formatted/anomalies/stationID and formatted/climatology/stationID
    elif anomalies:
        for subdir in ["anomalies", "climatology"]:
            if not os.path.exists(os.path.join(dataset.location, "formatted", subdir)):
                try:
                    os.mkdir(os.path.join(dataset.location, "formatted", subdir))                               
                except OSError:
                    # likely that another instance has made the directory, so just skip
                    pass
            if not os.path.exists(os.path.join(dataset.location, "formatted", subdir, station.id)):
                try:
                    os.mkdir(os.path.join(dataset.location, "formatted", subdir, station.id))
                except OSError:
                    # likely that another instance has made the directory, so just skip
                    pass

    return # check_directories

# #*********************************************
# def netcdf_write(filename, data, month=0):
#     """
#     Write the final netcdf file

#     :param str filename: output filename (full path)
#     :param array data: data to write
#     :param int month: which month to write - for netcdf variable name
#     """

#     outfile = ncdf.Dataset(filename, 'w', format='NETCDF4')

#     nclon = outfile.createDimension('longitude', len(box_centre_lons))
#     nclat = outfile.createDimension('latitude', len(box_centre_lats))
#     if data.shape[0] == 1: # climatology array
#         nctime = outfile.createDimension('time', 1)
#     else:
#         nctime = outfile.createDimension('time', len(REFERENCEYEARS))

#     ncbnd = outfile.createDimension('ncbnd', 2)

#     lonvar = outfile.createVariable('longitude', 'f4', ('longitude',))
#     latvar = outfile.createVariable('latitude', 'f4', ('latitude',))
#     timevar = outfile.createVariable('time', 'i4', ('time',))

#     lonBvar = outfile.createVariable('longitude_bounds', 'f4', ('longitude', 'ncbnd'))
#     latBvar = outfile.createVariable('latitude_bounds', 'f4', ('latitude', 'ncbnd'))

#     # do these need to be rolled?

#     lonvar[:] = box_centre_lons
#     latvar[:] = box_centre_lats
#     if data.shape[0] == 1: # climatology array
#         timevar[:] = 1960
#     else:
#         timevar[:] = REFERENCEYEARS

#     lonBvar[:] = np.array(list(zip(box_edge_lons[:-1], box_edge_lons[1:])))
#     latBvar[:] = np.array(list(zip(box_edge_lats[:-1], box_edge_lats[1:])))

#     if month == 0:
#         # then annual
#         datavar = outfile.createVariable("Ann", 'f4', ('time', 'latitude', 'longitude'))
#     else:
#         # then monthly = select correct name
#         name = calendar.month_abbr[month]
#         datavar = outfile.createVariable(name, 'f4', ('time', 'latitude', 'longitude'))
    
#     datavar.missing_value = HADEX_MDI
#     datavar[:] = data[:, 0, :, :]   # dummy array dimension
    
#     outfile.author = "Robert Dunn, robert.dunn@metoffice.gov.uk"
#     outfile.history = "Created by " + os.path.basename(__file__)
# #    outfile.title = index
#     outfile.file_created = dt.datetime.strftime(dt.datetime.now(), "%a %b %d, %H:%M %Y")
#     outfile.close()

#     return # netcdf_write

#*********************************************
def read_qc_flags(dataset, index, timescale):
    """
    Read the QC flag information - must be stored per index
    
    :param dataset object: dataset object
    :param str timescale: "MON" or "ANN" or ""
    :param str index: index to read
    """

    flags = {}
    try:
        with open(os.path.join(dataset.location, "{}.qc.{}.{}.txt".format(dataset.name, index, timescale)), "r", encoding="latin-1") as flag_file:
            
            for line in flag_file:
                
                # allow for zero length lines if necessary
                if len(line.strip()) > 0:
                    stn, flag = line.split()
                    flags[stn] = flag

    except IOError:
        flags = -1

    return flags # read_qc_flags

#*********************************************
def read_inventory(dataset, final=False, timescale="", index="", subdir="", qc_flags="", anomalies="None"):
    '''
    Read the Climpact metadata file as the station inventory

    :param dataset object: dataset object
    :param bool final: to read the reduced dataset
    :param str subdir: subdirectory of station data

    if not final, then the following also have meaning
    :param str timescale: "MON" or "ANN" or ""
    :param str index: index to read
    :param str qc_flags: QC flags to process (use non "" string as boolean)
    :param str anomalies: run code on anomalies or climatology rather than raw data
    '''


    if final:
        assert timescale != "" 
        assert index != ""
        
        if anomalies == "None":
            indata = np.genfromtxt(os.path.join(dataset.location, "{}.metadata.{}.{}.txt".format(dataset.name, index, timescale)), skip_header=1, dtype=(str), encoding="latin-1")
        else:
            indata = np.genfromtxt(os.path.join(dataset.location, "{}.metadata.{}.{}.{}.txt".format(dataset.name, index, timescale, anomalies)), skip_header=1, dtype=(str), encoding="latin-1")
    else:
        indata = np.genfromtxt(os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name)), skip_header=1, dtype=(str), encoding="latin-1")

    stn_qc_flags = read_qc_flags(dataset, index, timescale)

    stations = []

    # if no lines, then just skip
    if len(indata) != 0:
        
        # if only a single line, need to increase array dimension
        if len(indata.shape) == 1:
            indata = np.expand_dims(indata, axis=0)

        for l, line in enumerate(indata):

            this_station = Station(line[0].split(".")[0], float(line[1]), float(line[2]), os.path.join(dataset.location, subdir), dataset.name)
            this_station.qc = "-"
            # use dictionary to assign values
            if qc_flags != "" and stn_qc_flags != -1:
                # expecting to find some flags, and files exist, then use values
                this_station.qc = stn_qc_flags[this_station.id]


            stations += [this_station]

    return np.array(stations) # read_inventory

#*********************************************
def read_station_index(station, index, timescale, mask_early=True):
    """
    Read the CLIMPACT format index data for a specific station

    :param station station: station to process
    :param str index: index to process
    :param str timescale: timescale to process (MON/ANN)
    :param bool mask_early: remove data before the STARTYEAR (1901)
    """

    def metadata_check(station, index, infile):
        """
        Read the CLIMPACT format index data header for a specific station
        
        :param station station: station to process
        :param str index: index to process
        :param file infilee: file handle to work through
        """

        # Cross check of header to index and station metadata to ensure reading correct info!

        for ll, line in enumerate(infile):
            line = line.split(",")

            if ll == 1:
                # test station ID
                if line[1].strip().replace('"', '') != station.id:
                    return False
                    
            elif ll == 2:
                # test station latitude
                if not nearly_equal(float(line[1].strip().replace('"', '')), station.latitude, sig_fig=1):
                    return False
                    
            elif ll == 3:
                # test station longitude
                if not nearly_equal(float(line[1].strip().replace('"', '')), station.longitude, sig_fig=2):
                    return False
                    
            elif ll == 6:
                # test index abbreviation
                if "spi" in index.lower():
                    if line[1].strip().replace('"', '').lower() != "spi":
                        return False
                elif "spei" in index.lower():
                    if line[1].strip().replace('"', '').lower() != "spei":
                        return False
                elif "heatwave" in index.lower():
                    # can't test yet
                    pass
                elif line[1].strip().replace('"', '').lower() != index.lower():
                    return False

            elif ll > 7:
                # and escape
                break

        return True # metadata_check

    # read file, check it's correct, then process and return
    try:
        with open(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, index.lower(), timescale)), "r", encoding="latin-1") as infile:

            # take results from cross-check
            is_good = metadata_check(station, index, infile)

    except IOError:
        #return [], []
        # this index for this station doesn't exist.
        # shouldn't be needed
        input("Unexpected missing file: {}".format(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, index.lower(), timescale))))

            
    # need second open to ensure that the file has re-wound.
    with open(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, index.lower(), timescale)), "r", encoding="latin-1") as infile:

        if not is_good:
            # return empty fields
            return np.array([]), np.array([])

        # else we can go further
        indata = np.genfromtxt(infile, skip_header=7, delimiter=",", dtype=(str), missing_values="NaN", filling_values=HADEX_MDI, encoding="latin-1")

        if "heatwave" not in index.lower():
            if timescale == "MON":
                data = np.array([float(i) for i in indata[:, 1]])
                data = data.reshape(-1, 12)

                times = np.array([int(i.split("-")[0]) for i in indata[:, 0]])
                times = times.reshape(-1, 12)[:, 0] # just get the Januaries

            else:
                 if len(indata.shape) == 1:
                    data = np.array([indata[1]]).astype(float)
                    times = np.array([indata[0]]).astype(int)
                 else:
                    data = indata[:, 1].astype(float)
                    times = np.array([int(i) for i in indata[:, 0]])
        elif "heatwave" in index.lower():
            # then there are 5 columns which need dealing with separately
            #     and also contain NA as missing values
            if len(indata.shape) == 1:
                data = np.array([indata[1]])
                if "NA" in data[0]:
                    data[0] = HADEX_MDI
                else:
                    data = data.astype(float)
                times = np.array([indata[0]]).astype(int)
            else:
                data = indata[:, 1]
                data = np.array([d.strip() for d in data]) # in case of spaces as still strings
                bad_locs, = np.where(data == "NA")
                data[bad_locs] = HADEX_MDI
                data = data.astype(float)
                times = np.array([int(i) for i in indata[:, 0]])

            # TODO - handle sub-indices successfully, somehow! placeholder for code

    # close with open() as loop

    # manually remove any remaining nans and set to missing
    nan_loc = np.where(data != data)
    if len(nan_loc[0]) > 0:
        data[nan_loc] = HADEX_MDI


    # mask data
    bad_locs = np.ma.where(data == -9999) #  Thailand 61-90 issue - easiest to fix for all as this shouldn't be a value anyway
    data[bad_locs] = HADEX_MDI

    data = np.ma.masked_where(data == HADEX_MDI, data)

    data.fill_value = HADEX_MDI

    # account for data which is all present, and hence a non-array mask
    if type(data.mask) == np.bool_:
        if data.mask:
            data.mask = np.ones(data.shape)
        elif not data.mask:
            data.mask = np.zeros(data.shape)

    # some data starts before 1901
    if mask_early:
        if times[0] < STARTYEAR.year:

            # but if it also ends before 1901
            if times[-1] < STARTYEAR.year:
                # all data is pre-1901
                times = np.ma.array([])
                data = np.ma.array([])

            else:
                # data ends after start of dataset, so cut off the front
                loc, = np.where(times == STARTYEAR.year)

                times = times[loc[0]:]
                data = data[loc[0]:]

        # else all data is post 1900

    return times, data # read_station_index


#*******************************************
def map_2_points(point0, point1, angle_int=False):
    '''
    Copied from FORTRAN code in utility_2.5.2.f90
    Adapted to work on all array at once.

    :param array point0: input point (lat,lon) - gridbox centre, repeated to match length of station
    :param array point1: input point (lat,lon) - station
    :returns: separation, angle - separation in metres and angle in radians (-pi to pi)
    '''

    if len(point0) != len(point1):
        print("lengths do not match")
        raise ValueError

    lat0, lon0 = point0[:, 0], point0[:, 1]
    lat1, lon1 = point1[:, 0], point1[:, 1]

    R_earth = 6378206.4

    coslt1 = np.cos(np.radians(lat1))
    sinlt1 = np.sin(np.radians(lat1))
    coslt0 = np.cos(np.radians(lat0))
    sinlt0 = np.sin(np.radians(lat0))
    
    cosl0l1 = np.cos(np.radians(lon1 - lon0))
    sinl0l1 = np.sin(np.radians(lon1 - lon0))

    cosc = sinlt0 * sinlt1 + coslt0 * coslt1 * cosl0l1

    gt1 = np.where(cosc > 1.)
    if len(gt1) > 0: cosc[gt1] = 1.
    lt1 = np.where(cosc < -1.)
    if len(lt1) > 0: cosc[lt1] = -1.

    #if cosc > 1. : cosc=1.
    #if cosc < -1.: cosc=-1.

    sinc = np.sqrt(1.0 - cosc**2)
    

    cosaz = np.ones(len(point0))
    sinaz = np.zeros(len(point0))

    nonzero = np.where(np.abs(sinc) > 1.e-7)

    cosaz[nonzero] = ((coslt0[nonzero] * sinlt1[nonzero]) - \
                        (sinlt0[nonzero] * coslt1[nonzero] * cosl0l1[nonzero])) / sinc[nonzero]
    sinaz[nonzero] = sinl0l1[nonzero] * coslt1[nonzero] / sinc[nonzero]

    #if np.abs(sinc) > 1.e-7:
    #    cosaz=(coslt0 * sinlt1 - sinlt0 *c oslt1 * cosl0l1) / sinc #Azimuth
    #    sinaz = sinl0l1 * coslt1 / sinc
    #else:		#Its antipodal
    #    cosaz = 1.d0
    #    sinaz = 0.d0

    separation = np.arccos(cosc) * R_earth / 1000.
    angle = np.arctan2(sinaz, cosaz)
    if angle_int:
        angle = angle.astype(int)

    return separation.astype(int), angle #  map_2_points

#*******************************************
def get_land_sea_mask(boxlats, boxlons, floor=True, write=False):
    """
    Create a land-sea mask from a 15min version

    :param array boxlats: array of box centre latitudes
    :param array boxlons: array of box centre longitudes
    :param bool floor: take any grid box if it has some land in it, else 50%
    :param bool write: output a txt file
    """

    # read in data
    not_got = True
    while not_got:
        try:
            all_data = np.genfromtxt(os.path.join(INFILELOCS, "land_mask_15min.msk"), dtype=(float), encoding="latin-1")

            # if successful read
            if all_data.shape[0] == 1036800:
                not_got = False

        except ValueError:
            # likely to be from multiple reads - wait a bit and try again
            continue

        if not_got:
            # likely to be from multiple reads - wait a bit and try again
            import time
            time.sleep(5*(1 + np.random.random()))

    # (-90 --> 90, 0 --> 360)

    # extract lat/lon
    lats = np.unique(all_data[:, 0])
    lons = np.unique(all_data[:, 1])

    # and make grid at original resolution
    land_fraction = np.zeros((len(lats), len(lons)))
    for v, val in enumerate(all_data):
        lat_loc, = np.where(lats == val[0])
        lon_loc, = np.where(lons == val[1])
        land_fraction[lat_loc, lon_loc] = val[2]


    # make suitably sized boxes from supplied coordinates
    deltalat = np.abs(np.max(np.diff(boxlats)))
    deltalon = np.abs(np.max(np.diff(boxlons)))
    
    box_edge_lons = np.arange(0 - deltalon/2., 360 + deltalon/2., deltalon)
    box_edge_lats = np.arange(-90 - (deltalat)/2., 90 + deltalat, deltalat)

    lon_bbox = np.array(list(zip(box_edge_lons[:-1], box_edge_lons[1:])))
    lat_bbox = np.array(list(zip(box_edge_lats[:-1], box_edge_lats[1:])))

    new_land_fraction = np.zeros([len(boxlats), len(boxlons)])
    
    # spin through each lat/lon of the new grid
    for lt in range(new_land_fraction.shape[0]):
        for ln in range(new_land_fraction.shape[1]):

            # lower left inclusive, upper right exclusive boundaries
            lat_locs, = np.where(np.logical_and(lats >= lat_bbox[lt, 0], lats < lat_bbox[lt, 1]))

            if lon_bbox[ln, 0] < 0:
                # have to do either side of zero separately
                lon_locs1, = np.where(lons >= lon_bbox[ln, 0] + 360)
                lon_locs2, = np.where(lons < lon_bbox[ln, 1])
                lon_locs = np.append(lon_locs1, lon_locs2)

            else:
                lon_locs, = np.where(np.logical_and(lons >= lon_bbox[ln, 0], lons < lon_bbox[ln, 1]))


            # make grid of lats/lons
            meshlons, meshlats = np.meshgrid(lon_locs, lat_locs)

            new_land_fraction[lt, ln] = np.mean(land_fraction[meshlats, meshlons])

    
    # all unset
    lsm = np.zeros([len(boxlats), len(boxlons)], dtype=(bool))

    if floor:
        lsm[new_land_fraction == 0] = 1
    else:
        lsm[new_land_fraction < 50.] = 1


    if write:
        resolution = [deltalat, deltalon]

        outsuffix = "{:3.0f}x{:3.0f}".format(180./resolution[0], 360./resolution[1])

        outfilename = os.path.join(INFILELOCS, "land_mask_{}.msk".format(outsuffix))
        if os.path.exists(outfilename): os.remove(outfilename)

        outfile = open(outfilename, "w")
        for lt in range(new_land_fraction.shape[0]):
            for ln in range(new_land_fraction.shape[1]):

                outfile.write("{:8.3f} {:8.3f} {:5.1f} {:5.1f}\n".format(box_centre_lats[lt], box_centre_lons[ln], new_land_fraction[lt, ln], 100.0))

        outfile.close()          

    return lsm # get_land_sea_mask

#***************************************
def CompletenessCheckGrid(data, endtime, starttime, threshold=0.9):
    """
    Check whether the grid box has data for 90% of the time

    :param array data: data array
    :param int endtime: last year of data
    :param int starttime: first year of data
    :param flt threshold: fraction of data needing to be present
    :returns: np masked array of mask.
    """
    nyears, nlat, nlon = data.shape

    good_boxes = np.empty([nlat, nlon])
    # set mask to True for all boxes
    good_boxes[:, :] = True

    # run through each grid box
    for lat in range(nlat):
        for lon in range(nlon):
            # get the timeseries for this gridbox

            # it's a masked array, so fill all masked bits
            #   with the fillvalue
            gridbox = data[:, lat, lon]
            
            if len(gridbox.compressed()) > threshold * (endtime - starttime + 1):
                good_boxes[lat, lon] = False
            else:
                good_boxes[lat, lon] = True

    return np.ma.make_mask(good_boxes) # CompletenessCheckGrid

#***************************************
def MedianPairwiseSlopes(xdata, ydata, mdi):
    '''
    Calculate the median of the pairwise slopes - assumes no missing values

    :param array xdata: x array
    :param array ydata: y array
    :param float mdi: missing data indicator
    :returns: float of slope
    '''
    # run through all pairs, and obtain slope
    slopes = []
    for i in range(len(xdata)):
        for j in range(i + 1, len(xdata)):
            if ydata[i] > mdi and ydata[j] > mdi:
                slopes += [(ydata[j] - ydata[i]) / (xdata[j] - xdata[i])]

    mpw = np.median(np.array(slopes))

    # copied from median_pairwise.pro methodology (Mark McCarthy)
    slopes.sort()

    good_data = np.where(ydata != mdi)[0]

    n = len(ydata[good_data])

    dof = n*(n-1)//2
    w = math.sqrt(n*(n-1)*((2.*n)+5.)/18.)

    # 95% confidence range
    rank_upper = ((dof+1.96*w)/2.)+1
    rank_lower = ((dof-1.96*w)/2.)+1

    if rank_upper >= len(slopes): rank_upper = len(slopes)-1
    if rank_upper < 0: rank_upper = 0
    if rank_lower < 0: rank_lower = 0

    upper = slopes[int(rank_upper)]
    lower = slopes[int(rank_lower)]
    

    return mpw, lower, upper      # MedianPairwiseSlopes

#*********************************************
def mean_absolute_deviation(data, median=False):    
    ''' Calculate the MAD of the data '''
    
    if median:
        mad = np.ma.mean(np.ma.abs(data - np.ma.median(data)))
        
    else:        
        mad = np.ma.mean(np.ma.abs(data - np.ma.mean(data)))

    return mad # mean_absolute_deviation

#*********************************************
def dms2dd(degrees, minutes, seconds, direction=""):
    """
    D:M:S to decimal degrees

    Direction of N/E gives positive, S/W gives negative
    """

    if float(degrees) < 0:
        # direction keyword should enforce correct direction later if necessary
        dd = float(degrees) - float(minutes)/60. - float(seconds)/(60.*60.)
    else:
        dd = float(degrees) + float(minutes)/60. + float(seconds)/(60.*60.)

    # only apply if necessary - some come with negatives set
    if direction != "":
        if direction == 'W' or direction == 'S':
            if dd > 0:
                dd *= -1
        if direction == "E" or direction == "N":
            if dd < 0:
                dd *= -1

    return dd # dms2dd

#*********************************************
def dd2dms(deg):
    """
    Decimal degrees to D:M:S]
    """
    sign = int(np.sign(deg))
    d = int(deg)
    md = abs(deg - d) * 60
    m = int(md)
    sd = (md - m) * 60
    return sign, [d, m, sd]

#*********************************************
def append_timeseries(ts, times, data):
    """
    Appends timeseries in the station object - time and data values
    """
    
    ts.times = np.append(ts.times, times)
    ts.data = np.append(ts.data, data)
    
    return ts # append_timeseries

#*********************************************
def check_completeness_indata(st_var):
    """
    Go through each year and month to find number of days present in each
    And hence years/months where valid indices could be calculated

    :param obj st_var: station variable
    """

    # go through years and months
    present_days_yr = []
    good_months = []

    # extract each year
    for year in range(STARTYEAR.year, ENDYEAR.year):
        this_year, = np.where(np.logical_and(st_var.times >= dt.date(year, 1, 1), st_var.times < dt.date(year + 1, 1, 1)))

        present_days = []

        # extract each month
        for month in range(1, 13):

            if month == 12:
                this_month, = np.where(np.logical_and(st_var.times >= dt.date(year, month, 1), st_var.times < dt.date(year+1, 1, 1)))
            else:
                this_month, = np.where(np.logical_and(st_var.times >= dt.date(year, month, 1), st_var.times < dt.date(year, month+1, 1)))

            # days of data in each month
            present_days += [len(this_month)]

        # number of months which have sufficient data for indices to be calculated
        # just assume all months are short for this first cut
        good_months += [len(np.where(np.array(present_days) > 28 - MAX_MISSING_DAYS_PER_MONTH)[0])]

        # days of data in each year
        present_days_yr += [len(this_year)]

    # find where year constraints would be met
    good_yrs_days, = np.where(np.array(present_days_yr) > 365 - MAX_MISSING_DAYS_PER_YEAR)
    good_yrs_days_months, = np.where(np.logical_and(np.array(present_days_yr) > 365 - MAX_MISSING_DAYS_PER_YEAR, np.array(good_months) > MAX_MISSING_MONTHS_PER_YEAR))

    return good_yrs_days, good_yrs_days_months # check_completeness

#*********************************************
def parse_PubA(WMO_id, dataset):
    """In case stations are given without IDs"""

    def fix_coord(coord):
        """Split the coords and convert to decimal """

        d, m, s = coord.split()
        s, dirn = s[:-1], s[-1]
        
        return  dms2dd(int(d), int(m), int(s), dirn) # fix_coord    

    # tab delimited pubA
    infilename = os.path.join(INFILELOCS, config.get("Files", "pubA"))

    stnIDs, stnNames, stnLats, stnLons = [], [], [], []
    with open(infilename, "r", encoding="latin-1") as infile:
        # read header and discard
        infile.readline()

        # spin through file and extract info
        for line in infile:
            sline = line.split("\t")
            
            stnIDs += [sline[5]]
            stnNames += [sline[7]]
            stnLats += [sline[8]]
            stnLons += [sline[9]]

    # if there is a match
    if WMO_id in stnIDs:

        loc, = np.where(np.array(stnIDs) == WMO_id)

        name = stnNames[loc[0]]

        lat = fix_coord(stnLats[loc[0]])
        lon = fix_coord(stnLons[loc[0]])

        return Station(WMO_id, lat, lon, os.path.join(dataset.location, "raw"), dataset.name) # parse_PubA
    
    else:
        return -1 # parse_PubA

#*********************************************
def select_qc_passes(stations, qc_flags=""):
    """
    Subset an array of Station Objects to only retain those with no QC flags of specified types

    :param array stations: array of station objects
    :param str qc_flags: QC flags to process: W, B, A, N, C, R, F, E; can be given as single string, e.g WNR

    :returns: array of stations
    """

    out_stations = []

    # longhand spin through all entries
    for stn in stations:

        if stn.qc == "-":
            # no flags set
            out_stations += [stn]

        else:
            # spin through the flags of interest
            for flag in qc_flags:
                if flag in stn.qc:
                    # flag found in this station
                    # single failure somewhere, so may as well stop looking
                    break
            else:
                # no break in loop
                out_stations += [stn]
    return np.array(out_stations) # select_qc_passes

#*********************************************
def read_blacklist(diagnostics = False):
    """
    Read in a file of the blacklisted stations
    """
    blacklist = []

    try:
        with open("{}/blacklist.dat".format(os.path.expanduser(ANCILS_LOC)), "r") as infile:
            for line in infile:
                if diagnostics:
                    print(line)
                sline = line.split()
                blacklist += [Station(sline[0], float(sline[1]), float(sline[2]), "", sline[3])]

    except IOError:
        # no such file
        pass

    return np.array(blacklist)# read_blacklist

#*********************************************
def make_filenames(index="TX90p", grid="ADW", anomalies="None", extra="", month_index=""):
    """
    Build appropriate filenames for all the information

    :param str index: the index
    :param str grid: the grid
    :param str anomalies: the anomalies
    :param str extra: e.g. num, or numdls

    :returns: string of filename
    """
    attribs = ncdfp.read_global_attributes(os.path.join(INFILELOCS, "attributes.dat"))
    dataset_name = attribs["title"]

    baseP_string = "{}-{}".format(str(REF_START)[-2:], str(REF_END)[-2:])
    resol_string = "{}x{}deg".format(DELTALAT, DELTALON)

    if month_index == "":
        outfilename = "{}_{}_{}-{}_{}_{}_{}".format(dataset_name, index, STARTYEAR.year, ENDYEAR.year-1, grid, baseP_string, resol_string)
    else:
        outfilename = "{}_{}_{}-{}_{}_{}_{}_{}".format(dataset_name, index, STARTYEAR.year, ENDYEAR.year-1, grid, month_index, baseP_string, resol_string)

    if extra != "":
        outfilename = "{}_{}".format(outfilename, extra)
    if anomalies != "None":
        outfilename = "{}_{}".format(outfilename, anomalies)
    
    return "{}.nc".format(outfilename) # make_filenames


#*********************************************
def parse(string, fields):
    """
    Replacement for struct and parser in Py 3
    """

    cumulative = 0
    outlist = []

    for field in fields:
        outlist += [string[cumulative: cumulative + field]]
        cumulative += field

    return np.array(outlist) # parse

#*********************************************
def match_reference_period(BP):
    """
    Match the data Base Period with the selected Reference Period
    """
    def make_year(string):
        """convert YY into YYYY and integer"""
        if int(string) < 50:
            value = 2000 + int(string)
        else:
            value = 1900 + int(string)
        return value # make_year
    #***************************

    start, end = BP.split("-")

    match = False
    if (start == end) and (start == "00"):
        match = True

    else:
        start = make_year(start)
        end = make_year(end)

        if start == REF_START and end == REF_END:
            match = True
    
    return match # match_reference_period

#************************************************************************
def nearly_equal(a, b, sig_fig=2):
    """
    Returns it two numbers are nearly equal within sig_fig decimal places

http://stackoverflow.com/questions/558216/function-to-determine-if-two-numbers-are-nearly-equal-when-rounded-to-n-signific

    :param flt a: number 1
    :param flt b: number 2
    :param int sig_fig: number of decimal places to check agreement to

    :returns bool:    
    """

    return (a == b or 
             round(a*10**sig_fig) == round(b*10**sig_fig)
           ) # nearly_equal

#******************************************************************************************
#******************************************************************************************
CONFIG_FILE = "ancils/configuration.txt"

if not os.path.exists(os.path.expanduser(CONFIG_FILE)):
    print("Configuration file missing - {}".format(os.path.expanduser(CONFIG_FILE)))
    sys.exit

# read in configuration file
config = configparser.ConfigParser()
config.read(os.path.expanduser(CONFIG_FILE))

#************************************************************************
# set paths
ROOT = config.get("Paths", "root")
ANCILS_LOC = config.get("Paths", "ancillaries")
INPUT_DATA_LOC = config.get("Paths", "input_location")

# make all the other paths, and files
if not os.path.exists(ROOT): os.mkdir(ROOT)

# tag indata and files with reference period
ref_start = config.getint("Climpact", "reference_start")
ref_end = config.getint("Climpact", "reference_end")
deltalon = config.getfloat("Grid", "deltalon")
deltalat = config.getfloat("Grid", "deltalat")
baseP_string = "{}{}".format(str(ref_start)[-2:], str(ref_end)[-2:])
resol_string = "{}x{}".format(deltalat, deltalon)

# set up the paths
INROOT = "{}/indata_{}/".format(ROOT, baseP_string)
if not os.path.exists(INROOT):
    os.mkdir(INROOT)

OUTROOT = "{}/outdata/".format(ROOT)
if not os.path.exists(OUTROOT):
    os.mkdir(OUTROOT)

FINALROOT = "{}/finaldata/".format(ROOT)
if not os.path.exists(FINALROOT):
    os.mkdir(FINALROOT)

DLSROOT = "{}/dls/".format(ROOT)
if not os.path.exists(DLSROOT):
    os.mkdir(DLSROOT)
DLSLOCS = "{}/{}_{}".format(DLSROOT, baseP_string, resol_string)
if not os.path.exists(DLSLOCS):
    os.mkdir(DLSLOCS)

PLOTLOCS = "{}/plots/".format(ROOT)
if not os.path.exists(PLOTLOCS):
    os.mkdir(PLOTLOCS)
if not os.path.exists(os.path.join(PLOTLOCS, "DLS")):
    os.mkdir(os.path.join(PLOTLOCS, "DLS"))

INFILELOCS = "{}/files_{}/".format(ROOT, baseP_string)
if not os.path.exists(INFILELOCS):
    os.mkdir(INFILELOCS)

# copy over files
copy(ANCILS_LOC, INFILELOCS)

# set relative to current as checked out as part of the repository
CLIMPACT_LOCS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "climpact2")

# extract the climpact 2 code if not already done so - deprecated as now climpact2-master under svn version control.
#if not os.path.exists(CLIMPACT_LOCS):
#    import zipfile
#    # unzip the file
#    zip_ref = zipfile.ZipFile(os.path.join(os.path.dirname(os.path.abspath(__file__)), "climpact2-master.zip"), 'r')
#    zip_ref.extractall(os.path.dirname(os.path.abspath(__file__)))
#    zip_ref.close()
#    print("{} unzipped".format(os.path.join(os.path.dirname(os.path.abspath(__file__)), "climpact2-master.zip")))

#************************************************************************
PARAM_FILE = config.get("Files", "parameters")
with open(os.path.join(os.path.expanduser(ANCILS_LOC), PARAM_FILE), "r", encoding="latin-1") as pf:
    parameters = json.load(pf)

#************************************************************************
# default values for dataset
HADEX_MDI = config.getfloat("Values", "mdi")
STARTYEAR = dt.datetime(config.getint("Values", "dataset_start_year"), config.getint("Values", "dataset_start_month"), 1)
ENDYEAR = dt.datetime(config.getint("Values", "dataset_end_year"), config.getint("Values", "dataset_end_month"), 1)
REFERENCEYEARS = np.arange(int(STARTYEAR.year), int(ENDYEAR.year))

#************************************************************************
# grid box defaults
DELTALON = config.getfloat("Grid", "deltalon")
DELTALAT = config.getfloat("Grid", "deltalat")
DOLSM = config.getboolean("Grid", "doLSM")

# to match HadEX2
# # set up box edges
# box_edge_lons = np.arange(0 - DELTALON/2., 360 + DELTALON/2., DELTALON)
# box_edge_lats = np.arange(-90 - (DELTALAT)/2., 90 + DELTALAT, DELTALAT)
# box_centre_lons = box_edge_lons[1:] - DELTALON/2.
# box_centre_lats = box_edge_lats[1:] - DELTALAT/2.

# # ensure don't be unphysical
# box_edge_lats[0] = -90
# box_edge_lats[-1] = 90

# # and cope with final boxes being 1/2 height
# box_centre_lats[0] = box_centre_lats[0] + DELTALAT/2.
# box_centre_lats[-1] = box_centre_lats[-1] - DELTALAT/2.

# just a simple grid with no half-width boxes
box_edge_lons = np.arange(0, 360 + DELTALON, DELTALON)
box_edge_lats = np.arange(-90, 90 + DELTALAT, DELTALAT)
box_centre_lons = box_edge_lons[1:] - DELTALON/2.
box_centre_lats = box_edge_lats[1:] - DELTALAT/2.


#************************************************************************
# input data requirements
FINISH_AFTER = dt.date(config.getint("Data", "finish_after"), 1, 1) # dt.date(1981,1,1)
MODERN_PERIOD = dt.date(config.getint("Data", "modern_period"), 1, 1) # dt.date(1951,1,1)
MIN_N_YEARS = config.getint("Data", "minimum_number_of_years")
MIN_LENGTH_OF_RECORD = dt.timedelta(days=(365*MIN_N_YEARS))
MAX_MISSING_DAYS_PER_YEAR = config.getint("Data", "maximum_missing_days_per_year")
MAX_MISSING_DAYS_PER_MONTH = config.getint("Data", "maximum_missing_days_per_month")
MAX_MISSING_MONTHS_PER_YEAR = config.getint("Data", "maximum_missing_months_per_year")
MAX_MISSING_YEARS = config.getint("Data", "maximum_missing_years")

MIN_YEARS_OF_INDEX_DATA = MIN_N_YEARS

#************************************************************************
# dls settings
LAT_BANDS = parameters["dls_bands"]
DEFAULT_DLS = config.getint("DLS", "default")
MAX_DLS = config.getint("DLS", "maximum_dls")
FIX_ZERO = config.getboolean("DLS", "fix_zero")

#************************************************************************
# Climpact defaults
REMOVE_EXTRA = config.getboolean("Climpact", "remove_extra_output")
REF_START = config.getint("Climpact", "reference_start")
REF_END = config.getint("Climpact", "reference_end")
NCORES = config.getint("Climpact", "ncores")
# standards for HadEX3 indices
WSDI_N = config.getint("Climpact", "wsdi_n")
CSDI_N = config.getint("Climpact", "csdi_n")
TB_HDD = config.getint("Climpact", "tb_hdd")
TB_CDD = config.getint("Climpact", "tb_cdd")
TB_GDD = config.getint("Climpact", "tb_gdd")
RX_UD = config.getint("Climpact", "rx_ud")
RNNMM_UD = config.getint("Climpact", "rnnmm_ud")
TXTN_UD = config.getint("Climpact", "txtn_ud")
SPEI = config.getint("Climpact", "spei")

#************************************************************************
# CAM gridding defaults
CAM_COMPLETENESS = config.getint("CAMgrid", "cam_completeness") # for anomalies in CAM (different to anomalies in the workflow)
CLIM_START = dt.datetime(config.getint("CAMgrid", "cam_clim_start_year"), config.getint("CAMgrid", "cam_clim_start_month"), 1)
CLIM_END = dt.datetime(config.getint("CAMgrid", "cam_clim_end_year"), config.getint("CAMgrid", "cam_clim_end_month"), 1)
STATIONS_IN_BOX = config.getint("CAMgrid", "stations_in_box")

#************************************************************************
# ADW gridding defaults
STATIONS_IN_DLS = config.getint("ADWgrid", "stations_in_dls")
M = config.getint("ADWgrid", "m")

#************************************************************************
# Climatology and Anomalies
CLIM_COMPLETENESS = config.getfloat("Climatologies", "climatology_completeness") #  for calculating anomalies 

#************************************************************************
# QC settings
QC_PERCENTILES_COMPLETENESS = config.getfloat("QC", "percentile_completeness")

#************************************************************************
# GHCND subset settings
ghcnd_peterson = config.getboolean("GHCND", "peterson")
ghcnd_hcn = config.getboolean("GHCND", "hcn")
ghcnd_gsn = config.getboolean("GHCND", "gsn")

#************************************************************************
# Plot settings
TREND_START = config.getint("PLOTS", "trend_start")
TREND_END = config.getint("PLOTS", "trend_end")
TREND_COMPLETENESS = config.getfloat("PLOTS", "trend_completeness")
TREND_FINAL_YEAR = config.getint("PLOTS", "trend_final_year")
WATERMARK = config.getboolean("PLOTS", "watermark")
FONTSIZE = config.getint("PLOTS", "fontsize")
HADEX_LOC = config.get("PLOTS", "hadex_loc")
HADEX2_LOC = config.get("PLOTS", "hadex2_loc")
GHCNDEX_LOC = config.get("PLOTS", "ghcndex_loc")
ERA5_LOC = config.get("PLOTS", "era5_loc")

#************************************************************************
# process all indices and get their metadata

ALL_INDICES = parameters["indices"]

#************************************************************************
# extract index information
INDICES = {}
with open(os.path.join(CLIMPACT_LOCS, "ancillary/climate.indices.csv"), "r", encoding="latin-1") as infile:
    
    old_name = ""
    for lc, line in enumerate(infile):

        if lc == 0: continue
        
        name, units, lname, definition, team, annual, baseperiod = line.split("\t")

        # sort name given strange capitalization
        for case_name in ALL_INDICES:
            if case_name.lower() == name:
                name = case_name
                break

        # ETSCI indices need a little more handling
        if name.lower() in ["spei", "spi"]:
            # these names cover 4 periods each which are written out.  All are monthly (without annual counterpart)
            for period in [3, 6, 12, SPEI]:
                new_name = "{}month_{}".format(period, name)
                INDICES[new_name] = Index(new_name, lname, units, definition, team, True, baseperiod[:-1])

        elif name.lower() in ["csdid", "wsdid", "hddheatn", "cddcoldn", "gddgrown"]:
            # in all of these the names inherit from the settings
            if name == "CSDId":
                name = "CSDI{}".format(CSDI_N)
            elif name == "WSDId":
                name = "WSDI{}".format(WSDI_N)
            elif name == "HDDheatn":
                name = "HDDheat{}".format(TB_HDD)
            elif name == "CDDcoldn":
                name = "CDDcold{}".format(TB_CDD)
            elif name == "GDDgrown":
                name = "GDDgrow{}".format(TB_GDD)

            INDICES[name] = Index(name, lname, units, definition, team, False, baseperiod[:-1])
        
        elif name.lower() in ["txdtnd", "txbdtnbd", "rxdday"]:
            # in all of these the "d"s go missing on writing out
            if name == "TXdTNd":
                name = "TXTN"
            elif name == "TXbdTNbd":
                name = "TXbTNb"
            elif name == "RXdday":
                name = "RXday"

            INDICES[name] = Index(name, lname, units, definition, team, False, baseperiod[:-1])
            
        elif name.lower() == "hw":
            # the heatwave indices are more complicated.
            # output file formats are trickier too - see read_station_index()
            for new_name in ["ECF_heatwave", "EHF_heatwave", "TX90_heatwave", "TN90_heatwave"]:
                INDICES[new_name] = Index(new_name, lname, units, definition, team, False, baseperiod[:-1])

        else:
            INDICES[name] = Index(name, lname, units, definition, team, False, baseperiod[:-1])

        
        # if second entry - therefore monthly version
        if old_name == name:
            # have read in annual version - just now indicate it's monthly too
            try:
                index = INDICES[name]

                assert index.long_name == lname
                assert annual == "F" # get the monthly index
                
                index.monthly = True
                # fix definition
                def_list = index.definition.split(" ")[1:]
                def_list[0] = def_list[0].capitalize()
                index.definition = " ".join(def_list)
                INDICES[name] = index
                
            except IndexError:
                input("index process")

        old_name = name

# fix for missing ETR
INDICES["ETR"] = Index("ETR", "Extreme Temperature Range", "degrees_C", "Mean annual difference between annual TX and annual TN", "ETCCDI", True, baseperiod[:-1])
                
MONTHLY_INDICES = []

for name, index in INDICES.items():
    try:
        if index.monthly:
            MONTHLY_INDICES += [index.name]
    except KeyError:
        pass

# make all plot locations
for name, index in INDICES.items():

    if not os.path.exists(os.path.join(PLOTLOCS, index.name)):
        os.mkdir(os.path.join(PLOTLOCS, index.name))

#------------------------------------------------------------
# END
#------------------------------------------------------------
