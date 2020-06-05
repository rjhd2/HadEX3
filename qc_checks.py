#!/usr/local/sci/bin/python
#------------------------------------------------------------
#                    SVN Info
#$Rev::                                         $:  Revision of last commit
#$Author::                                      $:  Author of last commit
#$Date::                                        $:  Date of last commit
#------------------------------------------------------------
#  
#  Run QC tests on index files
#
#  Uses work done by Tanya Lippmann (UNSW) in developing these checks
#------------------------------------------------------------
"""
Runs QC checks on the index data.

Set flags for whole station for
W - World record values 
B - Base period average [not used down stream]
A - Annual in monthly 
N - Negative Values
C - Consistency between indices
R - Correlation
F - File exists
E - Empty file
V - Coverage over base period
M - Metadata

qc_checks.py invoked by typing::

  python qc_checks.py --index "TX90p" --diagnostics

Input arguments:

--index         Which ETCCDI index to use

--diagnostics   Output extra info (default False)
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')

import utils
import dls_cal as dls

TOLERANCE = 1.e-10 # for floating point errors on tests

# WORLD RECORD VALUES
WORLD_RECORDS_VALUES = {'Tx':56.7, 'Tn':-89.2, 'Rx1day':1825, 'PRCPTOT':26470}
WORLD_RECORDS_INDICES = ['TXx', 'TXn', 'TNx', 'TNn', 'Rx1day', 'PRCPTOT']
# http://wmo.asu.edu/#continental
 
# THESE INDICES SHOULD NOT HAVE NEGATIVE VALUES:
NO_NEGATIVES = set(['TX10p', 'TX90p', 'TN10p', 'TN90p', "R10mm", "R20mm", "SDII", "SU", "TR", "CSDI", "WSDI", "DTR", "ETR", "Rx1day", "Rx5day", "FD", "ID", "GSL", "Rnnmm", "CDD", "CWD", "PRCPTOT", "R99p", "R95p", "R99pTOT", "R95pTOT", "TNlt2", "TNltm2", "TNltm20", "WSDI3", "CSDI3", "TMge5", "TMlt5", "TMlt10", "TMgt10", "TXge30", "TXge35", "TXTN", "RXday", "TXbTNb"])

# THESE INDICES SHOULD AVERAGE 10% over the base period 1961-1990. EXCEPT TN50p & TX50p which would average 50%:
PERCENTILE_INDICES = set(['TX10p', 'TX90p', 'TN10p', 'TN90p', 'TXgt50p', 'R95p', 'R99p', 'R95pTOT', 'R99pTOT', "TXbTNb", "TXTN"])
PERCENTILE50_INDICES = set(['TXgt50p'])

# CONSISTENCY - note fix for second TXx call
GT_INDEX = set(["TXx", "TNx", "TXx2", "TXn", "R10mm", "R95p", "PRCPTOT", "Rx5day", "TXm", "TMm"])
LT_INDEX = {"TXx":"TXn", "TNx":"TNn", "TXx2":"TNx", "TXn":"TNn", "R10mm":"R20mm", "R95p":"R99p", "PRCPTOT":"R95p", "Rx5day":"Rx1day", "TXm" : "TMm", "TMm" : "TNm"}

# ANNUAL VALUE SHOULD BE ONE OF THE MONTHLY ONES
ANNUAL_IN_MONTHLY = set(["TNx", "TNn", "TXx", "TXn", "Rx1day", "Rx5day"]) 

# TEST IF STATIONS VERY CLOSE - DENSE NETWORKS OF RAIN GAUGES CAUSE ISSUES LATER ON
NEARBY = set(["Rx1day", "Rx5day", "R10mm", "R20mm", "SDII", "CDD", "CWD", "PRCPTOT", "R99p", "R95p", "R99pTOT", "R95pTOT"])
NEARBY_SEP = 5.0

# CORRELATION THRESHOLDS
NOT_CORRELATION = set(["FD", "ID", "SU", "TR", "R10mm", "R20mm", "TXge35", "TXge30", "TNlt2", "TNltm2", "TNltm20" ])
CORR = 0.99
SEP = 250.0

#*******************************************
def touch(fname, times=None):
    """ Emulates Linux touch utility """
    # https://stackoverflow.com/questions/1158076/implement-touch-using-python

    with open(fname, 'a'):
        os.utime(fname, times)

    return # touch

#*******************************************
def set_qc_attr(station, flag):
    """
    Sets the qc attribute - appending if already flag applied or overwriting no-flag marker 
    """

    if station.qc == "-":
        setattr(station, "qc", "{}".format(flag))
    else:
        # don't add many flag if test done on multiple months
        if flag not in station.qc:
            setattr(station, "qc", "{}{}".format(station.qc, flag))

    return # set_qc_attr

#*******************************************
def setup_outfile(index, all_datasets):
    """
    Clean up before new run - so that flags are current
    """
    # if the file exists from a previous run, then delete it
    
    # QC and Metadata info
    for name in ["QC", "Metadata"]:
        if os.path.exists(os.path.join(utils.INFILELOCS, "{}_fails_{}.txt".format(name, index))):
            os.remove(os.path.join(utils.INFILELOCS, "{}_fails_{}.txt".format(name, index)))

        touch(os.path.join(utils.INFILELOCS, "{}_fails_{}.txt".format(name, index)))

    # QC flag info
    if index in utils.MONTHLY_INDICES:
        timescales = ["MON", "ANN"]
    else:
        timescales = ["ANN"]

    for dataset in all_datasets:
        for timescale in timescales:
            if os.path.exists(os.path.join(dataset.location, "{}.qc.{}.{}.txt".format(dataset.name, index, timescale))):
                os.remove(os.path.join(dataset.location, "{}.qc.{}.{}.txt".format(dataset.name, index, timescale)))
                          
    return # setup_outfile

#*******************************************
def write_qc_flag_info(index, timescale, message, station, year):
    """
    Write the QC logfile
    """

    with open(os.path.join(utils.INFILELOCS, "QC_fails_{}.txt".format(index)), "a") as outfile:

        outfile.write("{} {} {} {} {}\n".format(timescale, station.source, station.id, message, year))

    return # write_qc_flag_info

#*******************************************
def write_metadata_flag_locations(index, timescale, string):
    """
    Write the metadata logfile - stations where lat/lon match
    """

    with open(os.path.join(utils.INFILELOCS, "Metadata_fails_{}.txt".format(index)), "a") as outfile:

        outfile.write("{} {}\n".format(timescale, string))

    return # write_qc_flag_locations

#*******************************************
def write_qc_flags(index, timescale, all_datasets, stations):
    """
    As QC runs through all source data in one go, need to split back into 
      different datasets before identifying the station and writing the files
    """

    # need to split up list into different datasets.
    
    stn_datasets = np.array([s.source for s in stations])
    dataset_names = np.array([d.name for d in all_datasets])
   
    for source in set(stn_datasets):
        
        dataset = all_datasets[dataset_names == source][0]
       
        stn_locs, = np.where(source == stn_datasets)

        # open and overwrite (as read in should keep QC info).
        with open(os.path.join(dataset.location, "{}.qc.{}.{}.txt".format(dataset.name, index, timescale)), "w") as outfile:

            for stn in stations[stn_locs]:
                outfile.write("{} {}\n".format(stn.id, stn.qc))

    return # write_qc_flags

#*******************************************
def read_datasets(all_datasets, timescale, index, qc_flags=" "):
    """
    Read in inventories from all the different contributing datasets

    :param list all_datasets: list of dataset objects
    :param str timescale: "MON" or "ANN"
    :param str index: index to process
    :param str qc_flags: always read in the QC flags (done always as using space " ") or select specific set

    :returns: array of stations
    """
    
    print("  Reading {} - {}".format(index, timescale))

    stations = np.array([])
    for dataset in all_datasets:

        try:
            ds_stations = utils.read_inventory(dataset, subdir="formatted/indices", final=True, timescale=timescale, index=index, qc_flags=qc_flags)
            good_stations = utils.select_qc_passes(ds_stations, qc_flags=qc_flags)
            stations = np.append(stations, good_stations)

            print("Adding {} ({} stations, {} BP), nstations = {}".format(dataset.name, len(good_stations), dataset.base_period, len(stations)))

        except IOError:
            # file missing
            print("No stations with data for {}".format(dataset.name))

    return stations # read_datasets

#*******************************************
def read_data(station, timescale, index, nyears, nmonths):
    """
    Read in the data for a specific station

    :param station station: station to be read
    :param str timescale: which timescale (MON/ANN)
    :param str index: index to be read
    :param int nyears: number of years - to define array
    :param int nmonths: number of months - to define array
    """


    # need to read in all the data to be able to process
    data = np.ma.zeros([nyears, nmonths])
    data[:] = -99.9
    data.fill_value = -99.9
    data.mask = np.ones(data.shape)
    
    times, indata = utils.read_station_index(station, index, timescale)
        
    match = np.in1d(utils.REFERENCEYEARS, times)
    match_back = np.in1d(times, utils.REFERENCEYEARS)
        
    if len(indata.shape) == 1:
        indata = indata.reshape([indata.shape[0], 1])
    
    data[match, :] = indata[match_back, :] # store all the info

    return data # read_data

#*********************************************
def check_lat_lon(locations, diagnostics=False):
    """
    Cross checks lat-lon pairs against all other stations

    :param array locations: array of lat, lon pairs
    :param bool diagnostics: extra output

    :returns: list of matches (list-indices)
    """

    matches = []

    for i, (lat, lon) in enumerate(locations):

        match, = np.where(np.logical_and(locations[:, 0] == lat, locations[:, 1] == lon))

        if len(match) == 1 and match[0] == i:
            # single hit
            continue
        elif match[0] == i:
            # a forward match, so no duplicates in the list
            matches += [match]

    return matches # check_lat_lon

#*********************************************
def file_exists_check(all_datasets, timescale, index, diagnostics=False):
    """
    Checks that file exists

    :param list all_datasets: list of dataset objects
    :param str timescale: "MON" or "ANN"
    :param str index: index to process
    :param bool diagnostics: extra output
    """
    print("File Exists - {}".format(timescale))
        
    nyears = len(utils.REFERENCEYEARS)
    if timescale == "MON":
        nmonths = 12
    else:
        nmonths = 1

    stations = read_datasets(all_datasets, timescale, index)    
    nstations = len(stations)

    # spin through each station
    nflags = 0
    for stn in stations:

        if not os.path.exists(os.path.join(stn.location, stn.id, "{}_{}_{}.csv".format(stn.id, index.lower(), timescale))):

            set_qc_attr(stn, "F")
            write_qc_flag_info(index, timescale, "File Exists Check", stn, "NULL")
            nflags += 1
            if diagnostics:
                print(stn, "F")

    write_qc_flags(index, timescale, all_datasets, stations)
    print("File Exists - {} flags".format(nflags))

    return # file_exists_check

#*********************************************
def empty_file_check(all_datasets, timescale, index, diagnostics=False):
    """
    Checks that file exists and has data

    :param list all_datasets: list of dataset objects
    :param str timescale: "MON" or "ANN"
    :param str index: index to process
    :param bool diagnostics: extra output
    """
    print("Empty File - {}".format(timescale))
        
    nyears = len(utils.REFERENCEYEARS)
    if timescale == "MON":
        nmonths = 12
    else:
        nmonths = 1

    stations = read_datasets(all_datasets, timescale, index)    
    nstations = len(stations)
    all_data = np.ma.zeros([len(stations), nyears, nmonths])

    nflags = 0
    # spin through each station
    for s, stn in enumerate(stations):

        indata = read_data(stn, timescale, index, nyears, nmonths)

        if len(indata.compressed()) == 0:
            set_qc_attr(stn, "E")
            write_qc_flag_info(index, timescale, "Empty File Check", stn, "NULL")
            nflags += 1
            if diagnostics:
                print(stn, "E")

        all_data[s, :, :] = indata
        
    write_qc_flags(index, timescale, all_datasets, stations)
    print("Empty File - {} flags".format(nflags))

    return stations, all_data # empty_file_check

#*********************************************
def metadata_exact_check(all_datasets, stations, timescale, index, diagnostics=False):
    """
    Calculates and then processes list-indices of stations which have exact matching lat/lon values
    
    :param list all_datasets: list of dataset objects
    :param array stations: array of station objects
    :param str timescale: "MON" or "ANN"
    :param str index: index to process
    :param bool diagnostics: extra output

    """
    print("Metadata - {}".format(timescale))

    # array of lats and lons for calculation of separations
    all_locations = np.array([[stn.latitude, stn.longitude] for stn in stations])
    names = [s.id for s in stations]

    #**************
    # Metadata Check
    matches = check_lat_lon(all_locations)   
    for match in matches:
        string = ""
        for m in match:
            stn = stations[m]
            string = string + "{}-{} lat {} lon {} ".format(stn.source, stn.id, stn.latitude, stn.longitude)

        print("{} {}".format(timescale, string))

        write_metadata_flag_locations(index, timescale, string) 

    # deconfounding done in wrapper

    return # metadata_exact_check

#*********************************************
def nearby_check(all_datasets, stations, timescale, index, diagnostics=False):
    """
    Calculates and then processes list-indices of stations which are closer than a threshold (5km)
 
    :param list all_datasets: list of dataset objects
    :param array stations: array of station objects
    :param str timescale: "MON" or "ANN"
    :param str index: index to process
    :param bool diagnostics: extra output

    """
    print("Nearby - {}".format(timescale))

    all_locations = np.array([[stn.latitude, stn.longitude] for stn in stations])
    names = [stn.id for stn in stations]

    # get the separations (km, radians)
    stn_separation, stn_angle = dls.get_separations(stations, all_locations)

    # convert to a masked array to enable removal of lower triangle
    stn_separation = np.ma.masked_array(stn_separation)
    stn_separation.mask = np.zeros(stn_separation.shape)

    tril = np.tril_indices(stn_separation.shape[0])
    stn_separation.mask[tril] = True
    # using != 0 as well ensures no overlap with metadata_exact_check()
    locsi, locsj = np.ma.where(np.logical_and(stn_separation > 0, stn_separation < NEARBY_SEP))

    # process the matches
    for i, j in zip(locsi, locsj):
        print("{},{}-{},{} sep = {}km".format(stations[i].source, names[i], stations[j].source, names[j], stn_separation[i, j]))
        
        # write metadata flag
        string = "{}-{} lat {} lon {} ".format(stations[i].source, stations[i].id, stations[i].latitude, stations[i].longitude)
        string = string + "{}-{} lat {} lon {} ".format(stations[j].source, stations[j].id, stations[j].latitude, stations[j].longitude)
        write_metadata_flag_locations(index, timescale, string) 

    # deconfounding done in wrapper

    return # nearby_check

#*********************************************
def select_metadata_match(all_datasets, stations, data, timescale, index, diagnostics=False):
    """
    Processes stations which have exact matching lat/lon values to select one to take forward
    
    :param list all_datasets: list of dataset objects
    :param list stations: list of station objects
    :param str timescale: "MON" or "ANN"
    :param str index: index to process
    :param bool diagnostics: extra output
    """

    nyears = len(utils.REFERENCEYEARS)
    if timescale == "MON":
        nmonths = 12
    else:
        nmonths = 1

    stn_ds = np.array([s.source for s in stations])
    stn_ids = np.array([s.id for s in stations])

    # use logfile as source
    with open(os.path.join(utils.INFILELOCS, "Metadata_fails_{}.txt".format(index)), "r") as infile:
        # spin through each line
        for line in infile:

            line = line.split()
            if line[0] != timescale:
                continue

            # extract the datasets and station IDs
            matches = line[1::5]
            ds, stn, indata = [], [], []

            # store in lists - along with the data.
            for m in matches:
                ds += [m.split("-")[0]]
                stn += ["-".join(m.split("-")[1:])] # some stations have "-" in the name.  No dataset has

                loc, = np.where(np.logical_and(stn_ds == ds[-1], stn_ids == stn[-1]))
                indata += [data[loc[0], :, :]]
            
            # and now do the tests
            keep = -1

            #**********
            # Check overlap
            """
            ACRE data will be very early period, so want to keep, even if later data from another source exists
            So check if there is overlap (generally useful anyway), but also if one source is ACRE, then keep
            """
            # ToDo


            #**********
            # Data Source
            """
            Sources from national archives maybe of higher quality, despite shorter record
              length.  Hence de-prioritise ECAD/SACAD/LACAD/GHCND/EOBS/LAOBS/SAOBS over other data
            """
            if len(np.unique(ds)) != 1:
                # then have more than one data source.
                country_loc, = np.where(np.in1d(ds, np.array(["ecad", "sacad", "lacad", "ghcnd", "eobs", "saobs", "laobs", "ghcndex"]), invert=True) == True)
                collection_loc, = np.where(np.in1d(ds, np.array(["ecad", "sacad", "lacad", "ghcnd", "eobs", "saobs", "laobs", "ghcndex"])) == True)
                
                if len(country_loc) == 1:
                    # one entry has source not from a collection - keep that one
                    keep = country_loc[0]
                    print("Keeping {}-{} - country".format(ds[country_loc[0]], stn[country_loc[0]]))
                elif len(country_loc) == 0:
                    # all entries have sources which are collections, but not a single collection
                    # filter later on length

                    if "ghcndex" in ds and "ghcnd" in ds:
                        # if both GHCNDEX and GHCND, then take GHCND
                        for d, s in zip(np.array(ds)[collection_loc], np.array(stn)[collection_loc]):
                            if d == "ghcndex":
                                loc, = np.where(np.logical_and(stn_ds == d, stn_ids == s))
                                this_stn = stations[loc][0]
                                set_qc_attr(this_stn, "M")                   
 
                        # need to reduce down the lists to test
                        indata = [indata[cl] for cl in collection_loc if ds[cl] != "ghcndex"]
                        stn = [stn[cl] for cl in collection_loc if ds[cl] != "ghcndex"]
                        ds = [ds[cl] for cl in collection_loc if ds[cl] != "ghcndex"]

                else:
                    # more than one entry can come from a single country - filter later on length
                    # AND/OR
                    # at least one entry have sources not from a collection
                    # flag any collection entries and move on
                    if len(collection_loc) != 0:
                        for d, s in zip(np.array(ds)[collection_loc], np.array(stn)[collection_loc]):

                            loc, = np.where(np.logical_and(stn_ds == d, stn_ids == s))
                            this_stn = stations[loc][0]
                            set_qc_attr(this_stn, "M")                    

                        # need to reduce down the lists to test
                        indata = [indata[cl] for cl in country_loc]
                        stn = [stn[cl] for cl in country_loc]
                        ds = [ds[cl] for cl in country_loc]

                    else:
                        # no collection entries
                        pass
            else:
                # then have only one data source.
                pass         
      

            #**********
            # Data Length
            if keep == -1:
                lengths = np.array([len(i.compressed()) for i in indata])
                locs, = np.where(lengths == max(lengths))

                if len(locs) == 1:
                    # one station is longest - use
                    keep = locs[0]
                    print("Keeping {}-{} - data length".format(ds[locs[0]], stn[locs[0]]))
                else:
                    # matches have same length - next test
                    pass

            if np.max(lengths) != 0:
                # then no data in any
                # both should have been flagged by empty file
                # but in this case flag both

                #**********
                # Runs latest
                if keep == -1:
                    latest = []
                    for i in indata:
                        latest += [np.where(i.mask == False)[0][-1]]
                    latest = np.array(latest)

                    locs, = np.where(latest == max(latest))

                    if len(locs) == 1:
                        # one station has data the latest - use
                        keep = locs[0]
                        print("Keeping {}-{} - runs latest".format(ds[locs[0]], stn[locs[0]]))
                    else:
                        # matches have data to same year - next test
                        pass

                #**********
                # Random selection of station as data availability the same
                if keep == -1:
                    # no test has selected a "winner"
                    keep = 0
                    print("Keeping {}-{} - random".format(ds[0], stn[0]))


                # And now assign flag values
                # remove "good" station from lists - so can place flags on others
                del ds[keep]
                del stn[keep]

            # run through remaining and set value
            for d, s in zip(ds, stn):

                loc, = np.where(np.logical_and(stn_ds == d, stn_ids == s))
                this_stn = stations[loc][0]
                set_qc_attr(this_stn, "M")
 
    return # select_metadata_match

#*********************************************
def metadata_checks(all_datasets, stations, data, timescale, index, nearby=True, diagnostics=False):
    """
    Runs two checks on metadata - one for exact match, one for a nearby match
    Allows for flexibility if necessary
    
    :param list all_datasets: list of dataset objects
    :param array stations: array of station objects
    :param array data: array of all data (nstations x nyears x nmonths)
    :param str index: index to process
    :param bool nearby: run the nearby check
    :param bool diagnostics: extra output
    """
 
    #**************
    # Metadata Checks
    metadata_exact_check(all_datasets, stations, timescale, index)

    #**************
    # Nearby Check
    if nearby:
        # Removed test, 19th November 2019
        nearby_check(all_datasets, stations, timescale, index)

    #**************
    # and now deconfound to sort flags
    print("Metadata - deconfound - {}".format(timescale))
    select_metadata_match(all_datasets, stations, data, timescale, index, diagnostics=diagnostics)
    write_qc_flags(index, timescale, all_datasets, stations)

    return # metadata_checks


#*********************************************
def world_record_check(all_datasets, stations, data, timescale, index, diagnostics=False):
    """
    Checks values against world records (Tx, Tn, P), using floating point tolerance

    :param list all_datasets: list of dataset objects
    :param array stations: array of station objects
    :param array data: array of all data (nstations x nyears x nmonths)
    :param str timescale: "MON" or "ANN"
    :param str index: index to process
    :param bool diagnostics: extra output

    """
    print("\n World Record - {}".format(timescale))

    
    # spin through each station
    nflags = 0
    for s, stn in enumerate(stations):

        if index in ["TNn", "TXn", "TNx", "TXx"]: 
            # check if too cold rather than too warm/wet
            locs = np.ma.where(np.logical_or(data[s] > WORLD_RECORDS_VALUES["Tx"]+TOLERANCE, np.logical_and(data[s] < WORLD_RECORDS_VALUES["Tn"]+TOLERANCE, data[s] != utils.HADEX_MDI)))
        else:
            locs = np.ma.where(data[s] > WORLD_RECORDS_VALUES[index]+TOLERANCE)
    
        if len(locs[0]) > 0:  
            set_qc_attr(stn, "W")
            write_qc_flag_info(index, timescale, "World Record Check", stn, utils.REFERENCEYEARS[locs[0]])
            nflags += 1
            if diagnostics:
                print(stn, "W")
            
    write_qc_flags(index, timescale, all_datasets, stations)
    print("World Record Check - {} flags".format(nflags))

    return # world_record_check

#*********************************************
def base_period_coverage_check(all_datasets, stations, data, timescale, index, diagnostics=False):
    """
    Checks that sufficient values exist over base period to allow thresholds to have been calculated

    :param list all_datasets: list of dataset objects
    :param array stations: array of station objects
    :param array data: array of all data (nstations x nyears x nmonths)
     :param str timescale: "MON" or "ANN"
    :param str index: index to process
    :param bool diagnostics: extra output
    """
    print("\n Base Period Coverage - {}".format(timescale))

    # spin through each station
    nflags = 0
    for s, stn in enumerate(stations):

        indata = data[s]

        base_period = indata[utils.REF_START - utils.STARTYEAR.year : utils.REF_END - utils.STARTYEAR.year + 1, :]

        assert base_period.shape[0] == utils.REF_END - utils.REF_START + 1

        n_present = np.ma.count(base_period, axis=0)

        insufficient, = np.where(n_present < int(utils.QC_PERCENTILES_COMPLETENESS * (utils.REF_END - utils.REF_START + 1)))

        # test to make sure that there would have been sufficient entries to calculate the percentiles!
        if len(insufficient) > 0:
            set_qc_attr(stn, "V")
            write_qc_flag_info(index, timescale, "Base Period Coverage Check", stn, insufficient)
            nflags += 1
            if diagnostics:
                print(stn, "V")

    write_qc_flags(index, timescale, all_datasets, stations)
    print("Base Period Coverage - {} flags".format(nflags))

    return # base_period_coverage_check

#*********************************************
def base_period_check(all_datasets, stations, data, timescale, index, diagnostics=False):
    """
    Checks that values over base period are within tolerated range

    :param list all_datasets: list of dataset objects
    :param array stations: array of station objects
    :param array data: array of all data (nstations x nyears x nmonths)
    :param str timescale: "MON" or "ANN"
    :param str index: index to process
    :param bool diagnostics: extra output
    """
    print("\n Base Period - {}".format(timescale))

    # spin through each station
    nflags = 0
    for s, stn in enumerate(stations):

        indata = data[s]

        base_period = indata[utils.REF_START - utils.STARTYEAR.year : utils.REF_END - utils.STARTYEAR.year + 1, :]

        assert base_period.shape[0] == utils.REF_END - utils.REF_START + 1

        average = np.ma.mean(base_period, axis=0)
         
        if index in PERCENTILE50_INDICES:
            EXPECTED = 50
            locs, = np.ma.where(np.logical_and(average > EXPECTED-2+TOLERANCE, average < EXPECTED+2-TOLERANCE))
        else:
            EXPECTED = 10
            locs, = np.ma.where(np.logical_and(average > EXPECTED-1+TOLERANCE, average < EXPECTED+1-TOLERANCE))

        if len(locs) > 0: 
            set_qc_attr(stn, "B")
            write_qc_flag_info(index, timescale, "Base Period Check", stn, locs)
            nflags += 1
            if diagnostics:
                print(stn, "B")

    write_qc_flags(index, timescale, all_datasets, stations)
    print("Base Period - {} flags".format(nflags))

    return # base_period_check

#*********************************************
def annual_in_monthly_check(all_datasets, astations, all_annual_data, index, upper=True, diagnostics=False):
    """
    Checks that annual values appear in monthly ones for appropriate indices
    As have to match across stations, reading in done in test

    :param list all_datasets: list of dataset objects
    :param array astations: array of station objects
    :param array all_annual_data: array of all data (nstations x nyears x nmonths)
    :param str index: index to process
    :param bool upper: whether annual value should be highest or lowest of monthly ones
    :param bool diagnostics: extra output

    """
    print("\n Annual in Monthly")
    nyears = len(utils.REFERENCEYEARS)

    mstations = read_datasets(all_datasets, "MON", index)
    # extract matching IDs
    m_ids = np.array([s.id for s in mstations])

    nflags = 0
    for astn, annual_station in enumerate(astations):
        
        # test if both annual and monthly exist.
        if os.path.exists(os.path.join(annual_station.location, annual_station.id, "{}_{}_{}.csv".format(annual_station.id, index.lower(), "ANN"))) \
                and\
                os.path.exists(os.path.join(annual_station.location, annual_station.id, "{}_{}_{}.csv".format(annual_station.id, index.lower(), "MON"))):

            annual_data = all_annual_data[astn].squeeze()
            monthly_data = read_data(annual_station, "MON", index, nyears, 12)

            good, = np.where(annual_data.mask == False)

            n_months = np.ma.count(monthly_data, axis=1)

            if not upper:
                extreme = np.ma.min(monthly_data, axis=1)

            else:
                extreme = np.ma.max(monthly_data, axis=1)

            sufficient, = np.where(n_months[good] == 12)

            # highlight where have annual value but not all 12 months.
            bad_locs, = np.where(annual_data[good][sufficient] != extreme[good][sufficient])

            if len(bad_locs) > 0:          
                set_qc_attr(annual_station, "A")
                set_qc_attr(mstations[m_ids == annual_station.id][0], "A")
                write_qc_flag_info(index, "ANN", "Annual in Monthly Check", annual_station, utils.REFERENCEYEARS[good[sufficient[bad_locs]]])
                nflags += 1
                if diagnostics:
                    print(stn, "A")

    # can't determine which is at fault
    write_qc_flags(index, "ANN", all_datasets, astations)
    write_qc_flags(index, "MON", all_datasets, mstations)
    print("Annual in Monthly - {} flags".format(nflags))

    return # annual_in_monthly_check

#*********************************************
# def monthly_vs_annual_check(all_datasets, index, diagnostics = False):
#     """
#     Checks that annual values only exist if all 12 monthly ones do
#     Missing day limitations means that may have all monthly values, but cumulatively
#       too few days to calculate annual
#     As matching between timescales, reading in done in test

#     :param list all_datasets: list of dataset objects
#     :param str index: index to process
#     :param bool diagnostics: extra output

#     """
#     print("\n Monthly vs Annual")
#     nyears = len(utils.REFERENCEYEARS)

#     mstations = read_datasets(all_datasets, "MON", index)    
#     astations = read_datasets(all_datasets, "ANN", index)

#     nflags = 0
#     for annual_station in astations:

#         annual_data = read_data(annual_station, "ANN", index, nyears, 1)[:, 0]
#         monthly_data = read_data(annual_station, "MON", index, nyears, 12)
        
#         n_months = np.ma.count(monthly_data, axis=1)

#         insufficient, = np.where(n_months < 12)

#         # highlight where have annual value but not all 12 months.
#         bad_locs, = np.where(annual_data.mask[insufficient] == False)

#         if len(bad_locs) > 0:          
#             write_qc_flag_info(index, "ANN", "Annual vs Monthly Check", annual_station, utils.REFERENCEYEARS[insufficient[bad_locs]])
#             nflags += 1

#     print("Monthly vs Annual - {} flags".format(nflags))
#     return # monthly_vs_annual_check

#*********************************************
def negative_value_check(all_datasets, stations, data, timescale, index, diagnostics=False):
    """
    Checks that all values for index are positive

    :param list all_datasets: list of dataset objects
    :param array stations: array of station objects
    :param array data: array of all data (nstations x nyears x nmonths)
    :param str timescale: "MON" or "ANN"
    :param str index: index to process
    :param bool diagnostics: extra output
    """
    print("\n Negative Value - {}".format(timescale))

    # spin through each station
    nflags = 0
    for s, stn in enumerate(stations):

        locs = np.ma.where(data[s] < 0.0-TOLERANCE)
    
        if len(locs[0]) > 0:
            set_qc_attr(stn, "N")
            write_qc_flag_info(index, timescale, "Negative Value Check", stn, utils.REFERENCEYEARS[locs[0]])
            nflags += 1
            if diagnostics:
                print(stn, "N")
            
    write_qc_flags(index, timescale, all_datasets, stations)
    print("Negative Value - {} flags".format(nflags))

    return # negative_value_check

#*********************************************
def consistency_check(all_datasets, gt_stations, all_gt_data, timescale, gt_index, lt_index, diagnostics=False):
    """
    Checks that values for paired indices are consistent - TXx > TNx etc
    Matching between indices - reading in done in test.

    :param list all_datasets: list of dataset objects
    :param array gt_stations: array of station objects
    :param array all_gt_data: array of all data (nstations x nyears x nmonths)
    :param str timescale: "MON" or "ANN"
    :param str gt_index: index with expected higher value to process
    :param str lt_index: index with expected lower valueto process
    :param bool diagnostics: extra output

    """
    print("\n Consistency between indices - {}".format(timescale))
    nyears = len(utils.REFERENCEYEARS)
    if timescale == "MON":
        nmonths = 12
    else:
        nmonths = 1

    lt_stations = read_datasets(all_datasets, timescale, lt_index)
    # extract matching IDs
    lt_ids = np.array([s.id for s in lt_stations])

    nstations = len(gt_stations)

    nflags = 0
    for gtstn, station in enumerate(gt_stations):


        # only run if files for both indices exist
        if (station.id in lt_ids) and \
              (os.path.exists(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, gt_index.lower(), timescale)))) and \
                (os.path.exists(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, lt_index.lower(), timescale)))):

            gt_data = all_gt_data[gtstn]
            lt_data = read_data(station, timescale, lt_index, nyears, nmonths)

            year, month = np.ma.where(lt_data > gt_data+TOLERANCE)
            
            if len(year) > 0:
                set_qc_attr(station, "C")
                set_qc_attr(lt_stations[lt_ids == station.id][0], "C")

                write_qc_flag_info(gt_index, timescale, "Consistency Check", station, utils.REFERENCEYEARS[year])
                nflags += 1
                if diagnostics:
                    print(stn, "C")

        else:
            if diagnostics:
                print("File doesn't exist for one of the indices")
                if not os.path.exists(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, gt_index.lower(), timescale))):
                    print(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, gt_index.lower(), timescale)))
                if not os.path.exists(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, lt_index.lower(), timescale))):
                    print(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, lt_index.lower(), timescale)))

    # can't determine which is the issue, flag both
    write_qc_flags(gt_index, timescale, all_datasets, gt_stations)
    write_qc_flags(lt_index, timescale, all_datasets, lt_stations)
    print("Consistency - {} flags".format(nflags))

    return # temporal_consistency_check

#*********************************************
def correlation_check(all_datasets, timescale, index, diagnostics=False):
    """
    Checks that annual value only exists if all monthly values are present

    :param list all_datasets: list of dataset objects
    :param str timescale: "MON" or "ANN"
    :param str index: index to process
    :param bool diagnostics: extra output

    """
    print("\n Correlation - {}".format(timescale))
    nyears = len(utils.REFERENCEYEARS)
    if timescale == "MON":
        nmonths = 13 # dls includes annual read in monthly in column 0
    else:
        nmonths = 1

    # two reads - need full list and reduced list
    all_stations = read_datasets(all_datasets, timescale, index)  # dummy string to ensure values are retained.
    stations = read_datasets(all_datasets, timescale, index, qc_flags="M")  

    all_locations = np.array([[stn.latitude, stn.longitude] for stn in stations])
    names = [stn.id for stn in stations]

    # read in all the station data
    all_data = dls.get_all_data(stations, index, timescale, nyears, nmonths)

    # get the separations (km, radians)
    stn_separation, stn_angle = dls.get_separations(stations, all_locations)

    for month in range(nmonths):
        print("month {}".format(month))
        if timescale == "MON" and month == 0:
            print("skipping annual as done separately")
            continue

        seps, cors = dls.separations_and_correlations(all_data[:, :, month], stn_separation, names, diagnostics=diagnostics, flatten=False)
        # use a 250km cutoff
        spurious_corrs_i, spurious_corrs_j = np.ma.where(np.logical_and(cors > CORR, seps > SEP))

        nflags = 0
        for i, j in zip(spurious_corrs_i, spurious_corrs_j):
            print("{},{}-{},{} r = {}, sep = {}km".format(stations[i].source, names[i], stations[j].source, names[j], cors[i, j], seps[i, j]))

            # flag both stations as can't be sure which is the correct one (easily!)
            set_qc_attr(stations[i], "R")
            set_qc_attr(stations[j], "R")           
            write_qc_flag_info(index, timescale, "Correlation Check", stations[i], month)
            write_qc_flag_info(index, timescale, "Correlation Check", stations[j], month)
            nflags += 2
            if diagnostics:
                print(stations[i], stations[j], "R")

        if False:
            # To see what whas going on with low-correlation, close-by stations
            spurious_corrs_i, spurious_corrs_j = np.ma.where(np.logical_and(cors < 0.9, seps < SEP))
            
            for i, j in zip(spurious_corrs_i, spurious_corrs_j):
                print("{},{}-{},{} r = {}, sep = {}km".format(stations[i].source, names[i], stations[j].source, names[j], cors[i, j], seps[i, j]))

        seps = np.reshape(seps, [-1]).compressed()
        cors = np.reshape(cors, [-1]).compressed()
        
        print("Correlation {} - {} flags".format(month, nflags))
 
    # match reduced list with complete list
    all_names = np.array(["{}-{}".format(s.source, s.id) for s in all_stations])
    names = np.array(["{}-{}".format(s.source, s.id) for s in stations])
    match = np.in1d(all_names, names)

    # update stations in complete list
    all_stations[match] = stations

    write_qc_flags(index, timescale, all_datasets, all_stations)

    return # correlation_check


#*********************************************
def main(index="TX90p", diagnostics=False):
    """
    The main QC function

    :param str index: which index to run
    :param bool diagnostics: extra verbose output
    """

    # move this up one level eventually?
    all_datasets = utils.get_input_datasets()

    # removes all files from previous runs
    setup_outfile(index, all_datasets)

    #**************
    # File Exists Check
    if index in utils.MONTHLY_INDICES:
        file_exists_check(all_datasets, "MON", index, diagnostics=diagnostics)
    # and always run on annual
    file_exists_check(all_datasets, "ANN", index, diagnostics=diagnostics)

    #**************
    # Empty Check - also reads in all the data
    if index in utils.MONTHLY_INDICES:
        monthly_stations, monthly_data = empty_file_check(all_datasets, "MON", index, diagnostics=diagnostics)
    # and always run on annual
    annual_stations, annual_data = empty_file_check(all_datasets, "ANN", index, diagnostics=diagnostics)
       
    #**************
    # World Record Check
    if index in WORLD_RECORDS_INDICES:
        if index in utils.MONTHLY_INDICES:
            world_record_check(all_datasets, monthly_stations, monthly_data, "MON", index, diagnostics=diagnostics)
        # and always run on annual
        world_record_check(all_datasets, annual_stations, annual_data, "ANN", index, diagnostics=diagnostics)
        
    #**************
    # Base Period Coverage heck
    if index in PERCENTILE_INDICES:
        if index in utils.MONTHLY_INDICES:
            base_period_coverage_check(all_datasets, monthly_stations, monthly_data, "MON", index, diagnostics=diagnostics)
        # and always run on annual
        base_period_coverage_check(all_datasets, annual_stations, annual_data, "ANN", index, diagnostics=diagnostics)

    #**************
    # Base Period Check
    if index in PERCENTILE_INDICES:
        if index in utils.MONTHLY_INDICES:
            base_period_check(all_datasets, monthly_stations, monthly_data, "MON", index, diagnostics=diagnostics)
        # and always run on annual
        base_period_check(all_datasets, annual_stations, annual_data, "ANN", index, diagnostics=diagnostics)
        
#    #**************
#    # Monthly vs Annual Check
#    if index in utils.MONTHLY_INDICES:
#        monthly_vs_annual_check(all_datasets, index, diagnostics=diagnostics)

    #**************
    # Annual in Monthly Check
    if index in ANNUAL_IN_MONTHLY:
        upper = True
        if index in ["TNn", "TXn"]: 
            upper = False     
        annual_in_monthly_check(all_datasets, annual_stations, annual_data, index, upper=upper, diagnostics=diagnostics)

    #**************
    # Negative Values Check
    if index in NO_NEGATIVES:
        if index in utils.MONTHLY_INDICES:
            negative_value_check(all_datasets, monthly_stations, monthly_data, "MON", index, diagnostics=diagnostics)
        # and always run on annual
        negative_value_check(all_datasets, annual_stations, annual_data, "ANN", index, diagnostics=diagnostics)

    #**************
    # Consistency Check
    if index in GT_INDEX:
        lt_index = LT_INDEX[index]
        # TXx the maximum for two other indices, and can't have duplicated keys, so work around
        if index == "TXx2":
            gt_index = "TXx"
        else:
            gt_index = index[:]

        if (gt_index in utils.MONTHLY_INDICES) and (lt_index in utils.MONTHLY_INDICES):
            consistency_check(all_datasets, monthly_stations, monthly_data, "MON", gt_index, lt_index, diagnostics=diagnostics)
        consistency_check(all_datasets, annual_stations, annual_data, "ANN", gt_index, lt_index, diagnostics=diagnostics)

    #**************
    # Metadata Checks
    #   Only cross check all tests
    if index in utils.MONTHLY_INDICES:
        metadata_checks(all_datasets, monthly_stations, monthly_data, "MON", index, nearby=False, diagnostics=diagnostics)
    # and always run on annual
    metadata_checks(all_datasets, annual_stations, annual_data, "ANN", index, nearby=False, diagnostics=diagnostics)

    #**************
    # Correlation Check [Duplicates] (and most intensive check)
    #   Do not run for threshold ID/FD/SU/TR as easily have lots of zeros or complete years with all/none set
    #     also R10mm/R20mm as also can have lots of zeros
    #   in tropics/high latitudes
    #   Also exclude stations where metadata cross checks have flagged (thins network out!)
    if index not in NOT_CORRELATION:
        if index in utils.MONTHLY_INDICES:
            correlation_check(all_datasets, "MON", index, diagnostics=diagnostics)
        correlation_check(all_datasets, "ANN", index, diagnostics=diagnostics)

    return # main
    

#************************************************************************
if __name__ == "__main__":

    pass
    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', dest='index', action='store', default="TX90p",
                        help='Which index to run')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default=False')

    args = parser.parse_args()
          
    main(index=args.index, diagnostics=args.diagnostics)


#*******************************************
# END
#*******************************************
