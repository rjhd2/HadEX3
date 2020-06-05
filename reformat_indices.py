#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 488                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2020-04-20 15:50:26 +0100 (Mon, 20 Apr 2020) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
'''
Reformats the raw files to CLIMPACT2 format

reformat_indices.py invoked by typing::

  python reformat_indices.py --indata "ecad" --index "TX90p" --unzip --convert  --diagnostics

Input arguments:

--indata        Which input data to use ("ecad, lacad, sacad, west_africa_pptn, west_africa_indices, south_america, arabia)

--index         Which ETCCDI index to use

--unzip         Unzip the ECAD files (default True)

--convert       Convert to CLIMPACT format (default True)

--diagnostics   Output extra info (default False)
'''
# NOTE - must run (unzip) PRCPTOT before R95p and R99p

import os
import calendar
import copy
import numpy as np

# RJHD utils
import utils
import convert_utils
import inventory_utils

#*********************************************
def unzip_ecad_family(dataset, index, diagnostics=False):
    """
    Unzip the ECAD downloaded files, and in some cases further process because of inconsistent naming

    :param datasetObj dataset: dataset object for metadata
    :param str index: ETCCDI index name to process
    :param bookl diagnostics: output diagnostic information
    """

    import zipfile
    import shutil

    #*********************
    # unzip the file - special bits for poorly labelled indices from ECAD

    if index in ["Rx1day", "Rx5day"]:
        index_name = "X".join(index.split("x"))
    else:
        index_name = index
            
    try:
        zip_ref = zipfile.ZipFile(os.path.join(dataset.location, "raw", "{}_index{}.zip".format(dataset.name.upper()[:-1], index_name)), 'r')
    except IOError:
        print("No file: {}".format(os.path.join(dataset.location, "raw", "{}_index{}.zip".format(dataset.name.upper()[:-1], index_name))))
        return


    # make the folder to hold the converted data
    if not os.path.exists(os.path.join(dataset.location, "raw", index)):
        os.mkdir(os.path.join(dataset.location, "raw", index))

    # extract the file and close the handler
    zip_ref.extractall(os.path.join(dataset.location, "raw", index))
    zip_ref.close()
    
    #*********************
    # move stations.txt metadata file to a unique filename                
    shutil.move(os.path.join(dataset.location, "raw", index, "stations.txt"), os.path.join(dataset.location, "{}_stations.txt".format(index)))

    #*********************
    # and rename files for poorly labelled indices from ECAD
    if index in ["Rx1day", "Rx5day"]:
        files = os.listdir(os.path.join(dataset.location, "raw", index))

        for infile in files:
            # to allow overwriting, only move if it has the [malformed] index_name
            if index_name in infile:

                root, stn = infile.split(index_name)
                shutil.move(os.path.join(dataset.location, "raw", index, infile), os.path.join(dataset.location, "raw", index, "{}{}{}".format(root, index, stn)))

    return # unzip_ecad_family

#*********************************************
def read_ecad_file(filename, dataset, index):
    """
    Read the ECAD index file - including metadata.

    :param str filename: full filename to process
    :param datasetObj dataset: dataset object for metadata
    :param str index: ETCCDI index name to process
    """

    if dataset.name == "ecad":
        OFFSET = 0
    elif dataset.name == "lacad":
        OFFSET = +1
    elif dataset.name == "sacad":
        OFFSET = 0

    try:
        with open(filename, "r", encoding="latin-1") as infile:

            # spin through each line for metadata
            for lc, line in enumerate(infile):

                if lc == 9 - 1 + OFFSET:
                    file_mdi = int(line.split("=")[-1].strip()[:-2])
                if lc == 11 - 1 + OFFSET:
                    station_name = line[74:115].strip()
                if lc == 12 - 1 + OFFSET:
                    lat_dms = line.split("ss")[-1].split(":")
                    lat_dir = "N"
                    if lat_dms[0][0] == "-":
                        lat_dir = "S"
                    station_lat = utils.dms2dd(lat_dms[0], lat_dms[1], lat_dms[2], direction=lat_dir)
                if lc == 13 - 1 + OFFSET:
                    lon_dms = line.split("ss")[-1].split(":")
                    lon_dir = "E"
                    if lon_dms[0][0] == "-":
                        lon_dir = "W"
                    station_lon = utils.dms2dd(lon_dms[0], lon_dms[1], lon_dms[2], direction=lon_dir)
                if lc == 14 - 1 + OFFSET:
                    station_elv = int(line.split("m")[-1])
                if lc == 16 - 1 + OFFSET:
                    file_scale = float(line.split(" with unit ")[-1].split(" ")[0])
                    break # finished getting metadata

        # get the rest of the raw data

        if index in ["CSDI", "WSDI", "R99pTOT", "R95pTOT", "SDII"] and dataset.name in ["ecad", "lacad", "sacad"]:
            # ECAD et al run these as monthly, ETCCDI treats as annual.
            all_data = np.genfromtxt(filename, skip_header=30 + OFFSET, encoding="latin-1")
            if len(all_data.shape) == 1:
                # expand array dimensions
                all_data = np.ma.expand_dims(all_data, axis=0) 

            # and just retain annual values
            all_data = all_data[:, :3]

        # rest should be standard
        elif index in utils.MONTHLY_INDICES:
            all_data = np.genfromtxt(filename, skip_header=30 + OFFSET, encoding="latin-1")
            if len(all_data.shape) == 1:
                # expand array dimensions
                all_data = np.ma.expand_dims(all_data, axis=0) 
        else:
            all_data = np.genfromtxt(filename, skip_header=21 + OFFSET, encoding="latin-1")
            if len(all_data.shape) == 1:
                # expand array dimensions
                all_data = np.ma.expand_dims(all_data, axis=0) 

    except IOError:
        print("{} missing".format(filename))
        raise IOError

    # and just return the annual and monthly info
    years = all_data[:, 1].astype(int)
    annual = all_data[:, 2]

    try:
        monthly = all_data[:, 9:]
    except IndexError:
        # doesn't exist or empty, so return all as missing
        monthly = np.zeros(annual.shape[0], 12)
        monthly[:] = utils.HADEX_MDI
        

    # return tuple of metadata, extra information and the data
    return (station_name, station_lat, station_lon, station_elv), file_mdi, file_scale, years, annual, monthly # read_ecad_file

#*********************************************
def read_hadex2_file(filename, dataset, index):
    """
    Read the HadEX2 index file - including metadata.

    :param str filename: full filename to process
    :param datasetObj dataset: dataset object for metadata
    :param str index: ETCCDI index name to process
    """

    indata = np.genfromtxt(filename, skip_header=1, encoding="latin-1")

    years = indata[:, 0]
    
    if indata.shape[1] == 2:
        # just annual data
        annual = indata[:, 1]
        monthly = np.zeros((annual.shape[0], 12))
        monthly[:] = utils.HADEX_MDI
        
    elif indata.shape[1] == 14:
        # annual and monthly data
        annual = indata[:, -1]
        monthly = indata[:, 1:13]     

    # return nulls for other output for match
    station_name, station_lat, station_lon, station_elv = "", "", "", "" # to match output from ECAD family
    file_mdi = utils.HADEX_MDI
    file_scale = 1.0

    # return tuple of metadata, extra information and the data
    return (station_name, station_lat, station_lon, station_elv), file_mdi, file_scale, years, annual, monthly  # read_hadex2_file

#*********************************************
def read_station(dataset, filename, station, index, doMonthly=False):
    """
    Read in the raw station metadata and data and convert into station object

    :param datasetObj dataset: dataset object for metadata
    :param str filename: location of raw data file
    :param stationObj station: station object to hold information
    :param str index: ETCCDI index being processed
    :param bool doMonthly: if a monthly index
    """
    
    # read in all the data, metadata and extra information
    if dataset.name in ["ecad", "lacad", "sacad"]:
        metadata, mdi, scale, years, annual, monthly = read_ecad_file(filename, dataset, index)
    elif dataset.name in ["hadex2", "ghcndex"]:
        metadata, mdi, scale, years, annual, monthly = read_hadex2_file(filename, dataset, index)

    # if no actual data then make empty arrays and exit
    if annual.shape[0] == 0:
        station.months = []
        station.years = []
        station.myears = []
        station.monthly = []
        station.annual = []
        return station, metadata # skip out early

    #*********************
    # process monthly indices
    # extra processing if monthly to convert columns into a single timeseries
    if doMonthly:
        # fake the months
        months = np.array([range(1, 13) for i in years]).reshape([-1])
        myears = np.array([[year for i in range(12)] for year in years]).reshape([-1])

        monthly = monthly.reshape([-1])
        
        # apply scale and fix missing values
        if scale != 1.0:
            monthly *= scale
            monthly[monthly == mdi*scale] = utils.HADEX_MDI

        # fix TX90p etc and R95pTOT etc from ECAD to be in correct units (ANNUAL)
        if (dataset.name in ["ecad", "lacad", "sacad"]) and \
                (index in ["TX90p", "TX10p", "TN90p", "TN10p"]):
            # need to convert from days to %
            ndays = np.array([calendar.monthrange(y, m)[1] for y, m in zip(myears, months)])
            good, = np.where(monthly != utils.HADEX_MDI)
            monthly[good] = 100. * (monthly[good] / ndays[good])

        # and store
        station.monthly = monthly
        station.myears = myears
        station.months = months

    else:
        # store empty values otherwise
        station.monthly = None
        station.myears = None
        station.months = None

    #*********************
    # process annual indices
    # apply scale and fix missing values
    if scale != 1.0:
        annual *= scale
        annual[annual == mdi*scale] = utils.HADEX_MDI

    # fix TX90p etc and R95pTOT etc from ECAD to be in correct units (ANNUAL)
    if (dataset.name in ["ecad", "lacad", "sacad"]) and \
            (index in ["TX90p", "TX10p", "TN90p", "TN10p"]):
        # need to convert from days to %
        ndays = np.array([366 if calendar.isleap(y) else 365 for y in years])
        good, = np.where(annual != utils.HADEX_MDI)
        annual[good] = 100. * (annual[good] / ndays[good])

    station.annual = annual
    station.years = years
 
    return station, metadata # read_station

#*********************************************
def convert_ecad_family(dataset, index, diagnostics=False):
    """
    Controller routine for ECAD index conversion
   
    :param datasetObj dataset: dataset to process (ECAD/LACAD/SACAD)
    :param str index: which ETCCDI index to process
    :param bool diagnostics: output diagnostic information
    """

    # read in the list of available stations for this index
    station_list = inventory_utils.read_ecad(dataset, index, diagnostics=diagnostics)

    # if no stations, then just exit (inventory_utils.read_ecad returns stations as [] if empty)
    if station_list == []:
        return

    # spin through each station (just metadata so far)
    for station in station_list:
        
        if diagnostics:
            print(station.id)

        if index in ["R95p", "R99p"]:
            """
            Need to fix inconsistency in index calculation for these two - see Donat et al, 2013 para 11.

            Have to re-read PRCPTOT and R??pTOT indices and recalculate.
            """

            #*********************
            # Read the PRCPTOT index

            # make a copy so that the dare are not stored in the station list.
            P_tot_station = copy.deepcopy(station)
            P_name = os.path.join(dataset.location, "raw", "PRCPTOT", "index{}{}.txt".format("PRCPTOT", P_tot_station.id))

            if index in utils.MONTHLY_INDICES:
                P_tot_station, P_metadata = read_station(dataset, P_name, P_tot_station, "PRCPTOT", doMonthly=True)
            else:
                P_tot_station, P_metadata = read_station(dataset, P_name, P_tot_station, "PRCPTOT")

            #*********************
            # Read the R??pTOT index

            # make a copy so that the dare are not stored in the station list.
            R_station = copy.deepcopy(station)
            R_name = os.path.join(dataset.location, "raw", "{}TOT".format(index), "index{}TOT{}.txt".format(index, R_station.id))

            if index in utils.MONTHLY_INDICES:
                R_station, R_metadata = read_station(dataset, R_name, R_station, "{}TOT".format(index), doMonthly=True)
            else:
                R_station, R_metadata = read_station(dataset, R_name, R_station, "{}TOT".format(index))

            #*********************
            # check all metadata matches!
            for r, p in zip(R_metadata, P_metadata):
                assert r == p

            #*********************
            # convert %age to actuals by multiplying by PRCPTOT and rescaling
            this_station = copy.deepcopy(R_station)
            if index in utils.MONTHLY_INDICES:
                R_station.monthly = np.ma.masked_where(R_station.monthly == utils.HADEX_MDI, R_station.monthly)
                P_tot_station.monthly = np.ma.masked_where(P_tot_station.monthly == utils.HADEX_MDI, P_tot_station.monthly)

                if len(R_station.monthly) == len(P_tot_station.monthly) and len(P_tot_station.monthly) != 0:
                    this_station.monthly = (R_station.monthly * P_tot_station.monthly) / 100.
                else:
                    this_station.monthly = []
                
                this_station.monthly.fill_value = utils.HADEX_MDI
                this_station.monthly = this_station.monthly.filled()

                this_station.myears = R_station.myears
                this_station.months = R_station.months

            if len(R_station.annual) == len(P_tot_station.annual) and len(P_tot_station.annual) != 0:
                R_station.annual = np.ma.masked_where(R_station.annual == utils.HADEX_MDI, R_station.annual)
                P_tot_station.annual = np.ma.masked_where(P_tot_station.annual == utils.HADEX_MDI, P_tot_station.annual)

                this_station.annual = (R_station.annual * P_tot_station.annual) / 100.

                this_station.annual.fill_value = utils.HADEX_MDI
                this_station.annual = this_station.annual.filled()
            else:
                this_station.annual = np.zeros(len(R_station.years))
                this_station.annual.fill(utils.HADEX_MDI)
            this_station.years = R_station.years

        else:
            # make a copy so that the dare are not stored in the station list.
            this_station = copy.deepcopy(station)
            filename = os.path.join(dataset.location, "raw", index, "index{}{}.txt".format(index, this_station.id))

            # extract station data
            if index in utils.MONTHLY_INDICES:
                this_station, metadata = read_station(dataset, filename, this_station, index, doMonthly=True)      
            else:
                this_station, metadata = read_station(dataset, filename, this_station, index)      
    
        # prepare to write the data - check all directories present
        utils.check_directories(dataset, this_station)

        # and write the data
        if index in utils.MONTHLY_INDICES:
            utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", this_station.id, "{}_{}_MON.csv".format(this_station.id, index.lower())), this_station, index, doMonthly=True)

        utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", this_station.id, "{}_{}_ANN.csv".format(this_station.id, index.lower())), this_station, index)

    return # convert_ecad_family

#*********************************************
def convert_hadex2(dataset, index, diagnostics=False):
    """
    Controller routine for HadEX2 index conversion for workshop stations
   
    :param datasetObj dataset: dataset to process (HadEX2)
    :param str index: which ETCCDI index to process
    :param bool diagnostics: output diagnostic information
    """

    # read in the list of available stations for this index
    station_list = inventory_utils.read_hadex2(dataset, index, diagnostics=diagnostics)

    # if no stations, then just exit (inventory_utils.read_hadex2 returns stations as [] if empty)
    if station_list == []:
        return

    # spin through each station (just metadata so far)
    for station in station_list:
        if diagnostics:
            print(station.id)
        
        filename = os.path.join(dataset.location, "raw", index, "{}_{}.txt".format(station.id, index))

        # read in station data           
        if index in utils.MONTHLY_INDICES:
            station, metadata = read_station(dataset, filename, station, index, doMonthly=True)      
        else:
            station, metadata = read_station(dataset, filename, station, index)      

        # prepare to write the data - check all directories present
        utils.check_directories(dataset, station)

        if index in utils.MONTHLY_INDICES:
            utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", station.id, "{}_{}_MON.csv".format(station.id, index.lower())), station, index, doMonthly=True)
        utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", station.id, "{}_{}_ANN.csv".format(station.id, index.lower())), station, index)

    return # convert_hadex2

#*********************************************
def convert_south_america(dataset, index, diagnostics=False):
    """
    Controller routine for South American index conversion for workshop stations (v20180813)
   
    :param datasetObj dataset: dataset to process (HadEX2)
    :param str index: which ETCCDI index to process
    :param bool diagnostics: output diagnostic information
    """

    COUNTRIES = ["Colombia", "Suriname", "Venezuela", "Peru"]
    
    # don't have inventory information, so have to create
    all_files = os.listdir(os.path.join(dataset.location, "raw"))

    # spin through each entry and find matches (capitalization issues)
    infiles = []
    for af in all_files:
        if index.lower() in af.lower() and af.split()[0] != "Registros":
            infiles += [os.path.join(dataset.location, "raw", af)]

    # set up storage for this metadata (INDEX_stations.txt doesn't clash with south_america_stations.txt from raw data)
    metadata_file = open(os.path.join(dataset.location, "{}_stations.txt".format(index)), "w")

    # now spin through each file
    for infilename in infiles:
        if diagnostics:
            print(infilename)

        country = infilename.split("/")[-1].split("_")[1]
        # only interested in those countries where don't have raw data
# removed 20191002 - want as many stations as possible
#        if country not in COUNTRIES:
#            continue
        
        all_indata = np.genfromtxt(infilename, delimiter=",", dtype="U16", encoding="latin-1")

        years = all_indata[1:, 0]
        stations = all_indata[0, 1:]

        # ignore blank lines or entries
        good_locs = np.where(years != "")
        years = years[good_locs].astype(int)
        
        indata = all_indata[1:, 1:][good_locs]
        indata[indata == ""] = "{:5.1f}".format(utils.HADEX_MDI)
        indata = indata.astype(float)

        # go through all stations, find metadata if available and write out if so.
        for s, raw_id in enumerate(stations):

            if diagnostics:
                print(country, raw_id)
            
            this_station = convert_utils.find_south_american_metadata(dataset, country, raw_id)

            if this_station != -1:

                this_station.years = years
                this_station.annual = indata[:, s]
                # store empty values for monthly data
                this_station.monthly = None
                this_station.myears = None
                this_station.months = None

                # update station id to avoid overwriting raw south_american data
                this_station.id = "{}i".format(this_station.id)

                # prepare to write the data - check all directories present
                utils.check_directories(dataset, this_station)

                # annual only data
                utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", this_station.id, "{}_{}_ANN.csv".format(this_station.id, index.lower())), this_station, index)

                # write metadata!
                metadata_file.write("{:10s} {:7.3f} {:7.3f}\n".format(this_station.id, this_station.latitude, this_station.longitude))

    metadata_file.close()

    return # convert_south_america

#*********************************************
def convert_west_africa_pptn(dataset, index, diagnostics=False):
    """
    Controller routine for West African index conversion from Theo Vischel
   
    :param datasetObj dataset: dataset to process (HadEX2)
    :param str index: which ETCCDI index to process
    :param bool diagnostics: output diagnostic information
    """
    
    # read inventory
    station_list = inventory_utils.read_west_africa_pptn(dataset, diagnostics=diagnostics)

    # set up storage for this metadata
    metadata_file = open(os.path.join(dataset.location, "{}_stations.txt".format(index)), "w")

    # spin through each station
    for this_station in station_list:

        # see if a file exists for this index (only some of the precip ones, really).
        try:
            if diagnostics:
                print(this_station.id)

            indata = np.genfromtxt(os.path.join(dataset.location, "raw", "{}_{}.csv".format(this_station.id, index.upper())), delimiter=",", skip_header=1, missing_values="NA", filling_values=utils.HADEX_MDI, dtype=float, encoding="latin-1")

            # extract the two columns
            years = indata[:, 0].astype(int)
            annual = indata[:, 1]

            # and write into the attributes
            this_station.years = years
            this_station.annual = annual
            # store empty values for monthly data
            this_station.monthly = None
            this_station.myears = None
            this_station.months = None

            # prepare to write the data - check all directories present
            utils.check_directories(dataset, this_station)
                
            # annual only data
            utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", this_station.id, "{}_{}_ANN.csv".format(this_station.id, index.lower())), this_station, index)

            # write metadata!
            metadata_file.write("{:10s} {:7.3f} {:7.3f}\n".format(this_station.id, this_station.latitude, this_station.longitude))

        except IOError:
            print("file missing {}".format(os.path.join(dataset.location, "raw", "{}_{}.csv".format(this_station.id, index.upper()))))

    metadata_file.close()
    return # convert_west_africa_pptn

#*********************************************
def convert_west_africa_indices(dataset, index, diagnostics=False):
    """
    Controller routine for West African index conversion from Aziz Barry
   
    :param datasetObj dataset: dataset to process (HadEX2)
    :param str index: which ETCCDI index to process
    :param bool diagnostics: output diagnostic information
    """
    
    # read inventory
    station_list = inventory_utils.read_west_africa_indices(dataset, diagnostics=diagnostics)

    # set up storage for this metadata
    metadata_file = open(os.path.join(dataset.location, "{}_stations.txt".format(index)), "w")

    # rename the indices as necessary
    if index in ["TNn", "TNx", "TXn", "TXx", "TR", "SU"]:
        rename = {"TNn" : "Tnn", "TNx" : "tnx", "TXn" : "txn", "TXx" : "txx", "SU" : "SU25", "TR" : "TR20"}
        index_name = rename[index]
    else:
        index_name = index    

    # spin through each station
    for this_station in station_list:

        infilename = ""
        if os.path.exists(os.path.join(dataset.location, "raw", "{}_{}.txt".format(this_station.id.upper(), index_name))):
            infilename = os.path.join(dataset.location, "raw", "{}_{}.txt".format(this_station.id.upper(), index_name))
        elif os.path.exists(os.path.join(dataset.location, "raw", "{}_{}.txt".format(this_station.id.lower(), index_name))):
            infilename = os.path.join(dataset.location, "raw", "{}_{}.txt".format(this_station.id.lower(), index_name))
        elif os.path.exists(os.path.join(dataset.location, "raw", "{}_{}.txt".format(this_station.id.capitalize(), index_name))):
            infilename = os.path.join(dataset.location, "raw", "{}_{}.txt".format(this_station.id.capitalize(), index_name))
        else:
            print("file for {} -- {} missing in {}".format(this_station.id, index_name, os.path.join(dataset.location, "raw")))


        # see if a file exists for this index (only some of the precip ones, really).
        try:
            if diagnostics:
                print(this_station.id)


            if sum(1 for line in open(infilename)) == 0:
                print("file empty {}".format(infilename))

            else:
                indata = np.genfromtxt(infilename, missing_values="NA", filling_values=utils.HADEX_MDI, dtype=float, encoding="latin-1")

                # extract the two columns
                years = indata[:, 0].astype(int)
                annual = indata[:, 1]

                # and write into the attributes
                this_station.years = years
                this_station.annual = annual
                # store empty values for monthly data
                this_station.monthly = None
                this_station.myears = None
                this_station.months = None

                # prepare to write the data - check all directories present
                utils.check_directories(dataset, this_station)

                # annual only data
                utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", this_station.id, "{}_{}_ANN.csv".format(this_station.id, index.lower())), this_station, index)

                # write metadata!
                metadata_file.write("{:10s} {:7.3f} {:7.3f}\n".format(this_station.id, this_station.latitude, this_station.longitude))

        except IOError:
            print("file missing {}".format(os.path.join(dataset.location, "raw", "{}_{}.csv".format(this_station.id, index.upper()))))

    metadata_file.close()
    return # convert_west_africa_indices

#*********************************************
def convert_arabia_6190(dataset, index, diagnostics=False):
    """
    Controller routine for Arabian index conversion from Markus Donat (long record subset, 1961-90)
   
    :param datasetObj dataset: dataset to process (HadEX2)
    :param str index: which ETCCDI index to process
    :param bool diagnostics: output diagnostic information
    """
    
    # read inventory
    station_list = inventory_utils.read_arabia_6190_index(dataset, index, diagnostics=diagnostics)

    # set up storage for this metadata
    metadata_file = open(os.path.join(dataset.location, "{}_stations.txt".format(index)), "w")

    # spin through each station
    for this_station in station_list:

        # see if a file exists for this index (only some of the precip ones, really).
        try:
            if diagnostics:
                print(this_station.id)

            indata = np.genfromtxt(os.path.join(dataset.location, "raw", index, "{}_{}.txt".format(this_station.id, index)), skip_header=1, missing_values=-99.9, filling_values=utils.HADEX_MDI, dtype=float, encoding="latin-1")


            # extract the two columns
            years = indata[:, 0].astype(int)
            annual = indata[:, -1]

            if index in utils.MONTHLY_INDICES:
                # not all monthly indices have monthly values for this source
                try:
                    # and write into the attributes
                    this_station.years = years
                    this_station.annual = indata[:, -1]

                    monthly = indata[:, 1: -1].reshape(-1)
                    months = np.array([range(1, 13) for i in years]).reshape([-1])
                    myears = np.array([[year for i in range(12)] for year in years]).reshape([-1])

                    if len(monthly) == 0:
                        monthly = np.ones(months.shape[0]) * utils.HADEX_MDI

                    # store values for monthly data
                    this_station.monthly = monthly
                    this_station.myears = myears
                    this_station.months = months 
                except IndexError:               
                    # and write into the attributes
                    this_station.years = years
                    this_station.annual = annual
                    # store empty values for monthly data
                    this_station.monthly = None
                    this_station.myears = None
                    this_station.months = None
            else:
                # and write into the attributes
                this_station.years = years
                this_station.annual = annual
                # store empty values for monthly data
                this_station.monthly = None
                this_station.myears = None
                this_station.months = None

            # prepare to write the data - check all directories present
            utils.check_directories(dataset, this_station)
                
            if index in utils.MONTHLY_INDICES:
                utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", this_station.id, "{}_{}_MON.csv".format(this_station.id, index.lower())), this_station, index, doMonthly=True)
            utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", this_station.id, "{}_{}_ANN.csv".format(this_station.id, index.lower())), this_station, index)

            # write metadata!
            metadata_file.write("{:10s} {:7.3f} {:7.3f}\n".format(this_station.id, this_station.latitude, this_station.longitude))

        except IOError:
            print("file missing {}".format(os.path.join(dataset.location, "raw", index, "{}_{}.txt".format(this_station.id, index))))

    metadata_file.close()
    return # convert_arabia_6190

#*********************************************
def convert_arabia_8110(dataset, index, diagnostics=False):
    """
    Controller routine for Arabian index conversion from Markus Donat (All stations, 1981-2010)
   
    :param datasetObj dataset: dataset to process (HadEX2)
    :param str index: which ETCCDI index to process
    :param bool diagnostics: output diagnostic information
    """
    
    # read inventory
    station_list = inventory_utils.read_arabia_8110_index(dataset, index, diagnostics=diagnostics)

    # set up storage for this metadata
    metadata_file = open(os.path.join(dataset.location, "{}_stations.txt".format(index)), "w")

    # spin through each station
    for this_station in station_list:

        # see if a file exists for this index (only some of the precip ones, really).
        try:
            if diagnostics:
                print(this_station.id)

            if index in ["TX90p", "TX10p", "TN90p", "TN10p"]:
                index_name = index.upper()
            elif index in ["FD", "ID", "SU", "TR", "Rx1day", "Rx5day"]:
                rename = {"FD" : "FD0", "ID" : "ID0", "SU" : "SU25", "TR" : "TR20", "Rx1day" : "RX1day", "Rx5day" : "RX5day"}
                index_name = rename[index]
            else:
                index_name = index

            indata = np.genfromtxt(os.path.join(dataset.location, "raw", index, "{}_{}.csv".format(this_station.id, index_name)), skip_header=1, missing_values=-99.9, filling_values=utils.HADEX_MDI, dtype=float, delimiter=',', encoding="latin-1")

            # extract the two columns
            years = indata[:, 0].astype(int)
            annual = indata[:, -1]

            if index in utils.MONTHLY_INDICES:
                # not all monthly indices have monthly values for this source
                try:
                    # and write into the attributes
                    this_station.years = years
                    this_station.annual = indata[:, -1]

                    monthly = indata[:, 1: -1].reshape(-1)
                    months = np.array([range(1, 13) for i in years]).reshape([-1])
                    myears = np.array([[year for i in range(12)] for year in years]).reshape([-1])

                    if len(monthly) == 0:
                        monthly = np.ones(months.shape[0]) * utils.HADEX_MDI

                    # store values for monthly data
                    this_station.monthly = monthly
                    this_station.myears = myears
                    this_station.months = months 
                except IndexError:               
                    # and write into the attributes
                    this_station.years = years
                    this_station.annual = annual
                    # store empty values for monthly data
                    this_station.monthly = None
                    this_station.myears = None
                    this_station.months = None

            else:
                # and write into the attributes
                this_station.years = years
                this_station.annual = annual
                # store empty values for monthly data
                this_station.monthly = None
                this_station.myears = None
                this_station.months = None

            # prepare to write the data - check all directories present
            utils.check_directories(dataset, this_station)
                
            if index in utils.MONTHLY_INDICES:
                utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", this_station.id, "{}_{}_MON.csv".format(this_station.id, index.lower())), this_station, index, doMonthly=True)
            utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", this_station.id, "{}_{}_ANN.csv".format(this_station.id, index.lower())), this_station, index)

            # write metadata!
            metadata_file.write("{:10s} {:7.3f} {:7.3f}\n".format(this_station.id, this_station.latitude, this_station.longitude))

        except IOError:
            print("file missing {}".format(os.path.join(dataset.location, "raw", index, "{}_{}.csv".format(this_station.id, index))))

    metadata_file.close()
    return # convert_arabia_8110

#*********************************************
def convert_south_africa(dataset, index, diagnostics=False):
    """
    Controller routine for SAWS index conversion from Andries Kruger
   
    :param datasetObj dataset: dataset to process (HadEX2)
    :param str index: which ETCCDI index to process
    :param bool diagnostics: output diagnostic information
    """

    # read inventory
    station_list = inventory_utils.read_south_africa(dataset, diagnostics=diagnostics)

    # set up storage for this metadata
    metadata_file = open(os.path.join(dataset.location, "{}_stations.txt".format(index)), "w")

    # spin through each station
    for this_station in station_list:

        try:
        # see if a file exists for this index (only some of the precip ones, really).
            if diagnostics:
                print(this_station.id)

            if index in ["CSDI", "DTR", "ETR", "FD", "GSL", "ID", "SU", "TN10p", "TN90p", "TNn", "TNx", "TR", "TX10p", "TX90p", "TXn", "TXx", "WSDI"]:
                
                # some indices have a non standard name
                if index in ["FD", "ID", "SU", "TR"]:
                    rename = {"FD" : "FD0", "ID" : "ID0", "SU" : "SU25", "TR" : "TR20"}
                    index_name = rename[index]
                elif index in ["TX90p", "TX10p", "TN90p", "TN10p"]:
                    rename = {"TX90p" : "TX90P", "TX10p" : "TX10P", "TN90p" : "TN90P", "TN10p" : "TN10P"}
                    index_name = rename[index]
                else:
                    index_name = index    

                # try first filename
                for dummy_id in [this_station.id, this_station.id.upper(), this_station.id.lower(), this_station.id.capitalize()]:
                    print(dummy_id)
                    infile = os.path.join(dataset.location, "raw", "{}_{}.csv".format(dummy_id, index_name))
                    if os.path.exists(infile):
                        break
                    # try second
                    infile = os.path.join(dataset.location, "raw", "{} daily temps_{}.csv".format(dummy_id, index_name))
                    if os.path.exists(infile):
                        break
                    # try third
                    infile = os.path.join(dataset.location, "raw", "{} daily temp_{}.csv".format(dummy_id, index_name))
                    if os.path.exists(infile):
                        break

                            

            elif index in ["CDD", "CWD", "PRCPTOT", "R10mm", "R20mm", "R95p", "R95pTOT", "R99p", "R99pTOT", "Rx1day", "Rx5day", "SDII"]:

                # some indices have a non standard name
                if index in ["Rx1day", "Rx5day"]:
                    index_name = "X".join(index.split("x"))
                else:
                    index_name = index
                for dummy_id in [this_station.id, this_station.id.upper(), this_station.id.lower()]:
                    print(dummy_id)
                    infile = os.path.join(dataset.location, "raw", "{}_{}.csv".format(dummy_id, index_name))
                    if os.path.exists(infile):
                        break
                    # try second
                    infile = os.path.join(dataset.location, "raw", "{} daily temps_{}.csv".format(dummy_id, index_name))
                    if os.path.exists(infile):
                        break
                    # try third
                    infile = os.path.join(dataset.location, "raw", "{} daily temp_{}.csv".format(dummy_id, index_name))
                    if os.path.exists(infile):
                        break

#                infile = os.path.join(dataset.location, "raw", "{}_{}.csv".format(this_station.id, index_name))

            else:
                infile = ""
        

            # escape if no infile [ET-SCI indices]
            if infile == "":
                continue

            if diagnostics:
                print(infile)
            #******************
            # and read the data
            indata = np.genfromtxt(infile, delimiter=",", skip_header=1, dtype=(float), encoding="latin-1")

            # extract the columns
            years = indata[:, 0].astype(int)

            if index in utils.MONTHLY_INDICES:
                # and write into the attributes
                this_station.years = years
                this_station.annual = indata[:, -1]

                monthly = indata[:, 1: -1].reshape(-1)
                months = np.array([range(1, 13) for i in years]).reshape([-1])
                myears = np.array([[year for i in range(12)] for year in years]).reshape([-1])

                if len(monthly) == 0:
                    monthly = np.ones(months.shape[0]) * utils.HADEX_MDI

                # store values for monthly data
                this_station.monthly = monthly
                this_station.myears = myears
                this_station.months = months 
            else:
                # and write into the attributes
                this_station.years = years
                this_station.annual = indata[:, 1]
                # store empty values for monthly data
                this_station.monthly = None
                this_station.myears = None
                this_station.months = None

            # and fix the IDs
            this_station.id = "_".join([chars.upper() for chars in this_station.id.split()])

            # prepare to write the data - check all directories present
            utils.check_directories(dataset, this_station)
            if index in utils.MONTHLY_INDICES:
                # monthly data
                utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", this_station.id, "{}_{}_MON.csv".format(this_station.id, index.lower())), this_station, index, doMonthly=True)

            # annual data
            utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", this_station.id, "{}_{}_ANN.csv".format(this_station.id, index.lower())), this_station, index)

            # write metadata!
            metadata_file.write("{:10s} {:7.3f} {:7.3f}\n".format(this_station.id, this_station.latitude, this_station.longitude))
    
        except IOError:
            print("file missing {}".format(os.path.join(dataset.location, "raw", "{}_{}.csv".format(this_station.id, index))))

    metadata_file.close()
    return # convert_south_africa

#*********************************************
def convert_ghcndex(dataset, index, diagnostics=False):
    """
    Controller routine for GHCNDEX index conversion for workshop stations
   
    :param datasetObj dataset: dataset to process (GHCNDEX)
    :param str index: which ETCCDI index to process
    :param bool diagnostics: output diagnostic information
    """

    # read in the list of available stations for this index
    station_list = inventory_utils.read_ghcndex(dataset, index, diagnostics=diagnostics)

    # if no stations, then just exit (inventory_utils.read_hadex2 returns stations as [] if empty)
    if station_list == []:
        return

    # spin through each station (just metadata so far)
    for station in station_list:
        if diagnostics:
            print(station.id)
        
        filename = os.path.join(dataset.location, "raw", index, "{}_{}.txt".format(station.id, index))

        # read in station data           
        if index in utils.MONTHLY_INDICES:
            station, metadata = read_station(dataset, filename, station, index, doMonthly=True)      
        else:
            station, metadata = read_station(dataset, filename, station, index)      

        # prepare to write the data - check all directories present
        utils.check_directories(dataset, station)

        if index in utils.MONTHLY_INDICES:
            utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", station.id, "{}_{}_MON.csv".format(station.id, index.lower())), station, index, doMonthly=True)
        utils.write_station_index(os.path.join(dataset.location, "formatted", "indices", station.id, "{}_{}_ANN.csv".format(station.id, index.lower())), station, index)

    return # convert_ghcndex

#*********************************************
def main(indata="ecad", index="TX90p", unzip=True, convert=True, diagnostics=False):
    """
    Calls the unzipping and conversion scripts.

    :param str indata: input dataset name
    :param str index: index to process
    :param bool unzip: select to turn on the unzipping
    :param bool convert: select to turn on the conversion
    :param bool diagnostics: output diagnostic information
    """

    # get all possible datasets
    all_datasets = utils.get_input_datasets()

    # and their names
    names = np.array([d.name for d in all_datasets])

    # if dataset selected and in the list of available, then run
    if indata in names and indata in ["ecad", "sacad", "lacad"]:
        
        if unzip:
            unzip_ecad_family(all_datasets[names == indata][0], index, diagnostics=diagnostics)

        if convert:
            convert_ecad_family(all_datasets[names == indata][0], index, diagnostics=diagnostics)

    elif indata in names and indata in ["hadex2"]:

        convert_hadex2(all_datasets[names == indata][0], index, diagnostics=diagnostics)
        
    elif indata in names and indata in ["south_america"]:

        convert_south_america(all_datasets[names == indata][0], index, diagnostics=diagnostics)
    
    elif indata in names and indata in ["west_africa_pptn"]:

        convert_west_africa_pptn(all_datasets[names == indata][0], index, diagnostics=diagnostics)

    elif indata in names and indata in ["west_africa_indices"]:

        convert_west_africa_indices(all_datasets[names == indata][0], index, diagnostics=diagnostics)

    elif indata in names and indata in ["arabia"]:
        
        if utils.REF_START == 1961 and utils.REF_END == 1990:
            convert_arabia_6190(all_datasets[names == indata][0], index, diagnostics=diagnostics)
        elif utils.REF_START == 1981 and utils.REF_END == 2010:
            convert_arabia_8110(all_datasets[names == indata][0], index, diagnostics=diagnostics)

    elif indata in names and indata in ["south_africa"]:

        convert_south_africa(all_datasets[names == indata][0], index, diagnostics=diagnostics)

    elif indata in names and indata in ["ghcndex"]:

        convert_ghcndex(all_datasets[names == indata][0], index, diagnostics=diagnostics)
        
    # fail gracefully
    else:
        print("data name not available: {}\n".format(indata))
        print("available data names: {}".format(" ".join(names)))

    return # main



#************************************************************************
#************************************************************************
if __name__ == "__main__":


    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--indata', dest='indata', action='store', default="ecad",
                        help='Which dataset to convert')
    parser.add_argument('--index', dest='index', action='store', default="TX90p",
                        help='Which index to run')
    parser.add_argument('--unzip', dest='unzip', action='store', default=True,
                        help='select unzipping')
    parser.add_argument('--convert', dest='convert', action='store', default=True,
                        help='elect conversion')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default = False')
 
    args = parser.parse_args()         

    if args.indata in ["ecad", "lacad", "sacad"]:

        main(indata=args.indata, index=args.index, unzip=args.unzip, \
                 convert=args.convert, diagnostics=args.diagnostics)

    elif args.indata in ["pacific", "malaysia", "brunei", "china", "myanmar", "vietnam", \
    	"philippines", "singapore", "thailand", "indonesia"]:
        
        print("nothing needed - already in Climpact2 format and in formatted/ directory")

    elif args.indata in ["hadex2", "west_africa_pptn", "west_africa_indices", "south_america", \
                         "arabia", "south_africa", "ghcndex"]:
        
        main(indata=args.indata, index=args.index, diagnostics=args.diagnostics)


    # else just exit again


#*******************************************
# END
#*******************************************
