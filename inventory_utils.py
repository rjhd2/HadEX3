#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 449                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2020-01-15 11:46:14 +0000 (Wed, 15 Jan 2020) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""
Inventory reading scripts to convert to standard metadata format - for each input datast type
"""
import os
import sys
import struct
import re
import numpy as np

# RJHD scripts
import utils

OBSERVABLES = ["TMAX", "TMIN", "PRCP"]

#*********************************************
def check_metadata(station_list, name, lat, lon):

    all_names = np.array([stn.id for stn in station_list])

    loc, = np.where(all_names == name)[0]

    test_match = station_list[loc]

    if np.abs(test_match.latitude - lat) > 0.1:
        print("Error - {} exists in list but metadata doesn't match".format(name))
        print("Latitude {}, {}".format(test_match.latitude, lat))
    if np.abs(test_match.longitude - lon) > 0.1:
        print("Error - {} exists in list but metadata doesn't match".format(name))
        print("Longitude {}, {}".format(test_match.longitude, lon))

    return # check_metadata

#*********************************************
def read_ghcnd(dataset, diagnostics=False):
    """
    Read the two GHCND inventory files.  Allows selection of GSN and HCN sub-networks
    Does some processing on record length to reduce number of stations.

    :param datasetObj dataset: dataset object holding metadata about the input dataset
    :param bool diagnostics: output diagnostic information
    """

    #********************
    def extract_hadex_ids():

        these_ids = []

        for index in utils.INDICES:
        
            hadex_list = "/project/hadobs2/hadex3/hadex2/station_list/HadEX2_station_list_all_{}.csv".format(index)

            try:
                with open(hadex_list, 'r') as infile:
                    for line in infile:
                        if re.search("Tom Peterson", line):
                            this_id = line.split(",")[8] # get the GHCND ID
                            if this_id[0] == '"':
                                this_id = this_id[1:]
                            if this_id[-1] == '"':
                                this_id = this_id[:-1]
                            these_ids += [this_id]

            except IOError:
                # likely don't have HadEX2 list for this index
                continue
        these_ids = np.unique(np.array(these_ids))

        return these_ids # extract_hadex_ids
    #********************

    hadex_ids = extract_hadex_ids()
    if diagnostics and utils.ghcnd_peterson:
        print("HadEX subset of GHCND stations: {}".format(len(hadex_ids)))

    names = []

    # process first metadata file
    fieldwidths = [11, 9, 10, 7, 3, 31, 4, 4, 6]

    # check GHCND stations to enable further selection by GCOS and HCN networks
    with open(os.path.join(dataset.location, "raw", "ghcnd-stations.txt"), "r", encoding="latin-1") as infile:
        
        for line in infile:
            
            fields = utils.parse(line, fieldwidths)

            ID = fields[0]

            if utils.ghcnd_peterson:
                # process HadEX list
                if ID in hadex_ids:
                    names += [ID]

            elif (utils.ghcnd_gsn and fields[6].strip() == "GSN") and (utils.ghcnd_hcn and fields[7].strip() == "HCN"):
                # if both
                names += [ID]

            elif utils.ghcnd_gsn and fields[6].strip() == "GSN":
                # if GSN
                names += [ID]

            elif utils.ghcnd_hcn and fields[7].strip() == "HCN":
                # if HCN
                names += [ID]

            elif not utils.ghcnd_gsn and not utils.ghcnd_hcn:
                # if none
                names += [ID]            

    if diagnostics:
        print("Initial list of GHCND stations: {} (GSN = {}, HCN = {}, Peterson = {})".format(len(names), utils.ghcnd_gsn, utils.ghcnd_hcn, utils.ghcnd_peterson))

    # now check if have parameters and sufficient length of record
    stations = []

    with open(os.path.join(dataset.location, "raw", "ghcnd-inventory.txt"), "r", encoding="latin-1") as infile:
        
        for line in infile:
            
            split = line.split()

            ID = split[0]
            latitude = float(split[1])
            longitude = float(split[2])
            parameter = split[3]
            end = int(split[5])
            length = int(split[5]) - int(split[4]) + 1
            # if preselected, has suitable data, length and end year
            if ID in names and \
                    parameter in OBSERVABLES and \
                    length >= utils.MIN_N_YEARS and \
                    end >= utils.FINISH_AFTER.year:

                if stations == []:
                    stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

                else:
                    if stations[-1].id != ID:
                        stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

        
    if diagnostics:
        print("Final list of {} stations: {}".format(dataset.name, len(stations)))

    return stations # read_ghcnd

#*********************************************
def read_spain(dataset, diagnostics=False):
    """
    Read the inventory files for Spain from Manola Brunet.

    :param datasetObj dataset: dataset object holding metadata about the input dataset
    :param bool diagnostics: output diagnostic information
    """

    stations = []

    try:
        with open(os.path.join(dataset.location, "raw", "estacion.txt"), "r", encoding="latin-1") as infile:

            for line in infile:

                split = line.split()

                if split == [] or split[0] == "Brunet,": break

                ID = split[0][:-4]
                latitude_dms = [i for i in split[1:4]]
                longitude_dms = [i for i in split[4:8]]

                latitude = utils.dms2dd(latitude_dms[0], latitude_dms[1], latitude_dms[2], direction="N")
                longitude = utils.dms2dd(longitude_dms[0], longitude_dms[1], longitude_dms[2], direction=longitude_dms[3])

                if stations == []:
                    stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

                else:
                    if stations[-1].id != ID:
                        stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

        if diagnostics:
            print("Final list of {} stations: {}".format(dataset.name, len(stations)))

    except IOError:
        print("No file: {}".format(os.path.join(dataset.location, "raw", "estacion.txt")))
        
    return stations # read_spain

#*********************************************
def read_russia(dataset, diagnostics=False):
    """
    Read the inventory files for Russia from Xuebin Zhang.

    :param datasetObj dataset: dataset object holding metadata about the input dataset
    :param bool diagnostics: output diagnostic information
    """

    stations = []

    try:
        with open(os.path.join(dataset.location, "raw", "stn_info.dat"), "r", encoding="latin-1") as infile:

            for line in infile:

                split = line.split()

                if split == []: break

                ID = split[0]
                longitude = float(split[2])
                latitude = float(split[1])

                if stations == []:
                    stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

                else:
                    if stations[-1].id != ID:
                        stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

        if diagnostics:
            print("Final list of {} stations: {}".format(dataset.name, len(stations)))

    except IOError:
        print("No file: {}".format(os.path.join(dataset.location, "raw", "stn_info.dat")))
        
    return stations # read_russia

#*********************************************
def read_acre(dataset, diagnostics=False):
    """
    Read the inventory files for ACRE from Rob Allan, reformatted at UNSW.

    :param datasetObj dataset: dataset object holding metadata about the input dataset
    :param bool diagnostics: output diagnostic information
    """

    stations = []

    try:
        with open(os.path.join(dataset.location, "raw", "acre_meta.csv"), "r", encoding="latin-1") as infile:

            for line in infile:

                split = line.split(",")

                if split == []: break

                ID = split[0].replace(" ", "").lower() # remove spaces and convert to lowercase
                longitude = float(split[2])
                latitude = float(split[1])

                if stations == []:
                    stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

                else:
                    if stations[-1].id != ID:
                        stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

        if diagnostics:
            print("Final list of {} stations: {}".format(dataset.name, len(stations)))

    except IOError:
        print("No file: {}".format(os.path.join(dataset.location, "raw", "acre_meta.csv")))
        
    return stations # read_acre

#*********************************************
def read_nz(dataset, diagnostics=False):
    """
    Read the inventory files for NZ from Jim Salinger.

    :param datasetObj dataset: dataset object holding metadata about the input dataset
    :param bool diagnostics: output diagnostic information
    """

    stations = []

    try:
        with open(os.path.join(dataset.location, "raw", "NZ_list_orig.txt"), "r", encoding="latin-1") as infile:

            for line in infile:

                split = line.split()

                if split == [] or split[0] == "Musselburgh,": break

                ID = split[0]
                latitude_dms = [i for i in split[1:4]]
                longitude_dms = [i for i in split[4:7]]

                latitude = utils.dms2dd(latitude_dms[0], latitude_dms[1], 0, direction=latitude_dms[2])
                longitude = utils.dms2dd(longitude_dms[0], longitude_dms[1], 0, direction=longitude_dms[2])

                if stations == []:
                    stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

                else:
                    if stations[-1].id != ID:
                        stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

        if diagnostics:
            print("Final list of {} stations: {}".format(dataset.name, len(stations)))

    except IOError:
        print("No file: {}".format(os.path.join(dataset.location, "raw", "NZ_list_orig.txt")))
        
    return stations # read_nz

#*********************************************
def read_honduras(dataset, diagnostics=False):
    """
    Read the inventory files for Honduras from Ernesto Salgado.

    :param datasetObj dataset: dataset object holding metadata about the input dataset
    :param bool diagnostics: output diagnostic information
    """
    import glob
    stations = []

    for filename in glob.glob(r'{}/Temperatura Minima*.csv'.format(os.path.join(dataset.location, "raw"))):

        try:
            with open(filename, "r", encoding="latin-1") as infile:
                for l, line in enumerate(infile):
                    if l == 2:
                        # process the single line of interest
                        line = line.split(",")

                        name = line[0]
                        # ID = line[7].split(" ")[-1] # if use WMO id
                        ID = filename.split("/")[-1].split("Abs. ")[-1][:-4] # if use filename - which is much easier
                        lat = line[10].split()
                        lon = line[12].split()
                        elev = line[15]

                        # convert the coordinates to decimal
                        latitude = utils.dms2dd(lat[0][:2], lat[1][:2], lat[2][:2], direction=lat[3])
                        longitude = utils.dms2dd(lon[0][:2], lon[1][:2], lon[2][:2], direction=lon[3])
                        
                        stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]
        except IndexError:
            print("Metadata missing in {}".format(filename))
        except IOError:
            print("No file: {}".format(filename))

    if diagnostics:
        print("Final list of {} stations: {}".format(dataset.name, len(stations)))

    return stations # read_honduras

#*********************************************
def read_ecad(dataset, index, diagnostics=False):
    """
    Read the ECAD family inventory files for the specified Index

    :param datasetObj dataset: dataset to process (ECAD/LACAD/SACAD)
    :param bool diagnostics: output diagnostic information
    :param str index: which ETCCDI index to process
    """

    stations = []

    # set up parseing info
    fieldwidths = [6, 41, 41, 10, 11, 4]
    
    try:
        with open(os.path.join(dataset.location, "{}_stations.txt".format(index)), "r", encoding="latin-1") as infile:

             for lc, line in enumerate(infile):

                if lc < 18 -1: continue
                if len(line.strip()) != sum(fieldwidths): continue #  just skip if malformed

                fields = utils.parse(line.strip(), fieldwidths)

                # need to remove trailing comma
                ID = "{:06d}".format(int(fields[0][:-1]))
                country = fields[2][:-1].strip()
                latitude_dms = fields[3][:-1].split(":")
                longitude_dms = fields[4][:-1].split(":")

                lat_dir = "N"
                if latitude_dms[0][0] == "-":
                    lat_dir = "S"
                lon_dir = "E"
                if longitude_dms[0][0] == "-":
                    lon_dir = "W"

                latitude = utils.dms2dd(latitude_dms[0], latitude_dms[1], latitude_dms[2], direction=lat_dir)
                longitude = utils.dms2dd(longitude_dms[0], longitude_dms[1], longitude_dms[2], direction=lon_dir)

                if dataset.name == "sacad" and country == "AUSTRALIA":
                    # now have data from BOM (Blair) so don't include these.
                    continue

                if stations == []:
                    stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw", index), dataset.name)]

                else:
                    if stations[-1].id != ID:
                        stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw", index), dataset.name)]

        if diagnostics:
            print("Final list of {} stations for {}: {}".format(index, dataset.name.upper(), len(stations)))

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "{}_stations.txt".format(index))))

    return stations # read_ecad

#*********************************************
def read_hadex2(dataset, index, diagnostics=False):
    """
    Read the HadEX2 inventory files for the specified Index

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    :param str index: which ETCCDI index to process
    """

    stations = []

    try:
        with open(os.path.join(dataset.location, "HadEX2_station_list_{}.csv".format(index)), "r", encoding="latin-1") as infile:

            for lc, line in enumerate(infile):

                sline = line.split(",")

                ID = sline[0][1:-1]
                latitude = float(sline[1])
                longitude = float(sline[2])

                stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw", index), dataset.name)]
        
        if diagnostics:
            print("Final list of {} stations for {}: {}".format(index, dataset.name.upper(), len(stations)))

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "HadEX2_station_list_{}.csv".format(index))))

    return stations # read_hadex2

#*********************************************
def read_generic_index(dataset, index, diagnostics=False):
    """
    Read the inventory files for the specified Index [v20180813]

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    :param str index: which ETCCDI index to process
    """

    stations = []

    try:
        with open(os.path.join(dataset.location, "{}_stations.txt".format(index)), "r", encoding="latin-1") as infile:

            for lc, line in enumerate(infile):

                sline = line.split()

                ID = sline[0]
                latitude = float(sline[1])
                longitude = float(sline[2])

                stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw", index), dataset.name)]
        
        if diagnostics:
            print("Final list of {} stations for {}: {}".format(index, dataset.name.upper(), len(stations)))

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "{}_stations.txt".format(index))))

    return stations # read_generic_index

#*********************************************
def read_generic_obs(dataset, diagnostics=False):
    """
    Read the South America & Brazil_SP inventory files

    :param datasetObj dataset: dataset to process 
    :param bool diagnostics: output diagnostic information
    """

    stations = []

    try:
        with open(os.path.join(dataset.location, "{}_stations.txt".format(dataset.name)), "r", encoding="latin-1") as infile:

            for lc, line in enumerate(infile):

                sline = line.split()

                ID = sline[0]
                latitude = float(sline[1])
                longitude = float(sline[2])

                stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]
        
        if diagnostics:
            print("Final list of stations for {}: {}".format(dataset.name.upper(), len(stations)))

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "{}_stations.txt".format(dataset.name))))

    return stations # read_generic_obs

#*********************************************
def read_australia(dataset, diagnostics=False):
    """
    Read the Australian inventory

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []

    fieldwidths = [7, 41, 9, 9, 8]

    try:
        with open(os.path.join(dataset.location, "raw", "Station list.txt"), "r", encoding="latin-1") as infile:

            for lc, line in enumerate(infile):

                fields = utils.parse(line, fieldwidths)

                ID = fields[0].strip()
                latitude = float(fields[2])
                longitude = float(fields[3])
                
                # check if ID exists in the subset of stations sent
                if os.path.exists(os.path.join(dataset.location, "raw", "hq17dtr{}".format(ID))):
                    stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]
                elif os.path.exists(os.path.join(dataset.location, "raw", "DC02D_Data_{}_46309839556073.txt".format(ID))):
                    stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]
        
        if diagnostics:
            print("Final list of stations for {}: {}".format(dataset.name.upper(), len(stations)))

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "raw", "Station list.txt")))

    return stations # read_australia

#*********************************************
def read_west_africa_pptn(dataset, diagnostics=False):
    """
    Read the West Africa inventory file 

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []

    try:
        with open(os.path.join(dataset.location, "raw", "stas_coords.csv"), "r", encoding="latin-1") as infile:

            for lc, line in enumerate(infile):

                if lc == 0:
                    # skip first line
                    continue

                sline = line.split(",")

                ID = sline[0]
                latitude = float(sline[1])
                longitude = float(sline[2])

                if "filled" in ID:
                    # skip these as not sure what's been done
                    continue
                stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]
        
        if diagnostics:
            print("Final list of stations for {}: {}".format(dataset.name.upper(), len(stations)))

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "raw", "stas_coords.csv")))

    return stations # read_west_africa_pptn

#*********************************************
def read_arabia_6190_index(dataset, index, diagnostics=False):
    """
    Read the Arabian inventory files for the specified index (for long records subset, re 1961-90)

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    :param str index: which ETCCDI index to process
    """

    stations = []

    try:
        metadata = np.genfromtxt(os.path.join(dataset.location, "raw", index, "HadEX2_station_list_Arab-Region_{}.csv".format(index)), delimiter=",", skip_header=1, dtype=str)

        IDs = metadata[:, 0]
        lats = metadata[:, 1].astype(float)
        lons = metadata[:, 2].astype(float)
        
        for s, stnid in enumerate(IDs):
            stations += [utils.Station(stnid, lats[s], lons[s], os.path.join(dataset.location, "raw", index), dataset.name)]
        
        if diagnostics:
            print("Final list of {} stations for {}: {}".format(index, dataset.name.upper(), len(stations)))

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "raw", index, "HadEX2_station_list_Arab-Region_{}.csv".format(index))))

    return stations # read_arabia_6190_index

#*********************************************
def read_arabia_8110_index(dataset, index, diagnostics=False):
    """
    Read the Arabian inventory files for the specified index (for all records, ref 1981-2010)

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    :param str index: which ETCCDI index to process
    """

    stations = []

    try:
        metadata = np.genfromtxt(os.path.join(dataset.location, "raw", "station-List.csv"), delimiter=",", dtype=str)

        IDs = metadata[:, 0]
        lats = metadata[:, 3].astype(float)
        lons = metadata[:, 4].astype(float)
        
        for s, stnid in enumerate(IDs):

            if stnid[:7] == "Morocco":
                final_id = stnid[:-8]
            else:
                final_id = stnid[:-4]

            stations += [utils.Station(final_id, lats[s], lons[s], os.path.join(dataset.location, "raw", index), dataset.name)]
        
        if diagnostics:
            print("Final list of {} stations for {}: {}".format(index, dataset.name.upper(), len(stations)))

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "raw", "station-List.csv")))

    return stations # read_arabia_8110_index

#*********************************************
def read_chile(dataset, diagnostics=False):
    """
    Read the Chilean inventory file

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []

    try:
        metadata = np.genfromtxt(os.path.join(dataset.location, "raw", "Climdex1961-2018_est.csv"), delimiter=",", skip_header=1, dtype=str, encoding="latin-1")
    
        lons = metadata[:, 0].astype(float)
        lats = metadata[:, 1].astype(float)
        ids = metadata[:, 2]
        elev = metadata[:, 3]
        names = metadata[:, 4]
        
        for n, name in enumerate(names):
            if name == "RAPA NUI":
                name = "IPASCUA"
                

            stations += [utils.Station(name, lats[n], lons[n], os.path.join(dataset.location, "raw"), dataset.name)]

        if diagnostics:
            print("Final list of stations for {}: {}".format(dataset.name.upper(), len(stations)))

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "raw", "Climdex1961-2018_est.csv")))

    return stations # read_chile

#*********************************************
def read_colombia(dataset, diagnostics=False):
    """
    Read the Colombian inventory file

    Two inventory files - one for stations provided using IDs, one for stations provided with names
    There is some overlap between the two.

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []
    try:
        metadata = np.genfromtxt(os.path.join(dataset.location, "raw", "COLOMBIA DailyPrecip Data_inventory.csv"), delimiter=",", skip_header=1, dtype=str)

        lons = metadata[:, 0].astype(float)
        lats = metadata[:, 1].astype(float)
        elev = metadata[:, 2]
        ids = metadata[:, 3]
        names = np.array([m.strip() for m in metadata[:, 4]])
        
        for s, sid in enumerate(ids):
            stations += [utils.Station(sid, lats[s], lons[s], os.path.join(dataset.location, "raw"), dataset.name)]

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "raw", "COLOMBIA DailyPrecip Data_inventory.csv")))


    try:

        with open(os.path.join(dataset.location, "raw", "COH2L_inventory.txt"), "r", encoding="latin-1") as infile:
            # mix of space and tab delimited, just split process, and join names.
            for line in infile:

                fields = line.split()
                
                ID = fields[0].strip()
                latitude = float(fields[1])
                longitude = float(fields[2])
                name = " ".join([f.strip() for f in fields[5:]])

                stations += [utils.Station(name, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "raw", "COH2L_inventory.txt")))

    if diagnostics:
        print("Final list of stations for {}: {}".format(dataset.name.upper(), len(stations)))

    return stations # read_colombia


#*********************************************
def read_canada(dataset, diagnostics=False):
    """
    Read the Canadian inventory file

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []
    all_ids = []

    for filename in ["Adj_Precipitation_Stations.csv", "Homog_Temperature_Stations.csv"]:
        try:
            metadata = np.genfromtxt(os.path.join(dataset.location, "raw", filename), delimiter=",", skip_header=4, dtype=str)

            lats = metadata[:, 7].astype(float)
            lons = metadata[:, 8].astype(float)
            ids = metadata[:, 2]
            elev = metadata[:, 9]
            names = metadata[:, 1]

            ids = np.array([i.strip() for i in ids])
            names = np.array([n.strip() for n in names])

            for s, stnid in enumerate(ids):
                if stnid not in all_ids:
                    stations += [utils.Station(stnid, lats[s], lons[s], os.path.join(dataset.location, "raw"), dataset.name)]
                    all_ids += [stnid]
                else:
                    print("ID {} already present in list".format(stnid))
                    check_metadata(stations, stnid, lats[s], lons[s])

            if diagnostics:
                print("Final list of stations for {}: {}".format(dataset.name.upper(), len(stations)))

        except IOError:
            print("No station list: {}".format(os.path.join(dataset.location, "raw", "")))

    return stations # read_canada

#*********************************************
def read_decade(dataset, diagnostics=False):
    """
    Read the DECADE files and make an inventory

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []
    all_names = []

    for parameter in ["TX", "PP"]:

        infilenames = os.listdir(os.path.join(dataset.location, "raw", parameter))

        for infilename in infilenames:
            
            with open(os.path.join(dataset.location, "raw", parameter, infilename), "r", encoding="latin-1") as infile:
                
                for lc, line in enumerate(infile):
                    line = line.strip()
                    if lc == 0:
                        name = line.split(":")[1].strip()
                    if lc == 4:
                        lat = float(line.split(":")[1])
                    if lc == 5:
                        lon = float(line.split(":")[1])

                    if line[:2] == "--":
                        break

                if name not in all_names:
                    stations += [utils.Station(name, lat, lon, os.path.join(dataset.location, "raw"), dataset.name)]
                    all_names += [name]
                else:
                    print("Name {} already present in list".format(name))
                    check_metadata(stations, name, lat, lon)
 
    return stations # read_decade

#*********************************************
def read_south_africa(dataset, diagnostics=False):
    """
    Read the SAWS files and make an inventory

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []
    all_names = []

    for parameter in ["temperatures", "precipitation"]:

        infilename = os.path.join(dataset.location, "raw", "south_africa_{}.txt".format(parameter))

        with open(os.path.join(dataset.location, "raw", parameter, infilename), "r", encoding="latin-1") as infile:

            for lc, line in enumerate(infile):
                line = line.split()
                
                lat = float(line[-2])
                lon = float(line[-1])

                if parameter == "temperatures":
                    name = " ".join(line[:-2])
                elif parameter == "precipitation":
                    if len(line[0]) == 7:
                        name = " ".join(line[2:-2])
                    else:
                        name = " ".join(line[1:-2])
                else:
                    print("unexpected")
                    sys.exit(1)
                        
                if name not in all_names:
                    stations += [utils.Station(name, lat, lon, os.path.join(dataset.location, "raw"), dataset.name)]
                    all_names += [name]
                else:
                    print("Name {} already present in list".format(name))
                    check_metadata(stations, name, lat, lon)

    return stations # read_south_africa

#*********************************************
def read_west_africa_indices(dataset, diagnostics=False):
    """
    Read the files from Aziz Barry and make an inventory

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []


    with open(os.path.join(dataset.location, "raw", "country_stations.csv"), "r", encoding="latin-1") as infile:

        for lc, line in enumerate(infile):
            if lc < 2:
                # skip header
                continue
                
            line = line.split(",")

            if line[0] == "":
                # entry for station
                name = line[1]
                longitude = float(line[2])
                latitude = float(line[3])
            
                stations += [utils.Station(name, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

    return stations # read_west_africa_indices

#*********************************************
def read_eobs_sources(dataset, parameter, diagnostics=False):
    """
    Read the EOBS-style station list

    :param datasetObj dataset: dataset to process (EOBS/LAOBS/SAOBS)
    :param str parameter: parameter to read (TX/TN/RR)
    :param bool diagnostics: output diagnostic information
    """

    stations = []

    # set up parseing info
    fieldwidths = [7, 41, 3, 10, 11, 5, 5, 9, 9, 6, 46]
    

    try:
        with open(os.path.join(dataset.location, "raw", dataset.name.upper(), "sources_{}.txt".format(parameter)), "r") as infile:

            for lc, line in enumerate(infile):

                if lc < 25 -1: continue

                fields = utils.parse(line, fieldwidths)

                # need to remove trailing comma
                ID = "{:06d}".format(int(fields[0][:-1]))
                name = fields[1][:-1].strip()
                country = fields[2][:-1].strip()
                latitude_dms = fields[3][:-1].split(":")
                longitude_dms = fields[4][:-1].split(":")
                elevation = fields[5][:-1]
                start = fields[7][:-1]
                end = fields[8][:-1]

                if dataset.name == "saobs" and country == "AU": # Australia
                    # now have data from BOM (Blair) so don't include these.
                    continue

                if name != "":

                    lat_dir = "N"
                    if latitude_dms[0][0] == "-":
                        lat_dir = "S"
                    lon_dir = "E"
                    if longitude_dms[0][0] == "-":
                        lon_dir = "W"

                    latitude = utils.dms2dd(latitude_dms[0], latitude_dms[1], latitude_dms[2], direction=lat_dir)
                    longitude = utils.dms2dd(longitude_dms[0], longitude_dms[1], longitude_dms[2], direction=lon_dir)

                    if stations == []:
                        station = utils.Station(ID, latitude, longitude, \
                                                    os.path.join(dataset.location, "raw", dataset.name.upper()), dataset.name)
                        station.elevation = int(elevation)
                        station.country = country
                        station.name = name
                        station.start = start
                        station.end = end

                        stations += [station]

                    else:
                        if stations[-1].id != ID:
                            station = utils.Station(ID, latitude, longitude, \
                                                        os.path.join(dataset.location, "raw", dataset.name.upper()), dataset.name)
                            station.elevation = int(elevation)
                            station.country = country
                            station.name = name
                            station.start = start
                            station.end = end

                            stations += [station]

        if diagnostics:
            print("List of stations for {} {}: {}".format(dataset.name.upper(), parameter, len(stations)))

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "raw", dataset.name.upper(), "sources_{}.txt".format(parameter))))

    return stations # read_eobs_sources

#*********************************************
def read_eobs(dataset, diagnostics=False):
    """
    Read the EOBS-style station lists for each parameter and merge IDs from each parameter together
    Need exact metadata matching

    :param datasetObj dataset: dataset to process (EOBS/LAOBS/SAOBS)
    :param bool diagnostics: output diagnostic information
    """

    stations_tx = np.array(read_eobs_sources(dataset, "tx", diagnostics=diagnostics))
    stations_tn = np.array(read_eobs_sources(dataset, "tn", diagnostics=diagnostics))
    stations_rr = np.array(read_eobs_sources(dataset, "rr", diagnostics=diagnostics))

    used_tn = np.zeros(len(stations_tn))
    used_rr = np.zeros(len(stations_rr))

    names_tx = np.array([s.name for s in stations_tx])
    names_tn = np.array([s.name for s in stations_tn])
    names_rr = np.array([s.name for s in stations_rr])

    final_stations = []
    lookup = {}

    # merge TN into TX
    for t, tx_stn in enumerate(stations_tx):
        
        locs, = np.where(names_tn == names_tx[t])

        if len(locs) == 0:
            # no TN matches TX
            lookup["{}".format(tx_stn.id)] = {"tx" : tx_stn.id, "tn" : -99, "rr" : -1}
            if diagnostics:
                print("tx {} no match in tn".format(tx_stn.id))

        elif len(locs) == 1:
            # single match - easy
            lookup["{}".format(tx_stn.id)] = {"tx" : tx_stn.id, "tn" : stations_tn[locs[0]].id, "rr" : -1}
            used_tn[locs[0]] = 1
            if diagnostics:
                print("tx {} single match in tn {}".format(tx_stn.id, stations_tn[locs[0]].id))

        else:
            # multiple matches - check metadata
            if diagnostics:
                print("tx {} multiple matches in tn".format(tx_stn.id))
                
            matched = False
            for l in locs:
                
                if tx_stn.latitude == stations_tn[l].latitude and \
                   tx_stn.longitude == stations_tn[l].longitude and \
                   tx_stn.elevation == stations_tn[l].elevation and \
                   tx_stn.start == stations_tn[l].start and \
                   tx_stn.end == stations_tn[l].end:
                    # all metadata matches
                    lookup["{}".format(tx_stn.id)] = {"tx" : tx_stn.id, "tn" : stations_tn[l].id, "rr" : -1}
                    used_tn[l] = 1
                    matched = True
                    if diagnostics:
                        print("tx {} single match in tn {}".format(tx_stn.id, stations_tn[l].id))
                    
            if not matched:
                lookup["{}".format(tx_stn.id)] = {"tx" : tx_stn.id, "tn" : -99, "rr" : -1}

        final_stations += [tx_stn]
                 
       
    # back fill unused TNs
    unused_tn = np.where(used_tn == 0)
    for stn in stations_tn[unused_tn]:
        final_stations += [stn]
        lookup["{}".format(stn.id)] = {"tx" : -99, "tn" : stn.id, "rr" : -1}
        if diagnostics:
            print("tn {} no match in tx".format(stn.id))
        
        
    # merge RR into Ts   
    names_t = np.array([s.name for s in final_stations])
   
    for t, t_stn in enumerate(final_stations):
        
        locs, = np.where(names_rr == names_t[t])

        if len(locs) == 0:
            # no TN matches TX
            lookup["{}".format(t_stn.id)]["rr"] = -99
            if diagnostics:
                print("t {} no match in rr".format(t_stn.id))

        elif len(locs) == 1:
            # single match - easy
            lookup["{}".format(t_stn.id)]["rr"] = stations_rr[locs[0]].id
            used_rr[locs[0]] = 1
            if diagnostics:
                print("t {} single match in rr {}".format(t_stn.id, stations_rr[locs[0]].id))

        else:
            # multiple matches - check metadata
            if diagnostics:
                print("t {} multiple matches in rr".format(t_stn.id))
                
            matched = False
            for l in locs:
                
                if t_stn.latitude == stations_rr[l].latitude and \
                   t_stn.longitude == stations_rr[l].longitude and \
                   t_stn.elevation == stations_rr[l].elevation and \
                   t_stn.start == stations_rr[l].start and \
                   t_stn.end == stations_rr[l].end:
                    # all metadata matches
                    lookup["{}".format(t_stn.id)]["rr"] = stations_rr[l].id
                    used_rr[l] = 1
                    matched = True
                    if diagnostics:
                        print("T {} single match in rr {}".format(t_stn.id, stations_rr[l].id))
                    
            if not matched:
                lookup["{}".format(t_stn.id)]["rr"] = -99

        # all stations already in final list - so don't need to add               
       
    # back fill unused RRs
    unused_rr = np.where(used_rr == 0)
    for stn in stations_rr[unused_rr]:
        final_stations += [stn]
        lookup["{}".format(stn.id)] = {"tx" : -99, "tn" : -99, "rr" : stn.id}
        if diagnostics:
            print("rr {} no match in t".format(stn.id))

    return final_stations, lookup # read_eobs

#*********************************************
def read_india(dataset, diagnostics=False):
    """
    Read the Indian files and make an inventory

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []
    all_names = []

    infilename = os.path.join(dataset.location, "raw", "India_station_data_Hadset_2.csv")

    with open(infilename, "r", encoding="latin-1") as infile:

        for lc, line in enumerate(infile):
            if lc < 3:
                # skip header
                continue

            line = line.split(",")

            name = line[3]
            lat = float(line[5])
            lon = float(line[6])

            if name not in all_names:
                stations += [utils.Station(name, lat, lon, os.path.join(dataset.location, "raw"), dataset.name)]
                all_names += [name]
            else:
                print("Name {} already present in list".format(name))
                check_metadata(stations, name, lat, lon)

    infilename = os.path.join(dataset.location, "raw", "LAT_LON_96INDIANSTATION.csv")

    with open(infilename, "r", encoding="latin-1") as infile:

        for lc, line in enumerate(infile):
            if lc < 1:
                # skip header
                continue

            line = line.split(",")

            name = line[1]
            if name == "#N/A":
                name = line[2]
            lat = float(line[4])
            lon = float(line[5])

            if name not in all_names:
                stations += [utils.Station(name, lat, lon, os.path.join(dataset.location, "raw"), dataset.name)]
                all_names += [name]
            else:
                print("Name {} already present in list".format(name))
                check_metadata(stations, name, lat, lon)


    return stations # read_india

#*********************************************
def read_japan(dataset, diagnostics=False):
    """
    Read the Japanese file and make an inventory

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []
    all_names = []

    infilename = os.path.join(dataset.location, "raw", "Readme_Japan_station_data.txt")

    try:
        with open(infilename, "r", encoding="latin-1") as infile:

            for lc, line in enumerate(infile):

                if lc < 12:
                    continue
                elif lc >= 22:
                    continue

                sline = line.split()

                name = sline[0]

                latdd = float(sline[2][:-3])
                latms = float(sline[3].split("'")[0])
                londd = float(sline[4][:-3])
                lonms = float(sline[5].split("'")[0])

                latitude = utils.dms2dd(latdd, latms, 0, direction="N")
                longitude = utils.dms2dd(londd, lonms, 0, direction="E")

                stations += [utils.Station(name, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

            if diagnostics:
                print("Final list of {} stations: {}".format(dataset.name, len(stations)))

    except IOError:
        print("No file: {}".format(os.path.join(dataset.location, "raw", "Readme_Japan_station_data.txt")))
        
    return stations # read_japan

#*********************************************
def read_mexico(dataset, diagnostics=False):
    """
    Read the Mexican file and make an inventory

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []
    all_names = []

    infilename = os.path.join(dataset.location, "raw", "Catalogo_completo.csv")

    try:
        with open(infilename, "r", encoding="latin-1") as infile:

            for lc, line in enumerate(infile):

                if lc < 1:
                    continue

                sline = line.split(",")

                name = "s{}rcl".format(sline[0][-5:])

                latdd = float(sline[13])
                latmm = float(sline[14])
                latss = float(sline[15])
                londd = float(sline[16])
                lonmm = float(sline[17])
                lonss = float(sline[18])

                latitude = utils.dms2dd(latdd, latmm, latss, direction="N")
                longitude = utils.dms2dd(londd, lonmm, latss, direction="W")

                stations += [utils.Station(name, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

            if diagnostics:
                print("Final list of {} stations: {}".format(dataset.name, len(stations)))

    except IOError:
        print("No file: {}".format(os.path.join(dataset.location, "raw", "Catalogo_completo.csv")))
        
    return stations # read_mexico

#*********************************************
def read_iran(dataset, diagnostics=False):
    """
    Read the Iranian file and make an inventory

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []

    try:
        metadata = np.genfromtxt(os.path.join(dataset.location, "raw", "list.csv"), delimiter=",", skip_header=1, dtype=str, encoding="latin-1")

        for info in metadata:

            name = info[0]
            longitude = info[2].astype(float)
            latitude = info[3].astype(float)

            stations += [utils.Station("1{}".format(name), latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]


    except IOError:
        print("No file: {}".format(os.path.join(dataset.location, "raw", "list.csv")))
    
    return stations # read_iran

#*********************************************
def read_brazil(dataset, diagnostics=False):
    """
    Read the Brazilian file and make an inventory

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """

    stations = []

    try:
        metadata = np.genfromtxt(os.path.join(dataset.location, "brazil_inventory.txt"), skip_header=1, dtype=(str), encoding="latin-1")

        for info in metadata:

            name = info[0].split(".")[0]
            latitude = info[1].astype(float)
            longitude = info[2].astype(float)
           
            stations += [utils.Station(name, latitude, longitude, os.path.join(dataset.location, "raw"), dataset.name)]

    except IOError:
        print("No file: {}".format(os.path.join(dataset.location, "brazil_inventory.txt")))
 
    return stations # read_brazil

#*********************************************
def read_ghcndex(dataset, index, diagnostics=False):
    """
    Read the GHCNDEX inventory files for the specified Index

    :param datasetObj dataset: dataset to process (GHCNDEX)
    :param bool diagnostics: output diagnostic information
    :param str index: which ETCCDI index to process
    """

    stations = []

    try:
        with open(os.path.join(dataset.location, "GHCNDEX_station_list_{}.csv".format(index)), "r", encoding="latin-1") as infile:

            for lc, line in enumerate(infile):

                sline = line.split()

                ID = sline[0]
                latitude = float(sline[1])
                longitude = float(sline[2])

                stations += [utils.Station(ID, latitude, longitude, os.path.join(dataset.location, "raw", index), dataset.name)]
        
        if diagnostics:
            print("Final list of {} stations for {}: {}".format(index, dataset.name.upper(), len(stations)))

    except IOError:
        print("No station list: {}".format(os.path.join(dataset.location, "GHCNDEX_station_list_{}.csv".format(index))))

    return stations # read_ghcndex
#*******************************************
# END
#*******************************************
