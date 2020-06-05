#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 462                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2020-02-04 09:46:10 +0000 (Tue, 04 Feb 2020) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
'''
Checks which stations actually have useful data, for each index.

select_final_stations.py invoked by typing::

  python select_final_stations.py --index "TX90p" --diagnostics

Input arguments:

--index         Which ETCCDI index to use

--diagnostics   Output extra info (default False)
'''
import os
import datetime as dt
import numpy as np

# RJHD
import utils

PRECIP_INDICES = ["PRCPTOT", "Rx1day", "Rx5day", "R10mm", "R20mm", "CDD", "CWD", "R95p", "R95pTOT", "R99p", "R99pTOT"]

#*********************************************
def assess_station(station, index, timescale, diagnostics=False):
    '''
    Quickly read the index files to determine whether (a) they exist and (b) they contain suitable data.

    Return True if both

    :param stationObj station: station object (for metadata)
    :param str index: which index to process
    :param str timescale: working on monthly or annual index
    :param bool diagnostics: output diagnostic information
    '''
    # check if it is mentioned in a blacklist
    blacklisted = utils.read_blacklist(diagnostics=diagnostics)

    for bstn in blacklisted:
        if station.id == bstn.id and station.source == bstn.source:
            if station.latitude == bstn.latitude and station.longitude == bstn.longitude:
                return False
            else:
                input("Unexpected issue - blacklisted and matched station have different metadata")

    # presume data
    has_data = True

    try:
        # Check if file is malformed in some way
        with open(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, index.lower(), timescale)), "r", encoding="latin-1") as infile:
            first_line = infile.readline()

            if first_line.split(",")[0] != '"Description: "':
                # malformed file
                if diagnostics:
                    print("station {} for index {} has malformed file".format(station.id, index))
                has_data = False
                return has_data # read_index_files
    except IOError:
        if diagnostics:
            print("File for station {} index {} doesn't exist".format(station.id, index))
        # this index for this station doesn't exist - will return False
        has_data = False

    try:
        # find the raw index file and read in the data
        with open(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, index.lower(), timescale)), "r", encoding="latin-1") as infile:
            indata = np.genfromtxt(infile, skip_header=7, delimiter=",", dtype=(str))
 
    except IOError:
        if diagnostics:
            print("File for station {} index {} doesn't exist".format(station.id, index))
        # this index for this station doesn't exist - will return False
        has_data = False
    except ValueError:
        # likely malformed file
        has_data = False

    #*********************************************
    # if successfully read the file, then indata is set
    if has_data:

        # if there is at least one line of data
        if indata.shape[0] != 0: 

            if len(indata.shape) == 1:
                # if only one line of data, match array size
                indata = np.ma.expand_dims(indata, axis=0) 

            # convert times to datetimes
            if timescale == "ANN":
                try:
                    times = np.array([dt.date(int(time), 1, 1) for time in indata[:, 0]])
                except ValueError:
                    # malformed date:
                    has_data = False
                    if diagnostics:
                        print("File for station {} index {} has malformed date".format(station.id, index))

            elif timescale == "MON":
                try:
                    times = np.array([dt.date(int(time[:4]), int(time[5:7]), 1) for time in indata[:, 0]])
                except ValueError:
                    # malformed date:
                    has_data = False
                    if diagnostics:
                        print("File for station {} index {} has malformed date".format(station.id, index))

        else:
            has_data = False
            if diagnostics:
                print("File for station {} index {} has no data".format(station.id, index))

    #*********************************************
    # if sensible time formats, then continue
    if has_data:

        if "heatwave" not in index.lower():
            # convert date to floats and mask where missing
            data = indata[:, 1].astype(float)           
        else:
            # heatwave data contain "NA" as missing value
            data = indata[:, 1]
            data = np.array([d.strip() for d in data]) # in case of spaces as still strings
            bad_locs, = np.where(data == "NA")
            data[bad_locs] = utils.HADEX_MDI
            data = data.astype(float)

        data = np.ma.masked_where(data == utils.HADEX_MDI, data)
        times = np.ma.array(times, mask=data.mask)

        # if sufficient data, return True
        if timescale == "MON" and len(data.compressed()) >= 12*utils.MIN_YEARS_OF_INDEX_DATA:
            pass # makes logic easier
        elif timescale == "MON":
            has_data = False
            if diagnostics:
                print("station {} for index {} has insufficient data {}/{} (monthly)".format(station.id, index, len(data.compressed()), 12*utils.MIN_YEARS_OF_INDEX_DATA))

        # if sufficient data, return True
        if timescale == "ANN" and len(data.compressed()) >= utils.MIN_YEARS_OF_INDEX_DATA:
            pass
        elif timescale == "ANN":
            has_data = False
            if diagnostics:
                print("station {} for index {} has insufficient data {}/{} (annual)".format(station.id, index, len(data.compressed()), utils.MIN_YEARS_OF_INDEX_DATA))


    #*********************************************
    # if have data to play with - check when it ends   
    if has_data:
        # so it could contribute to the dataset
        if len(times.compressed()) > 0:
            # if data ends after start of dataset (1901)
            if times.compressed()[-1].year >= utils.STARTYEAR.year:
                pass
            else:
                has_data = False
                if diagnostics:
                    print("station {} for index {} only has data before {}".format(station.id, index, utils.STARTYEAR.year))

        # ensure that it has data in the recent period
        # this excludes all of e.g. ACRE
        # testing in September 2019
        if len(times.compressed()) > 0:
            # if data ends after most recent period - for studying more recent extremes
            if times.compressed()[-1].year >= utils.FINISH_AFTER.year:
                pass
            else:
                has_data = False
                if diagnostics:
                    print("station {} for index {} only has data before {}".format(station.id, index, utils.FINISH_AFTER.year))

    #*********************************************
    # If have data to play with - check most recent period in more depth
# commented out 29-July 2019
#    if has_data:

#        locs, = np.where(times > utils.MODERN_PERIOD)        
#        if len(times[locs].compressed()) > 1:

#            # check for large data gaps - if monthly or annual
#            if np.max(np.diff(times[locs].compressed())).days <= (utils.MAX_MISSING_YEARS + 1) * 366:
#                # The +1 as difference between timestamps already 1 year for ANN or 1 month for MON
#                pass
#            else:
#                has_data = False
#                if diagnostics:
#                    print("station {} for index {} has too large a data gap ({} days)".format(station.id, index, np.max(np.ma.diff(times[locs].compressed())).days))
        
    #*********************************************
    # more detailed checks for some input data sources and index combinations
# commented out 14 November - see how many this ends up with!
    # if has_data:
    #     if (("ecad" in station.location) and (index in PRECIP_INDICES)) or \
    #        (("ghcndex" in station.location) and (index in PRECIP_INDICES)):

    #         # lots of stations in ECAD precip indices - overwhelm distribution *and* processing time for gridding.

    #         # if sufficient data, return True
    #         if timescale == "MON" and len(data.compressed()) >= 1.5 * (12 * utils.MIN_YEARS_OF_INDEX_DATA):
    #             pass # makes logic easier
    #         elif timescale == "MON":
    #             has_data = False
    #             if diagnostics:
    #                 print("extra checks for {} in {}".format(index, station.location))
    #                 print("station {} for index {} has insufficient data {}/{} (monthly)".format(station.id, index, len(data.compressed()), 1.5*12*utils.MIN_YEARS_OF_INDEX_DATA))

    #         # if sufficient data, return True
    #         if timescale == "ANN" and len(data.compressed()) >= 1.5 * utils.MIN_YEARS_OF_INDEX_DATA:
    #             pass
    #         elif timescale == "ANN":
    #             has_data = False
    #             if diagnostics:
    #                 print("extra checks for {} in {}".format(index, station.location))
    #                 print("station {} for index {} has insufficient data {}/{} (annual)".format(station.id, index, len(data.compressed()), 1.5*utils.MIN_YEARS_OF_INDEX_DATA))


    #         if len(times.compressed()) > 0:
    #             # only keep series that go very recently
    #             if times.compressed()[-1].year >= utils.MODERN_PERIOD.year + 50: # 2001, not used as of October 2019 as all need to be in or after 2009
    #                 pass
    #             else:
    #                 has_data = False
    #                 if diagnostics:
    #                     print("extra checks for {} in {}".format(index, station.location))
    #                     print("station {} for index {} only has data after {}".format(station.id, index, utils.MODERN_PERIOD.year+50))

    #         locs, = np.where(times > utils.MODERN_PERIOD)        
    #         if len(times[locs].compressed()) > 1:

    #             # check for large data gaps - if monthly or annual
    #             if np.max(np.diff(times[locs].compressed())).days <= ((utils.MAX_MISSING_YEARS)/2. + 1) * 366:
    #                 # The +1 as difference between timestamps already 1 year for ANN or 1 month for MON
    #                 pass
    #             else:
    #                 has_data = False
    #                 if diagnostics:
    #                     print("extra checks for {} in {}".format(index, station.location))
    #                     print("station {} for index {} has too large a data gap ({} days)".format(station.id, index, np.max(np.ma.diff(times[locs].compressed())).days))

    #*********************************************
#    else:
#        has_data = False
#        if diagnostics:
#            print("station {} for index {} has no (suitable) data".format(station.id, index))

#    print(has_data)
    return has_data # assess_station


#*********************************************
def main(index="TX90p", diagnostics=False):
    """
    For all datasets, finds stations that exist for given index (and appropriate timescales)
    Checks for presence of data and write final station listing
    
    :param str index: which index to process
    :param bool diagnostics: output diagnostic information
    """
    
    # check if need to do monthly ones
    if index in utils.MONTHLY_INDICES:
        timescales = ["ANN", "MON"]
    else:
        timescales = ["ANN"]

    # read in all datasets
    all_datasets = utils.get_input_datasets()

    # for appropriate number of timescales
    for ts in timescales:
        print("{}".format(ts))

        # spin through each dataset
        for d, dataset in enumerate(all_datasets):

            dataset_stations = utils.read_inventory(dataset, subdir="formatted/indices")

            if diagnostics:
                print("{} - {}".format(dataset.name, index))

            final_inventory = []

            # check each station
            for stn in dataset_stations:

                if diagnostics:
                    print("{} - {}".format(dataset.name, stn.id))

                if assess_station(stn, index, ts, diagnostics=diagnostics):
                    final_inventory += [stn]
                    if diagnostics:
                        print("{}\n".format(len(final_inventory)))
                else:
                    if diagnostics:
                        print("\n")

            # then write everything out.
            utils.write_climpact_inventory_header(os.path.join(dataset.location, "{}.metadata.{}.{}.txt".format(dataset.name, index, ts)))

            for stn in final_inventory:
                utils.write_climpact_inventory(os.path.join(dataset.location, "{}.metadata.{}.{}.txt".format(dataset.name, index, ts)), stn)
                
            print("{} - {} stations".format(dataset.name, len(final_inventory)))

    return # main

#************************************************************************
#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', dest='index', action='store', default="TX90p",
                        help='Which index to run')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default = False')
 
    args = parser.parse_args()         

    main(index=args.index, diagnostics=args.diagnostics)

#------------------------------------------------------------
# END
#------------------------------------------------------------
