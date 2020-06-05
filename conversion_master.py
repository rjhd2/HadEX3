#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 444                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2020-01-13 16:13:06 +0000 (Mon, 13 Jan 2020) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""
Selects correct conversion script to run for specified dataset
merge_metadata.py invoked by typing::

  python conversion_master.py --indata "ecad" --diagnostics

Input arguments:

--indata        Which input data to use (ecad, sacad, lacad)

--diagnostics   Output extra info (default False)
"""
import os
import copy
import subprocess
import datetime as dt
import numpy as np

# RJHD scripts
import utils
import convert_utils 
import inventory_utils 

#*********************************************
def read_station(name, loc, station, diagnostics=False, lookup=None):
    """
    Reads the station data into Timeseries objects

    :param str name: indata name
    :param str loc: indata location
    :param stationObj station: station object
    :param book diagnostics: output diagnostic information
    :param dict lookup: lookup of station ID across parameters
    """
   
    # set up blank timeseries
    tmax = utils.Timeseries("TMAX", np.array([]), np.array([]))
    tmin = utils.Timeseries("TMIN", np.array([]), np.array([]))
    prcp = utils.Timeseries("PRCP", np.array([]), np.array([]))

    # read in using correct subroutine
    if name in ["ghcnd"]:
        tmax, tmin, prcp = convert_utils.read_ghcnd_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["spain", "nz", "south_america", "mexico", "brazil_sp"]:
        tmax, tmin, prcp = convert_utils.read_etccdi_format_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["russia"]:
        if loc.split("/")[-1][:5] == "v2018":
            # first submission was already in a nicer format (2020-01)
            tmax, tmin, prcp = convert_utils.read_etccdi_format_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics, extension="dat")
        elif loc.split("/")[-1][:5] == "v2020":
            # second needs sorting
            tmax, tmin, prcp = convert_utils.read_russia(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["iran"]:
        tmax, tmin, prcp = convert_utils.read_etccdi_format_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics, skip_header=1)
    elif name in ["acre"]:
        tmax, tmin, prcp = convert_utils.read_acre_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["honduras"]:
        tmax, tmin, prcp = convert_utils.read_honduras_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["australia"]:
        tmax, tmin, prcp = convert_utils.read_australia_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["chile"]:
        tmax, tmin, prcp = convert_utils.read_chile_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["colombia"]:
        tmax, tmin, prcp = convert_utils.read_colombia_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["canada"]:
        tmax, tmin, prcp = convert_utils.read_canada_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["decade"]:
        tmax, tmin, prcp = convert_utils.read_decade_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["india"]:
        tmax, tmin, prcp = convert_utils.read_india_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["japan"]:
        tmax, tmin, prcp = convert_utils.read_japan_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["brazil"]:
        tmax, tmin, prcp = convert_utils.read_brazil_station(loc, station, tmax, tmin, prcp, diagnostics=diagnostics)
    elif name in ["eobs", "laobs", "saobs"]:
        tmax, tmin, prcp = convert_utils.read_eobs_station(loc, station, tmax, tmin, prcp, lookup, diagnostics=diagnostics)

    # store with the station object
    station.tmax = tmax
    station.tmin = tmin
    station.prcp = prcp

    return station # read_station

#*********************************************
def fill_missing_months(station, diagnostics=False):
    """
    Checks the station timeseries and inputs missing months

    :param stationObj station: station object
    :param book diagnostics: output diagnostic information
    """

    # if at least one data field is present
    if len(station.prcp.times) > 0 or len(station.tmax.times) > 0 or len(station.tmin.times) > 0:
    
        start = np.min([station.prcp.times[0], station.tmax.times[0], station.tmin.times[0]])

        end = np.max([station.prcp.times[-1], station.tmax.times[-1], station.tmin.times[-1]])

        difference = end - start

        # super set of days and data
        all_days = np.array([start + dt.timedelta(days=d) for d in range(difference.days + 1)])

        for parameter in convert_utils.OBSERVABLES:
            all_data = np.ones(all_days.shape[0]) * utils.HADEX_MDI

            st_var = getattr(station, parameter.lower())

            # match times and fill
            # matching is quicker with arrays of integers, so extract timedeltas
            all_days_diffs = np.array([(ad - start).days for ad in all_days])
            st_var_diffs = np.array([(ad - start).days for ad in st_var.times])

            match = np.in1d(all_days_diffs, st_var_diffs)
            match_back = np.in1d(st_var_diffs, all_days_diffs)
            all_data[match] = st_var.data[match_back]

            # old way using the dt.datetimes - very slow in Py3
    #        match = np.in1d(all_days, st_var.times)
    #        match_back = np.in1d(st_var.times, all_days)
    #        all_data[match] = st_var.data[match_back]
            # (match_back needed in case of malformed dates - esp in ACRE)

            st_var.times = all_days
            st_var.data = np.copy(all_data)

            setattr(station, parameter, st_var)
        
    return station # fill_missing_months

#*********************************************
def process_dataset(dataset, select_on_completeness=False, diagnostics=False):
    """
    Read in the specified dataset inventory.  Check each station for suitability for HadEX3

    :param datasetObj dataset: dataset object holding metadata about the input dataset
    :param bool select_on_completeness: Assess dataset completeness (default - False)
    :param bool diagnostics: output diagnostic information

    """

    # read in the dataset inventory metadata
    if dataset.name == "ghcnd":
        dataset_stations = inventory_utils.read_ghcnd(dataset, diagnostics=diagnostics)
    if dataset.name == "spain":
        dataset_stations = inventory_utils.read_spain(dataset, diagnostics=diagnostics)
    if dataset.name == "nz":
        dataset_stations = inventory_utils.read_nz(dataset, diagnostics=diagnostics)
    if dataset.name == "honduras":
        dataset_stations = inventory_utils.read_honduras(dataset, diagnostics=diagnostics)
    if dataset.name == "russia":
        dataset_stations = inventory_utils.read_russia(dataset, diagnostics=diagnostics)
    if dataset.name == "acre":
        dataset_stations = inventory_utils.read_acre(dataset, diagnostics=diagnostics)
    if dataset.name == "australia":
        dataset_stations = inventory_utils.read_australia(dataset, diagnostics=diagnostics)
    if dataset.name == "south_america":
        # need to pre-process these data as spread across multiple files and no metadata.
        convert_utils.pre_process_south_america(dataset, diagnostics=diagnostics)
        dataset_stations = inventory_utils.read_generic_obs(dataset, diagnostics=diagnostics)
    if dataset.name == "colombia":
        dataset_stations = inventory_utils.read_colombia(dataset, diagnostics=diagnostics)
    if dataset.name == "chile":
        dataset_stations = inventory_utils.read_chile(dataset, diagnostics=diagnostics)
    if dataset.name == "canada":
        dataset_stations = inventory_utils.read_canada(dataset, diagnostics=diagnostics)
    if dataset.name == "decade":
        dataset_stations = inventory_utils.read_decade(dataset, diagnostics=diagnostics)
    if dataset.name == "india":
        dataset_stations = inventory_utils.read_india(dataset, diagnostics=diagnostics)
    if dataset.name == "japan":
        dataset_stations = inventory_utils.read_japan(dataset, diagnostics=diagnostics)
    if dataset.name == "mexico":
        dataset_stations = inventory_utils.read_mexico(dataset, diagnostics=diagnostics)
    if dataset.name == "iran":
        dataset_stations = inventory_utils.read_iran(dataset, diagnostics=diagnostics)
    if dataset.name == "brazil":
        dataset_stations = inventory_utils.read_brazil(dataset, diagnostics=diagnostics)
    if dataset.name == "eobs":
        dataset_stations, lookup = inventory_utils.read_eobs(dataset, diagnostics=diagnostics)
    if dataset.name == "laobs":
        dataset_stations, lookup = inventory_utils.read_eobs(dataset, diagnostics=diagnostics)
    if dataset.name == "saobs":
        dataset_stations, lookup = inventory_utils.read_eobs(dataset, diagnostics=diagnostics)
    if dataset.name == "brazil_sp":
        # need to pre-process these data as spread across two file type sets and no metadata.
        convert_utils.pre_process_brazil_sp(dataset, diagnostics=diagnostics)
        dataset_stations = inventory_utils.read_generic_obs(dataset, diagnostics=diagnostics)

    # check that output directory exists
    if os.path.exists(os.path.join(dataset.location, "formatted")):
        # remove any previous processing - and files
        subprocess.check_call(["rm", "-rf", os.path.join(dataset.location, "formatted")])
#        shutil.rmtree(os.path.join(dataset.location, "formatted"), ignore_errors = True)
        os.mkdir(os.path.join(dataset.location, "formatted"))
    else:
        os.mkdir(os.path.join(dataset.location, "formatted"))

    # write headers on inventory list
    utils.write_climpact_inventory_header(os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name)))

    # loop through all stations
    for station in dataset_stations:

        if diagnostics:
            print(station)

        # some might be missing in the inventory lists
        if station.latitude != station.latitude:
            # lat not number so just skip
            continue
        if station.longitude != station.longitude:
            # lon not number so just skip
            continue
            
        # need to copy so that not all the data are stored in the dataset_stations list
        this_station = copy.deepcopy(station)

        if dataset.name in ["eobs", "laobs", "saobs"]:
            read_station(dataset.name, dataset.location, this_station, diagnostics=diagnostics, lookup=lookup)
        else:
            read_station(dataset.name, dataset.location, this_station, diagnostics=diagnostics)

        # Climpact can't cope with missing data, so fill with MDI
        this_station = fill_missing_months(this_station)

        # need to select stations which have suitable data!
        good_end = False
        good_length = False
        good_completeness = False

        # as long as one of the parameters has enough data, then use all.
        for parameter in convert_utils.OBSERVABLES:
 
            st_var = getattr(this_station, parameter.lower())

            # if there is data
            if len(st_var.times) > 0: 
                # extract information
                start = st_var.times[0]
                end = st_var.times[-1]
                nobs = len(st_var.data)
                diffs = np.diff(st_var.times) # first differences, to find large gaps.

                # derived values
                length = end - start + dt.timedelta(days=1)
                completeness = nobs / float(length.days)

                #****
                # Commented out on 29-8-2018 - want early period data, even if it stops short of 1951
                #****
                # data series has to end in more recent period
                # if end >= utils.FINISH_AFTER:
                good_end = True

                # need at least a few (20) of years of data (non-contiguous)
                if length >= utils.MIN_LENGTH_OF_RECORD:
                    good_length = True

                # select on completeness within the record (too many missing months/years aren't helpful)
                #    optional check - off by default
                if select_on_completeness:
                    good_yrs_days, good_yrs_days_months = utils.check_completeness_indata(st_var)
                    if len(good_yrs_days) > (utils.MIN_LENGTH_OF_RECORD.days / 365.) or \
                            len(good_yrs_days_months) > (utils.MIN_LENGTH_OF_RECORD.days / 365.):
                        good_completeness = True
                else:
                    good_completeness = True

            if diagnostics:
                if select_on_completeness:
                    print("{}:  end: {} length: {} complete: {}".format(parameter, good_end, good_length, good_completeness))
                else:
                    print("{}:  end: {} length: {}".format(parameter, good_end, good_length))

        # ensure no spaces in names
        if " " in this_station.id:
            this_station.id = "".join(this_station.id.split(" "))


        # only write if all three tests matched at some point in the record
        #   across all variables
        if good_end and good_length and good_completeness:               
            if diagnostics: 
                print(station.id)

            utils.write_climpact_inventory(os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name)), this_station)

            utils.write_climpact_format(os.path.join(dataset.location, "formatted"), this_station)
        else:
            if diagnostics: 
                print("skipping {}".format(station.id))

    return # process_dataset

#*********************************************
def main(indata="ghcnd", diagnostics=False):
    """
    Extract relevant dataset from the command line switchs

    :param str indata: input dataset to process
    :param bool diagnostics: output diagnostic information
    """

    # get all possible datasets
    all_datasets = utils.get_input_datasets()

    # and their names
    names = np.array([d.name for d in all_datasets])

    # if dataset selected and in the list of available, then run
    if indata in names:        
        process_dataset(all_datasets[names == indata][0], diagnostics=diagnostics)

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
    parser.add_argument('--indata', dest='indata', action='store', default="ghcnd",
                        help='Which dataset to convert')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default = False')
 
    args = parser.parse_args()         

    if args.indata in ["ghcnd", "spain", "nz", "honduras", "russia", "acre", "south_america", \
                           "australia", "colombia", "chile", "canada", "decade", "eobs", "laobs", \
                       "saobs", "india", "japan", "mexico", "iran", "brazil", "brazil_sp"]:
        # extra check to ensure this isn't run with index data
        main(indata=args.indata, diagnostics=args.diagnostics)
    
#    elif args.indata in []:
#        # these data are already in correct format for Climpact2
#        if diagnostics:
#            print("{} indata already in correct format for Climpact2".format(args.indata))

#*******************************************
# END
#*******************************************
