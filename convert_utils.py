##!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 448                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2020-01-14 10:01:43 +0000 (Tue, 14 Jan 2020) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""
Conversion utility scripts - for each input datast type
"""
import os
import datetime as dt
import calendar
import glob
import pandas
import numpy as np

from six.moves import input

# RJHD scripts
import utils

OBSERVABLES = ["TMAX", "TMIN", "PRCP"]

#*********************************************
def F_to_C(T):
    """ Convert from Fahrenheit to Celsius"""

    # don't translate missing data values!
    goods, = np.where(T != utils.HADEX_MDI)

    T[goods] = ((T[goods]-32.)*(5./9.))

    return T # F_to_C


#*********************************************
def remove_trace_and_cumulative(pcp):
    """
    Remove codes for trace and cumulative amounts from precipitation data

    """

    traces, = np.where(pcp == -3)
    cumulative, = np.where(pcp == -4)

    if len(traces) > 0:
        pcp[traces] = 0
    if len(cumulative) > 0:
        pcp[cumulative] = 0

    return pcp # remove_trace_and_cumulative

#*********************************************
def check_duplicate_times(times, diagnostics=False):
    """
    Helper routine to check for duplicated timestamps in input data

    :param array times: np array of date objects
    :param bool diagnostics: output diagnostic information
    """

#    https://stackoverflow.com/questions/30003068/get-a-list-of-all-indices-of-repeated-elements-in-a-numpy-array
    sort_order = np.argsort(times)
    sorted_times = times[sort_order]
    values, start, count = np.unique(sorted_times, return_counts=True, return_index=True)

    # sets of indices
    locs = np.split(sort_order, start[1:])

    # filter them with respect to their size, keeping only items occurring more than once
    values = values[count > 1]
    # python 3 filter() is now an iterable, so force to be a list is a quick fix 14-2-2019 RJHD
    locs = list(filter(lambda x: x.size > 1, locs))

    # print out if required
    if diagnostics:
        for v in values:
            print("Duplicate times at {}".format(v.isoformat()))

    return locs # check_duplicate_times

#*********************************************
def process_duplicates(times, tx, tn, pcp, diagnostics=False):
    """
    Identify duplicates and then remove these from the data

    :param array times: array of date time values
    :param array tx: array of Tmax values
    :param array tn: array of Tmin values
    :param array pcp: array of P values
    :param bool diagnostics: verbose output.

    :returns: times, tx, tn, pcp
    """

    non_unique = check_duplicate_times(times, diagnostics=diagnostics)

    if len(non_unique) > 0:
        mask = np.ones(times.shape[0])

        for entry in non_unique:
            for pos in entry:
                mask[pos] = 0

        mask = [mask.astype("bool")]

        times = times[mask]
        tx = tx[mask]
        tn = tn[mask]
        pcp = pcp[mask]

        assert len(times) == len(tx)
        assert len(np.unique(times)) == len(times)

    return times, tx, tn, pcp # process_duplicates

#*********************************************
def remove_bad_dates(bad_dates, times, tx, tn, pcp):
    """
    Remove entries from malformed dates

    :param array bad_dates: locations of malformed dates
    :param array times: array of date time values
    :param array tx: array of Tmax values
    :param array tn: array of Tmin values
    :param array pcp: array of P values
   
    :returns: times, tx, tn, pcp
    """

    mask = np.ones(tx.shape[0])
    for bd in bad_dates:
        mask[bd] = 0

    mask = [mask.astype("bool")]

    # times already missing these values
    tx = tx[mask]
    tn = tn[mask]
    pcp = pcp[mask]

    assert len(times) == len(tx)

    return times, tx, tn, pcp # remove_bad_dates

#*********************************************
def fill_missing_timeseries(tmax, tmin, prcp):
    '''
    A whole parameter missing - fill with missing data

    If just different lengths, then processed in parent routine
    '''

    # if mismatch in T data or P data need to fill with missing values
    if len(tmax.times) == 0 and len(prcp.times) != 0:            
        tmax.times = prcp.times
        tmax.data = np.zeros(len(prcp.data))
        tmax.data[:] = utils.HADEX_MDI
    elif len(tmax.times) == 0 and len(tmin.times) != 0:            
        tmax.times = tmin.times
        tmax.data = np.zeros(len(tmin.data))
        tmax.data[:] = utils.HADEX_MDI

    if len(tmin.times) == 0 and len(prcp.times) != 0:            
        tmin.times = prcp.times
        tmin.data = np.zeros(len(prcp.data))
        tmin.data[:] = utils.HADEX_MDI
    elif len(tmin.times) == 0 and len(tmax.times) != 0:            
        tmin.times = tmax.times
        tmin.data = np.zeros(len(tmax.data))
        tmin.data[:] = utils.HADEX_MDI

    if len(prcp.times) == 0 and len(tmax.times) != 0:            
        prcp.times = tmax.times
        prcp.data = np.zeros(len(tmax.data))
        prcp.data[:] = utils.HADEX_MDI
    elif len(prcp.times) == 0 and len(tmin.times) != 0:            
        prcp.times = tmin.times
        prcp.data = np.zeros(len(tmin.data))
        prcp.data[:] = utils.HADEX_MDI

    # if different lengths of record for Tx/Tn and P, then not a problem
    #  matched on writing out 

    return tmax, tmin, prcp # fill_missing_timeseries


#*********************************************
def find_south_american_metadata(dataset, country, stn_id):
    """
    Have three places to look for South American Metadata
    Country file
    Combined (CRC) file
    Pub A

    """


    # try three files for metadata
    have_metadata = False
    this_station = -1

    # first try metadata for country (but don't look if CRC as that's only in the second option)
    if country != "CRC":
        if os.path.exists(os.path.join(dataset.location, "raw", "Metadatos_estaciones_{}.csv".format(country))):
            metadata = np.genfromtxt(os.path.join(dataset.location, "raw", "Metadatos_estaciones_{}.csv".format(country)), skip_header=4, dtype=str, delimiter=",", encoding="latin-1")

            names = metadata[:, 0]
            ids = metadata[:, 1]
            lons = metadata[:, 3].astype(float)
            lats = metadata[:, 4].astype(float)

            if stn_id in ids:
                loc, = np.where(ids == stn_id)[0]
                this_station = utils.Station("{}_{}".format(country, stn_id), lats[loc], lons[loc], os.path.join(dataset.location, "raw"), dataset.name)
                have_metadata = True

    # then try more general file
    if not have_metadata and os.path.exists(os.path.join(dataset.location, "raw", "Metadatos_estaciones_CRC.csv")):
        metadata = np.genfromtxt(os.path.join(dataset.location, "raw", "Metadatos_estaciones_CRC.csv"), skip_header=1, dtype=str, delimiter=",", encoding="latin-1")

        countries = metadata[:, 0]
        names = metadata[:, 2]
        ids = metadata[:, 1]
        lons = metadata[:, 6].astype(float)
        lats = metadata[:, 5].astype(float)

        # test both the ID and country (not all are WMO ids)
        loc, = np.where(np.logical_and(ids == stn_id, country == countries))

        if len(loc) == 1:
            this_station = utils.Station("{}_{}".format(country, stn_id), lats[loc[0]], lons[loc[0]], os.path.join(dataset.location, "raw"), dataset.name)
            have_metadata = True

    # finally try Pub A (and fix ID)
    if not have_metadata:
        this_station = utils.parse_PubA(stn_id, dataset)
        if this_station != -1:
            this_station.id = "{}_{}".format(country, stn_id)
            have_metadata = True

    return this_station # find_south_american_metadata

#*********************************************
def read_ghcnd_station(loc, station, tmax, tmin, prcp, diagnostics=False):
    """
    Reads the GHCND station data into timeseries objects

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information (unused)
    """

    # set up parseing info
    fieldwidths = [11, 4, 2, 4]
    for i in range(31):
        fieldwidths += [5, 3] # value and flag
    fieldwidths = tuple(fieldwidths)


    try:
        # read as pandas dataframe
        df = pandas.read_fwf(os.path.join(loc, "raw", "{}.dly".format(station.id)), widths=fieldwidths, header=None)
        # convert to numpy
        indata = df.to_numpy()        

        # using old code in line-by-line for the moment (03-2019)
        for fields in indata:

            if fields[0] != station.id:
                print("file ID number doesn't match station: {} != {}".format(fields[0], station.id))
                break

            # get the useful data
            year = int(fields[1])
            month = int(fields[2])
            parameter = fields[3]

            # extract just the measurements from the value,flag pairs, convert to correct units
            month_data = np.array(fields[4:]).reshape(-1, 2)[:, 0].astype(int)/10.
            month_data = month_data.reshape(-1)
            month_data[month_data == -9999/10.] = utils.HADEX_MDI

            # now to make the timeseries
            times = []
            for d in range(calendar.monthrange(year, month)[1]):
                times += [dt.date(year, month, d+1)]
            times = np.array(times)

            month_data = month_data[:len(times)]

            # append the timeseries
            if parameter == "TMAX":
                tmax = utils.append_timeseries(tmax, times, month_data)
            elif parameter == "TMIN":
                tmin = utils.append_timeseries(tmin, times, month_data)                   
            elif parameter == "PRCP":
                prcp = utils.append_timeseries(prcp, times, month_data)

    except IOError:
        print("File doesn't exist {}".format(os.path.join(loc, "raw", "{}.dly".format(station.id))))
        # file doesn't exist
        # return empty lists

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_ghcnd_station

#*********************************************
def read_honduras_station(loc, station, tmax, tmin, prcp, diagnostics=False):
    """
    Reads the Honduras station data into timeseries objects

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information (unused)
    """
    import glob
    #***********************
    def process_data(line):
        """
        Process each line of .csv file, extracting date and data to return 2 arrays
        """
        
        # extract data
        if "-" in line[1]:
            # years/months stored in single columns "YYYY-MM"
            year = int(line[1][:4])
            month = int(line[1][-2:])
            raw_data = [utils.HADEX_MDI if e.strip() == "" else float(e) for e in line[2:-4]]
        else:
            # years/months stored in two columns
            year = int(line[1])
            month = int(line[2])
            raw_data = [utils.HADEX_MDI if e.strip() == "" else float(e) for e in line[3:-4]]

        # spin through entries
        times = []
        data = []
        day = 0
        for d in raw_data:
            day += 1

            # cope with odd length months
            try:
                times += [dt.date(year, month, day)]
                data += [d]
            
            except ValueError:
                # this day doesn't exist - skip remaining entries in the line (e.g. Feb 30th)
                break

        return np.array(times), np.array(data) # process_data

    # do minimum temperatures first, and then maximum.  No precip
    try:
        with open("{}/Temperatura Minima Abs. {}.csv".format(os.path.join(loc, "raw"), station.id), "r", encoding="latin-1") as minfile:
            for l, line in enumerate(minfile):                       
                if l >= 8 and line[0] != ",":
                    times, data = process_data(line.split(","))
                    tmin = utils.append_timeseries(tmin, times, data)
    except IOError:
        print("No file: {}/Temperatura Minima Abs. {}.csv".format(os.path.join(loc, "raw"), station.id))

    # maximum temperature file
    try:
        with open('{}/Temperatura Maxima Abs. {}.csv'.format(os.path.join(loc, "raw"), station.id), "r", encoding="latin-1") as maxfile:
            for l, line in enumerate(maxfile):                       
                if l >= 8 and line[0] != ",":
                    times, data = process_data(line.split(","))
                    tmax = utils.append_timeseries(tmax, times, data)           
    except IOError:
        print("No file: {}/Temperatura Maxima Abs. {}.csv".format(os.path.join(loc, "raw"), station.id))

    # as no preciptation data
    prcp.times = tmax.times
    prcp.data = np.zeros(len(tmax.data))
    prcp.data[:] = utils.HADEX_MDI

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_honduras_station

#*********************************************
def read_etccdi_format_station(loc, station, tmax, tmin, prcp, diagnostics=False, delimiter="", extension="txt", skip_header=0):
    """
    Reads the station data into timeseries objects (which are already in ETCCDI format)
    ETCDDI format: YYYY MM DD PP TX TN

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information
    :param str delimiter: separator if not space or tab
    :param str extension: file extension, txt (default), dat
    :param int skip_header: header lines to skip
    """
   
    # go through each line in the file
    #   just skip invalid lines
    try:
        all_data = np.genfromtxt(os.path.join(loc, "raw", "{}.{}".format(station.id, extension)), invalid_raise=False, delimiter=delimiter, encoding="latin-1", skip_header=skip_header)

    except OSError:
        print("File not found: {}".format(os.path.join(loc, "raw", "{}.{}".format(station.id, extension))))
        return tmax, tmin, prcp
        
    if all_data.shape == (0,):
        print("File empty: {}".format(os.path.join(loc, "raw", "{}.{}".format(station.id, extension))))
        return tmax, tmin, prcp
        
    # extract information into arrays
    years = all_data[:, 0]
    month = all_data[:, 1]
    day = all_data[:, 2]
    pcp = all_data[:, 3]
    tx = all_data[:, 4]
    tn = all_data[:, 5]

    #*********************
    # make times array
    times = []
    bad_dates = []
    for y, year in enumerate(years):
        try:
            times += [dt.date(int(year), int(month[y]), int(day[y]))]
        except ValueError:
            # impossible_date
            bad_dates += [y]
            if diagnostics:
                print("Station {}: bad date at {}-{}-{}".format(station.id, int(year), int(month[y]), int(day[y])))

    times = np.array(times)

    #*********************
    # remove bad dates
    if len(bad_dates) > 0:
        times, tx, tn, pcp = remove_bad_dates(bad_dates, times, tx, tn, pcp)

    #*********************
    # check for duplicates
    times, tx, tn, pcp = process_duplicates(times, tx, tn, pcp, diagnostics=diagnostics)

        
    #*********************
    # remove trace & cumulative precip
    pcp = remove_trace_and_cumulative(pcp)

    #*********************
    # append the timeseries
    tmax.data = tx
    tmax.times = times
    tmin.data = tn
    tmin.times = times
    prcp.data = pcp
    prcp.times = times

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_etccdi_format_station

#*********************************************
def read_acre_station(loc, station, tmax, tmin, prcp, diagnostics=False, extension="txt"):
    """
    Reads the station data into timeseries objects
    Some are in ETCCDI format, others swap PP and TX/TN around, varying delimiters

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information
    :param str extension: file extension, txt (default), dat
    """

    def fix_mdi(indata, testvalue):
        
        locs, = np.where(indata == testvalue)

        indata[locs] = utils.HADEX_MDI

        return indata # fix_mdi

    def swap(a, b):

        temp = np.copy(b)
        b = np.copy(a)
        a = temp

        return a, b # swap
        

    #*********************
    #*********************
    import qc_checks # to get at world record values
   
    with open(os.path.join(loc, "raw", "{}.{}".format(station.id, extension)), "r", encoding="latin-1") as infile:

        line = infile.readline()
            
        # find delimiter
        if len(line.split(" ")) > 1:
            delimiter = " "
        elif len(line.split(",")) > 1:
            delimiter = ","

    all_data = np.genfromtxt(os.path.join(loc, "raw", "{}.{}".format(station.id, extension)), invalid_raise=False, delimiter=delimiter, encoding="latin-1")

    #*********************
    # Presume a system, fail if not.  Going for Tx, Tn, P, with T in Fahrenheit
    years = all_data[:, 0]
    month = all_data[:, 1]
    day = all_data[:, 2]
    tx = all_data[:, 3]
    tn = all_data[:, 4]
    pcp = all_data[:, 5]

    #*********************
    # fix any poorly coded MDIs (-999 rather than -99.9)
    tx = fix_mdi(tx, 10*utils.HADEX_MDI)
    tn = fix_mdi(tn, 10*utils.HADEX_MDI)
    pcp = fix_mdi(pcp, 10*utils.HADEX_MDI)

    # mask to improve later
    tx = np.ma.masked_where(tx == utils.HADEX_MDI, tx)
    tn = np.ma.masked_where(tn == utils.HADEX_MDI, tn)
    pcp = np.ma.masked_where(pcp == utils.HADEX_MDI, pcp)

    #*********************
    # use TX/TN order to identify where columns are in the wrong order
    bad_ts, = np.ma.where(tn > tx)

    if len(bad_ts) > 0:
        # either TX and TN swapped OR P in the wrong column

        bad_tns, = np.where(tn != utils.HADEX_MDI)
        bad_txs, = np.where(tx != utils.HADEX_MDI)
        bad_ps, = np.where(pcp != utils.HADEX_MDI)

        if len(bad_tns) == len(bad_txs):
            # likely to be swap for TX/TN
            tx, tn = swap(tx, tn)
        else:     
            # test percentage difference and against PCP
            t_diff = abs(float(len(bad_tns) - len(bad_txs))/len(bad_txs))
            p_diff = abs(float(len(bad_ps) - len(bad_txs))/len(bad_txs))
            
            if t_diff < 0.05 and p_diff > 0.5:
                # small difference between Ts, larger between Tx/Tn and P
                # so correct order, just different amounts of missing data between Tx and Tn (5%)
                pass
            else:            
                input("Ts to fix")

    #*********************
    # use world record check value to spot where likely to be swapped between F and C
    bad_units, = np.where(tx > qc_checks.WORLD_RECORDS_VALUES["Tx"])

    if len(bad_units) > 0:
        if diagnostics:
            print("converting to Fahrenheit")
        # presume Fahrenheit
        tx = F_to_C(tx)
        tn = F_to_C(tn)
    
       
    #*********************
    # make times array
    times = []
    bad_dates = []
    for y, year in enumerate(years):
        try:
            times += [dt.date(int(year), int(month[y]), int(day[y]))]
        except ValueError:
            # impossible_date
            bad_dates += [y]
            if diagnostics:
                print("Station {}: bad date at {}-{}-{}".format(station.id, int(year), int(month[y]), int(day[y])))

    times = np.array(times)

    #*********************
    # remove bad dates
    if len(bad_dates) > 0:
        times, tx, tn, pcp = remove_bad_dates(bad_dates, times, tx, tn, pcp)
        
    #*********************
    # check for duplicates
    times, tx, tn, pcp = process_duplicates(times, tx, tn, pcp, diagnostics=diagnostics)
        
    #*********************
    # remove trace & cumulative precip
    pcp = remove_trace_and_cumulative(pcp)

    #*********************
    # append the timeseries
    tmax.data = tx
    tmax.times = times
    tmin.data = tn
    tmin.times = times
    prcp.data = pcp
    prcp.times = times

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_acre_station

#*********************************************
def pre_process_south_america(dataset, diagnostics=False):
    """
    The South American data are in large files, 1-3 per country, in linear order

    Need to split out into individual station files and create an inventory

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """
    #***************************
    def fix(value):
        """
        Fix the read in values to remove the missing data indicator and convert to floats
        """
        if value.strip() == "\\N":
            value = utils.HADEX_MDI
        else:
            value = float(value)

        return value # fix

    #*********************
    #*********************
    # remove all files from before
    for filename in glob.glob(r'{}/*.txt'.format(os.path.join(dataset.location, "raw"))):
        os.remove(filename)


    # spin through parent files
    for filename in sorted(glob.glob(r'{}/Registros*.csv'.format(os.path.join(dataset.location, "raw")))):
        print(filename)
        
        # extract country
        country = filename.split("/")[-1].split(" ")[2].split("_")[0]
        if country == "(SMN)": country = filename.split("/")[-1].split(" ")[1]

        # read each line at a time
        with open(filename, "r", encoding="latin-1") as infile:
            for line in infile:
                if diagnostics:
                    print(line)

                # get the delimiter
                delimiter = ","
                if len(line.split(delimiter)) == 1:
                    delimiter = ";"

                # skip the header
                sline = line.split(delimiter)
                if sline[0] == "omm_id": continue

                # extract the data
                stn_id = int(sline[0])
                date = dt.datetime.strptime(sline[1], "%d/%m/%Y")

                tx, tn, pcp = sline[2], sline[3], sline[4]
                    
                tx = fix(tx)
                tn = fix(tn)
                pcp = fix(pcp)

                # make and write the file
                if not os.path.exists(os.path.join(dataset.location, "raw", "{}_{}.txt".format(country, stn_id))):
                    # if it doesn't exist, then create
                    outfile = open(os.path.join(dataset.location, "raw", "{}_{}.txt".format(country, stn_id)), "w")

                else:
                    # if it exists, then append
                    outfile = open(os.path.join(dataset.location, "raw", "{}_{}.txt".format(country, stn_id)), "a")

                outfile.write("{:-8d}{:-8d}{:-8d}{:-8.1f}{:-8.1f}{:-8.1f}\n".format(date.year, date.month, date.day, pcp, tx, tn))
                outfile.close()
                               
    # make inventory
    metadata_file = open(os.path.join(dataset.location, "{}_stations.txt".format(dataset.name)), "w")

    for filename in sorted(glob.glob(r'{}/*.txt'.format(os.path.join(dataset.location, "raw")))):

        # extract the ID from the filename
        file_root = filename.split("/")[-1].split(".")[0]
        country = file_root.split("_")[0]
        stn_id = file_root.split("_")[1]
        
        this_station = find_south_american_metadata(dataset, country, stn_id)

#        print this_station
        if this_station != -1:
            metadata_file.write("{:10s} {:7.3f} {:7.3f}\n".format(this_station.id, this_station.latitude, this_station.longitude))

    metadata_file.close()

    return # pre_process_south_america

#*********************************************
def read_australia_station(loc, station, tmax, tmin, prcp, diagnostics=False):
    """
    Reads the Australian station data into timeseries objects

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information (unused)
    """

    def fix_t(indata):
        
        locs, = np.where(indata != -999)

        indata[locs] = indata[locs]/10.

        locs, = np.where(indata == -999)

        indata[locs] = utils.HADEX_MDI        

        return indata # fix_t
    
    #*********************
    # Temperature data
    try:
        # extract all T data
        all_T_data = np.genfromtxt(os.path.join(loc, "raw", "hq17dtr{}".format(station.id)), delimiter="", dtype=str, encoding="latin-1")

        # pull out combined ID and date field, and split
        ID_date = all_T_data[:, 0]

        years, months, days = [], [], []
        for entry in ID_date:
            ID = entry[:6]
            years += [int(entry[6:10])]
            months += [int(entry[10:12])]
            days += [int(entry[12:14])]

        # make times array
        T_times = []
        bad_dates = []
        for y, year in enumerate(years):
            try:
                T_times += [dt.date(int(year), int(months[y]), int(days[y]))]
            except ValueError:
                # impossible_date
                bad_dates += [y]
                if diagnostics:
                    print("Station {}: bad date at {}-{}-{}".format(station.id, int(year), int(months[y]), int(days[y])))
        T_times = np.array(T_times)

        # has spaces for 2019 - which causes issues later on
        good_times = np.where(T_times < dt.date(2019, 1, 1))
        T_times = T_times[good_times]

        # extract the temperature fields
        tx = all_T_data[good_times, 1].astype(float)
        tn = all_T_data[good_times, 2].astype(float)

        # sort the x10 multiplier
        tx = fix_t(tx.squeeze())
        tn = fix_t(tn.squeeze())
        
        # add to series
        tmax = utils.append_timeseries(tmax, T_times, tx)
        tmin = utils.append_timeseries(tmin, T_times, tn)

    except IOError:
        print("File doesn't exist {}".format(os.path.join(loc, "raw", "hq17dtr{}".format(station.id))))
        # file doesn't exist
        # return empty lists
        pass


    #*********************
    # Precipitation data
    try:
        # extract all T data
        all_P_data = np.genfromtxt(os.path.join(loc, "raw", "DC02D_Data_{}_46309839556073.txt".format(station.id)), delimiter=",", dtype=str, skip_header=1, encoding="latin-1")

        # pull out ID and date fields, and split
        ID = all_P_data[:, 1]
        date = all_P_data[:, 2]

        years, months, days = [], [], []
        for entry in date:
            dd, mm, yy = entry.split("/")
            years += [int(yy)]
            months += [int(mm)]
            days += [int(dd)]

        # make times array - longhand to allow for checking
        P_times = []
        bad_dates = []
        for y, year in enumerate(years):
            try:
                P_times += [dt.date(int(year), int(months[y]), int(days[y]))]
            except ValueError:
                # impossible_date
                bad_dates += [y]
                if diagnostics:
                    print("Station {}: bad date at {}-{}-{}".format(station.id, int(year), int(months[y]), int(days[y])))
        P_times = np.array(P_times)

        # has spaces for 2019 - which causes issues later on
        good_times, = np.where(P_times < dt.date(2019, 1, 1))
        P_times = P_times[good_times]

        # extract the precipitation fields
        pptn = all_P_data[good_times, 3]
        bads = np.where(pptn == "      ")
        pptn[bads] = utils.HADEX_MDI
        pptn = pptn.astype(float)

        p_qc = all_P_data[good_times, 4]
        ndays = all_P_data[good_times, 5]
        bads = np.where(ndays == "  ")
        ndays[bads] = 0
        ndays = ndays.astype(int)

        accumul_days = all_P_data[good_times, 6]
        bads = np.where(accumul_days == "  ")
        accumul_days[bads] = 0
        accumul_days = accumul_days.astype(int)

        # sort the quality info       
        bads = np.where(np.logical_and(p_qc != "Y", p_qc != "N"))
        pptn[bads] = utils.HADEX_MDI
        # not just !=, as these aspects not set if no rain.
        bads = np.where(np.logical_or(ndays > 1, accumul_days > 1))
        pptn[bads] = utils.HADEX_MDI

        # add to series
        prcp = utils.append_timeseries(prcp, P_times, pptn)

    except IOError:
        print("File doesn't exist {}".format(os.path.join(loc, "raw", "DC02D_Data_{}_46309839556073.txt".format(station.id))))
        # file doesn't exist
        # return empty lists
        
    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_australia_station

#*********************************************
def read_chile_station(loc, station, tmax, tmin, prcp, diagnostics=False):
    """
    Reads the Chilean station data into timeseries objects

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information (unused)
    """

    #****************
    def read_file(infile):
        '''
        Read in file and make into single column if required
        '''

        indata = np.genfromtxt(infile, dtype=(float), delimiter=" ", encoding="latin-1")

        if len(indata.shape) == 2:
            # have 2-d array
            indata = indata.reshape(-1)

        return indata # read_file
    #****************
    
    # sort out names
    name = "_".join([s.capitalize() for s in station.id.split("_")])

    for suffix in ["tmax", "tmin", "rr"]:

        # try to read in one file
        infilename = os.path.join(station.location, "{}_{}.txt".format(name, suffix))
        if os.path.exists(infilename):
            homog = False
            obs_data = read_file(infilename)

        else:
            # else try for homogenised file
            infilename = os.path.join(station.location, "{}_{}h.txt".format(name, suffix))
            homog = True
            if os.path.exists(infilename): 
                obs_data = read_file(infilename)

        # set up start/end date of records
        end = dt.date(2019, 1, 1)
        if name not in ["Calama", "Ipascua"]:
            start = dt.date(1961, 1, 1)
        elif name == "Calama":
            start = dt.date(1971, 1, 1)
        elif name == "Ipascua":
            start = dt.date(1981, 1, 1)

        # make list of dates appropriate for location
        dt_dates = [start]
        while True:
            dt_dates += [dt_dates[-1] + dt.timedelta(days=1)]
            if dt_dates[-1] == end:
                dt_dates = dt_dates[:-1]
                break

        # no faffing with Feb 29ths required
        if name not in ["Calama", "Ipascua"]:
            final_dates = np.array(dt_dates)
            final_obs = np.array(obs_data)
 
        # exclude 29th for non leap years
        elif name in ["Calama", "Ipascua"]:
            final_obs = []
            counter = 0
            for date in dt_dates:
                if date.month == 2 and date.day == 29:
                    final_obs += [utils.HADEX_MDI]
                else:
                    final_obs += [obs_data[counter]]
                    counter += 1

            final_dates = np.array(dt_dates) 
            final_obs = np.array(final_obs) 
            
        # check we've got an array of length that we expect!
        assert final_dates.shape[0] == final_obs.shape[0]

        # and store
        if suffix == "tmax":
            tmax.times = final_dates
            tmax.data = final_obs
        if suffix == "tmin":
            tmin.times = final_dates
            tmin.data = final_obs
        if suffix == "rr":
            prcp.times = final_dates
            prcp.data = final_obs

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_chile_station

#*********************************************
def read_colombia_station(loc, station, tmax, tmin, prcp, diagnostics=False):
    """
    Reads the Colombian station data into timeseries objects
     Only precipitation

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information (unused)
    """

    # test which type of station do we have!
    try:
        test = int(station.id)
        isID = True
    except ValueError:
        isID = False


    if isID:
        # blocks of years, rows are days, columns are months
        # best to do line by line and fill in
        infilename = os.path.join(station.location, "COLOMBIA DailyPrecip Data_{}.csv".format(station.id))

        with open(infilename, "r", encoding="latin-1") as infile:
            for lc, line in enumerate(infile):
                if lc == 0:
                    # skip header
                    pass
                else:
                    split = line.split(",")
                    year = int(split[2])
                    day = int(split[3])

                    # if missing month, then have empty entry
                    for s, entry in enumerate(split):
                        if entry == "" or entry == "\n":
                            split[s] = utils.HADEX_MDI
                    

                    # use first real entry to find start of data series
                    if lc == 1:
                        start = dt.date(year, 1, 1)
                        end = dt.date(2019, 1, 1)
                        
                        # make list of dates
                        dt_dates = [start]
                        while True:
                            dt_dates += [dt_dates[-1] + dt.timedelta(days=1)]
                            if dt_dates[-1] == end:
                                dt_dates = np.array(dt_dates[:-1])
                                data = np.zeros(dt_dates.shape[0])
                                data[:] = utils.HADEX_MDI
                                break

                    # work through columns, finding date and storing value
                    for month in range(1, 13):
                        try:
                            loc, = np.where(dt_dates == dt.date(year, month, day))
                            data[loc] = float(split[3 + month])
                        except ValueError:
                            print("{}/{}/{} isn't a valid date".format(year, month, day))
        times = dt_dates       
                        
    else:
        infiles = glob.glob(r'{}/COH2L*(*)-{}.csv'.format(station.location, station.id))

        if len(infiles) == 1:
            infilename = infiles[0]
            # open until first entry is a date (presume 1/1/x)
            with open(infilename, "r", encoding="latin-1") as infile:
                for lc, line in enumerate(infile):
                    if line[0] == "1":
                        break

            # use this to get header and skip appropriately
            indata = np.genfromtxt(infilename, dtype=(str), delimiter=",", skip_header=lc, encoding="latin-1")

            # convert to date time
            times = [dt.datetime.strptime(d_t, "%m/%d/%Y %H:%M") for d_t in indata[:, 0]]
            data = indata[:, 1]

        else:
            # single entry, set to missing
            times = np.array([dt.datetime(2000, 1, 1), dt.datetime(2001, 1, 1)])
            data = np.ma.array([utils.HADEX_MDI, utils.HADEX_MDI], mask=[True, True])

    
    prcp.times = times
    prcp.data = data
    
    # fake up other variables
    tmin.times = prcp.times
    tmin.data = np.zeros(len(prcp.data))
    tmin.data[:] = utils.HADEX_MDI

    tmax.times = prcp.times
    tmax.data = np.zeros(len(prcp.data))
    tmax.data[:] = utils.HADEX_MDI

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_colombia_station

#*********************************************
def read_canada_station(loc, station, tmax, tmin, prcp, diagnostics=False):
    """
    Reads the Canadian station data into timeseries objects

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information (unused)
    """
   
    #*********************
    # Temperature data

    # set up fieldwidths
    fieldwidths = [6, 3]
    for i in range(31):
        fieldwidths += [7, 1] # value and flag
    fieldwidths = tuple(fieldwidths)

    for parameter in ["x", "n"]:
        # extract all T data
        try:
            df = pandas.read_fwf(os.path.join(loc, "raw", "d{}{}.txt".format(parameter, station.id)), widths=fieldwidths, header=4)
            # convert to numpy
            indata = df.to_numpy()
            
            for fields in indata:
                year = int(fields[0])
                month = int(fields[1])

                # make times
                times = []
                n_days = calendar.monthrange(year, month)[1]
                for d in range(n_days):
                    times += [dt.date(year, month, d+1)]
                times = np.array(times)

                # get data (first column is a flag)
                data = np.ones(times.shape[0])
                data[:] = [float(entry) for entry in fields[2: 2+(2*n_days): 2]]
                data[data == -9999.9] = utils.HADEX_MDI

                # and store
                if parameter == "x":
                    tmax = utils.append_timeseries(tmax, times, data)
                if parameter == "n":
                    tmin = utils.append_timeseries(tmin, times, data)

        except IOError:
            print("File doesn't exist {}".format(os.path.join(loc, "raw", "d{}{}.txt".format(parameter, station.id))))
            # file doesn't exist
            # return empty lists
            pass


    #*********************
    # Precipitation data
    
    fieldwidths = [5, 3]
    for i in range(31):
        fieldwidths += [8, 1] # value and flag
    fieldwidths = tuple(fieldwidths)

    try:
        df = pandas.read_fwf(os.path.join(loc, "raw", "dt{}.txt".format(station.id)), widths=fieldwidths, header=1)
        # convert to numpy
        indata = df.to_numpy()
        
        for fields in indata:

            year = int(fields[0])
            month = int(fields[1])

            # make times
            times = []
            n_days = calendar.monthrange(year, month)[1]
            for d in range(n_days):
                times += [dt.date(year, month, d+1)]
            times = np.array(times)

            # get data (first column is a flag)
            data = np.ones(times.shape[0])
            data[:] = [float(entry) for entry in fields[2: 2+(2*n_days): 2]]
            data[data == -9999.99] = utils.HADEX_MDI

            prcp = utils.append_timeseries(prcp, times, data)

    except IOError:
        print("File doesn't exist {}".format(os.path.join(loc, "raw", "dt{}.txt".format(station.id))))
        # file doesn't exist
        # return empty lists
        pass

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_canada_station

#*********************************************
def read_decade_station(loc, station, tmax, tmin, prcp, diagnostics=False):
    """
    Reads the station data from the DECADE project into timeseries objects

    Some meta data mis-matches highlighted during development so cross checking on read

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information (unused)
    """
    import pandas as pd
    for parameter in ["TX", "TN", "PP"]:

        try:
            # get metadata for cross check
            with open(os.path.join(loc, "raw", parameter, "{}.dat".format(station.id)), "r", encoding="latin-1") as infile:
                for lc, line in enumerate(infile):
                    line = line.strip()
                    if lc == 4:
                        lat = float(line.split(":")[1])
                    if lc == 5:
                        lon = float(line.split(":")[1])

            # # now read the rest of the data
            # indata = np.genfromtxt(os.path.join(loc, "raw", parameter, "{}.dat".format(station.id)), skip_header=14, dtype=str, encoding="latin-1")

            # # extract and convert
            # years = indata[:, 1].astype(int)
            # months = indata[:, 2].astype(int)
            # days = indata[:, 3].astype(int)
            # data = indata[:, 4].astype(float)
            # data[data == -999.9] = utils.HADEX_MDI

            # now read the rest of the data
            DF = pd.read_csv(os.path.join(loc, "raw", parameter, "{}.dat".format(station.id)), skiprows=12, delim_whitespace=True)
            years = DF["Year"].to_numpy().astype(int)[1:]
            months = DF["Mo"].to_numpy().astype(int)[1:]
            days = DF["Da"].to_numpy().astype(int)[1:]
            if parameter in ["TX", "TN"]:
                data = DF["Param=3"].to_numpy().astype(float)[1:]
            elif parameter == "PP":
                data = DF["Param=1"].to_numpy().astype(float)[1:]
            data[data == -999.9] = utils.HADEX_MDI

            # extract Gross Inhomogeneity Frequency
            inhomog_rate = DF["gro.inh"].to_numpy().astype(float)[1:]
            # mask anything which isn't low
            bad_locs, = np.where(inhomog_rate > 1)
            data[bad_locs] = utils.HADEX_MDI
            
            # make times array
            times = np.array([dt.datetime(year, months[y], days[y]) for y, year in enumerate(years)])
            
            # if metadata match sufficiently, store
            if np.abs(station.latitude - lat) < 0.1 and np.abs(station.longitude - lon) < 0.1:

                if parameter == "TX":
                    tmax = utils.append_timeseries(tmax, times, data)
                if parameter == "TN":
                    tmin = utils.append_timeseries(tmin, times, data)
                if parameter == "PP":
                    prcp = utils.append_timeseries(prcp, times, data)
                
            else:
                print("Mismatch between metadata from inventory and in file for station {} parameter {}".format(station.id, parameter))

        except IOError:
            print("File doesn't exist {}".format(os.path.join(loc, "raw", parameter, "{}.dat".format(station.id))))
            # file doesn't exist
            # return empty lists
            pass

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_decade_station

#*********************************************
def read_eobs_station(loc, station, tmax, tmin, prcp, lookup, diagnostics=False):
    """
    Reads the station data from EOBS/LAOBS/SAOBS into timeseries objects

    Some meta data mis-matches highlighted during development so cross checking on read

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information (unused)
    """

    if "laobs" in station.location:
        HEADER = {"tx" : 20, "tn" : 20, "rr": 20}
    else:
        HEADER = {"tx" : 19, "tn" : 19, "rr": 19}
    INSCALE = 0.1

    for parameter in ["tx", "tn", "rr"]:

        if lookup[station.id][parameter] == -1:
            input("unexpected untested ID")
        elif lookup[station.id][parameter] == -99:
            print("No station ID for this key {} in {}".format(station.id, parameter))
            continue
        
        try:
            skip = HEADER[parameter]
#            if "eobs" in station.location and parameter == "rr":
#                skip += 1
            
            indata = np.genfromtxt(os.path.join(station.location, "{}_SOUID{}.txt".format(parameter.upper(), lookup[station.id][parameter])), delimiter=",", skip_header=skip, dtype=(str))
            
            # if only a single line
            if len(indata.shape) == 1:
                indata = np.expand_dims(indata, 0)

            # all data
            dates = indata[:, 2]
            value = indata[:, 3].astype(float)
            qc = indata[:, 4].astype(int)
            
            # fix scale & missing data
            MDI = -9999
            present, = np.where(value != MDI)
            value[present] = value[present] * INSCALE
            not_present, = np.where(value == MDI)
            value[not_present] = utils.HADEX_MDI

            # remove malformed dates
            good, = np.where(dates != "00000000")
            dates = dates[good]
            value = value[good]
            qc = qc[good]
            
            # convert dates into something useful
            dt_dates = np.array([dt.datetime(int(d[:4]), int(d[4:6]), int(d[6:])) for d in dates])
            
            # 0 = valid, 1 = suspect, 9 = missing
            value = np.ma.masked_where(qc != 0, value)
            times = np.ma.masked_where(qc != 0, dt_dates)
            value.fill_value = utils.HADEX_MDI

            if parameter == "tx":
                tmax = utils.append_timeseries(tmax, times, value)
            if parameter == "tn":
                tmin = utils.append_timeseries(tmin, times, value)
            if parameter == "rr":
                prcp = utils.append_timeseries(prcp, times, value)
                
        except IOError:
            print("File doesn't exist {}".format(os.path.join(station.location, "{}_SOUID{}.txt".format(parameter.upper(), lookup[station.id][parameter]))))
            # file doesn't exist
            # return empty lists
            pass

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_eobs_station

#*********************************************
def read_india_station(loc, station, tmax, tmin, prcp, diagnostics=False):
    """
    Reads the station data into timeseries objects 

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information
    """
   
    # all data from all stations in two files - one for P, one for Ts

    # just do longhand - there's time to go through each line.

    def id_match(line, stnid):

        match = False
        
        if line == stnid:
            match = True
        elif line.split("/")[0] == stnid:
            match = True
        elif line.split("-")[0] == stnid:
            match = True
        elif line == "CALCUTTA" and stnid == "KOLKATTA":
            match = True

        return match


    for parameter in ["T", "P"]:

        try:
            if parameter == "T":
                filename = "HADSET3_TEMP_40STN.PRT"

                if diagnostics:
                    print("temperatures")
                
                dates, txs, tns = [], [], []
                with open(os.path.join(loc, "raw", filename), encoding="latin-1") as infile:

                    read = False

                    for ll, line in enumerate(infile):

                        line = line.split()
                        # skip headings
                        if line[0] == "YEAR":
                            continue

                        # start or end section to read read
                        if len(line) == 5:
                            if id_match(line[2], station.id):
                                read = True
                            else:
                                # next section - turn reading off
                                read = False
                            continue

                        # extract data
                        if read:
                            dummy, year, month, day, tx, tn = line
                            try:
                                dates += [dt.datetime(int(year), int(month), int(day))]
                                txs += [float(tx)]
                                tns += [float(tn)]

                            except ValueError:
                                # unphysical dates
                                pass

                # and append timeseries
                tmax.data = np.array(txs)
                tmax.times = np.array(dates)
                tmin.data = np.array(tns)
                tmin.times = np.array(dates)
        
                # deal with missing values
                bads, = np.where(tmax.data == 99.9)
                if len(bads) > 0:
                    tmax.data[bads] = utils.HADEX_MDI
                bads, = np.where(tmin.data == 99.9)
                if len(bads) > 0:
                    tmin.data[bads] = utils.HADEX_MDI
                
            elif parameter == "P":
                filename = "DAILYRFDATA_INDIA.PRT"
                
                if diagnostics:
                    print("precipitation")
                
                dates, pps = [], []
                with open(os.path.join(loc, "raw", filename), encoding="latin-1") as infile:
                    read = False
                    for ll, line in enumerate(infile):
                        # skip header
                        if ll == 0: 
                            continue

                        line = line.split()

                        # start or end section to read read
                        if len(line) > 4:
                            if id_match(line[1], station.id):
                                read = True
                            else:
                                read = False
                            continue

                        # extract data
                        if read:
                            year, month, day, pp = line
                            try:
                                dates += [dt.datetime(int(year), int(month), int(day))]
                                pps += [float(pp)]
                            except ValueError:
                                # unphysical dates
                                pass

  
                # and append timeseries
                prcp.data = np.array(pps)
                prcp.times = np.array(dates)

                # deal with missing values
                bads, = np.where(prcp.data == -1.1)
                if len(bads) > 0:
                    prcp.data[bads] = utils.HADEX_MDI
                        
        except IOError:
            print("Files don't exist in {}".format(loc))
            # file doesn't exist
            # return empty lists
            pass

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_india_station


#*********************************************
def read_japan_station(loc, station, tmax, tmin, prcp, diagnostics=False):
    """
    Reads the station data into timeseries objects 

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information
    """

    infiles = os.listdir(os.path.join(loc, "raw", station.id))

    # temporary storage to sort all later
    txs = np.array([])
    tns = np.array([])
    pps = np.array([])
    dts = np.array([])


    for infilename in infiles:
        
        indata = np.genfromtxt(os.path.join(loc, "raw", station.id, infilename), skip_header=6, dtype=str, encoding="latin-1", delimiter=",")

        date = indata[:, 0]
        tx = indata[:, 1]
        txf = indata[:, 2]
        tn = indata[:, 4]
        tnf = indata[:, 5]
        pp = indata[:, 7]
        ppf = indata[:, 9]

        datetimes = [dt.datetime.strptime(d, "%Y/%m/%d") for d in date]

        # only retain good data
        bad = np.where(txf != "8")
        tx[bad] = utils.HADEX_MDI
        bad = np.where(tnf != "8")
        tn[bad] = utils.HADEX_MDI
        bad = np.where(ppf != "8")
        pp[bad] = utils.HADEX_MDI

        # append
        txs = np.append(txs, tx.astype(float))
        tns = np.append(tns, tn.astype(float))
        pps = np.append(pps, pp.astype(float))
        dts = np.append(dts, datetimes)
        
    # sort into ascending time order
    order = np.argsort(dts)
    dts = dts[order]
    txs = txs[order]
    tns = tns[order]
    pps = pps[order]
    
    # append to timeseries
    tmax = utils.append_timeseries(tmax, dts, txs)
    tmin = utils.append_timeseries(tmin, dts, tns)
    prcp = utils.append_timeseries(prcp, dts, pps)

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_japan_station

#*********************************************
def read_brazil_station(loc, station, tmax, tmin, prcp, diagnostics=False):
    """
    Reads the station data from Brazil into timeseries objects

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information (unused)
    """

    def fix_missing(data):

        for missing in ("", "-"):

            locs = np.argwhere(data == missing)
            if len(locs) > 0:
                data[locs] = str(utils.HADEX_MDI)

        return data

    infilenames = glob.glob(r'{}/DADOS_DIARIOS_*_{}_*.CSV'.format(station.location, station.id))

    if len(infilenames) == 1:

        indata = np.genfromtxt(infilenames[0], delimiter=";", skip_header=9, dtype=(str), encoding="latin-1")

        date = indata[:, 0]
        # sort dates
        times = [dt.datetime.strptime(d, "%d/%m/%Y") for d in date]

        for param in (1, 2, 3):

            data = indata[:, param]

            if len(data) == len(data[data == ""]):
                # entire series is missing
                
                data = np.ones(len(data)) * utils.HADEX_MDI

            else:
                # replace commas
                data = np.array([d.replace(",", ".") for d in data])
                
                # fix blank entries
                data = fix_missing(data)

                data = data.astype(float)

            # issues with data type and rounding of MDI
            data[data < -90.] = utils.HADEX_MDI

            # save timeseries
            if param == 1:
                prcp = utils.append_timeseries(prcp, times, data)
            if param == 2:
                tmax = utils.append_timeseries(tmax, times, data)
            if param == 3:
                tmin = utils.append_timeseries(tmin, times, data)
        
    elif len(infilenames) == 0:
        input("file unexpectedly missing")
    elif len(infilenames) > 1:
        input("multiple files matching")

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_brazil_station

#*********************************************
def pre_process_brazil_sp(dataset, diagnostics=False):
    """
    The Brazilian Sao Paulo data are in two file types:
    Agua Funda has 3 files, one for each variable
    A spreadsheet for other stations, one variable per column

    Need to split out into individual station files and create an inventory

    :param datasetObj dataset: dataset to process (HadEX2)
    :param bool diagnostics: output diagnostic information
    """
    #*********************
    def adjust_dates(param, dates, alldates):

        pad_param = np.ones(alldates.shape) * utils.HADEX_MDI

        match = np.in1d(alldates, dates)
        match_back = np.in1d(dates, alldates)

        pad_param[match] = param[match_back]

        return pad_param # adjust_dates

    import locale
    
    # make inventory
    metadata_file = open(os.path.join(dataset.location, "{}_stations.txt".format(dataset.name)), "w")

    # 1) Agua Funda

    for filename in sorted(glob.glob(r'{}/Agua Fun*'.format(os.path.join(dataset.location, "raw")))):
        print(filename)
        
        if filename[-7:] == "ain.txt":
            # Rainfall
            precip = []
            p_dates = []
            # read each line at a time
            with open(filename, "r", encoding="latin-1") as infile:
                for lc, line in enumerate(infile):

                    sline = line.split()
                    if sline[0] == "DATA":
                        continue

                    precip += [float(sline[-1].strip())]
                    p_dates += [dt.datetime.strptime(sline[0], "%d/%m/%Y")]

            p_dates = np.array(p_dates)
            precip = np.array(precip)
 
        elif filename[-7:] == "max.txt":
            # Tmax
            tmax = []
            x_dates = []
            # read each line at a time
            with open(filename, "r", encoding="latin-1") as infile:
                for lc, line in enumerate(infile):

                    sline = line.split(",")
                    if sline[0] == "\"DATA\"":
                        continue

                    tmax += [float(sline[-1].strip()[1:-1])]
                    x_dates += [dt.datetime.strptime(sline[0].split()[0][1:], "%Y-%m-%d")]
        
            x_dates = np.array(x_dates)
            tmax = np.array(tmax)
 
        elif filename[-7:] == "min.txt":
            # Tmin
            tmin = []
            n_dates = []
            # read each line at a time
            with open(filename, "r", encoding="latin-1") as infile:
                for lc, line in enumerate(infile):

                    sline = line.split(",")
                    if sline[0] == "\"DATA\"":
                        continue

                    tmin += [float(sline[-1].strip()[1:-1])]
                    n_dates += [dt.datetime.strptime(sline[0].split()[0][1:], "%Y-%m-%d")]

            n_dates = np.array(n_dates)
            tmin = np.array(tmin)
 
    # and write out
    start = np.min([n_dates[0], x_dates[0], p_dates[0]])
    end = np.max([n_dates[-1], x_dates[-1], p_dates[-1]])
    difference = end - start
    all_days = np.array([start + dt.timedelta(days=d) for d in range(difference.days + 1)])
    
    tmax = adjust_dates(tmax, x_dates, all_days)
    tmin = adjust_dates(tmin, n_dates, all_days)
    precip = adjust_dates(precip, p_dates, all_days)

    # ETCCDI format
    with open(os.path.join(dataset.location, "raw", "AguaFunda.txt"), "w") as outfile:
        for dc, date in enumerate(all_days):
            outfile.write("{:-8d}{:-8d}{:-8d}{:-8.1f}{:-8.1f}{:-8.1f}\n".format(date.year, date.month, date.day, precip[dc], tmax[dc], tmin[dc]))

    # from Jose Marengo email
    metadata_file.write("{:10s} {:7.3f} {:7.3f}\n".format("AguaFunda", -23.65, -46.61))


    # 2) big csv file.

    # Portuguese names in dates
    locale.setlocale(locale.LC_TIME, "pt_BR")

    filename = "SP-DADOS DIARIOS-1961-2019-Teste.csv"
    dates = []
    all_data = []
    with open(os.path.join(dataset.location, "raw", filename), "r", encoding="latin-1") as infile:
        for lc, line in enumerate(infile):
            sline = line.split("|")

            # process lines appropriately
            if sline[0] == "" or sline[0].split()[0] == "Instituto" or sline[0] == "Altitude (m)":
                continue
            elif sline[0] == "Latitude":
                latitudes = sline[1:]
                latitudes = [utils.dms2dd(lat[0:2], lat[2:4], lat[4:6], direction="S") for lat in latitudes]
            elif sline[0] == "Longitude":
                longitudes = sline[1:]
                longitudes = [utils.dms2dd(lon[0:2], lon[2:4], lon[4:6], direction=lon[6]) for lon in longitudes]
            elif sline[0] == "cod/nome":
                ids = sline[1:]
                ids = [name.split()[0] for name in ids]

            else:
                dates += [dt.datetime.strptime(sline[0], "%d-%b-%Y")]
                data = np.array(sline[1:])
                data[data == "NULL"] = str(utils.HADEX_MDI)
                data[data == "NULL\n"] = str(utils.HADEX_MDI)
                all_data += [data]

    # convert to arrays to use np tools
    all_data = np.array(all_data).astype(float)
    ids = np.array(ids)

    # now extract data in 3 column batches
    stats = np.unique(ids)
    for stn in stats:
        locs, = np.where(ids == stn)

        prcp = all_data[:, locs[0]]
        tx = all_data[:, locs[1]]
        tn = all_data[:, locs[2]]
        
        # ETCCDI format for data
        with open(os.path.join(dataset.location, "raw", "{}.txt".format(stn)), "w") as outfile:
            for dc, date in enumerate(dates):
                outfile.write("{:-8d}{:-8d}{:-8d}{:-8.1f}{:-8.1f}{:-8.1f}\n".format(date.year, date.month, date.day, prcp[dc], tx[dc], tn[dc]))

        # and the metadata
        metadata_file.write("{:10s} {:7.3f} {:7.3f}\n".format(stn, latitudes[locs[0]], longitudes[locs[0]]))    
    metadata_file.close()

    return # pre_process_brazil_sp

#*********************************************
def read_russia(loc, station, tmax, tmin, prcp, diagnostics=False):
    """
    Reads the updated Russian station data into timeseries objects

    :param str loc: location of station data
    :param stationObj station: station Object
    :param TimeseriesObj tmax: Tmax timeseries
    :param TimeseriesObj tmin: Tmin timeseries
    :param TimeseriesObj prcp: Precip timeseries
    :param bool diagnostics: output diagnostic information (unused)
    """

    # set up parseing info
    fieldwidths = [6, 5, 3, 3, 6, 2, 6, 2, 6, 2, 2]
    fieldwidths = tuple(fieldwidths)

    try:
        # read as pandas dataframe
        df = pandas.read_fwf(os.path.join(loc, "raw", "{}".format(station.id)), widths=fieldwidths, header=None)
        # convert to numpy
        indata = df.to_numpy()    

        # blank entries are "nan"
        indata[indata != indata] = utils.HADEX_MDI

        years = indata[:, 1].astype(int)
        months = indata[:, 2].astype(int)
        days = indata[:, 3].astype(int)

        # make times array            
        times = []
        bad_dates = []
        for y, year in enumerate(years):
            try:
                times += [dt.date(int(year), int(months[y]), int(days[y]))]
            except ValueError:
                # impossible_date
                bad_dates += [y]
                if diagnostics:
                    print("Station {}: bad date at {}-{}-{}".format(station.id, int(year), int(months[y]), int(days[y])))

        times = np.array(times)

        tn = indata[:, 4]
        tx = indata[:, 6]
        pp = indata[:, 8]
        
        #*********************
        # remove bad dates
        if len(bad_dates) > 0:
            times, tx, tn, pp = remove_bad_dates(bad_dates, times, tx, tn, pp)

        tmax = utils.append_timeseries(tmax, times, tx)
        tmin = utils.append_timeseries(tmin, times, tn)
        prcp = utils.append_timeseries(prcp, times, pp)

    except IOError:
            print("File doesn't exist {}".format(os.path.join(loc, "raw", "".format(parameter, station.id))))
            # file doesn't exist
            # return empty lists
            pass

    tmax, tmin, prcp = fill_missing_timeseries(tmax, tmin, prcp)

    return tmax, tmin, prcp # read_russia


#*******************************************
# END
#*******************************************
