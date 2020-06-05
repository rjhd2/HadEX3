#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 460                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2020-01-31 12:27:09 +0000 (Fri, 31 Jan 2020) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
'''
Runs the gridding algorithm with the correct switches

Uses generates a land-sea mask from a 15min one to the box size specified
in the configuration file.

gridding.py invoked by typing::

  python gridding.py --grid "ADW" --index "TX90p" --hadex2_adw --month --qc_flags "" --anomalies --diagnostics

Input arguments:

--grid          Which gridding algorithm to use - ADW (Angular Distance Weighting) or CAM (Climate Anomaly Method)

--index         Which ETCCDI index to use

--hadex2_adw    Use form of ADW applied in HadEX2 (slight difference to formula)

--month         Month to run on.  Memory limits make it better to run each month separately.  0 - annual, 1-12 months

--qc_flags      Which QC flags to use

--anomalies     Uses anomalies for ADW

--diagnostics   Output extra info (default False)
'''
import os
import gc
import sys
import numpy as np

# RJHD utilities
import netcdf_procs as ncdfp
import utils

#*********************************************
#*******************************************
def read_data(station, index, timescale, nyears, nmonths):
    """
    Read in the data for a specific station

    :param station station: station to be read
    :param str index: index to be read
    :param str timescale: which timescale (MON/ANN)
    :param int nyears: number of years - to define array
    :param int nmonths: number of months - to define array
    """


    # need to read in all the data to be able to process
    data = np.ma.zeros([nyears, nmonths])
    data[:] = utils.HADEX_MDI
    data._FillValue = utils.HADEX_MDI
    
    times, indata = utils.read_station_index(station, index, timescale)
        
    # for climatology, which is single year
    if nyears == 1:
        data[:] = indata[:]

    else:
        # match the existing years and place in correct location in array
        match = np.in1d(utils.REFERENCEYEARS, times)
        match_back = np.in1d(times, utils.REFERENCEYEARS)

        if len(indata.shape) == 1:
            indata = indata.reshape([indata.shape[0], 1])

        data[match, :] = indata[match_back, :] # store all the info

    return data # read_data

#*******************************************
def get_all_data(station, index, timescale, nyears, month_index):
    """
    Read in the data for a specific station - and merge annual and monthly if appropriate

    :param station station: station to be read
    :param str index: index to be read
    :param str timescale: which timescale (MON/ANN)
    :param int nyears: number of years - to define array
    :param int month_index: which month to read   

    """


    if month_index == 0:
        # just read in the annual data - even if monthly index
        all_data = read_data(station, index, "ANN", nyears, 1)

    elif timescale == "MON" and month_index != 0:
        # if monthly index and monthly data, read in all data
        all_monthly_data = read_data(station, index, timescale, nyears, 12)
        # and select month in question
        all_data = all_monthly_data[:, month_index - 1] # offset by 1 as not a 13 column array

    # and apply the mask
    all_data = np.ma.masked_where(all_data == utils.HADEX_MDI, all_data)

    if len(all_data.shape) == 1: 
        all_data = np.ma.expand_dims(all_data, axis=1)


    return all_data # get_all_data

#*******************************************
def make_square(vector):
    """
    Take a 1-D mask array and make square
    with a numpy_and type of arangement, so mask False only where both x and y would be
    """

    mesh1, mesh2 = np.meshgrid(vector, vector)
    combination = mesh1 + mesh2

    new = np.ones(mesh1.shape)

    locs = np.where(combination == 0)
    new[locs] = 0

    return new.astype(bool) # make_square


#*******************************************
def set_up_grids(nyears, nmonths):
    """
    Set up empty masked arrays for the grids

    :param int nyears: number of years
    :param int nmonths: number of months
    """

    GridData = np.ma.zeros([nyears, nmonths, len(utils.box_centre_lats), len(utils.box_centre_lons)]) # data
    GridStations = np.ma.zeros([nyears, nmonths, len(utils.box_centre_lats), len(utils.box_centre_lons)]) # station count (in box)
    GridDLSStations = np.ma.zeros([nyears, nmonths, len(utils.box_centre_lats), len(utils.box_centre_lons)]) # station count (in DLS)
    
    GridData.fill = utils.HADEX_MDI
    GridStations.fill = utils.HADEX_MDI
    GridDLSStations.fill = utils.HADEX_MDI
    
    GridData.mask = np.ones(GridData.shape)
    GridStations.mask = np.ones(GridStations.shape)
    GridDLSStations.mask = np.ones(GridDLSStations.shape)

    return GridData, GridStations, GridDLSStations # set_up_grids

#*******************************************
def cam(all_datasets, index, timescale, nyears, qc_flags="", month_index=0, diagnostics=False, anomalies="None"):
    """
    Climate anomaly method gridding

    :param array all_datasets: array of dataset objects
    :param str index: which index to run
    :param str timescale: which timescale (MON/ANN)
    :param int nyears: number of years - to define array
    :param str qc_flags: which QC flags to process W, B, A, N, C, R, F, E, V, M
    :param int month_index: which month to read 
    :param bool diagnostics: output diagnostic information
    :param str anomalies: run code on anomalies or climatology rather than raw data

    """
    # change to do one month at a time to parallelise a little
    nmonths = 1
    print("Running Month {}".format(month_index))

    GridData, GridStations, dummy = set_up_grids(nyears, nmonths)

    # get the stations which actually have data
    # spin through all datasets
    stations = np.array([])
    for dataset in all_datasets:

        try:
            # choose appropriate subdirectory.
            if anomalies == "None":
                subdir = "formatted/indices"
            elif anomalies == "anomalies":
                subdir = "formatted/anomalies"
            elif anomalies == "climatology":
                subdir = "formatted/climatology"

            ds_stations = utils.read_inventory(dataset, subdir=subdir, final=True, timescale=timescale, index=index, anomalies=anomalies, qc_flags=qc_flags)

            good_stations = utils.select_qc_passes(ds_stations, qc_flags=qc_flags)
            
            stations = np.append(stations, good_stations)

            print("Adding {}, nstations = {}".format(dataset.name, len(stations)))

        except IOError:
            # file missing
            print("No stations with data for {}".format(dataset.name))

    # which lat and lon sequence does the station sit in.
    #   As will be using box centres, need to send list with same length as box_centres
    assign_stations_to_grid_boxes(stations, utils.box_edge_lats[1:], utils.box_edge_lons[1:])

    lon_sequence = np.array([stn.box_lon_sequence for stn in stations])
    lat_sequence = np.array([stn.box_lat_sequence for stn in stations])


    #*********************
    # run through each grid box
    for tlats, latitude in enumerate(utils.box_centre_lats):

        print(str(tlats)+"/"+str(len(utils.box_centre_lats)), latitude)

        lat_index, = np.where(lat_sequence == tlats)
        if len(lat_index) == 0:
            # no stations at this latitude so don't bother going any further
            continue  

        for tlons, longitude in enumerate(utils.box_centre_lons):

            lon_index, = np.where(lon_sequence == tlons)
            if len(lon_index) == 0:
                # no stations so don't bother
                continue

            # get the common stations to both selections (this lat and this lon sequence)
            lat_lon_match = np.intersect1d(lat_index, lon_index)
            if len(lat_lon_match) < 0:
                # no stations so don't bother
                continue

            # have at least one station in this grid box
            # go through each grid box

            box_data = np.ma.zeros([nyears, nmonths, len(lat_lon_match)])
            box_data.mask = np.ones(box_data.shape)

            for month in range(nmonths):

                print(len(lat_lon_match))
                for s, li in enumerate(lat_lon_match):
                    print(stations[li])

                    # read in the stations - use same routine as for ADW
                    data = get_all_data(stations[li], index, timescale, nyears, month_index)

                    if anomalies == "climatology":
                        # just read in to store
                        box_data[:, :, s] = data

                    else:
                        # calculate the anomalies

                        # back calculate times
                        good_times = utils.REFERENCEYEARS[data.mask[:, 0] == False]

                        # if no data then skip
                        if len(good_times) == 0: continue

                        #*********************
                        # anomalise
                        clim_years = np.where((utils.REFERENCEYEARS >= utils.CLIM_START.year) & (utils.REFERENCEYEARS < utils.CLIM_END.year), True, False) 

                        clim_data = data[clim_years]

                        #*********************
                        # check sufficient data points
                        completeness = np.ma.count(clim_data, axis=0)
                        locs, = np.where(completeness >= utils.CAM_COMPLETENESS)
                        if len(locs) == 0: continue

                        # single month at a time
                        climatology = np.ma.mean(clim_data, axis=0)

                        stn_anomalies = data - climatology

                        box_data[:, :, s] = stn_anomalies

                # done all stations in the box, take the mean
                GridData[:, month, tlats, tlons] = np.ma.mean(box_data, axis=-1)[:, 0]
                GridStations[:, month, tlats, tlons] = np.ma.count(box_data, axis=-1)[:, 0]
                # need at least N stations (default=3), so mask

                insufficient_stations = np.ma.where(GridStations[:, month, tlats, tlons] < utils.STATIONS_IN_BOX)
                GridData.mask[insufficient_stations, month, tlats, tlons] = True
                GridStations.mask[insufficient_stations, month, tlats, tlons] = True
    return GridData, GridStations # cam


#*******************************************
def adw(all_datasets, index, timescale, nyears, qc_flags="", month_index=0, diagnostics=False, hadex2_adw=False, anomalies="None"):
    """
    Angular Distance Weighting

    :param array all_datasets: array of dataset objects
    :param str index: which index to run
    :param str timescale: which timescale (MON/ANN)
    :param int nyears: number of years - to define array
    :param str qc_flags: which QC flags to process W, B, A, N, C, R, F, E
    :param int month_index: which month to read 
    :param bool diagnostics: output diagnostic information
    :param bool hadex2_adw: use the HadEX2 (erroneous) ADW method
    :param str anomalies: run code on anomalies or climatology rather than raw data

    """


    # http://journals.ametsoc.org/doi/pdf/10.1175/1520-0442%282000%29013%3C2217%3ARTCSTC%3E2.0.CO%3B2

    # change to do one month at a time to parallelise a little
    nmonths = 1
    loopwise = False
    print("Running Month {}".format(month_index))

    def calculate_cosine_term(stns_in_dls, box_angle, station_angle, hadex2_adw=False):
        """
        Helper routine to calculate the cosine term for the weighting function
        """

        if hadex2_adw:
            box_angle_k = station_angle[stns_in_dls][:, stns_in_dls]
            box_angle_i = np.tile(box_angle[stns_in_dls], (stns_in_dls.shape[0], 1))
        else:
            # do change as RJHD thinks - November 2017
            #   So that it is the angles between the stations and the box centre
            #   rather than the station and other stations.
            box_angle_i = np.tile(box_angle[stns_in_dls], (stns_in_dls.shape[0], 1))
            box_angle_k = box_angle_i.T
                            
        cosines = np.cos(box_angle_k - box_angle_i)
        box_angle_k = 0
        box_angle_i = 0

        return cosines # calculate_cosine_term

    def calculate_weighting_term(distance_weight, stns_in_dls):
        """
        Helper routine to calculate the weighting term
        """

        # tile the distance_weight array so replicated for all sid
        dist_weight_array = np.tile(distance_weight, (len(stns_in_dls), 1))
        dist_weight_array = np.ma.swapaxes(dist_weight_array, 0, 1)
                
        # set diagonal to zero (k != l)
        diag = np.arange(dist_weight_array.shape[-1])
        dist_weight_array[diag, diag] = 0.0 

        return dist_weight_array # calculate_weighting_term

    def calculate_top(dist_weight_array, cosines, nyears, mask):
        """
        Helper routine to calculate the top part of the weighting function
        """

        top_part = dist_weight_array * (1.0 - cosines)

        # now repeat this nyears times
        top_part = np.tile(top_part, (nyears, 1, 1))
        # as doing the sum, can just set masked elements to zero
        top_part[mask == True] = 0

        return np.sum(top_part, axis=-2) # calculate_top

    def calculate_bottom(dist_weight_array, nyears, mask):
        """
        Helper routine to calculate the bottom part of the weighting function
        """

        bottom_part = np.tile(dist_weight_array, (nyears, 1, 1))
        # as doing the sum, can just set masked elements to zero
        bottom_part[mask == True] = 0
        
        return np.sum(bottom_part, axis=-2) # calculate_bottom

    def calculate_adw(distance_weight, mask, top, bottom):
        """
        Helper routine to calculate the angular distance weights
        """

        distance_weight = np.tile(distance_weight, (nyears, 1))
        distance_weight[mask == True] = 0

        return distance_weight * (1 + (top/bottom)) # calculate_adw

    def calculate_separations_and_angles(stations, station_locs):
        """
        Calculate the station-station separation and bearing arrays
        """
        separation = np.zeros((stations.shape[0], stations.shape[0]))
        angle = np.copy(separation)

        for s, stn in enumerate(stations):
            this_stn = np.empty([len(stations), 2])
            this_stn[:, 0] = stn.latitude
            this_stn[:, 1] = stn.longitude
            separation[s, :], angle[s, :] = utils.map_2_points(this_stn, station_locs)

        return separation, angle


    #*******************************************
    
    # set up the grids
    GridData, GridStations, GridDLSStations = set_up_grids(nyears, nmonths)

    # get the DLS
    raw_dls = np.genfromtxt(os.path.join(utils.DLSLOCS, "dls_{}.txt".format(index)), dtype=(float), skip_header=4)

    dls_lat = raw_dls[:, 0]
    dls = raw_dls[:, 1:]

    # get the stations which actually have data
    # spin through all datasets
    stations = np.array([])
    for dataset in all_datasets:

        try:
            # choose appropriate subdirectory.
            if anomalies == "None":
                subdir = "formatted/indices"
            elif anomalies == "anomalies":
                subdir = "formatted/anomalies"
            elif anomalies == "climatology":
                subdir = "formatted/climatology"

            ds_stations = utils.read_inventory(dataset, subdir=subdir, final=True, timescale=timescale, index=index, anomalies=anomalies, qc_flags=qc_flags)

            good_stations = utils.select_qc_passes(ds_stations, qc_flags=qc_flags)
            
            stations = np.append(stations, good_stations)

            print("Adding {}, nstations = {}".format(dataset.name, len(stations)))

        except IOError:
            # file missing
            print("No stations with data for {}".format(dataset.name))

    # may have no stations for particular ETSCI combinations
    if len(stations) == 0:
        if diagnostics:
            print("No stations for {} - {}".format(index, timescale))
        return GridData, GridStations, GridDLSStations # adw

    station_locs = np.array([[stn.latitude, stn.longitude] for stn in stations])
    
    # get the distance and bearing arrays
    station_separation, station_angle = calculate_separations_and_angles(stations, station_locs)

    #*********************
    # read in all the data in one step for all the stations
    all_station_data = np.ma.zeros([nyears, nmonths, len(stations)])
    all_station_data.mask = np.ones(np.shape(all_station_data)) # mask everything
    latitudes = np.zeros(len(stations))
    longitudes = np.zeros(len(stations))
    

    # big and slow read loop
    for s, stat in enumerate(stations):
        data = get_all_data(stat, index, timescale, nyears, month_index)
        all_station_data[:, :, s] = data # store all the info
        all_station_data.mask[:, :, s] = data.mask # store all the info
        latitudes[s] = stat.latitude
        longitudes[s] = stat.longitude



    #*********************
    # run through each grid box
    for tlats, latitude in enumerate(utils.box_centre_lats):

        print(str(tlats)+"/"+str(len(utils.box_centre_lats)), latitude)

        for tlons, longitude in enumerate(utils.box_centre_lons):

            # distance of this box centre to all stations
            this_box = np.empty([len(stations), 2])
            this_box[:, 0] = latitude
            this_box[:, 1] = longitude
            box_separation, box_angle = utils.map_2_points(this_box, station_locs)

            # find those stations close enough to contribute
            # need to adjust if doing all months so that can read in all, but restrict to relevant ones if necessary
            stns_in_dls, = np.where(box_separation <= np.max(dls[tlats]))
            stns_in_dls_separations = box_separation[stns_in_dls]

            if len(stns_in_dls) < utils.STATIONS_IN_DLS: 
                # none of the months have DLS such that sufficient stations are included
                #   skip to next box
                if diagnostics: print("skipping lat {}, lon {} - no stations in range (max DLS = {})".format(latitude, longitude, np.max(dls[tlats])))
                continue

# REMOVED FOR SINGLE READ
            # set up blank array to store all station data that can contribute
#            stations_contrib_to_box_data = np.ma.zeros([nyears, nmonths, len(stns_in_dls)])
#            stations_contrib_to_box_data.mask = np.ones(np.shape(stations_contrib_to_box_data)) # mask everything
# REMOVED FOR SINGLE READ

            # this is for the stations actually within the grid box!
            stations_in_box = np.zeros([nyears, nmonths])                    

            print(" nstats {}".format(len(stns_in_dls)))

            # get the stations contributing to the box (i.e. within a DLS)
            stations_contrib_to_box_data = np.ma.copy(all_station_data[:, :, stns_in_dls])
            stations_in_box = np.ma.count(all_station_data[:, :, stns_in_dls], axis=2)

            # or get the stations located in the box
            # lat_locs, = np.where(np.logical_and(utils.box_edge_lats[tlats] < latitudes, latitudes <= utils.box_edge_lats[tlats+1]))
            # lon_locs, = np.where(np.logical_and(utils.box_edge_lons[tlons] < longitudes, longitudes <= utils.box_edge_lons[tlons+1]))
            # # station matches both latitude and longitude constraints
            # both_lat_and_lon = np.in1d(lat_locs, lon_locs)
            # in_box_locs = lat_locs[both_lat_and_lon]
            # if len(in_box_locs) > 0:
            #     stations_in_box = np.ma.count(all_station_data[:, :, in_box_locs], axis=2)
            
# REMOVED FOR SINGLE READ
#             # read in all the stations - and do this once
#             for s, li in enumerate(stns_in_dls):
# #                if diagnostics: 
# #                    print(stations[li], stns_in_dls_separations[s])

#                 # read in the station - matching done in subroutine
#                 data = get_all_data(stations[li], index, timescale, nyears, month_index) 

#                 stations_contrib_to_box_data[:, :, s] = data # store all the info
#                 stations_contrib_to_box_data.mask[:, :, s] = data.mask # store all the info

#                 # and number of stations in the box
#                 # need to subtract if not present at any year or with smaller DLS
#                 if (utils.box_edge_lats[tlats] < stations[li].latitude <= utils.box_edge_lats[tlats+1]) \
#                         and (utils.box_edge_lons[tlons] < stations[li].longitude <= utils.box_edge_lons[tlons+1]):
#                     stations_in_box[data.mask == False] += 1
# REMOVED FOR SINGLE READ

            # go through each month label - not used as parallelised instead
            for month in range(nmonths):
                
                if diagnostics:
                    print("latitude {} ({}), longitude {} ({}), dls {}, nstations {}".format(latitude, tlats, longitude, tlons, dls[tlats][month_index], len(np.where(box_separation[stns_in_dls] < dls[tlats][month_index])[0])))

                # get weights
                distance_weight = np.exp(utils.M * -box_separation[stns_in_dls]/dls[tlats][month_index]) # for each contributing station

                # filter out stations too far away for this month DLS
                sep_locs, = np.where(stns_in_dls_separations > dls[tlats][month_index])
                stations_contrib_to_box_data[:, month, sep_locs] = 0
                stations_contrib_to_box_data.mask[:, month, sep_locs] = True 

                if loopwise:
                    pass
#                     # this is the original longhand version.  
#                     for year in range(nyears):

#                         if np.ma.count(stations_contrib_to_box_data[:, month], axis=1)[year] < utils.STATIONS_IN_DLS:
#                             # insufficient stations - don't bother
#                             continue

#                         # which stations do contribute
#                         this_year_mask = -stations_contrib_to_box_data.mask[year, month]

#                         if diagnostics:
#                             # testing long-hand looping
#                             # using Caesar et al terminology - http://onlinelibrary.wiley.com/doi/10.1029/2005JD006280/pdf
#                             w_is=[]
#                             for i in stns_in_dls[this_year_mask]:
#                                 w_i = np.exp(utils.M * -box_separation[i]/dls[tlats][month_index])

#                                 tops=[]
#                                 bottoms=[]

#                                 for k in stns_in_dls[this_year_mask]: # all other stations
#                                     if i != k:

#                                         w_k = np.exp(utils.M * -box_separation[k]/dls[tlats][month_index])
                                        
#                                         bottoms += [w_k]
#                                         if hadex2_adw:
#                                             tops += [w_k * (1.0 - np.cos(station_angle[k,i] - box_angle[i]))]
#                                         else:
#                                             tops += [w_k * (1.0 - np.cos(box_angle[k] - box_angle[i]))]
                                
#                                 w_is += [w_i * (1 + np.sum(tops)/np.sum(bottoms))]
# #                            print("weights", w_is/sum(w_is))

#                         # tile the distance_weight array so replicated for all sid
#                         dist_weight_array = np.tile(distance_weight[this_year_mask], (np.ma.count(stations_contrib_to_box_data[year, month]), 1))
#                         dist_weight_array = np.swapaxes(dist_weight_array, 0, 1)

#                         # set diagonal to zero (k != l)
#                         diag = np.arange(dist_weight_array.shape[0])
#                         dist_weight_array[diag, diag] = 0.0

#                         # make a mesh of these so can subtract.
#                         box_angle_i = np.tile(box_angle[stns_in_dls][this_year_mask],(stns_in_dls[this_year_mask].shape[0],1))

#                         if hadex2_adw:
#                             box_angle_k = station_angle[stns_in_dls[this_year_mask]][:,stns_in_dls[this_year_mask]]
#                         else:
#                             box_angle_k = box_angle_i.T

#                         # get array for top
#                         top_sum = dist_weight_array * (1.0 - np.cos(box_angle_k - box_angle_i))
                        
#                         top = np.sum(top_sum, axis=0)
#                         bottom = np.sum(dist_weight_array, axis=0)

#                         # un-normalised weights
#                         angular_distance_weight = distance_weight[this_year_mask] * (1 + top / bottom)

#                         final_weights = angular_distance_weight/np.ma.sum(angular_distance_weight)
                        
# #                        if diagnostics: print("weights", final_weights)  # should match print line above)

#                         GridData[year, month, tlats, tlons] = np.ma.sum(final_weights * stations_contrib_to_box_data[year, month].compressed())

#                         GridStations[year, month, tlats, tlons] = stations_in_box[year, month]
# #                        if diagnostics: 
# #                            print(GridData[year, month, tlats, tlons])
# #                            raw_input("stop {}".format(hadex2_adw))

                else: # not loopwise
                    
                    """
                    W_k = weight_k * (1 + a_k)
                    
                    a_k = top/bottom
                    
                    bottom = sum_1_nstations(w_k)
                    
                    top = sum_1_nstations(w_k * (1 - cos(theta_k - theta_l))),  k != l
                    """
                    
                    # aim to remove this loop (longer than months loop, so saves more time?)
                    # if sufficient stations
                    insufficient_station_count = np.ma.count(stations_contrib_to_box_data[:, month, :], axis=1)
                    
                    if max(insufficient_station_count) < utils.STATIONS_IN_DLS:
                        # no year has sufficient stations
                        if diagnostics: print("skipping lat {}, lon {}, month {} - no stations in range (DLS = {})".format(latitude, longitude, month+1, dls[tlats]))
                        continue

                    this_month_mask = stations_contrib_to_box_data.mask[:, month, :]
                    mask = np.array([make_square(m) for m in this_month_mask])


                    # calculate angular part
                    cosines = calculate_cosine_term(stns_in_dls, box_angle, station_angle, hadex2_adw=hadex2_adw)
                    
                    # calculate weight part
                    dist_weight_array = calculate_weighting_term(distance_weight, stns_in_dls)
                    
                    # cosine and dist_weight the same each year - doesn't change
                    # calculate the top without years, clear memory, then expand and mask
                    top = calculate_top(dist_weight_array, cosines, nyears, mask)
                    cosines = 0

                    bottom = calculate_bottom(dist_weight_array, nyears, mask)
                    dist_weight_array = 0
                    
                    angular_distance_weight = calculate_adw(distance_weight, this_month_mask, top, bottom) 
                    top = 0
                    bottom = 0
                    
                    normalisation = np.ma.sum(angular_distance_weight, axis=1)
                    final_weights = angular_distance_weight/normalisation[:, None]

                    normalisation = 0
                    angular_distance_weight = 0

#                    if diagnostics: print(final_weights)

                    GridData[:, month, tlats, tlons] = np.ma.sum(final_weights * stations_contrib_to_box_data[:, month, :], axis=1)
                    GridData.mask[insufficient_station_count < utils.STATIONS_IN_DLS, month, tlats, tlons] = True
                    GridStations[:, month, tlats, tlons] = stations_in_box[:, month]
                    GridDLSStations[:, month, tlats, tlons] = insufficient_station_count

#                    if diagnostics: 
#                        print(np.max(insufficient_station_count))
#                        print(GridData[:, month, tlats, tlons])
#                        raw_input("stop {}".format(hadex2_adw))

                gc.collect()
                sys.stdout.flush()
    return GridData, GridStations, GridDLSStations # adw

#*******************************************
def rsm():
    
    return # rsm

#*******************************************
def fdm():
    
    return # fdm


#*******************************************
def assign_stations_to_grid_boxes(stations, lats, lons):
    '''
    assign stations to relevant grid box

    this is a relatively slow step, so only want to do this once.
    '''

    # need to fix longitudes to be ~0-360
    longitudes = np.array([s.longitude for s in stations])
    longitudes[longitudes < utils.box_edge_lons[0]] = longitudes[longitudes < utils.box_edge_lons[0]] + 360

    for s, station in enumerate(stations):

        # https://stackoverflow.com/questions/2236906/first-python-list-index-greater-than-x
        station.box_lon_sequence = next(i for i, v in enumerate(lons) if v > longitudes[s])
        station.box_lat_sequence = next(i for i, v in enumerate(lats) if v > station.latitude)

    return

#*********************************************
def main(grid="ADW", index="TX90p", month_index=0, diagnostics=False, hadex2_adw=False, qc_flags="", anomalies="None"):
    """
    :param str grid: gridding type ADW/CAM
    :param str index: which index to run
    :param int month_index: which month to apply (0 = Annual, 1-12 for months)    
    :param str qc_flags: which QC flags to process W, B, A, N, C, R
    :param bool diagnostics: output diagnostic information
    :param str anomalies: run code on anomalies or climatology rather than raw data
    """

    # ensure correct timescale is selected
    if args.index in utils.MONTHLY_INDICES:
        if month_index == 0:
            timescale = "ANN"
        else:
            timescale = "MON"
    else:
        if month_index == 0:
            timescale = "ANN"
        else:
            print("Monthly requested for annual-only index.\n Exiting")
            return
    

    # move this up one level eventually
    all_datasets = utils.get_input_datasets()

    # set up the data arrays
    if anomalies == "climatology":
        nyears = 1
    else:
        nyears = len(utils.REFERENCEYEARS)

    if grid == "CAM":
        GridData, GridStations = cam(all_datasets, index, timescale, nyears, qc_flags=qc_flags, month_index=month_index, diagnostics=diagnostics, anomalies=anomalies)

    elif grid == "ADW":
        GridData, GridStations, GridDLSStations = adw(all_datasets, index, timescale, nyears, qc_flags=qc_flags, month_index=month_index, diagnostics=diagnostics, anomalies=anomalies, hadex2_adw=hadex2_adw)


    if utils.DOLSM:
        nmonths = 1
        # apply LSM
        lsm = utils.get_land_sea_mask(utils.box_centre_lats, utils.box_centre_lons, floor=False) # not taking only purely non-land boxes.  Have to have sufficient amount of land!
        # resize to match
        lsm_sized = np.tile(np.tile(lsm, (1, 1, 1, 1)), (nyears, nmonths, 1, 1))

        GridData.mask = np.logical_or(GridData.mask, lsm_sized)
        GridStations.mask = np.logical_or(GridStations.mask, lsm_sized)
        if grid == "ADW":
            GridDLSStations.mask = np.logical_or(GridDLSStations.mask, lsm_sized)

    # correct fill_value
    GridData.fill_value = utils.HADEX_MDI
    GridStations.fill_value = utils.HADEX_MDI

    # append appropriate name to filename if anomalies or climatology
      
    filename = utils.make_filenames(index=index, grid=grid, anomalies=anomalies, month_index=month_index)
    
    ncdfp.netcdf_write(os.path.join(utils.OUTROOT, filename), index, GridData.filled(), utils.REFERENCEYEARS, utils.box_centre_lats, utils.box_centre_lons, single_month=month_index)
    
    filename = utils.make_filenames(index=index, grid=grid, anomalies=anomalies, extra="num", month_index=month_index)
    ncdfp.netcdf_write(os.path.join(utils.OUTROOT, filename), index, GridStations.filled(), utils.REFERENCEYEARS, utils.box_centre_lats, utils.box_centre_lons, single_month=month_index, station_count=True)

    if grid == "ADW":
        
        filename = utils.make_filenames(index=index, grid=grid, anomalies=anomalies, extra="numdls", month_index=month_index)
        ncdfp.netcdf_write(os.path.join(utils.OUTROOT, filename), index, GridDLSStations.filled(), utils.REFERENCEYEARS, utils.box_centre_lats, utils.box_centre_lons, single_month=month_index, station_count=True)



    return # main


#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--month', dest='month', action='store', default=0,
                        help='Month index to use 0 = Annual, 1-12 = months', type=int)
    parser.add_argument('--index', dest='index', action='store', default="TX90p",
                        help='Which index to run')
    parser.add_argument('--grid', dest='grid', action='store', default="ADW",
                        help='gridding routine to use - ADW or CAM')
    parser.add_argument('--hadex2_adw', dest='hadex2_adw', action='store_true', default=False,
                        help='Use ADW method as in HadEX2 (angles), default=False')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default=False')
    parser.add_argument('--qc_flags', dest='qc_flags', action='store', default="",
                        help='Which QC flags to use when filtering stations, default=""')
    parser.add_argument('--anomalies', dest='anomalies', action='store', default="None",
                        help='To use the anomalies or climatology rather than actuals, default="None"') 

    args = parser.parse_args()

    # if only annual index and not annual selected, then exit.
    if args.index not in utils.MONTHLY_INDICES:
        if args.month != 0:
            print("specifying month for annual only index: {}".format(args.index))
            sys.exit()

    # if month not possible, then exit.
    if args.month > 12:
        print("impossible month: {}".format(args.month))
        sys.exit()
            

    main(grid=args.grid, index=args.index, hadex2_adw=args.hadex2_adw, month_index=args.month, diagnostics=args.diagnostics, qc_flags=args.qc_flags, anomalies=args.anomalies)


#------------------------------------------------------------
# END
#------------------------------------------------------------
