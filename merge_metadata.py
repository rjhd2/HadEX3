#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 449                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2020-01-15 11:46:14 +0000 (Wed, 15 Jan 2020) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
'''
Merges index level station listings to single metadata file

merge_metadata.py invoked by typing::

  python merge_metadata.py --indata "ecad" --diagnostics

Input arguments:

--indata        Which input data to use (ecad, sacad, lacad)

--diagnostics   Output extra info (default False)
'''
import os
import re
import numpy as np

# RJHD utils
import utils
import inventory_utils

#*********************************************
def main(indata="ecad", diagnostics=False):
    """
    Find all metadata files for given dataset across all indices and merge together

    :param str indata: input dataset name
    :param bool diagnostics: output diagnostic information
    """

    # get all possible datasets
    all_datasets = utils.get_input_datasets()

    # and their names
    names = np.array([d.name for d in all_datasets])

    # if dataset selected and in the list of available, then run
    if indata in names:
        dataset = all_datasets[names == indata][0]

        all_stations = []
        all_names = []

        # spin through indices
        for index in utils.ALL_INDICES:
            
            if diagnostics:
                print("processing {}".format(index))

            # read in info
            if dataset.name in ["ecad", "sacad", "lacad"]:
                index_stations = inventory_utils.read_ecad(dataset, index, diagnostics=diagnostics)
            elif dataset.name in ["hadex2"]:
                index_stations = inventory_utils.read_hadex2(dataset, index, diagnostics=diagnostics)
            elif dataset.name in ["south_america"]:
                index_stations = inventory_utils.read_generic_index(dataset, index, diagnostics=diagnostics)
            elif dataset.name in ["west_africa_pptn"]:
                index_stations = inventory_utils.read_generic_index(dataset, index, diagnostics=diagnostics)
            elif dataset.name in ["west_africa_indices"]:
                index_stations = inventory_utils.read_generic_index(dataset, index, diagnostics=diagnostics)
            elif dataset.name in ["arabia"]:
                if utils.REF_START == 1961 and utils.REF_END == 1990:
                    index_stations = inventory_utils.read_arabia_6190_index(dataset, index, diagnostics=diagnostics)
                elif utils.REF_START == 1981 and utils.REF_END == 2010:
                    index_stations = inventory_utils.read_arabia_8110_index(dataset, index, diagnostics=diagnostics)

            elif dataset.name in ["south_africa"]:
                index_stations = inventory_utils.read_generic_index(dataset, index, diagnostics=diagnostics)
            elif dataset.name in ["ghcndex"]:
                index_stations = inventory_utils.read_ghcndex(dataset, index, diagnostics=diagnostics)

            # if no station metadata, then move on to next one
            if index_stations == []:
                continue
            
            # extract names
            station_names_for_index = [stn.id for stn in index_stations]

            # check if new
            for n, name in enumerate(station_names_for_index):
                if name in all_names:
                    # station exists already
                    loc = all_names.index(name)
                    # check all works out
                    try:
                        assert index_stations[n].latitude == all_stations[loc].latitude
                        assert index_stations[n].longitude == all_stations[loc].longitude
                    except AssertionError:
                        if (index_stations[n].latitude - all_stations[loc].latitude) > 0.1:
                            print("Station {} has mismatch in latitude: {} != {}".format(name, index_stations[n].latitude, all_stations[loc].latitude))
#                            sys.exit(1)
                        if (index_stations[n].longitude - all_stations[loc].longitude) > 0.1:
                            print("Station {} has mismatch in longitude: {} != {}".format(name, index_stations[n].longitude, all_stations[loc].longitude))
#                            sys.exit(1)

                else:
                    all_names += [name]
                    all_stations += [index_stations[n]]

        # sort alphabetically
        sort_order = np.argsort(np.array(all_names))

        all_stations = np.array(all_stations)[sort_order]

        # for some data sources there are index and raw data supplied
        #   For these, raw data produces the file, so for the index data, cross check  
        if os.path.exists(os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name))) and dataset.name in ["south_america", "south_africa"]:

            subset_stations = np.array([])

            # cross check all stations
            for this_station in all_stations:
                keep = True

                with open(os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name))) as infile:
                    for line in infile:
                        if re.search(this_station.id, line):
                            keep = False
                            if diagnostics:
                                print("{} {} already in metadata file".format(dataset.name, this_station.id))
                if keep:
                    subset_stations = np.append(subset_stations, this_station)

            all_station = subset_stations

        else:
            # write out the header
            utils.write_climpact_inventory_header(os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name)))

        # always prefer raw data over index data
        for this_station in all_stations:
            # write out the metadata for this station and index
            utils.write_climpact_inventory(os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name)), this_station)

    return # main

#************************************************************************
#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--indata', dest='indata', action='store', default="ecad",
                        help='Which dataset to convert')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default = False')

    args = parser.parse_args()         

    if args.indata in ["ecad", "lacad", "sacad", "hadex2", "south_america", "west_africa_pptn", \
                       "west_africa_indices", "arabia", "south_africa", "ghcndex"]:

        main(indata=args.indata, diagnostics=args.diagnostics)

    else:
        print("Nothing needs be done")
    # else just exit again


#*******************************************
# END
#*******************************************
