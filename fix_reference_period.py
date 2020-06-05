#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 177                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2018-09-26 15:46:03 +0100 (Wed, 26 Sep 2018) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""
Removes or adjusts files where the reference period wasn't correct.

Needed for e.g. ACRE to get Climpact to run.  Removing files shouldn't affect station selection

fix_reference_period.py invoked by typing::

  python fix_reference_period.py --indata "ecad" --diagnostics

Input arguments:

--indata        Which input data to use (ECAD the only sensible option)

--diagnostics   Output extra info (default False)
"""
import os
import numpy as np

# RJHD
import utils

# I think this is exhaustive :)
PERCENTILE_INDICES = ["TX90p", "TX10p", "TN90p", "TN10p", "WSDI", "CSDI", "R95p", "R99p", "R95pTOT", "R99pTOT"]

#*********************************************
def process_dataset(dataset, diagnostics=False):
    """
    Read in the specified dataset inventory.  Check each station for suitability for HadEX3

    :param datasetObj dataset: dataset object holding metadata about the input dataset
    :param bool diagnostics: output diagnostic information
    """

    # read in the dataset inventory metadata
    dataset_stations = utils.read_inventory(dataset, subdir="formatted/indices")

    # spin through stations
    for station in dataset_stations:

        # spin through indices
        for index in PERCENTILE_INDICES:

            # select appropriate timescales
            if index in utils.MONTHLY_INDICES:
                timescales = ["ANN", "MON"]
            else:
                timescales = ["ANN"]

            # and spin through those
            for timescale in timescales:
                if os.path.exists(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, index.lower(), timescale))):
                    os.remove(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, index.lower(), timescale)))

                    if diagnostics:
                        print("Removing {}".format(os.path.join(station.location, station.id, "{}_{}_{}.csv".format(station.id, index.lower(), timescale))))

    return # process_dataset

#*********************************************
def main(indata="acre", diagnostics=False):
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

        dataset = all_datasets[names == indata][0]

        if dataset.name in ["acre"]:
            # need to run ACRE as base period set separately
            process_dataset(dataset, diagnostics=diagnostics)

        elif dataset.base_period == "00-00":
            # the input dataset is raw observations, run by climpact, so no need
            #    as will have matched HadEX3 base period
            pass
        elif utils.match_reference_period(dataset.base_period):
            # the input datasets's reference period matches that of HadEX3 version
            pass
        else:
            # base period of input dataset doesn't match HadEX3, so remove
            #    appropriate indices
            process_dataset(dataset, diagnostics=diagnostics)

    return # main

#************************************************************************
#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--indata', dest='indata', action='store', default="acre",
                        help='Which dataset to process')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default = False')
 
    args = parser.parse_args()         

    if args.indata in ["ghcnd", "ecad", "sacad", "lacad", "spain", "pacific", "nz", "honduras", "hadex2", \
                       "russia", "acre", "south_america", "west_africa_pptn", "west_africa_indices", \
                       "australia", "arabia", "chile", "colombia", "canada", "decade", \
                       "south_africa", "eobs", "laobs", "saobs", "myanmar", "malaysia", "brunei", "india", \
                       "china", "japan", "mexico", "iran", "brazil", "vietnam", "philippines", "singapore", \
                       "brazil_sp", "ghcndex", "thailand", "indonesia"]:
        main(indata=args.indata, diagnostics=args.diagnostics)
    else:
        print("Adjustments not needed for {}".format(args.indata))

#------------------------------------------------------------
# END
#------------------------------------------------------------
