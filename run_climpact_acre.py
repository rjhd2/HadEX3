#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 332                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-04-16 16:06:06 +0100 (Tue, 16 Apr 2019) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""

Runs the Rscript climpact2 with the correct switches.
Setup for ACRE (and other early data) is to use a temporary metadata file
to ensure the base period covers the period of record
These should never be used for the percentile indices (fix_reference_period.py addresses this)

"""
import os
import subprocess
import shutil
import numpy as np

import utils


#*********************************************
def main(indata="acre", diagnostics=False):
    """
    Call the R package climpact2 with appropriate settings to calculate the indices

    :param str indata: name of dataset to process
    :param bool diagnostics: output diagnostic information
    """


    # get all possible datasets
    all_datasets = utils.get_input_datasets()

    # and their names
    names = np.array([d.name for d in all_datasets])

    # select the matching one
    if indata in names:
        dataset = all_datasets[names == indata][0]

        # check that there are stations to process for this dataset
        stations = utils.read_inventory(dataset)
        if len(stations) != 0:

            for station in stations:

                # read the station data
                infile = os.path.join(dataset.location, "formatted", "{}.txt".format(station.id))
                indata = np.genfromtxt(infile)

                # get the first year and last year
                ref_start = int(indata[0][0])
                ref_end = int(indata[-1][0])

                # write a temporary inventory file for just this station
                utils.write_climpact_inventory_header(os.path.join(dataset.location, "{}_temp.metadata.txt".format(dataset.name)))
                utils.write_climpact_inventory(os.path.join(dataset.location, "{}_temp.metadata.txt".format(dataset.name)), station)

                try:
                    with utils.cd(utils.CLIMPACT_LOCS):
                        # call the R process - which should automatically do everything and make suitable files etc
                        #  runs in subfolder with context manager, so returning to parent once done.

                        # ACRE (and others?) have stations that do not overlap the reference period.
                        #   Means that the QC process throws them out if insufficient overlap between data and reference period

                        print(" ".join(["Rscript", "climpact2.batch.stations.r", os.path.join(dataset.location, "formatted"), os.path.join(dataset.location, "{}_temp.metadata.txt".format(dataset.name)), str(ref_start), str(ref_end), str(utils.NCORES)]))

                        subprocess.check_call(["Rscript", "climpact2.batch.stations.r", os.path.join(dataset.location, "formatted"), os.path.join(dataset.location, "{}_temp.metadata.txt".format(dataset.name)), str(ref_start), str(ref_end), str(utils.NCORES)])

                except subprocess.CalledProcessError:
                    # handle errors in the called executable
                    raise Exception

                except OSError:
                    # executable not found
                    print("Cannot find Rscript")
                    raise OSError

                # remove temporary metadata file
                os.remove(os.path.join(dataset.location, "{}_temp.metadata.txt".format(dataset.name)))

        # fail gracefully
        else:
            print("No stations available in {}".format(indata))
            print("  Climpact2 not run")

        # remove plots, qc, thres and trend folders (save space)
        if utils.REMOVE_EXTRA:
            for subdir in ["plots", "qc", "thres", "trend"]:
                shutil.rmtree(os.path.join(dataset.location, "formatted", subdir))


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
    parser.add_argument('--indata', dest='indata', action='store', default="acre",
                        help='Which dataset to convert')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default = False')

    args = parser.parse_args()         

    if args.indata in ["acre"]:

        main(indata=args.indata, diagnostics=args.diagnostics)


#*******************************************
# END
#*******************************************
