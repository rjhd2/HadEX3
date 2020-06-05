#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 398                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-08-28 14:25:51 +0100 (Wed, 28 Aug 2019) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""

Runs the Rscript climpact2 with the correct switches.

"""
import os
import subprocess
import shutil
import numpy as np

import utils


#*********************************************
def main(indata="ghcnd", diagnostics=False):
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

        '''
        Process call structure

        climpact2.batch.stations.r ./sample_data/ ./sample_data/climpact2.sample.batch.metadata.txt 1971 2000 4
        '''

        # check that there are stations to process for this dataset
        stations = utils.read_inventory(dataset)
        if len(stations) != 0:

            try:
                with utils.cd(utils.CLIMPACT_LOCS):
                    # call the R process - which should automatically do everything and make suitable files etc
                    #  runs in subfolder with context manager, so returning to parent once done.

                    # ACRE (and others?) have stations that do not overlap the reference period.
                    #   Means that the QC process throws them out if insufficient overlap between data and reference period

                    if dataset.name == "acre":
                        ref_start = 1901
                        ref_end = 1930
                    else:
                        ref_start = utils.REF_START
                        ref_end = utils.REF_END


                    print(" ".join(["Rscript", "climpact2.batch.stations.r", os.path.join(dataset.location, "formatted"), os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name)), str(ref_start), str(ref_end), str(utils.NCORES)]))

                    subprocess.check_call(["Rscript", "climpact2.batch.stations.r", os.path.join(dataset.location, "formatted"), os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name)), str(ref_start), str(ref_end), str(utils.NCORES)])

            except subprocess.CalledProcessError:
                # handle errors in the called executable
                raise Exception
                 
            except OSError:
                # executable not found
                print("Cannot find Rscript")
                raise OSError

        # fail gracefully
        else:
            print("No stations available in {}".format(indata))
            print("  Climpact2 not run")

        # remove plots, qc, thres and trend folders (save space)
        if utils.REMOVE_EXTRA:
            for subdir in ["plots", "qc", "thres", "trend"]:
                try:
                    shutil.rmtree(os.path.join(dataset.location, "formatted", subdir))
                except FileNotFoundError:
                    print("{} doesn't exist".format(os.path.join(dataset.location, "formatted", subdir)))


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

    if args.indata in ["ghcnd", "spain", "nz", "honduras", "russia", "south_america", "australia", "chile", "colombia", \
                       "canada", "decade", "eobs", "laobs", "saobs", "india", "japan", "mexico", "iran", "brazil", \
                       "brazil_sp"]:

        main(indata=args.indata, diagnostics=args.diagnostics)


#*******************************************
# END
#*******************************************
