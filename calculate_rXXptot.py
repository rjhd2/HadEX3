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
Calculates the RXXpTOT indices if files not present.

calculate_rXXptot.py invoked by typing::

  python calcuate_rXXptot.py --indata "ghcndex" --diagnostics

Input arguments:

--indata        Which input data to use ("west_africa_pptn", "ghcndex", "south_africa")

--diagnostics   Output extra info (default False)
"""
import os
import numpy as np

# RJHD
import utils

PARTNERS = {"R95pTOT" : "R95p", "R99pTOT" : "R99p"}

#*********************************************
def main(indata="ghcndex", index="R95pTOT", diagnostics=False):
    """
    Read PRCPTOT and other indices and write out


    """

    # check if need to do monthly ones
    if index in utils.MONTHLY_INDICES:
        timescales = ["ANN", "MON"]
    else:
        timescales = ["ANN"]
    
    # get all possible datasets
    all_datasets = utils.get_input_datasets()
    # and their names
    names = np.array([d.name for d in all_datasets])

    # if dataset selected and in the list of available, then run
    if indata in names:
        dataset = all_datasets[names == indata][0]

        dataset_stations = utils.read_inventory(dataset, subdir="formatted/indices")

        # check each station
        for stn in dataset_stations:

            if diagnostics:
                print("{} - {}".format(dataset.name, stn.id))

            # for appropriate number of timescales
            for ts in timescales:

                if os.path.exists(os.path.join(stn.location, stn.id, "{}_{}_{}.csv".format(stn.id, PARTNERS[index].lower(), ts))) and os.path.exists(os.path.join(stn.location, stn.id, "{}_{}_{}.csv".format(stn.id, "PRCPTOT".lower(), ts))):
                    rtimes, rXXp = utils.read_station_index(stn, PARTNERS[index].lower(), ts)
                    ptimes, prcptot = utils.read_station_index(stn, "PRCPTOT".lower(), ts)

                    match = np.in1d(rtimes, ptimes)
                    match_b = np.in1d(ptimes, rtimes)

                    if len(match) != 0 and len(match_b) != 0:

                        rXXptot = (100 * rXXp) / prcptot
                        rXXptot_times = rtimes[match]

                        if ts == "MON":
                            myears = []
                            months = []
                            for y in rXXptot_times:
                                for m in range(1, 13):
                                    myears += [y]
                                    months += [m]                    

                            stn.monthly = rXXptot.filled().reshape(-1)
                            stn.myears = myears
                            stn.months = months 
                            path = os.path.join(dataset.location, "formatted", "indices", stn.id, "{}_{}_MON.csv".format(stn.id, index.lower()))
                            if not os.path.exists(path):
                                utils.write_station_index(path, stn, index, doMonthly=True)

                        else:            
                            stn.years = rXXptot_times
                            stn.annual = rXXptot.filled()
                            path = os.path.join(dataset.location, "formatted", "indices", stn.id, "{}_{}_ANN.csv".format(stn.id, index.lower()))
                            if not os.path.exists(path):
                                utils.write_station_index(path, stn, index)
        
    return # main


#************************************************************************
#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--indata', dest='indata', action='store', default="ghcndex",
                        help='Which dataset to process')
    parser.add_argument('--index', dest='index', action='store', default="R95pTOT",
                        help='Which index to run')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default = False')
 
    args = parser.parse_args()         

    if args.index not in ["R95pTOT", "R99pTOT"]:
        print("Calculation not needed for {}".format(args.index))

    else:
        if args.indata in ["west_africa_pptn", "ghcndex", "south_africa"]:
            main(indata=args.indata, index=args.index, diagnostics=args.diagnostics)
            
        else:
            print("Calculation not needed for {}".format(args.indata))

#------------------------------------------------------------
# END
#------------------------------------------------------------
