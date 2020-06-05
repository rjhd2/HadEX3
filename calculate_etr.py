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
Calculates the ETR index if files not present.

calculate_etr.py invoked by typing::

  python calcuate_etr.py --indata "ghcnd" --diagnostics

Input arguments:

--indata        Which input data to use 

--diagnostics   Output extra info (default False)
"""
import os
import numpy as np

# RJHD
import utils

#*********************************************
def main(indata="ghcnd", diagnostics=False):
    """
    Read TXn and TNn and write out ETR as difference


    """
    index = "ETR"

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

                if os.path.exists(os.path.join(stn.location, stn.id, "{}_{}_{}.csv".format(stn.id, "txx", ts))) and os.path.exists(os.path.join(stn.location, stn.id, "{}_{}_{}.csv".format(stn.id, "tnn", ts))):
                    xtimes, txx = utils.read_station_index(stn, "TXx", ts)
                    ntimes, tnn = utils.read_station_index(stn, "TNn", ts)

                    match = np.in1d(xtimes, ntimes)
                    match_b = np.in1d(ntimes, xtimes)

                    if len(match) != 0 and len(match_b) != 0:

                        etr = txx[match]-tnn[match_b]
                        etr_times = xtimes[match]

                        if ts == "MON":
                            myears = []
                            months = []
                            for y in etr_times:
                                for m in range(1, 13):
                                    myears += [y]
                                    months += [m]                    

                            stn.monthly = etr.filled().reshape(-1)
                            stn.myears = myears
                            stn.months = months 
                            path = os.path.join(dataset.location, "formatted", "indices", stn.id, "{}_{}_MON.csv".format(stn.id, index.lower()))
                            if not os.path.exists(path):
                                utils.write_station_index(path, stn, "ETR", doMonthly=True)

                        else:            
                            stn.years = etr_times
                            stn.annual = etr.filled()
                            path = os.path.join(dataset.location, "formatted", "indices", stn.id, "{}_{}_ANN.csv".format(stn.id, index.lower()))
                            if not os.path.exists(path):
                                utils.write_station_index(path, stn, "ETR")

    return # main


#************************************************************************
#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--indata', dest='indata', action='store', default="ghcnd",
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
        print("Calculation not needed for {}".format(args.indata))

#------------------------------------------------------------
# END
#------------------------------------------------------------
