#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 439                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2019-11-26 13:35:46 +0000 (Tue, 26 Nov 2019) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
"""
Merge individual months together

merge_files.py invoked by typing::

  python merge_files.py --grid "ADW" --index "TX90p" --hadex2_adw --anomalies --diagnostics

Input arguments:

--grid          Which gridding algorithm to use - ADW (Angular Distance Weighting) or CAM (Climate Anomaly Method)

--index         Which ETCCDI index to use

--hadex2_adw    Use form of ADW applied in HadEX2 (slight difference to formula)

--anomalies     Uses anomalies for ADW

--diagnostics   Output extra info (default False)

"""
import os
import numpy as np

# RJHD utilities
import netcdf_procs as ncdfp
import utils


#*********************************************
def main(index="TX90p", grid="ADW", diagnostics=False, hadex2_adw=False, anomalies="None"):
    """
    Merge the monthly files together, or, if just annual index, add extra metadata

    :param str grid: gridding type ADW/CAM
    :param str index: which index to run
    :param str anomalies: run code on anomalies or climatology rather than raw data
    :param bool diagnostics: output diagnostic information

    """
    if grid == "ADW":
        suffixes = ["", "num", "numdls"]
    else:
        suffixes = ["", "num"]

    # both the data and the counts
    for suffix in suffixes:
        station_count = False
        if suffix in ["_num", "_numdls"]:
            station_count = True

        # from when testing and comparing to HadEX2's ADW routine
        if hadex2_adw:
            suffix = "{}_H2ADW".format(suffix)

        if index in utils.MONTHLY_INDICES:
            nmonths = 13
        else:
            nmonths = 1

        # spin through the months
        for month_index in range(nmonths): # annual = 0

            infilename = os.path.join(utils.OUTROOT, utils.make_filenames(index=index, grid=grid, anomalies=anomalies, extra=suffix, month_index=month_index))

            if os.path.exists(infilename):
                indata, lons, lats, times = ncdfp.netcdf_read(infilename, month_index)
                try:
                    # store data
                    data[:, month_index, :, :] = indata
                except NameError:
                    # if first pass through, set up the array.
                    data = np.zeros((times.shape[0], nmonths, lats.shape[0], lons.shape[0]))
                    data[:, month_index, :, :] = indata[:]
                # delete indata to ensure it isn't reused
                indata = 0
                print(infilename)
            else:
                print("File {} doesn't exist".format(infilename))

        # create filename and write.
        outfilename = os.path.join(utils.FINALROOT, utils.make_filenames(index=index, grid=grid, anomalies=anomalies, extra=suffix, month_index=""))
        ncdfp.netcdf_write(outfilename, index, data, times, lats, lons, station_count=station_count)

        # delete data array to ensure try/except clause works on next pass
        del data

    return # main

#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--index', dest='index', action='store', default="TX90p",
                        help='Which index to read')
    parser.add_argument('--grid', dest='grid', action='store', default="ADW",
                        help='gridding routine output to read - ADW or CAM')
    parser.add_argument('--hadex2_adw', dest='hadex2_adw', action='store_true', default=False,
                        help='Read in files using ADW method as in HadEX2 (angles), default=False')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default=False')
    parser.add_argument('--anomalies', dest='anomalies', action='store', default="None",
                        help='To use the anomalies or climatology rather than actuals, default="None"') 

    args = parser.parse_args()

    main(grid=args.grid, index=args.index, hadex2_adw=args.hadex2_adw, diagnostics=args.diagnostics, anomalies=args.anomalies)

#------------------------------------------------------------
# END
#------------------------------------------------------------
