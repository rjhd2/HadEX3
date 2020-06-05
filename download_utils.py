#!/usr/local/sci/bin/python2.7
#------------------------------------------------------------
#                    SVN Info
#$Rev:: 493                                           $:  Revision of last commit
#$Author:: rdunn                                      $:  Author of last commit
#$Date:: 2020-04-29 16:53:01 +0100 (Wed, 29 Apr 2020) $:  Date of last commit
#------------------------------------------------------------
# START
#------------------------------------------------------------
'''
Downloads the necessary input files from a remote location

download.py invoked by typing::

  python download.py --indata "ecad" --index "TX90p" --diagnostics

Input arguments:

--indata        Which input data to use (ECAD the only sensible option)

--index         Which ETCCDI index to use

--diagnostics   Output extra info (default False)
'''

import urllib.request, urllib.error, urllib.parse
import os
import shutil
import glob
import subprocess
import time
import tarfile
import zipfile
from socket import error as SocketError
import errno
import numpy as np

# RJHD utils
import utils

HOSTS = {"ghcnd" : "ftp.ncdc.noaa.gov", \
#             "ecad" : "https://www.ecad.eu", \
             "ecad" : "https://knmi-ecad-assets-prd.s3.amazonaws.com", \
             "lacad" : "http://lacad.ciifen.org", \
             "sacad" : "http://sacad.database.bmkg.go.id", \
#             "eobs" : "https://www.ecad.eu", \
             "eobs" : "https://knmi-ecad-assets-prd.s3.amazonaws.com", \
             "laobs" : "http://lacad.ciifen.org", \
             "saobs" : "http://sacad.database.bmkg.go.id"}

REMOTE_LOC = {"ghcnd" : "pub/data/ghcn/daily", \
                  "ecad" : "download/millennium/data", \
                  "lacad" : "download/millennium/data", \
                  "sacad" : "download/millennium/data"}

#*********************************************
def doftp(host, remote_loc, filename, local_loc, diagnostics=False):
    """
    Interface to Met Office doftp command

    :param str host: FTP hostename
    :param str remote_loc: remote directory
    :param str filename: remote filename
    :param str local_loc: local directory
    :param bool diagnostics: extra output
    """

    try:
        if diagnostics:
            print("doftp -host {} -cwd {} -get {}={}/{}".format(host, remote_loc, filename, local_loc, filename))
        subprocess.check_call(["doftp", '-host', host, '-cwd', remote_loc, '-get', "{}={}/{}".format(filename, local_loc, filename)])

        if diagnostics:
            print("     Downloaded")
                   
    # handle the error
    except subprocess.CalledProcessError as e:
        print(e.output)
        print("waiting 10 sec and trying again")
        import time
        time.sleep(10) # wait 10 seconds and onto next while loop

    except OSError:
        # executable not found
        print("Issue with doftp")
        raise OSError("doftp not found")

    return # doftp
                         
#*********************************************
def http(host, remote_loc, filename, local_loc, diagnostics=False):
    """
    Using urllib2 http access

    :param str host: HTTP hostename
    :param str remote_loc: remote directory
    :param str filename: remote filename
    :param str local_loc: local directory
    :param bool diagnostics: extra output
    """

    success = 0
    try:
        if diagnostics:
            print("{}/{}/{}".format(host, remote_loc, filename))

        remote_file = urllib.request.urlopen("{}/{}/{}".format(host, remote_loc, filename))
        
#        if filename.split(".")[-1] == "zip":
#            # need to write out a "bytes" file, not a string one
        local_file = open(os.path.join(local_loc, filename), "wb")
#        else:
#            local_file = open(os.path.join(local_loc, filename), "w")
    
        local_file.write(remote_file.read())
        local_file.close()
        if diagnostics:
            print("     Downloaded")
        success = 1

    # handle the error
    except urllib.error.HTTPError as e:
        print(e.code)
        print(e.read())

        if e.code == 404:
            # server doesn't have the file, so waiting won't help
            print("File {} Not Found - no further attempts".format(filename))
            success = 1
            
    except urllib.error.URLError as e:
        print("Some sort of urllib2 error for {}".format("{}/{}/{}".format(host, remote_loc, filename)))
        if e.reason == "Not Found":
            # server doesn't have the file, so waiting won't help
            print("File {} Not Found - no further attempts".format(filename))
            success = 1

    except SocketError as e:
        if e.errno != errno.ECONNRESET:
            raise # Not the error we are looking for
        
        # likely that too many calls at once.  
        # Sleep for a while and try again.  Handled in caller

    return success # http

#*********************************************
def local_copy(source, destination, extension="", diagnostics=False):
    """
    Perform local copy from networked storage to working area 
 
    :param str source: source directory
    :param str destination: destination directory
    :param str extension: optional filename extension
    """
    
    for filename in glob.glob(r'{}*{}'.format(os.path.expanduser(source), extension)):

        if not os.path.exists(os.path.join(destination, filename.split("/")[-1])):

            shutil.copy(filename, destination)        

            if diagnostics:
                print(filename)
        else:
            if diagnostics:
                print(" exists".format(filename))
    
        # force update of timestamps
        os.utime(os.path.join(destination, filename.split("/")[-1]), None)

    return # local_copy

#*********************************************
def get_spain(dataset, diagnostics=False):
    """ data from Manola Brunet """

    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

    return # get_spain

#*********************************************
def get_russia(dataset, diagnostics=False):
    """ data from Xuebin Zhang """

    orig_file_loc = "{}/{}/{}/Russia_noqc/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

    return # get_russia

#*********************************************
def get_acre(dataset, diagnostics=False):
    """ data from Rob Allan """

    orig_file_loc = "{}/{}/{}/climdex_fmt/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

    return # get_acre

#*********************************************
def get_honduras(dataset, diagnostics=False):
    """ data from Ernesto Salgado"""

    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    # get .csv files only - converted from .xls on 14/08/2018
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension=".csv", diagnostics=diagnostics)

    return # get_honduras

#*********************************************
def get_south_america(dataset, diagnostics=False):
    """ data from Maria de los Milagros Skansi"""

    # 6th November - copying both raw and index data to merge later

    # get the raw observation set
    orig_file_loc = "{}/{}/{}/observations/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    # v20181023 - raw data - is already in .csv
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension=".csv", diagnostics=diagnostics)

    # get the indices set
    orig_file_loc = "{}/{}/{}/indices/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    # indices - get .csv files only - converted from .xls on 24/08/2018 (v20180813)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension=".csv", diagnostics=diagnostics)

    return # get_south_america

#*********************************************
def get_nz(dataset, diagnostics=False):
    """ data from Jim Salinger """

    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

    return # get_nz

#*********************************************
def get_hadex2(dataset, index="TX90p", diagnostics=False):
    """ 
    Get data from HadEX + HadEX2 workshops
    Index based, so a more shuffling required
    """

    import re # regular expressions to check entries.

    orig_file_loc = "{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name)
    filename = "{}/station_list/HadEX2_station_list_all_{}.csv".format(orig_file_loc, index)
    outfile = open(os.path.join(dataset.location, "HadEX2_station_list_{}.csv".format(index)), "w")

    # if data doesn't exist for index, exit
    if not os.path.exists(filename):
        if diagnostics:
            print("no file listing for {} {}".format(dataset.name, index))
        return

    # ensure destination exists
    if not os.path.exists(os.path.join(dataset.location, "raw", index)):
        os.mkdir(os.path.join(dataset.location, "raw", index))

    # find lines in file with "workshop" in them
    with open(filename, "r", encoding="utf-8") as infile:
        for line in infile:

            if re.search("workshop", line) or re.search("West_Indian_ocean", line): 
                sline = line.split(",")

                # and copy the file
                try:
                    shutil.copy(os.path.join(orig_file_loc, dataset.version, index, "{}_{}.txt".format(sline[0][1:-1], index)), os.path.join(dataset.location, "raw", index))      
                    if diagnostics:
                        print(os.path.join(orig_file_loc, dataset.version, index, "{}_{}.txt".format(sline[0][1:-1], index)))
                    
                    # write the relevant metadata to the file
                    outfile.write(line)

                except IOError:
                    # source file doesn't exist.  Fail nicely
                    if diagnostics:
                        print("{} doesn't exist".format(os.path.join(orig_file_loc, dataset.version, index, "{}_{}.txt".format(sline[0][1:-1], index))))
                    
    outfile.close()
                                
    return # get_hadex2

#*********************************************
def get_pacific(dataset, diagnostics=False):
    """ data from Simon McGree """

    def fix_coords(coord):
        """Process coordinate entry in file"""
        
        coord = coord.strip()
        coord, direction = coord.split("\xb0")
        coord = float(coord)

        # ensure correct direction
        if direction in ["S", "W"]:
            coord *= -1

        return coord # fix_coords


    # data were calculated by Climpact2, so already in right format

    orig_dir_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    destination = os.path.join(dataset.location, "formatted", "indices")

    try:
        shutil.rmtree(destination)
    except OSError:
        # already removed
        pass

    try:
        shutil.copytree(os.path.join(orig_dir_loc, "indices"), destination)
        # ensure update of timestamps
        for root, diry, files in os.walk(destination):
            for fname in files:
                os.utime(os.path.join(root, fname), None)

        if diagnostics:
            print("copied {} to {}".format(os.path.join(orig_dir_loc, "indices"), destination))

    except OSError:
        # likely already copied
        pass


    # and translate the metadata
    metadata_file = os.path.join(dataset.location, "pacific.metadata.txt")
    if not os.path.exists(metadata_file):
        
        utils.write_climpact_inventory_header(metadata_file)

        with open(os.path.join(orig_dir_loc, "trendToIndices", "Station_list.csv"), "r", encoding="latin-1") as infile:

            for line in infile:
                line = line.split(",")
                if line[0] == "Station": continue

                assert len(line) == 3
                name, lon, lat = line

                lon = fix_coords(lon)
                lat = fix_coords(lat)               

                station = utils.Station(name, lat, lon, dataset.location, dataset.name)
                utils.write_climpact_inventory(metadata_file, station)

    return # get_pacific

#*********************************************
def get_ecad(dataset, index="TX90p", diagnostics=False, download=False):
    """ Get ECAD family (ECAD/LACAD/SACAD) from remote sites """

    # Correct index name from ECAD to HadEX version
    if index in ["Rx1day", "Rx5day"]:
        index_name = "X".join(index.split("x"))
    else:
        index_name = index
   
    if not os.path.exists(os.path.join(dataset.location, "raw", "{}_index{}.zip".format(dataset.name.upper()[:-1], index_name))):
        
        if download:
            success, counter = 0, 0
            while success == 0:
                success = http(HOSTS[dataset.name], REMOTE_LOC[dataset.name], "{}_index{}.zip".format(dataset.name.upper()[:-1], index_name), os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

                # allow processes to finish up *or* server to become free
                time.sleep(60)

                counter += 1
                if counter > 5:
                    print("too many download attempts {} {}".format(dataset.name, index))
                    raise RuntimeError
        else:
            # pull from local storage
            orig_file_loc = "{}/{}/{}/{}".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version, "{}_index{}.zip".format(dataset.name.upper()[:-1], index_name))
            local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)
            
    else:
        if diagnostics:
            print(" exists".format(os.path.join(dataset.location, "raw", "{}_index{}.zip".format(dataset.name.upper()[:-1], index_name))))

    return # get_ecad

#*********************************************
def get_ghcnd(dataset, diagnostics=False, download=False):
    """ 
    Get GHCND from remote site or locally as remote download can fail due to size

    Downloads the tar.gz file, unpacks this and then renames the resulting directory
    """

    if utils.ghcnd_gsn or utils.ghcnd_peterson:
        # these override
        tar_archive = "ghcnd_all.tar.gz"
    elif utils.ghcnd_hcn:
        tar_archive = "ghcnd_hcn.tar.gz"

    # get the main tar file
    if not os.path.exists(os.path.join(dataset.location, tar_archive)):

        if download:
            print("downloading GHCND tar file")
            doftp(HOSTS[dataset.name], REMOTE_LOC[dataset.name], tar_archive, dataset.location, diagnostics=diagnostics)

        else:
            print("Copying GHCND tar file")
            orig_file_loc = "{}/{}/{}/{}".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version, tar_archive)
            local_copy(orig_file_loc, dataset.location, diagnostics=diagnostics)
    else:
        print("Archive already present")
            

    # "raw" created during setup, so if empty, delete,
    #  extract the files, and rename directory
    if len(os.listdir(os.path.join(dataset.location, "raw"))) == 0:
        os.rmdir(os.path.join(dataset.location, "raw"))

        print("untar-ing archive")
        tar = tarfile.open(os.path.join(dataset.location, tar_archive))
        tar.extractall(path=os.path.join(dataset.location))
        tar.close()
        # extracts to directory called "ghcnd_all" - so just rename this to "raw"
        os.rename(os.path.join(dataset.location, tar_archive.split(".")[0]), os.path.join(dataset.location, "raw"))

        # ensure update of timestamps
        for root, diry, files in os.walk(os.path.join(dataset.location, "raw")):
            for fname in files:
                os.utime(os.path.join(root, fname), None)

        if diagnostics: print("   complete")
    else:
        if diagnostics: print("{} not empty - GHCND not extracted".format(os.path.join(dataset.location, "raw")))

    # get the inventory and metadata
    if not os.path.exists(os.path.join(dataset.location, "raw", "ghcnd-inventory.txt")):
        if download:
            doftp(HOSTS[dataset.name], REMOTE_LOC[dataset.name], "ghcnd-inventory.txt", os.path.join(dataset.location, "raw"), diagnostics=diagnostics)
        else:
            orig_file_loc = "{}/{}/{}/{}".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version, "ghcnd-inventory.txt")
            local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)
    if not os.path.exists(os.path.join(dataset.location, "raw", "ghcnd-stations.txt")):
        if download:
            doftp(HOSTS[dataset.name], REMOTE_LOC[dataset.name], "ghcnd-stations.txt", os.path.join(dataset.location, "raw"), diagnostics=diagnostics)
        else:
            orig_file_loc = "{}/{}/{}/{}".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version, "ghcnd-stations.txt")
            local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

    return # get_ghcnd

#*********************************************
def get_west_africa_pptn(dataset, diagnostics=False):
    """
    Data from Theo Vischel, Guillaume Chagnaud, Youssouph Sané, Gérémy Panthou, Francis Nkrumah.
    """
    orig_file_loc = "{}/{}/{}/1900-2016/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension=".csv", diagnostics=diagnostics)

    # and get the metadata file
    local_copy("{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version), os.path.join(dataset.location, "raw"), extension=".csv", diagnostics=diagnostics)

    return # get_west_africa_pptn

#*********************************************
def get_australia(dataset, diagnostics=False):
    """ Data from Blair Trewin """

    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

    return # get_australia

#*********************************************
def get_arabia_6190(dataset, index="TX90p", diagnostics=False):
    """ Data from Arabian workshop via Markus Donat (long record stations only, ref 1961-90)"""

    orig_file_loc = "{}/{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version, index)
    
    # ensure destination exists (extra subdirectory)
    if not os.path.exists(os.path.join(dataset.location, "raw", index)):
        os.mkdir(os.path.join(dataset.location, "raw", index))

    local_copy(orig_file_loc, os.path.join(dataset.location, "raw", index), diagnostics=diagnostics)

    #inventory
    local_copy("{}/{}/{}/station-List.csv".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version), os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

    return # get_arabia_6190

#*********************************************
def get_arabia_8110(dataset, index="TX90p", diagnostics=False):
    """ Data from Arabian workshop via Markus Donat (all stations, ref 1981-2010)"""

    orig_file_loc = "{}/{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version, index)
    
    # ensure destination exists (extra subdirectory)
    if not os.path.exists(os.path.join(dataset.location, "raw", index)):
        os.mkdir(os.path.join(dataset.location, "raw", index))

    local_copy(orig_file_loc, os.path.join(dataset.location, "raw", index), diagnostics=diagnostics)

    #inventory
    local_copy("{}/{}/{}/station-List.csv".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version), os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

    return # get_arabia_8110

#*********************************************
def get_chile(dataset, diagnostics=False):
    """ data from Claudia Villarroel """

    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

    return # get_chile

#*********************************************
def get_colombia(dataset, diagnostics=False):
    """ data from Jose Daniel Pabon Caicedo """

    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="csv", diagnostics=diagnostics)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="txt", diagnostics=diagnostics)

    return # get_colombia

#*********************************************
def get_canada(dataset, diagnostics=False):
    """ data from Lucie Vincent and Vincent Cheng """

    # data
    orig_file_loc = "{}/{}/{}/data/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="txt", diagnostics=diagnostics)
    # inventory
    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="csv", diagnostics=diagnostics)

    return # get_colombia

#*********************************************
def get_decade(dataset, diagnostics=False):
    """ data from DECADE project - Hunziker & Bronnimann """

    # data
    for dat in ["TN", "TX", "PP"]:
        orig_file_loc = "{}/{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version, dat)
        if not os.path.exists(os.path.join(dataset.location, "raw", dat)):
            os.mkdir(os.path.join(dataset.location, "raw", dat))
        local_copy(orig_file_loc, os.path.join(dataset.location, "raw", dat), extension="dat", diagnostics=diagnostics)

    return # get_decade

#*********************************************
def get_south_africa(dataset, diagnostics=False):
    """ data from SAWS - Andries Kruger """
    
    orig_file_loc = "{}/{}/{}/temperatures/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="csv", diagnostics=diagnostics)
    orig_file_loc = "{}/{}/{}/precipitation/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="csv", diagnostics=diagnostics)

    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="txt", diagnostics=diagnostics)

    return # get_south_africa

#*********************************************
def get_west_africa_indices(dataset, diagnostics=False):
    """ data from Aziz Barry """
    
    orig_file_loc = "{}/{}/{}/indices/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="txt", diagnostics=diagnostics)
 
    orig_file_loc = "{}/{}/{}/country_stations.csv".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

    return # get_west_africa_indices

#*********************************************
def get_eobs_style(dataset, data_name, variable_name, diagnostics=False, download=False):
    '''get E-OBS style data'''    

    if not os.path.exists(os.path.join(dataset.location, "raw", "{}_nonblend_{}.zip".format(data_name, variable_name))):
        
        if download:
            if dataset.name in ["laobs", "saobs"]:
                http(HOSTS[dataset.name], "utils/downloadfile.php?file=download", "{}_nonblend_{}.zip".format(data_name, variable_name), os.path.join(dataset.location, "raw"), diagnostics=diagnostics)
            elif dataset.name == "eobs":
                http(HOSTS[dataset.name], "download", "{}_nonblend_{}.zip".format(data_name, variable_name), os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

        else:
            # copy from local storage
            orig_file_loc = "{}/{}/{}/{}".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version, "{}_nonblend_{}.zip".format(data_name, variable_name))
            local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)

    else:
        if diagnostics:
            print("{} exists".format(os.path.join(dataset.location, "raw", "{}_nonblend_{}.zip".format(data_name, variable_name))))

    if download:
        http(HOSTS[dataset.name], "download", "stations_{}.txt".format(variable_name), \
             os.path.join(dataset.location, "raw"), diagnostics=diagnostics)
    else:
        orig_file_loc = "{}/{}/{}/{}".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version, "stations_{}.txt".format(variable_name))
        local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), diagnostics=diagnostics)


    zip_ref = zipfile.ZipFile(os.path.join(dataset.location, "raw", "{}_nonblend_{}.zip".format(data_name, variable_name)), 'r')

    # make the folder to hold the converted data
    if not os.path.exists(os.path.join(dataset.location, "raw", dataset.name.upper())):
        os.mkdir(os.path.join(dataset.location, "raw", dataset.name.upper()))
        
    zip_ref.extractall(os.path.join(dataset.location, "raw", dataset.name.upper()))
    zip_ref.close()

    # move stations.txt metadata file to a unique filename   
    if os.path.exists(os.path.join(dataset.location, "raw", "stations_{}.txt".format(variable_name))):
        shutil.move(os.path.join(dataset.location, "raw", "stations_{}.txt".format(variable_name)), \
                    os.path.join(dataset.location, "raw/{}".format(dataset.name.upper()), "stations_{}.txt".format(variable_name)))

    if os.path.exists(os.path.join(dataset.location, "raw/{}".format(dataset.name.upper()), "sources.txt".format(variable_name))):
        shutil.move(os.path.join(dataset.location, "raw/{}".format(dataset.name.upper()), "sources.txt"), \
                    os.path.join(dataset.location, "raw/{}".format(dataset.name.upper()), "sources_{}.txt".format(variable_name)))

    return # get_eobs_style

#*********************************************
def get_eobs(dataset, diagnostics=False, download=False):
    '''get ECAD raw obs data (EOBS)'''

    get_eobs_style(dataset, "ECA", "tx", diagnostics=diagnostics, download=download)
    if diagnostics: print("ECA TX done")
    get_eobs_style(dataset, "ECA", "tn", diagnostics=diagnostics, download=download)
    if diagnostics: print("ECA TN done")
    get_eobs_style(dataset, "ECA", "rr", diagnostics=diagnostics, download=download)
    if diagnostics: print("ECA PP done")

    return # get_laobs

#*********************************************
def get_laobs(dataset, diagnostics=False, download=False):
    '''get LACAD raw obs data'''

    get_eobs_style(dataset, "LACA", "tx", diagnostics=diagnostics download=download)
    if diagnostics: print("LACA TX done")
    get_eobs_style(dataset, "LACA", "tn", diagnostics=diagnostics download=download)
    if diagnostics: print("LACA TN done")
    get_eobs_style(dataset, "LACA", "rr", diagnostics=diagnostics download=download)
    if diagnostics: print("LACA PP done")

    return # get_laobs

#*********************************************
def get_saobs(dataset, diagnostics=False, download=False):
    '''get SACAD raw obs data'''

    get_eobs_style(dataset, "SACA", "tx", diagnostics=diagnostics download=download)
    if diagnostics: print("SACA TX done")
    get_eobs_style(dataset, "SACA", "tn", diagnostics=diagnostics download=download)
    if diagnostics: print("SACA TN done")
    get_eobs_style(dataset, "SACA", "rr", diagnostics=diagnostics download=download)
    if diagnostics: print("SACA PP done")

    return # get_saobs

#*********************************************
def get_singapore_workshop(dataset, diagnostics=False):
    """ data from Singapore Workshop (except Indonesia) """
    
    # data were calculated by Climpact2, so already in right format

    orig_dir_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    destination = os.path.join(dataset.location, "formatted", "indices")

    try:
        shutil.rmtree(destination)
    except OSError:
        # already removed
        pass

    try:
        shutil.copytree(os.path.join(orig_dir_loc, "indices"), destination)
        # ensure update of timestamps
        for root, diry, files in os.walk(destination):
            for fname in files:
                os.utime(os.path.join(root, fname), None)

        if diagnostics:
            print("copied {} to {}".format(os.path.join(orig_dir_loc, "indices"), destination))

    except OSError:
        # likely already copied
        pass


    # no inventory list as yet (1/5/2019) - extract from files
    metadata_file = os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name))
    if not os.path.exists(metadata_file):
        
        utils.write_climpact_inventory_header(metadata_file)

        # spin through each station and read one file for the information
        for dirname in os.listdir(destination):
            
            # singapore stations have a little more information, so filter
            if dataset.name == "singapore":
               
                if dirname.split("_")[-1] == "".join(dataset.base_period.split("-")) and \
                   dirname.split("_")[-2] == "cpDLY":
                    # daily changepoints and the correct reference period
                    pass
                elif dirname.split("_")[-1] == "".join(dataset.base_period.split("-")) and \
                   dirname.split("_")[-2] == "cpMLY" and dirname.split("_")[0] == "S07":
                    # monthly changepoints and the correct reference period
                    # only have monthly for this station
                    pass
                else:
                    continue
                    
            stationfiles = os.listdir(os.path.join(destination, dirname))

            if len(stationfiles) != 0:
                                    
                with open(os.path.join(destination, dirname, stationfiles[0]), "r", encoding="latin-1") as infile:

                    for line in infile:

                        line = line.split()

                        if line[0] == '"Latitude:':
                            lat = float(line[1].split(",")[1][1:-1])
                        if line[0] == '"Longitude:':
                            lon = float(line[1].split(",")[1][1:-1])
                        if line[0] == "time,":
                            break

                    station = utils.Station(dirname, lat, lon, dataset.location, dataset.name)
                    utils.write_climpact_inventory(metadata_file, station)

                    if diagnostics:
                        print("{}: {}, {}".format(dirname, lat, lon))

    return # get_singapore_workshop

#*********************************************
def get_indonesia(dataset, diagnostics=False):
    """ data from Indonesia from Singapore workshop """
    
    # data were calculated by Climpact2, so already in right format
    #   but directory structure different to other sources from this workshop

    orig_dir_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    destination = os.path.join(dataset.location, "formatted", "indices")

    try:
        shutil.rmtree(destination)
    except OSError:
        # already removed
        pass
    os.mkdir(destination)

    # will have to walk each station to find data
    for path, dirs, files in os.walk(orig_dir_loc):
        if "indices" in path:
            iscsv = False
            for thisfile in files:
                # use single file for the station name
                if thisfile[-5:] == "N.csv":
                    stn_name = glob.glob("{}/*dtr_ANN*".format(path))[0].split("/")[-1].split("_dtr_")[0]
                    iscsv = True
                    break
            if iscsv:
                shutil.copytree(path, os.path.join(destination, stn_name))
                if diagnostics:
                    print("copied {} to {}".format(path, os.path.join(destination, stn_name)))
                   

    # now need to remove spaces
    for path, dirs, files in os.walk(os.path.join(destination)):
        if " " in path:
            for this_file in files:
                os.rename(os.path.join(path, this_file), os.path.join(path, this_file.replace(" ", "_")))
            os.rename(path, path.replace(" ", "_"))
            os.listdir(path.replace(" ", "_"))
                         

    # no inventory list as yet - extract from files
    metadata_file = os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name))
    if not os.path.exists(metadata_file):
        
        utils.write_climpact_inventory_header(metadata_file)

        # spin through each station and read one file for the information
        for dirname in os.listdir(destination):
            
                  
            stationfiles = os.listdir(os.path.join(destination, dirname))
                    
            with open(os.path.join(destination, dirname, stationfiles[0]), "r", encoding="latin-1") as infile:

                for line in infile:

                    line = line.split()

                    if line[0] == '"Latitude:':
                        lat = float(line[1].split(",")[1][1:-1])
                    if line[0] == '"Longitude:':
                        lon = float(line[1].split(",")[1][1:-1])
                    if line[0] == "time,":
                        break

                station = utils.Station(dirname, lat, lon, dataset.location, dataset.name)
                utils.write_climpact_inventory(metadata_file, station)
                
                if diagnostics:
                    print("{}: {}, {}".format(dirname, lat, lon))

    return # get_indonesia

#*********************************************
def get_china(dataset, diagnostics=False):
    """ data from Ying Sun for China """
    
    # data were calculated by Climpact2, so already in right format

    orig_dir_loc = "{}/{}/{}/319_all/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)

    destination = os.path.join(dataset.location, "formatted", "indices")

    try:
        shutil.rmtree(destination)
    except OSError:
        # already removed
        pass

    try:
        shutil.copytree(orig_dir_loc, destination)
        # ensure update of timestamps
        for root, diry, files in os.walk(destination):
            for fname in files:
                os.utime(os.path.join(root, fname), None)

        if diagnostics:
            print("copied {} to {}".format(os.path.join(orig_dir_loc, "indices"), destination))

    except OSError:
        # likely already copied
        pass

    # no inventory list as yet (22/5/2019) - extract from files
    metadata_file = os.path.join(dataset.location, "{}.metadata.txt".format(dataset.name))
    if not os.path.exists(metadata_file):
        
        utils.write_climpact_inventory_header(metadata_file)

        # spin through each station and read one file for the information
        for dirname in os.listdir(destination):

            stationfiles = os.listdir(os.path.join(destination, dirname))

            with open(os.path.join(destination, dirname, stationfiles[0]), "r", encoding="latin-1") as infile:

                for line in infile:

                    line = line.split()

                    if line[0] == '"Latitude:':
                        lat = float(line[1].split(",")[1][1:-1])
                    if line[0] == '"Longitude:':
                        lon = float(line[1].split(",")[1][1:-1])
                    if line[0] == "time,":
                        break

                station = utils.Station(dirname, lat, lon, dataset.location, dataset.name)
                utils.write_climpact_inventory(metadata_file, station)
  
    return # get_china

#*********************************************
def get_india(dataset, diagnostics=False):
    """ data from Arvind Srivastava """
    
    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="PRT", diagnostics=diagnostics)

    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="csv", diagnostics=diagnostics)
 
    return # get_india

#*********************************************
def get_japan(dataset, diagnostics=False):
    """ data from Hisayuki Kubota """
    
    orig_dir_loc = "{}/{}/{}/Japan_Station_data/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    destination = os.path.join(dataset.location, "raw")

    try:
        shutil.rmtree(destination)
    except OSError:
        # already removed
        pass

    try:
        shutil.copytree(orig_dir_loc, destination)
        # ensure update of timestamps
        for root, diry, files in os.walk(destination):
            for fname in files:
                os.utime(os.path.join(root, fname), None)

        if diagnostics:
            print("copied {} to {}".format(orig_dir_loc, destination))

    except OSError:
        # likely already copied
        pass

    # and the readme (inventory) file
    orig_file_loc = "{}/{}/{}/Japan_Station_data/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="txt", diagnostics=diagnostics)
 
    return # get_japan

#*********************************************
def get_mexico(dataset, diagnostics=False):
    """ data from Jorge Vazquez-Aguirre """
    
    orig_dir_loc = "{}/{}/{}/mxrclapr19/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    destination = os.path.join(dataset.location, "raw")

    subdirs = os.listdir(orig_dir_loc)

    for subd in subdirs:
        if diagnostics:
            print(subd)

        local_copy("{}/".format(os.path.join(orig_dir_loc, subd)), os.path.join(dataset.location, "raw"), extension="txt", diagnostics=diagnostics)
        

    # and the readme (inventory) file
    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="csv", diagnostics=diagnostics)
 
    return # get_mexico

#*********************************************
def get_iran(dataset, diagnostics=False):
    """ data from Fatemeh Rahimzadeh """
    
    orig_dir_loc = "{}/{}/{}/extrem/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    destination = os.path.join(dataset.location, "raw")

    local_copy(orig_dir_loc, os.path.join(dataset.location, "raw"), extension="txt", diagnostics=diagnostics)

    # and the readme (inventory) file
    orig_file_loc = "{}/{}/{}/extrem/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    local_copy(orig_file_loc, os.path.join(dataset.location, "raw"), extension="csv", diagnostics=diagnostics)
 
    return # get_iran

#*********************************************
def get_brazil(dataset, diagnostics=False):
    """ data from Andrea Ramos """
    
    orig_dir_loc = "{}/{}/{}/TMAX_TMIN_PREC_1961-2019/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    destination = os.path.join(dataset.location, "raw")

    local_copy(orig_dir_loc, os.path.join(dataset.location, "raw"), extension="CSV", diagnostics=diagnostics)
        
    # there isn't an inventory file, so create one
    # no inventory list as yet (22/5/2019) - extract from files
    metadata_file = os.path.join(dataset.location, "{}_inventory.txt".format(dataset.name))
    if not os.path.exists(metadata_file):
        
        utils.write_climpact_inventory_header(metadata_file)

        # spin through each station and read one file for the information
        for filename in os.listdir(destination):

            if filename.split("_")[0] == "DADOS":
                with open(os.path.join(destination, filename), "r", encoding="latin-1") as infile:

                    for line in infile:

                        line = line.split(";")

                        if line[0] == 'LATITUDE:':
                            lat = float(line[1].replace(",", "."))
                        if line[0] == 'LONGITUDE:':
                            lon = float(line[1].replace(",", "."))
                        if line[0] == "ALTITUDE:":
                            break

                    stnid = filename.split("_")[3]
                    station = utils.Station(stnid, lat, lon, dataset.location, dataset.name)
                    utils.write_climpact_inventory(metadata_file, station)

    return # get_brazil

#*********************************************
def get_brazil_sp(dataset, diagnostics=False):
    """ data from Jose Marengo (Sao Paulo) """
    
    orig_dir_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    destination = os.path.join(dataset.location, "raw")

    local_copy(orig_dir_loc, os.path.join(dataset.location, "raw"), extension="txt", diagnostics=diagnostics)
    local_copy(orig_dir_loc, os.path.join(dataset.location, "raw"), extension="csv", diagnostics=diagnostics)

    # two sets of files - 3 files for Agua Funda, 1 for a number of other stations

    return # get_brazil_sp

#*********************************************
def get_ghcndex(dataset, index="TX90p", diagnostics=False):
    """ 
    Get data from HadEX + HadEX2 workshops
    Index based, so a more shuffling required
    """

    import re # regular expressions to check entries.

    orig_file_loc = "{}/{}/{}/".format(utils.INPUT_DATA_LOCATION, dataset.name, dataset.version)
    filename = "{}/ghcnd_stations_{}_40+datayrs.txt".format(orig_file_loc, index)
    outfile = open(os.path.join(dataset.location, "GHCNDEX_station_list_{}.csv".format(index)), "w")

    # ensure destination exists
    if not os.path.exists(os.path.join(dataset.location, "raw", index)):
        os.mkdir(os.path.join(dataset.location, "raw", index))

    try:
        with open(filename, "r", encoding="utf-8") as infile:
            for line in infile:

                sline = line.split()

                if line[:3] != "ASN":
                    # Australia only - skip the rest
                    continue

                # and copy the file
                try:
                    shutil.copy(os.path.join(orig_file_loc, "stn-indices", index, "{}_{}.txt".format(sline[0], index)), os.path.join(dataset.location, "raw", index))      
                    if diagnostics:
                        print(os.path.join(orig_file_loc, "stn-indices", index, "{}_{}.txt".format(sline[0], index)))

                    # write the relevant metadata to the file
                    outfile.write(line)

                except IOError:
                    # source file doesn't exist.  Fail nicely
                    if diagnostics:
                        print("{} doesn't exist".format(os.path.join(orig_file_loc, "stn-indices", index, "{}_{}.txt".format(sline[0], index))))
    except IOError:
        if diagnostics:
            print("{} doesn't exist".format(filename))

    outfile.close()
                                
    return # get_ghcndex

#*********************************************
def main(indata="ghcnd", index="TX90p", aggregated=False, diagnostics=False):
    """
    Extract relevant dataset from the command line switches

    :param str indata: input dataset to process
    :param str index: ETCCDI index name to process
    :param bool diagnostics: output diagnostic information
    """

    # get all possible datasets
    all_datasets = utils.get_input_datasets()

    # and their names
    names = np.array([d.name for d in all_datasets])

    # if dataset selected and in the list of available, then run
    if indata in names:

        # index data which isn't split by index - so make up a set
        if aggregated:
            if indata == "south_america":
                get_south_america(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "west_africa_pptn":
                get_west_africa_pptn(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "west_africa_indices":
                get_west_africa_indices(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "pacific":
                get_pacific(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "south_africa":
                get_south_africa(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata in ["malaysia", "brunei", "myanmar", "vietnam", "philippines", "singapore", "thailand"]:
                get_singapore_workshop(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "china":
                get_china(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "indonesia":
                get_indonesia(all_datasets[names == indata][0], diagnostics=diagnostics)

        else:
            # raw obs
            if indata == "ghcnd":
                get_ghcnd(all_datasets[names == indata][0], diagnostics=diagnostics, download=False)
            if indata == "australia":
                get_australia(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "spain":
                get_spain(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "honduras":
                get_honduras(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "nz":
                get_nz(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "russia":
                get_russia(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "acre":
                get_acre(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "chile":
                get_chile(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "colombia":
                get_colombia(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "canada":
                get_canada(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "decade":
                get_decade(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "india":
                get_india(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "japan":
                get_japan(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "mexico":
                get_mexico(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "iran":
                get_iran(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "brazil":
                get_brazil(all_datasets[names == indata][0], diagnostics=diagnostics)
            if indata == "brazil_sp":
                get_brazil_sp(all_datasets[names == indata][0], diagnostics=diagnostics)

            if utils.REF_START == 1981 and utils.REF_END == 2010:
                # only download these if needed - running 1981-2010 Base Period
                if indata == "eobs":
                    get_eobs(all_datasets[names == indata][0], diagnostics=diagnostics, download=False)
                if indata == "laobs":
                    get_laobs(all_datasets[names == indata][0], diagnostics=diagnostics, download=False)
                if indata == "saobs":
                    get_saobs(all_datasets[names == indata][0], diagnostics=diagnostics, download=False)
            else:
                if indata in ["eobs", "laobs", "saobs"] and diagnostics:
                    print("Not downloading {} as reference period doesn't match".format(indata))

            # index data
            if utils.REF_START == 1961 and utils.REF_END == 1990:
                # only download these if needed - running 1961-1990 Base Period
                if indata == "ecad":
                    get_ecad(all_datasets[names == indata][0], index=index, diagnostics=diagnostics, download=False)
                if indata == "lacad":
                    get_ecad(all_datasets[names == indata][0], index=index, diagnostics=diagnostics, download=False)
                if indata == "sacad":
                    get_ecad(all_datasets[names == indata][0], index=index, diagnostics=diagnostics, download=False)
            else:
                if indata in ["ecad", "lacad", "sacad"] and diagnostics:
                    print("Not downloading {} as reference period doesn't match".format(indata))

            if indata == "hadex2":
                get_hadex2(all_datasets[names == indata][0], index=index, diagnostics=diagnostics)
            if indata == "ghcndex":
                get_ghcndex(all_datasets[names == indata][0], index=index, diagnostics=diagnostics)
            if indata == "arabia":
                if utils.REF_START == 1961 and utils.REF_END == 1990:
                    get_arabia_6190(all_datasets[names == indata][0], index=index, diagnostics=diagnostics)
                elif utils.REF_START == 1981 and utils.REF_END == 2010:
                    get_arabia_8110(all_datasets[names == indata][0], index=index, diagnostics=diagnostics)
  

    # fail gracefully
    else:
        print("data name not available: {}\n".format(indata))
        print("available data names: {}".format(" ".join(names)))

    return # main

#************************************************************************
if __name__ == "__main__":

    import argparse

    # set up keyword arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--indata', dest='indata', action='store', default="ghcnd",
                        help='Which dataset to attempt download')
    parser.add_argument('--diagnostics', dest='diagnostics', action='store_true', default=False,
                        help='Output diagnostic information, default=False')
    parser.add_argument('--index', dest='index', action='store', default="TX90p",
                        help='Which index to download (ECAD-family only), default=TX90p')
    parser.add_argument('--aggregated', dest='aggregated', action='store_true', default=False,
                        help='Run the aggregated index downloads')
 
    args = parser.parse_args()         

    if args.indata in ["ghcnd", "ecad", "sacad", "lacad", "spain", "pacific", "nz", "honduras", "hadex2", \
                       "russia", "acre", "south_america", "west_africa_pptn", "west_africa_indices", \
                       "australia", "arabia", "chile", "colombia", "canada", "decade", \
                       "south_africa", "eobs", "laobs", "saobs", "malaysia", "india", "brunei", "china", \
                       "myanmar", "japan", "mexico", "iran", "brazil", "vietnam", "philippines", \
                       "singapore", "indondesia", "brazil_sp", "ghcndex", "thailand", "indonesia"]:
        main(indata=args.indata, diagnostics=args.diagnostics, index=args.index, aggregated=args.aggregated)
    else:
        if args.diagnostics:
            print("download not required for {}".format(args.indata))

#*******************************************
# END
#*******************************************
