.. HadEX3 documentation master file, created by
   sphinx-quickstart on Thu May 23 11:35:12 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HadEX3 documentation!
==================================

Introduction
------------

This set of code is generates a set of files and plots for the HadEX3
dataset, created in preparation for the IPCC AR6 report.  A wide range
of submitted data are ingested, reformatted, processed into ETCCDI
indices and then gridded.

Running EEA
-----------

1. Check out the code
   The code is stored in my own svn repository [one of those legacy
   things].  This is at file:///home/h05/rdunn/svn/hadex3/trunk

2. Edit configuration files
   The configuation and other source files for this code are stored
   in::

   ~hadobs/hadex/input_files_v3

   1. configuration.txt
      This file contains the main settings to run the code.  Adjust
      the runyear, the start and end years/months as required.  All
      other settings should be able to be left alone::

        [Paths]
        root=/scratch/rdunn/hadex3/
        ancillaries=~hadobs/hadex/input_files_v3/
        [Files]
        parameters=parameters.json
        pubA=Pub9volA160415x.flatfile
        [Values]
        mdi=-99.9
        dataset_start_year=1901
        dataset_start_month=1
        dataset_end_year=2019
        dataset_end_month=1
        [Climatologies]
        climatology_completeness=0.85
        [DLS]
        default=200
        maximum_dls=3000
        fix_zero=True
        [Grid]
        deltalon=1.875
        deltalat=1.25
        doLSM=True
        [CAMgrid]
        cam_completeness=10
        cam_clim_start_year=1961
        cam_clim_start_month=1
        cam_clim_end_year=1991
        cam_clim_end_month=1
        stations_in_box=1
        [ADWgrid]
        stations_in_dls=3
        m=4
        [Data]
        finish_after=1981
        modern_period=1951
        minimum_number_of_years=20
        maximum_missing_days_per_year=15
        maximum_missing_days_per_month=3
        maximum_missing_months_per_year=1
        maximum_missing_years=2
        [Climpact]
        remove_extra_output=True
        reference_start=1961
        reference_end=1990
        ncores=1
        wsdi_n=3
        csdi_n=3
        tb_hdd=18
        tb_cdd=18
        tb_gdd=10
        rx_ud=3
        rnnmm_ud=30
        txtn_ud=7
        spei=24
        [QC]
        percentile_completeness=0.7
        [GHCND]
        hcn=True
        gsn=True

      
   2. parameters.json
      This contains the list of the latitude bands in one step and the
      set of indices in another.  

   3. attributes.dat
      This contains a set of global attributes for the final netCDF
      files

   4. input_datasets.dat
      This stores the list of input data collections for HadEX3, along
      with their version and the reference period used if they are
      already in index form.

   5. land_mask_15min.msk
      A land-sea mask at 15min (0.25 degree) resolution used to build
      one for the grid size used.

Using the Rose Suite
--------------------

It will be easiest to use the Rose suite mi-au039 to run all this code
for you.  However, this does depend on having write access to the svn
repository where the code is kept, as some of the settings are hard
coded (at the moment!).

1. Check out the Rose suite *mi-au039*
2. Update the input files as outlined above
3. Run the following from inside your roses/mi-au039 directory::
   rose suite-run 

Should take about 48-72 hours.  The multiple indices lends to a
parallised system, and so I do not recommend running this by hand.


Scripts
-------


Download
^^^^^^^^
.. automodule:: download_utils
   :members: main

Reformat
^^^^^^^^
.. automodule:: reformat_indices
   :members: main

Conversion
^^^^^^^^^^
.. automodule:: conversion_master
   :members: main

Merge Metadata
^^^^^^^^^^^^^^
.. automodule:: merge_metadata
   :members: main

Miscellaneous Tasks
^^^^^^^^^^^^^^^^^^^
.. automodule:: fix_reference_period
   :members: main
.. automodule:: calculate_etr
   :members: main
.. automodule:: calculate_rXXptot
   :members: main

Station Selection
^^^^^^^^^^^^^^^^^
.. automodule:: select_final_stations
   :members: main

Quality Control
^^^^^^^^^^^^^^^
.. automodule:: qc_checks
   :members: main

DLS Calculation
^^^^^^^^^^^^^^^
.. automodule:: dls_cal
   :members: main

Gridding
^^^^^^^^
.. automodule:: gridding
   :members: main

Merging Files
^^^^^^^^^^^^^
.. automodule:: merge_files
   :members: main

Plotting
^^^^^^^^
.. automodule:: plot_timeseries
   :members: main

.. automodule:: plot_timeseries_iris
   :members: main

.. automodule:: plot_stations
   :members: main

.. automodule:: plot_trend_maps
   :members: main

.. automodule:: plot_seasonal_trend_maps
   :members: main

.. automodule:: plot_temporal_difference
   :members: main

.. automodule:: plot_maps_interRefPeriod
   :members: main

.. automodule:: plot_timeseries_interRefPeriod
   :members: main


Contents:

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
