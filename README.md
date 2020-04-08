HadEX3

Code for the creation of the gridded HadEX3 dataset.

Most of this code was developed during 2018-20 at the UK Met Office Hadley Centre, and takes inputs
of raw daily station observations or station based indices.  These are then combined using the
Angular Distance Weighting method of Shepard (1968) to produce final gridded netCDF files.

Alternative gridding and combination schemes have been coded up in these files but were not used in
the initial release.

To reference this material please use:

Dunn, Alexander, Donat,et al, 2020, Development of an updated global land in-situ-based dataset of temperature
and precipitation extremes: HadEX3, JGR-Atmospheres (submitted)

Where indices have been calculated, we have used the Climpact2 scripts from the UNSW, available
here:
https://github.com/ARCCSS-extremes/climpact2


The use of these scripts is unsupported, and may be done so under the licence detailed. 
Documentation is available through the Sphinx auto-documentation for each script.


