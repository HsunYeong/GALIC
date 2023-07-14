#*******************************************************************************
# This file is part of the GALIC code developed by D. Yurin and V. Springel.
#
# Copyright (c) 2014
# Denis Yurin (denis.yurin@h-its.org) 
# Volker Springel (volker.springel@h-its.org)
#*******************************************************************************
#!/bin/bash            # this line only there to enable syntax highlighting in this file


#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1
#OUTPUT_IN_DOUBLEPRECISION # snapshot files will be written in double precision


#--------------------------------------- Output/Input options
HAVE_HDF5                     # needed when HDF5 I/O support is desired


#DEBUG_ENABLE_FPU_EXCEPTIONS   #enables floating point exceptions

#--------------------------------------- Special behaviour
#RADIAL_WEIGHTING_IN_DENSITY_RESPONSE

##---------------------------- Modifications
#VER_1_1   #enables version GALIC 1.1 with velocity dispertions patch
#VAR_1_1_KPARAMETER_MOD # changes vstr = k*vphi (experimental!)
#VER_1_1_GNUPLOT_LOG

