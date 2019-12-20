# enconding: utf-8
##########################################
#*************************************************************************
#*  Copyright (C) 2016 by Shiwei Zhao                                     *
#*  swzhao@scut.edu.cn                                                    *
#*  Programming Language: Python 2                                        *
#*  This program is free software; it is licensed under the terms of the  *
#*  GNU General Public License v2 or later. See file LICENSE for details. *
#*************************************************************************/

##############################################################################################
##This script was written for the visulization of 3D histograms in the following two
##publications. If you find it useful and use it in your research, please consider citing 
##one of the publications.
##[1] Zhao, S, and Zhou, X.(2017). Effects of Particle Asphericity on the Macro-
##and Micro-Mechanical Behaviors of Granular Assemblies. Granular Matter 19 (2):38.
##[2] Zhao, S, Evans, TM, and Zhou X.(2018). Three-Dimensional Voronoi
##Analysis of Monodisperse Ellipsoids during Triaxial Shear. Powder Technology 323:323-36.
##########################################################################
###example
###Note: the following scripts can be grouped into another individual file, and then you need to
###import the module ahead by "from vtk3Dhistogram import *".
from vtk3Dhistogram import *
myvtk = VTK3Dhist()
file_input = './test.dat'
file_output = './contactnormal'
myvtk.VTK3Dhistogram_basic(file_input,file_output)

