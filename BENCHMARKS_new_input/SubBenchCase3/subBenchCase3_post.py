#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 14:12:08 2019

@author: abauville
"""
import numpy as np
import h5py
import matplotlib.pyplot as plt

from mdoodz import read

filename = '../Output/Output00002.gzip.h5'
#read.printFileTree(filename)
grid = read.Grid(filename)
time = read.Time(filename)
topo = read.Topo(filename)
compo = read.get_viz_grid_data(filename,'compo')


## Plot 
## =======================
plt.clf()
plt.subplot(211)
plt.plot(topo.height)

plt.subplot(212)
plt.imshow(compo,origin='lower',extent=(grid.xmin,grid.xmax,grid.zmin,grid.zmax))



