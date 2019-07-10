#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 14:34:37 2019

@author: abauville
"""
import numpy as np
from mdoodz.model import Model
from mdoodz.scaling import Scaling

from mdoodz.utils import use_numba, maybe_numba



class TopoChain():
    """
    This class contains data related to topography. 
    Initialization arguments (for default values see __init__ below):
        fact: refinement factor between the grid and the topography vector
        TopoLevel: value at which the z coordinates are initialized
        phase: number at which the phase vector is initialized
        over_alloc_fac: Nb_part_max/Nb_part
        

    The class contains the following attributes:
        Nb_part: number of points defining the topography vectors
        Nb_part_max: number used for allocating the topography vectors
        x: vector containing the x coordinates of the topography
        z: vector containing the z coordinates of the topography
        phase: phase number of the topography nodes
        
        
    """
    def __init__(self,
                 model,
                 fact=24,
                 Topo_level=0.0,
                 over_alloc_fac=2.0,
                 phase=0):
        
        
        # Calculate Nb_part using algorithm from SetTopoChainHorizontalCoords
        dxGrid = model.get_dx()
        dx = dxGrid/fact;
        Nb_part = int(round(model.Nx*fact-(fact-1) - 2))
        
        
        self.Nb_part = Nb_part
        self.Nb_part_max = int(round(over_alloc_fac*self.Nb_part))
        self.x = np.zeros(Nb_part)
        self.z = np.ones(Nb_part) * Topo_level
        self.Vx = np.ones(Nb_part)
        self.Vz = np.ones(Nb_part)
        self.phase = np.zeros(Nb_part,dtype=int)
        
        # self._isScaled = model._isScaled
        
        # Assign topo_chain_x using algorithm from SetTopoChainHorizontalCoords
        for k in range(self.Nb_part):
            self.x[k]     = model.xmin + k*dx + dx;
        
        
    # def scale(self,scaling):
    #     if not isinstance(scaling,Scaling):
    #         raise TypeError("'scale' must be an instance of Scaling")
    #     if self._isScaled == True:
    #         raise ValueError("Trying to scale an already isScaled Polygon")
    #     self.x /= scaling.L
    #     self.z /= scaling.L
    #     self._isScaled = True
        
    # def unscale(self,scaling):
    #     if not isinstance(scaling,Scaling):
    #         raise TypeError("'scale' must be an instance of Scaling")
            
    #     if self._isScaled == False:
    #         raise ValueError("Trying to unscale a non-isScaled Polygon")
    #     self.x *= scaling.L
    #     self.z *= scaling.L 
    #     self._isScaled = False

