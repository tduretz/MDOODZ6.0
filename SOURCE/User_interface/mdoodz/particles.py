#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 14:34:37 2019

@author: abauville
"""
import numpy as np
import mdoodz.geometry as geometry
from mdoodz.model import Model
from mdoodz.utils import use_numba, maybe_numba

class Particles():
    def __init__(self,
                 model, 
                 topo_chain=None,
                 Nx_part = 4,
                 Nz_part = 4,
                 min_part_cell=16,
                 over_alloc_fac=4.1,
                 d = 1.0,
                 phi = 0.0,
                 X = 0.0,
                 T = 1.0,
                 phase = 0):
        
        if not isinstance(model,Model):
            raise ValueError("model must be an instance of mdoodz.model.Model")
        
        # Checks
        if model.free_surf == 1:
            if topo_chain is None:
                raise ValueError("model.free_surf==1 therefore, a non-empty Topo_chain object must be passed.")
        else:
            if topo_chain is not None:
                raise ValueError("model.free_surf==0 but a non-empty Topo_chain object has been passed. Check your inputs for consistency")
        
        if min_part_cell<Nx_part*Nz_part:
            raise ValueError("min_part_cell should be lower than Nx_part*Nz_part")
        
        
        self.Nx_part = Nx_part
        self.Nz_part = Nz_part
        self.min_part_cell = min_part_cell
        
        self._xmin = model.xmin
        self._xmax = model.xmax
        self._zmin = model.zmin
        self._zmax = model.zmax
        
        
        # Assign Nb_part, x, z
        self.x = None
        self.z = None
        self.Nb_part = 0
        if model.free_surf==1:
            self._init_xz_free_surf(model,topo_chain) 
        else:            
            self._init_xz(model,topo_chain)             
        
        self.Nb_part_max = self.Nb_part * over_alloc_fac
        
        
#        self.Nb_part = (part_per_cell*(model.Nx-1))*(part_per_cell*(model.Nz-1))
#        
#        self.x, self.z = np.meshgrid(np.linspace(model.xmin,model.xmax,(model.Nx-1)*part_per_cell),
#                                     np.linspace(model.zmin,model.zmax,(model.Nz-1)*part_per_cell))
#        self.x      = self.x.flatten()
#        self.z      = self.z.flatten()
        self.P      = np.ones(self.Nb_part)
        self.Vx     = -1.0*self.x*model.EpsBG;
        self.Vz     =  1.0*self.z*model.EpsBG;
        self.d      = d*np.ones(self.Nb_part)                           # same grain size everywhere
        self.phi    = phi*np.ones(self.Nb_part)                            # zero porosity everywhere
        self.X      = X*np.ones(self.Nb_part)                              # X set to 0
        self.T      = T*np.ones(self.Nb_part)                           # same temperature everywhere
    
        self.phase  = phase*np.ones(self.Nb_part,dtype=int)
  
        self.isScaled = model.isScaled

    #void init_xz( markers *particles, grid *mesh, params model, surface topo, scale scaling )
    def _init_xz(self,model,topo_chain):        
    # Set the original particle layout throughout the mesh.
    # The particlers are set with regular spacing.   
    # The difference with the original C function is that topo_chain is used instead of topo
    
    # add: two pass: first pass assigns Nbpart, second one assigns x and z
    
        @maybe_numba(use_numba)
        def main_loop(model_Nx, model_Nz, model_Dx, model_Dz, 
                      model_xmin, model_zmin,
                      part_Nx, part_Nz): 
    
            dx = model_Dx
            dz = model_Dz
            # Compute the spacing between self:
            dx_part = dx/(part_Nx);
            dz_part = dz/(part_Nz);
               
            
            Nb_part = (model_Nx-1)*(model_Nz-1)*part_Nx*part_Nz
            part_x = np.zeros(Nb_part)
            part_z = np.zeros(Nb_part)
            iP = 0             
            for j in range((model_Nx-1)*part_Nx):
                coord_x = model_xmin + dx_part * (j+0.5)
                for i in range((model_Nz-1)*part_Nz):
                    part_x[iP]  = coord_x;
                    part_z[iP]  = model_zmin + dz_part * (i+0.5)# + dz_*(i + 0.5)  
                    iP+=1
                # for i
            # for j
            return (part_x, part_z)
        # end function main_loop
        self.x, self.z = main_loop(
                      model.Nx, model.Nz, model.get_dx(), model.get_dz(), 
                      model.xmin, model.zmin,
                      self.Nx_part, self.Nz_part)
        self.Nb_part = self.x.size
    # end def init_xz
    
    
    def _init_xz_free_surf(self,model,topo_chain): 
    
    # Set the original particle layout throughout the mesh.
    # The particlers are set with regular spacing.   
    # The difference with the original C function is that topo_chain is used instead of topo
    
    # add: two pass: first pass assigns Nbpart, second one assigns x and z
        @maybe_numba(use_numba)
        def main_loop(model_Nx, model_Nz, model_Dx, model_Dz, 
                      model_xmin, model_zmin,
                      part_Nx, part_Nz, 
                      topo_chain_x, topo_chain_z):  

            part_x = np.zeros(1)
            part_z = np.zeros(1)
            dx = model_Dx
            dz = model_Dz
            # Compute the spacing between self:
            dx_part = dx/(part_Nx);
            dz_part = dz/(part_Nz);
    
            for iPass in range(2): # pass 1 determines self.Nbpart and assigns x and z; pass 2 fills x and z      
                iP = 0
                for j in range((model_Nx-1)*part_Nx):
                    coord_x = model_xmin + dx_part * (j+0.5)
    
                    I = np.argmin(np.abs(topo_chain_x-coord_x))
                    if (topo_chain_x[I]>coord_x):
                        I=I-1
                    xLoc = coord_x-topo_chain_x[I]
                    b = topo_chain_z[I]
                    a = (topo_chain_z[I+1]-topo_chain_z[I])/(topo_chain_x[I+1]-topo_chain_x[I])
                    h = b + a*xLoc;
                    for i in range((model_Nz-1)*part_Nz):
                        coord_z = model_zmin + dz_part * (i+0.5)
                        
                        if ( coord_z < h ) :
                            if iPass==1:
                                part_x[iP]  = coord_x;
                                part_z[iP]  = coord_z;
                            iP+=1
                        else:
                            break
                    # end for i
                # end for j
            
                Nb_part = iP
                if iPass==0:
                    part_x = np.zeros(Nb_part)
                    part_z = np.zeros(Nb_part)
            # end for iPass
            return (part_x, part_z)
        # end function Local

        self.x, self.z = main_loop(
                      model.Nx, model.Nz, model.get_dx(), model.get_dz(), 
                      model.xmin, model.zmin,
                      self.Nx_part, self.Nz_part, 
                      topo_chain.x, topo_chain.z)
        self.Nb_part = self.x.size
    # end def init_xz
    
    
    def assign_phase_from_geometry(
        self,
        geometryObject,
        atol = 1e-6,
        ):
        print(type(geometryObject) )
        if not isinstance(geometryObject,(geometry.Polygon,geometry.BinaryNoise)):
            raise ValueError("'geometryObject' must be an instance of 'mdoodz.geometry.Polygon' or 'mdoodz.geometry.BinaryNoise'")
        geometryObject.assign_phase_to_particles(self,atol=atol)
        
    
    
# end class Particles
