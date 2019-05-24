#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 14:34:37 2019

@author: abauville
"""
import numpy as np
import json
import time
import Geometry

# Global variables
s1Type = np.int32
s2Type = np.float64
s3Type = np.float64
s4Type = np.str

# =============================================================================
# Define classes
class Model():
    def __init__(self,
            Nx=2,
            Nz=2,
            Nt=1,
            xmin=-1.0,
            xmax= 1.0,
            zmin=-1.0,
            zmax= 1.0,
            dt=1.0,
            Courant=0.5,
            isScaled=False):
        
        self.Nx = Nx
        self.Nz = Nz
        self.Nt = Nt
        self.xmin = xmin
        self.xmax = xmax
        self.zmin = zmin
        self.zmax = zmax
        self.EpsBG = 0.0
        self.dt = dt
        self.Courant = Courant
        
        # For breakpoint
        self.dt0    = 1.0
        self.time   = 1.0
        self.L0     = 1.0
        
        # Switches
        self.free_surf = 1
        self.iselastic = 1
        self.fstrain = 1
        self.rec_T_P_x_z = 1
        self.isthermal = 1
        
        # Multi-value switches
        self.eqn_state = 1
    
    
        self.isScaled = isScaled
        
    
    def getDx(self):
        return (self.xmax-self.xmin)/(self.Nx-1)
    def getDz(self):
        return (self.zmax-self.zmin)/(self.Nz-1)
        
    def scale(self,scaling):
        if type(scaling)!=Scaling:
            raise TypeError("'scale' must be an instance of Input.Scaling")
        if self.isScaled == True:
            raise ValueError("Trying to scale an already isScaled Polygon")
        self.xmin /= scaling.L
        self.xmax /= scaling.L
        self.zmin /= scaling.L
        self.zmax /= scaling.L        
        self.isScaled = True
        
    def unscale(self,scaling):
        if type(scaling)!=Scaling:
            raise TypeError("'scale' must be an instance of Input.Scaling")
            
        if self.isScaled == False:
            raise ValueError("Trying to unscale a non-isScaled Polygon")
        self.xmin *= scaling.L
        self.xmax *= scaling.L
        self.zmin *= scaling.L
        self.zmax *= scaling.L    
        self.isScaled = False
        
        
class Scaling():
    def __init__(self,L = 1.0, T = 1.0, eta = 1.0, V = 1.0):
        self.L = L
        self.T = T
        self.eta = eta
        self.V = V

class Topo_chain():
    """
    This class contains data related to topography. 
    Initialization arguments (for default values see __init__ below):
        fact: refinement factor between the grid and the topography vector
        TopoLevel: value at which the z coordinates are initialized
        phase: number at which the phase vector is initialized
        overAllocationFact: Nb_part_max/Nb_part
        

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
                 overAllocationFact=2.0,
                 phase=0):
        
        
        # Calculate Nb_part using algorithm from SetTopoChainHorizontalCoords
        dxGrid = model.getDx()
        dx = dxGrid/fact;
        Nb_part = model.Nx*fact-(fact-1) - 2;
        
        
        self.Nb_part = Nb_part
        self.Nb_part_max = int(round(overAllocationFact*self.Nb_part))
        self.x = np.zeros(Nb_part)
        self.z = np.ones(Nb_part) * Topo_level
        self.Vx = np.ones(Nb_part)
        self.Vz = np.ones(Nb_part)
        self.phase = np.zeros(Nb_part,dtype=int)
        
        self.isScaled = model.isScaled
        
        # Assign topo_chain_x using algorithm from SetTopoChainHorizontalCoords
        for k in range(self.Nb_part):
            self.x[k]     = model.xmin + k*dx + dx;
        
        
    def scale(self,scaling):
        if type(scaling)!=Scaling:
            raise TypeError("'scale' must be an instance of Input.Scaling")
        if self.isScaled == True:
            raise ValueError("Trying to scale an already isScaled Polygon")
        self.x /= scaling.L
        self.z /= scaling.L
        self.isScaled = True
        
    def unscale(self,scaling):
        if type(scaling)!=Scaling:
            raise TypeError("'scale' must be an instance of Input.Scaling")
            
        if self.isScaled == False:
            raise ValueError("Trying to unscale a non-isScaled Polygon")
        self.x *= scaling.L
        self.z *= scaling.L 
        self.isScaled = False

        
class Particles():
    def __init__(self,
                 model, 
                 topo_chain="None",
                 Nx_part = 4,
                 Nz_part = 4,
                 d = 1.0,
                 phi = 0.0,
                 X = 0.0,
                 T = 1.0,
                 phase = 0):
        
        
        # Checks
        if model.free_surf == 1:
            
            if topo_chain=="None":
                raise ValueError("model.free_surf==1 therefore, a non-empty Topo_chain object must be passed.")
        else:
            if topo_chain!="None":
                raise ValueError("model.free_surf==0 but a non-empty Topo_chain object has been passed. Check your inputs for consistency")
        
        self.Nx_part = Nx_part
        self.Nz_part = Nz_part
        self.PutPartInBox(model,topo_chain)
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

    #void PutPartInBox( markers *particles, grid *mesh, params model, surface topo, scale scaling )
    def PutPartInBox(self,model,topo_chain):        
    # Set the original particle layout throughout the mesh.
    # The particlers are set with regular spacing.   
    # The difference with the original C function is that topo_chain is used instead of topo
    
    # add: two pass: first pass assigns Nbpart, second one assigns x and z
    
    
        dx = model.getDx()
        dz = model.getDz()
        # Compute the spacing between self:
        dx_part = dx/(self.Nx_part);
        dz_part = dz/(self.Nz_part);
           
        
        if (model.free_surf==1):

            for iPass in range(2): # pass 1 determines self.Nbpart and assigns x and z; pass 2 fills x and z      
                iP = 0
                for j in range((model.Nx-1)*self.Nx_part):
                    coord_x = model.xmin + dx_part * (j+0.5)

                    I = np.argmin(np.abs(topo_chain.x-coord_x))
                    if (topo_chain.x[I]>coord_x):
                        I=I-1
                    xLoc = coord_x-topo_chain.x[I]
                    b = topo_chain.z[I]
                    a = (topo_chain.z[I+1]-topo_chain.z[I])/(topo_chain.x[I+1]-topo_chain.x[I])
                    h = b + a*xLoc;
                    for i in range((model.Nz-1)*self.Nz_part):
                        coord_z = model.zmin + dz_part * (i+0.5)
                        
                        
    
                        
                        if ( coord_z < h ) :
                            if iPass==1:
                                self.x[iP]  = coord_x;
                                self.z[iP]  = coord_z;
                            iP+=1
                        else:
                            break
                    # end for i
                # end for j
            
                self.Nb_part = iP
                if iPass==0:
                    self.x = np.zeros(self.Nb_part)
                    self.z = np.zeros(self.Nb_part)
            # end for iPass
            


        else: # if model.free_surf
            self.Nb_part = (model.Nx-1)*(model.Nz-1)*self.Nx_part*self.Nz_part
            self.x = np.zeros(self.Nb_part)
            self.z = np.zeros(self.Nb_part)
            iP = 0             
            for j in range((model.Nx-1)*self.Nx_part):
                coord_x = model.xmin + dx_part * (j+0.5)
                for i in range((model.Nz-1)*self.Nz_part):
                    self.x[iP]  = coord_x;
                    self.z[iP]  = model.zmin + dz_part * (i+0.5)# + dz_*(i + 0.5)  
                    iP+=1
                # for i
            # for j
        
        # end if model.free_surf
    # end def PutPartInBox
    
    
    def assignPhaseFromPolygon(
        self,
        polygon,
        atol = 1e-6,
        ):
        if not isinstance(polygon,Geometry.Polygon):
            raise ValueError("'Polygon' must be an instance of 'Geometry.Polygon'")
        polygon.assignPhaseToParticles(self,atol=atol)
    
    
# end class Particles
    
    
    
    


# File writing functions 
# =============================================================================
def writeJsonFile(topo_chain,
                  filename="input.json"):
    
    myJsonFile = dict(TopoChain = vars(topo_chain))   
    with open(filename, "w") as write_file:
        json.dump(myJsonFile, write_file , indent=4, sort_keys=True, separators=(',', ': '), ensure_ascii=False, cls=NumpyEncoder)

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


      
def writeIniParticles(model,topo_chain, particles, filename = "IniParticles.dat"):

    if model.isScaled:
        raise ValueError("the model should not be scaled to write the input file. (This error raising should disappear when development is mature)")
        
    s1 = 4 # sizeof(int)
    s2 = 8 # sizeof(double)
    s3 = 8 # sizeof(DoodzFP)
    s4 = 1 # sizeof(char)
    
    
#    print("Writing file...")
    
    with open(filename, "w") as file:
        myArray = np.array([s1,s2,s3,s4],dtype=np.int8)
        myArray.tofile(file)
        #fwrite( &s1, 1, 1, file);
        #fwrite( &s2, 1, 1, file);
        #fwrite( &s3, 1, 1, file);
        #fwrite( &s4, 1, 1, file);    
        
        if model.free_surf == 1:
            myArray = np.array([topo_chain.Nb_part],dtype=s1Type)   
            myArray.tofile(file)
            #    fwrite( &topo_chain->Nb_part,  s1,                   1, file );
                
            topo_chain.x.astype(s3Type).tofile(file)
            topo_chain.z.astype(s3Type).tofile(file)
#            topo_chain.Vx.astype(s3Type).tofile(file)
#            topo_chain.Vz.astype(s3Type).tofile(file)
            #    fwrite( topo_chain->x,         s3, topo_chain->Nb_part, file );
            #    fwrite( topo_chain->z,         s3, topo_chain->Nb_part, file );
            #    fwrite( topo_chain->Vx,        s3, topo_chain->Nb_part, file );
            #    fwrite( topo_chain->Vz,        s3, topo_chain->Nb_part, file );
    
    
    
    
        myArray = np.array([particles.Nb_part],dtype=s1Type)   
        myArray.tofile(file)
        #fwrite( &particles->Nb_part,  s1, 1, file);
        particles.x.astype(s3Type).tofile(file)
        particles.z.astype(s3Type).tofile(file)
#        particles.P.astype(s3Type).tofile(file)
#        particles.Vx.astype(s3Type).tofile(file)
#        particles.Vz.astype(s3Type).tofile(file)
#        particles.phi.astype(s3Type).tofile(file)
#        particles.X.astype(s3Type).tofile(file)
        particles.phase.astype(s1Type).tofile(file)
        #fwrite( particles->x,     s3, particles->Nb_part, file);
        #fwrite( particles->z,     s3, particles->Nb_part, file);
        #fwrite( particles->P,     s3, particles->Nb_part, file);
        #fwrite( particles->Vx,    s3, particles->Nb_part, file);
        #fwrite( particles->Vz,    s3, particles->Nb_part, file);
        #fwrite( particles->phi,   s3, particles->Nb_part, file);
        #fwrite( particles->X  ,   s3, particles->Nb_part, file);
        #fwrite( particles->phase, s1, particles->Nb_part, file);
        #
    
    
#    print("file writing took %.2f s\n" % (time.time()-tic))
    
# =============================================================================
