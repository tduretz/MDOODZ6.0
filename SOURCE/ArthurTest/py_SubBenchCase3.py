#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 11:24:22 2019

@author: abauville
"""
import numpy as np
import matplotlib.pyplot as plt
import time

# import mdoodz classes
from mdoodz import Material, Model, Particles, Scaling, TopoChain

# import mdoodz modules (modules contain classes and functions)
from mdoodz import write, geometry

fastPlotting = True
if fastPlotting:
    try:
        import pandas as pd
        import datashader as ds    
    except:
        fastPlotting = False

# Scaling 
# ==========================================
scaling     = Scaling(eta = 1e4,
                      L   = 1e-1,
                      V   = 1e2,
                      T   = 20)


# Model
# ==========================================

model       = Model(Nx      = round(300 *1.0),
                    Nz      = round(250 *1.0),
                    Nt      = 1,
                    xmin    =  0.000000e0,
                    zmin    = -11.50000000e-2,
                    xmax    =  54.000000e-2,
                    zmax    =  0.8000000e-2,
                    dt      = 2.5e+1,
                    Courant = 0.5,
                    free_surf=1
                    )

# Polygons
# ==========================================
OcPlate = geometry.Box([.20 , .54 , -.012 , .0], phase = 0)
Slab    = geometry.Box([.20 , .14 , -.012 , .0], phase = 0)
Slab.rotate(xc=0.2,angle=34.0*np.pi/180.0)


# Scale everything
# ==========================================
#model.scale(scaling)
#Slab.scale(scaling)
#OcPlate.scale(scaling)


# Topo chain
# ==========================================
topo_chain  = TopoChain(model,Topo_level=0.000,fact=24)

# Particles
# ==========================================
print("Initializing particles: ", end = ''); tic = time.time()
Tbg  = 773.15/scaling.T;                         # reference temperature
gsbg = 2e-3/scaling.L;                           # reference grain size
H    = 1.2e-2/scaling.L;                         # plate thickness

particles   = Particles(model, 
                              topo_chain=topo_chain,
                              d = gsbg,
                              phi = 0.0,
                              X = 0.0,
                              T = Tbg,
                              phase = 1)

print("%.1f s" % (time.time()-tic))

# ==== Assign phase from polygons
print("Assigning phase to particles: ", end = ''); tic = time.time()
particles.assign_phase_from_polygon(Slab)
particles.assign_phase_from_polygon(OcPlate)
print("%.1f s" % (time.time()-tic))


## Plot
## ==========================================
print("Plotting: ", end = ''); tic = time.time();
plt.clf()    
if fastPlotting:
    df = pd.DataFrame(np.array([particles.x,particles.z,particles.phase]).T,columns=('x','y','phase'))
    refineFac = 2
    cvs = ds.Canvas(plot_width=model.Nx*refineFac, plot_height=model.Nz*refineFac,
                       x_range=(model.xmin,model.xmax), y_range=(model.zmin,model.zmax),
                       x_axis_type='linear', y_axis_type='linear')
    agg = cvs.points(df, 'x', 'y', ds.mean('phase'))
    plt.imshow(agg.variable,extent=[model.xmin,model.xmax,model.zmin,model.zmax],origin='lower')

else:    
    plt.scatter(particles.x,particles.z,c=particles.phase)

cb=plt.colorbar()

plt.fill(OcPlate.x,OcPlate.z,faceColor='none',edgeColor='r')
plt.fill(Slab   .x,Slab   .z,faceColor='none',edgeColor='r')

plt.axis("equal")
plt.xlim([model.xmin,model.xmax])

print("%.1f s" % (time.time()-tic))


# Write files
# ==========================================
#print("Writing particle file: ", end = ''); tic = time.time()
#Write.writeIniParticles(model,topo_chain, particles)
#print("%.1f s" % (time.time()-tic))
