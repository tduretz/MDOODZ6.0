#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 11:24:22 2019

@author: abauville
"""

import sys
sys.path.insert(0, '/Users/tduretz/REPO_GIT/MDOODZ6.0/SOURCE/User_interface')
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
model       = Model(#**** RESTART ****#
                    istep = 2,
                    irestart = 0,
                    
                    #**** OUTPUT FILES ****#
                    writer = 1,
                    writer_step = 2,
                    
                    #**** INPUT FILE ****#
                    input_file = "setup_r501x501_fsp40.bin",
                    
                   
                    #**** SPACE-TIME ****#
                    Nx      = 300,
                    Nz      = 250,
                    Nt      = 2,
                    xmin    =  0.000000e0,
                    zmin    = -11.50000000e-2,
                    xmax    =  54.000000e-2,
                    zmax    =  0.8000000e-2,
                    dt          = 2.5e+1,
                    Courant     = 0.5,
                    
                    #**** SWITCHES ****#
                    penalty         = 1e14,
                    eta_avg         = 1,
                    cpc             =-1,
                    surf_remesh     = 0,
                    abs_tol_div     = 1e-11,
                    rel_tol_div     = 1e-11,
                    decoupled_solve = 1,
                    ismechanical    = 1,
                    RK              = 4,
                    free_surf       = 1,
                    free_surf_stab  = 0.15,
                    
                    #**** SETUP DEPENDANT ****#
                    EpsBG           = 0, #Background strain rate
                    user0           = 0, #Activate open BC EAST side
                    user1           = 0, 
                    user2           = 1, #Activate slab dip 
                    user3           = 0,
                    
                    #**** GRAVITY ****#
                    gx = 0.0000,
                    gz = -9.81,
                    
                    #**** MIN/MAX VISCOSITY ****#
                    mineta   = 1.0e0,
                    maxeta   = 1.0e6,
                    
                    #**** BOUNDARY CONDITIONS ****#
                    BC_setup_type = 1
                    )

# Polygons
# ==========================================
ocPlate = geometry.Box([.20 , .54 , -.012 , .0], phase = 0)
slab    = geometry.Box([.20 , .14 , -.012 , .0], phase = 0)
slab.rotate(xc=0.2,angle=34.0*np.pi/180.0)


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
particles.assign_phase_from_geometry(slab)
particles.assign_phase_from_geometry(ocPlate)
print("%.1f s" % (time.time()-tic))



# Material
# ==========================================
slab_mat = Material(model,
                    ID   = 0,
                    rho  = 1495.00,
                    mu   = 1e10,
                    Cv   = 1050,
                    k    = 2.3,
                    Qr   = 1.5e-6,
                    C    = 1e70,
                    phi  = 30,
                    Slim = 500e9,
                    alp  = 10.0e-6,
                    bet  = 1e-11,
                    drho = 0,
                    cstv = 1,
                    pwlv = 0,
                    linv = 0,
                    gbsv = 0,
                    expv = 0,
                    gsel = 0,
                    eta0 = 3.5e5,
                    npwl = 3.3,
                    Qpwl = 186.5e3)

mantle_mat = Material(model,
                      ID   = 1,
                      rho  = 1415.00,
                      mu   = 1e10,
                      Cv   = 1050,
                      k    = 2.3,
                      Qr   = 1.5e-6,
                      C    = 1e70,
                      phi  = 30,
                      Slim = 500e9,
                      alp  = 10.0e-6,
                      bet  = 1e-11,
                      drho = 0,
                      cstv = 1,
                      pwlv = 0,
                      linv = 0,
                      gbsv = 0,
                      expv = 0,
                      gsel = 0,
                      eta0 = 32,
                      npwl = 3.3,
                      Qpwl = 186.5e3)


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

plt.fill(ocPlate.x,ocPlate.z,faceColor='none',edgeColor='r')
plt.fill(slab   .x,slab   .z,faceColor='none',edgeColor='r')

plt.axis("equal")
plt.xlim([model.xmin,model.xmax])

print("%.1f s" % (time.time()-tic))


# Write files
# ==========================================
write.text_file(model, scaling, particles, filename='Input.txt')
write.ini_particles_file(model, particles, topo_chain,filename='Input.dat')