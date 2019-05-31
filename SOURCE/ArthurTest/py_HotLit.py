#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 17:20:14 2019

@author: abauville
"""

#void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {
    
#int np;

import numpy as np
import matplotlib.pyplot as plt
import time

# import mdoodz classes
from mdoodz import Material, Model, Particles, Scaling, TopoChain

# import mdoodz modules (modules contain classes and functions)
from mdoodz import write, geometry

try:
    from numba import jit #, prange
    useNumba = True
except:
    print("Warning: the package numba was not found. just-in-time compilation is off")
    useNumba = False
    
def maybe_numba(useNumba):
    return jit(nopython=True) if useNumba else lambda x:x


print("Start!!")
wallTimeTic = time.time()
fastPlotting = True
if fastPlotting:
    try:
        import pandas as pd
        import datashader as ds    
        print("fastPlotting is on baby!")
    except:
        fastPlotting = False
        print("fastPlotting is off... too bad")




# Scaling 
# ==========================================
scaling     = Scaling()


# Model
# ==========================================
print("Initializing model: ", end = ''); tic = time.time()
model       = Model(Nx      = int(2000/1),#int(400),
                    Nz      = int(2000/1),#int(220),
                    Nt      = 1000,
                    xmin    = -100.000000e3,
                    zmin    = -95.000000e3,
                    xmax    =  100.000000e3,
                    zmax    = 5.000000e3 ,
                    dt      = 1.2623e+11,
                    Courant = 0.3,
                    penalty = 1e19,
                    eta_avg = 1,
                    cpc     = 1,
                    surf_remesh = 1,
                    abs_tol_div = 1e-14,
                    rel_tol_div = 1e-14,
                    DefectCorrectionForm = 1,
                    
                    EpsBG           = 5.e-16, #Background strain rate
                    user0           = 2.0e3,
                    user1           = 0,
                    user2           = 0,
                    user3           = 0
                    )


print("%.1f s" % (time.time()-tic))

# Particles
# ==========================================
print("Initializing particles: ", end = ''); tic = time.time()
Tbg  = 773.15/scaling.T;                         # reference temperature
gsbg = 2e-3/scaling.L;                           # reference grain size
H    = 1.2e-2/scaling.L;                         # plate thickness
model.free_surf=0
particles   = Particles(model, 
#                        topo_chain=topo_chain,
                        d = gsbg,
                        phi = 0.0,
                        X = 0.0,
                        T = Tbg,
                        phase = 1)
print("%.1f s" % (time.time()-tic))


print("Jpoh part 1: ", end = ''); tic = time.time()
# Define dimensions                          
Lx      = model.xmax - model.xmin
Lz      = model.zmax - model.zmin
gsbg    = 2.0e-3/scaling.L;                            # reference grain size
HLit    = 98.e3/scaling.L;                           # lithosphere thickness
HCrust  = 30.e3/scaling.L;                            # crust thickness
Huc     = 2.5e3/scaling.L;                             # Upper crust
Hlc     = HCrust - Huc;                            # Lower crust
#
#Tsurf   = 0.0/scaling.T, Tpart;                   # surface temperature
#Tmant   = (1330.0+zeroC)/scaling.T;                   # adiabatic mantle temperature
rad=model.user0/scaling.L 
la= 1.0*rad
sa = 1.0*rad
theta=(0.0)*np.pi/180;             # Dimensions for granite seed

costheta = np.cos(theta)
sintheta = np.sin(theta)

#printf("Rad: \t %.10f\n", rad);
spacing = 2.0e3/scaling.L;
lc     = 14.0e3/scaling.L;
#printf("Model paramters: \t %.4f \t %.4f \t %.4f \t %.4f\n", model.xmax, model.xmin, model.zmax, model.zmin);
#printf("Model extents: \t %.4f \t %.4f\n", Lx, Lz);


# JPoh Addition:===
#* Perturbation switches *#
rand_pert = 0;
layering = 0;
passive_marker = 1;
weak_seed = 1;
if rad==0.0:
    weak_seed=0

#* Creation of crustal seed markers *#
seed_number = 15;                               # number of circles
zc_1 = -8.0e3/scaling.L;
zc_2 = -20.0e3/scaling.L;
xc_freq = Lx / seed_number;
#X, Z, Xn, Zn, X2, Z2, Xn2, Zn2;
#int i;

# Crustal weak seed */
seed_x = -10.0e3/scaling.L;
seed_z = -8.5e3/scaling.L;
#seed_X, seed_Z, seed_Xn, seed_Zn;

# Vectors to store values
xc_Layer1 = np.zeros(seed_number)
xc_Layer2 = np.zeros(seed_number)
zc_Layer1 = np.zeros(seed_number)
zc_Layer2 = np.zeros(seed_number)

#srand((unsigned) time(NULL));
for i in range(seed_number):#( i = 0; i <= (seed_number); i++) {
    xc_Layer1[i] = (model.xmin + 5.e3/scaling.L) + (i * xc_freq);
    xc_Layer2[i] = (model.xmin + 10.e3/scaling.L) + (i * xc_freq);
        
    zc_Layer1[i] = zc_1;
    zc_Layer2[i] = zc_2;
        



#* Creating noise layer at the Moho *#
hi_num = 2.0;
low_num = 1.0;
line_z = HCrust;                            # Change value here to place your noise layer
dx_line = Lx/(model.Nx - 1);
rand_array = np.zeros(model.Nx)#[model.Nx];
rand_array2 = np.zeros(model.Nx)
rand_line = np.zeros(model.Nx)
#m, c, y, y1, y2, x1, x2;
#m_uc, c_uc, y_uc, y1_uc, y2_uc, x1_uc, x2_uc;

#srand((unsigned) time(NULL));
for i in range(model.Nx):
#for ( i = 0; i <= (model.Nx-1) ; i++) {
    # Generate noise layer for crust
    rand_array[i] = ((np.random.rand() * hi_num - low_num)) * (1.e3/scaling.L);
#    rand_array[i] = ((((double)rand() / (RAND_MAX)) * hi_num - low_num)) * (1.e3/scaling.L);
    rand_array[i] *= 2.;                           # Amplitude of error
    rand_array[i] += line_z;
    
    # Noise layer for 2nd random line
    rand_array2[i] = ((np.random.rand() * hi_num - low_num)) * (1.e3/scaling.L);
    rand_array2[i] *= 1.;                           # Amplitude of error
    rand_array2[i] += zc_2;
    #rand_array[i] = rand_array[i]/scaling.L;
    rand_line[i] = model.xmin;
    if (i > 0) :
        rand_line[i] = rand_line[i-1] + dx_line;

    

#}
        
        
print("%.1f s" % (time.time()-tic))     
        
        
        

#--------------------------------------------------#

@maybe_numba(useNumba)
def particles_assignPhase(particles_x,particles_z,model_Nx):
    Nb_part = particles_x.size
    Phase = np.zeros(Nb_part)
    # Loop over particles
    for ip in range(Nb_part):
        px = particles_x[ip]
        pz = particles_z[ip]
        
        # POTATO_SARAH: ==
        if (pz<-(HLit)): phase = 3;
        if (pz>-(HLit)): phase = 2;
        if (pz>-(HCrust)): phase = 1;
        
        # Initiate noise layer between Moho and crust*/
        if (rand_pert == 1) :
            for i in range(model_Nx):
                if (px> rand_line[i] and px<= rand_line[i+1]):
                    y1 = -(rand_array[i]);
                    y2 = -(rand_array[i+1]);
                    x1 = rand_line[i];
                    x2 = rand_line[i+1];
                    
                    m = (y1 - y2) / (x1 - x2);               # Calculate gradient of line
                    c = y1 - (m * x1);                       # Calculation of the y-intercept
                    y = m * px + c;            # Determine where the marker would be on the line
                    
                    if (pz>=y):
                        phase = 1;
                                    
                    elif (pz<y and pz>=-(HCrust)):
                        phase = 2;
                    
                # end if
            # end for i in model.Nx
        # end if rand_pert
    
        # Add upper crust layer
        if (pz>-(Huc)): phase = 0;
    
        
        # Passive markers
        if (passive_marker == 1) :
            if (phase == 1) :
                pz_mod = (pz/spacing)%2
                if pz_mod<1.0 and pz_mod>=0.0:
                    phase = 4            
        
        # Crustal weak seed
        if (weak_seed == 1) :
            seed_X = (px-seed_x)/rad;
            seed_Z = (pz-seed_z)/rad;
            seed_Xn = seed_X*costheta + seed_Z*sintheta;
            seed_Zn = seed_X*sintheta - seed_Z*costheta;
            
            if ( seed_Xn*seed_Xn + seed_Zn*seed_Zn - 1.0 < 0.0) :
                phase = 5;
    
        Phase[ip] = phase
    
    # end loop over particles
    return Phase
#end function particles_assignPhase()  

print("Particles: assign phase: ", end = ''); tic = time.time()    
particles.phase = particles_assignPhase(particles.x,particles.z,model.Nx)
#particles.phase = particles_assignPhase(X,Z,Nx)
print("%.1f s" % (time.time()-tic))

## Plot
## ==========================================
print("Plotting: ", end = ''); tic = time.time();
plt.clf()    
if fastPlotting:
    df = pd.DataFrame(np.array([particles.x,particles.z,particles.phase]).T,columns=('x','y','phase'))
    refineFac = 1
    cvs = ds.Canvas(plot_width=model.Nx*refineFac, plot_height=model.Nz*refineFac,
                       x_range=(model.xmin,model.xmax), y_range=(model.zmin,model.zmax),
                       x_axis_type='linear', y_axis_type='linear')
    agg = cvs.points(df, 'x', 'y', ds.mean('phase'))
    plt.imshow(agg.variable,extent=[model.xmin,model.xmax,model.zmin,model.zmax],origin='lower')

else:    
    plt.scatter(particles.x,particles.z,c=particles.phase)

cb=plt.colorbar()

#plt.fill(OcPlate.x,OcPlate.z,faceColor='none',edgeColor='r')
#plt.fill(Slab   .x,Slab   .z,faceColor='none',edgeColor='r')

plt.axis("equal")
plt.xlim([model.xmin,model.xmax])

print("%.1f s" % (time.time()-tic))


print("Total wall time: %.1f s" % (time.time()-wallTimeTic))  