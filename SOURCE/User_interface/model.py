#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 11:30:07 2019

@author: abauville
"""
from scaling import Scaling

class Model():
    def __init__(self,
         # Simulation start/restart from Breakpoint        
                 istep=0,
                 irestart=0,
         # Output
                 writer=0,
                 writer_step=1,
                 writer_markers=0,
                 writer_debug=0,
                 writer_energies=0,
         # Input
                 input_file="blah.bin",
         # Domain size       
                 Nx=10,
                 Nz=10,
                 Nt=1,
                 xmin=-1.0,
                 xmax= 1.0,
                 zmin=-1.0,
                 zmax= 1.0,
                 dt=0.0,
                 Courant=0.5,
         # Numerics
                 penalty=1.0e10,
                 abs_tol_div=1.0e-14,
                 rel_tol_div=1.0e-5,
                 auto_penalty=0.0,
                 decoupled_solve=1,
                 
         # Switches
                 ismechanical=1,
                 dt_constant=0,
                 RK=4,
                 isperiodic_x=0,
                 ispureshear_ALE=0,
                 isinertial=0,
                 iselastic=0,
                 isthermal=0,
                 line_search=0,
                 free_surf=0,
                 free_surf_stab=0,
                 eqn_state=0,
                 thermal_eq=0,
                 cooling_time=1000e6*365.25*3600*24,#/scaling->t
                 subgrid_diff=0,
                 shear_heat=1,
                 adiab_heat=0,
                 isPl_soft=0,
                 surf_processes=0,
                 surf_remesh=1,
                 cpc=1,
                 advection=1,
                 loc_iter=1,
                 therm_pert=0,
                 fstrain=0,
                 cut_noise=0,
                 accu=1.0e13,
                 rheo_on_cells=0,
                 DefectCorrectionForm=0,
                 HsOnly=0,
                 HomoFields=0,
                 rec_T_P_x_z=0,
                 rm_break=1,
                 eta_VP=0.0,
                 topografix=0,
                 aniso=0,
                 
         # Setup Dependent
                 EpsBG=0.0,
                 PrBG=0.0,
                 
         # Surface processes 
                 surf_diff=0.0,
                 surf_ised1=0.0,
                 surf_ised2=0.0,
                 surf_sedirate=0.0,
                 surf_baselev=0.0,
                 
         # Initial thermal perturbation
                 therm_pert_x0=0.0,
                 therm_pert_z0=0.0,
                 therm_pert_rad=0.0,
                 therm_pert_dT=0.0,
                 
         # For rheological database reasons...
                 force_act_vol_ast=0,
                 act_vol_dis_ast=0.0,
                 act_vol_dif_ast=0.0,
                 
         #  Model user's delights
                 user0=0.0,
                 user1=0.0,
                 user2=0.0,
                 user3=0.0,
                 user4=0.0,
                 user5=0.0,
                 user6=0.0,
                 user7=0.0,
                 user8=0.0,
                 
         # Derived quantities
                 # dx, dz, dt0, p_avg: defined in C
                 eta_avg=0, # 0 : arithmetic mean
                 
         # Gravity
                 gx=0.0,
                 gz=0.0,
                 
     # DEFORMATION MAP PARAMETERS
         # Create deformation maps or not (default no)
                 def_maps=0,
         # Resolution
                 nT=11,
                 nE=11,
                 nd=11,
         # Temperature, strain rate, grain size MIN/MAX & pressure
                 Tmin=100.0,
                 Tmax=1000.0,
                 Emin=-30.0,
                 Emax=-4.0,
                 dmin=-7.0,
                 dmax=-2.0,
                 Pn=5.0e8,
                 
                 
         # Nonlinear iteration parameters
                 Newton=0,
                 nit_max=1,
                 tol_u=5.0e-6,
                 tol_p=5.0e-6,
                 mineta=1.0e18,
                 maxeta=1.0e24,
         # Direct solver parameters
                 lsolver=0
                 ):
        
        
        
        
# Simulation start/restart from Breakpoint        
        self.istep = istep
        self.irestart = irestart
# Output
        self.writer = writer
        self.writer_step = writer_step
        self.writer_markers = writer_markers
        self.writer_debug = writer_debug
        self.writer_energies = writer_energies
# Input
        self.input_file = input_file
# Domain size       
        self.Nx = Nx
        self.Nz = Nz
        self.Nt = Nt
        self.xmin = xmin
        self.xmax = xmax
        self.zmin = zmin
        self.zmax = zmax
        self.dt = dt
        self.Courant = Courant
# Numerics
        self.penalty = penalty
        self.abs_tol_div = abs_tol_div 
        self.rel_tol_div = rel_tol_div
        self.auto_penalty = auto_penalty
        self.decoupled_solve = decoupled_solve
         
# Switches
        self.ismechanical = ismechanical
        self.dt_constant = dt_constant
        self.RK = RK
        self.isperiodic_x = isperiodic_x
        self.ispureshear_ALE = ispureshear_ALE
        self.isinertial = isinertial
        self.iselastic = iselastic
        self.isthermal = isthermal
        self.line_search = line_search
        self.free_surf = free_surf
        self.free_surf_stab = free_surf_stab
        self.eqn_state = eqn_state
        self.thermal_eq = thermal_eq
        self.cooling_time = cooling_time
        self.subgrid_diff = subgrid_diff
        self.shear_heat = shear_heat
        self.adiab_heat = adiab_heat
        self.isPl_soft = isPl_soft
        self.surf_processes = surf_processes
        self.surf_remesh = surf_remesh
        self.cpc = cpc
        self.advection = advection
        self.loc_iter = loc_iter
        self.therm_pert = therm_pert
        self.fstrain = fstrain
        self.cut_noise = cut_noise
        self.accu = accu
        self.rheo_on_cells = rheo_on_cells
        self.DefectCorrectionForm = DefectCorrectionForm
        self.HsOnly = HsOnly
        self.HomoFields = HomoFields
        self.rec_T_P_x_z = rec_T_P_x_z
        self.rm_break = rm_break
        self.eta_VP = eta_VP
        self.topografix = topografix
        self.aniso = aniso
         
# Setup Dependent
        self.EpsBG = EpsBG
        self.PrBG = PrBG
         
# Surface processes 
        self.surf_diff = surf_diff
        self.surf_ised1 = surf_ised1
        self.surf_ised2 = surf_ised2
        self.surf_sedirate = surf_sedirate
        self.surf_baselev = surf_baselev
        
# Initial thermal perturbation
        self.therm_pert_x0 = therm_pert_x0
        self.therm_pert_z0 = therm_pert_z0
        self.therm_pert_rad = therm_pert_rad
        self.therm_pert_dT = therm_pert_dT
         
# For rheological database reasons...
        self.force_act_vol_ast = force_act_vol_ast
        self.act_vol_dis_ast = act_vol_dis_ast
        self.act_vol_dif_ast = act_vol_dif_ast
        
 #  Model user's delights
        self.user0 = user0
        self.user1 = user1
        self.user2 = user2
        self.user3 = user3
        self.user4 = user4
        self.user5 = user5
        self.user6 = user6
        self.user7 = user7
        self.user8 = user8
         
# Derived quantities
        # dx, dz, dt0, p_avg: defined in C
        self.eta_avg = eta_avg, # 0 : arithmetic mean
         
# Gravity
        self.gx = gx,
        self.gz = gz,
         
# DEFORMATION MAP PARAMETERS
# Create deformation maps or not (default no)
        self.def_maps = def_maps
# Resolution
        self.nT = nT
        self.nE = nE
        self.nd = nd
# Temperature, strain rate, grain size MIN/MAX & pressure
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Emin = Emin
        self.Emax = Emax
        self.dmin = dmin
        self.dmax = dmax
        self.Pn = Pn
        
        
# Nonlinear iteration parameters
        self.Newton = Newton
        self.nit_max = nit_max
        self.tol_u = tol_u
        self.tol_p = tol_p
        self.mineta = mineta
        self.maxeta = maxeta
# Direct solver parameters
        self.lsolver = lsolver
    
        self.isScaled=False
    
    def get_dx(self):
        return (self.xmax-self.xmin)/(self.Nx-1)
        
    def get_dz(self):
        return (self.zmax-self.zmin)/(self.Nz-1)
        
    def scale(self,scaling):
        if not isinstance(scaling,Scaling):
            raise TypeError("'scaling' must be an instance of Scaling")
        self.xmin /= scaling.L
        self.xmax /= scaling.L
        self.zmin /= scaling.L
        self.zmax /= scaling.L
