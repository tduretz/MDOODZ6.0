// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2019  MDOODZ Developper team
//
// This file is part of MDOODZ.
//
// MDOODZ is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MDOODZ is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MDOODZ.  If not, see <http://www.gnu.org/licenses/>.
// =========================================================================

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "time.h"
#include "cholmod.h"
#include "header_MDOODZ.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

#ifdef _VG_
#define printf(...) printf("")
#endif


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int main( int nargs, char *args[] ) {
    
	int          istep, irestart, writer = 0, writer_step;
	char         *fin_name, *PartFileName;
	params       model;
	grid         mesh;
    Nparams      Nmodel;
	markers      particles, topo_chain, topo_chain_ini;
	mat_prop     materials;
	clock_t      t_omp, t_omp_step;
	scale        scaling;
    SparseMat    Stokes, Jacob;
    DirectSolver CholmodSolver;
    surface      topo, topo_ini;
    int          Nmax_picard, nstag, nstagmax=3, safe=0;
    SparseMat    StokesA, StokesB, StokesC, StokesD;
    SparseMat    JacobA,  JacobB,  JacobC,  JacobD;
    int          Nx, Nz, Ncx, Ncz;
    
    // Initialise some of the knobs
    model.Newton          = 1;
    model.lsolver         = 1;
    if ( model.Newton == 1 ) model.lsolver         = 1;
    
    // Initialise integrated quantities
    mesh.W    = 0.0; // Work
    mesh.Ut   = 0.0; // heat
    mesh.Ue   = 0.0; // elastic energy
    
#ifdef _NEW_INPUT_
    // Input file name
    if ( nargs < 3 ) {
        printf( "NEW INPUT: You should enter the setup file and the initial particle file names as command line arguments.\nExiting...\n" );
        exit(1);
    }
    else {
        asprintf(&fin_name,"%s", args[1]);
        asprintf(&PartFileName,"%s", args[2]);
    }
#else
    // Input file name
    if ( nargs < 2 ) {
        printf( "OLD INPUT: You should (at least) enter the setup file name as a command line argument.\nExiting...\n" );
        exit(1);
    }
    else {
        asprintf(&fin_name,"%s", args[1]);
    }
#endif
    
    printf("\n********************************************************\n");
    printf("************ Starting MDOODZ 6.0 simulation ************\n");
    printf("********************************************************\n");
    
    // Read input data
    ReadInputFile( fin_name, &istep, &irestart, &writer, &writer_step, &model, &scaling, &materials, &particles, &Nmodel);
    
    printf("*************************************\n");
    printf("****** Allocate and initialise ******\n");
    printf("*************************************\n");
    
    // Multi resolution allocation pointers
    GridAlloc( &mesh, &model );
    
    // Initialise grid coordinates
    SetGridCoordinates( &mesh, &model, model.Nx, model.Nz );
    
    // Initial solution fields (Fine mesh)
    InitialiseSolutionFields( &mesh, &model );
    
    // Get grid indices
    GridIndices( &mesh );
    
    Nx = mesh.Nx; Nz = mesh.Nz; Ncx = Nx-1; Ncz = Nz-1;
    
    printf("*************************************\n");
    printf("******* Initialize particles ********\n");
    printf("*************************************\n");
    
    // Allocate particle fields
    PartAlloc( &particles, &model );
    
    // Allocate marker chain
    if ( model.free_surf == 1 ) AllocateMarkerChain( &topo,     &topo_chain,     model );
    if ( model.free_surf == 1 ) AllocateMarkerChain( &topo_ini, &topo_chain_ini, model );
    
    // Set new particle distribution
    if ( irestart == 0 ) {
        
        model.step = 0;
        model.time = 0.0;
        
        // Initialise particle fields
        PartInit( &particles, &model );
        
#ifdef _NEW_INPUT_
            // Initial grid tags
            model.BC_setup_type = 1; // eventually it should be set from the input file
            SetBCs_new( &mesh, &model, scaling , &particles, &materials );
            
            LoadIniParticles( PartFileName, &particles, &mesh, &topo_chain, &topo_chain_ini, &model, scaling );
            
            if ( model.free_surf == 1 ) {
                // Project topography on vertices
                ProjectTopography( &topo, &topo_chain, model, mesh, scaling, mesh.xg_coord, 0 );
                
                // Marker chain polynomial fit
                MarkerChainPolyFit( &topo, &topo_chain, model, mesh );
                
                // Call cell flagging routine for free surface calculations
                CellFlagging( &mesh, model, topo, scaling );
            }
            
#else
            // Initial grid tags
            SetBCs( &mesh, &model, scaling , &particles, &materials );
            if ( model.free_surf == 1 ) {
                
                // Define the horizontal position of the surface marker chain
                SetTopoChainHorizontalCoords( &topo,     &topo_chain,     model, mesh, scaling );
                SetTopoChainHorizontalCoords( &topo_ini, &topo_chain_ini, model, mesh, scaling );
                
                // Define the vertical position of the surface marker chain
                BuildInitialTopography( &topo,     &topo_chain,     model, mesh, scaling );
                BuildInitialTopography( &topo_ini, &topo_chain_ini, model, mesh, scaling );
                
                // Project topography on vertices
                ProjectTopography( &topo, &topo_chain, model, mesh, scaling, mesh.xg_coord, 0 );
                
                // Marker chain polynomial fit
                MarkerChainPolyFit( &topo, &topo_chain, model, mesh );
                
                // Call cell flagging routine for free surface calculations
                CellFlagging( &mesh, model, topo, scaling );
            }
            // Set particles coordinates
            PutPartInBox( &particles, &mesh, model, topo, scaling );
            
            // Set phases on particles
            SetParticles( &particles, scaling, model, &materials );
            if ( model.free_surf == 1 ) CleanUpSurfaceParticles( &particles, &mesh, topo, scaling ); /////////!!!!!!!!!
            
#endif
        
        if ( model.free_surf == 1 ) CleanUpSurfaceParticles( &particles, &mesh, topo, scaling ); /////////!!!!!!!!!
        
        Interp_P2N ( particles, materials.eta0,  &mesh, mesh.eta_s, mesh.xg_coord,  mesh.zg_coord, 0, 0, &model );
        Interp_P2C ( particles, materials.eta0,  &mesh, mesh.eta_n, mesh.xg_coord,  mesh.zg_coord, 0, 0 );
        
        Interp_P2C ( particles, particles.T,    &mesh, mesh.T,     mesh.xg_coord,  mesh.zg_coord, 1, 0 );
        Interp_P2C ( particles, materials.Cv,   &mesh, mesh.Cv,       mesh.xg_coord,  mesh.zg_coord, 0, 0 );
        Interp_P2U ( particles, materials.k_eff,    &mesh, mesh.kz,       mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1,     mesh.Nz,   0, mesh.BCv.type , &model);
        Interp_P2U ( particles, materials.k_eff,    &mesh, mesh.kx,       mesh.xg_coord,  mesh.zvx_coord, mesh.Nx,       mesh.Nz+1, 0, mesh.BCu.type , &model);
        
        if ( model.eqn_state > 0) {
            Interp_P2N ( particles, particles.rho,  &mesh, mesh.rho_s, mesh.xg_coord,  mesh.zg_coord, 1, 0, &model );
            Interp_P2C ( particles, particles.rho,  &mesh, mesh.rho_n, mesh.xg_coord,  mesh.zg_coord, 1, 0 );
        }
        else {
            Interp_P2N ( particles, materials.rho,  &mesh, mesh.rho_s, mesh.xg_coord,  mesh.zg_coord, 0, 0, &model );
            Interp_P2C ( particles, materials.rho,  &mesh, mesh.rho_n, mesh.xg_coord,  mesh.zg_coord, 0, 0 );
        }
        
        printf("*************************************\n");
        printf("****** Initialize temperature *******\n");
        printf("*************************************\n");
        
#ifdef _NEW_INPUT_
        SetBCs_new( &mesh, &model, scaling , &particles, &materials );
#else
        SetBCs( &mesh, &model, scaling , &particles, &materials );
#endif
        
        // Get energy and related material parameters from particles
        Interp_P2C ( particles, materials.Cv,   &mesh, mesh.Cv,   mesh.xg_coord, mesh.zg_coord,  0, 0 );
        Interp_P2C ( particles, materials.Qr,   &mesh, mesh.Qr,   mesh.xg_coord, mesh.zg_coord,  0, 0 );
        Interp_P2U ( particles, materials.k_eff,    &mesh, mesh.kz,   mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz,   0, mesh.BCv.type , &model);
        Interp_P2U ( particles, materials.k_eff,    &mesh, mesh.kx,   mesh.xg_coord, mesh.zvx_coord,  mesh.Nx, mesh.Nz+1,   0, mesh.BCu.type , &model);
        Interp_P2C ( particles, particles.T, &mesh, mesh.T, mesh.xg_coord, mesh.zg_coord,  1, 0 );
        
#ifdef _NEW_INPUT_
        SetBCs_new( &mesh, &model, scaling , &particles, &materials );
#else
        SetBCs( &mesh, &model, scaling , &particles, &materials );
#endif
        if ( model.thermal_eq == 1 ) ThermalSteps( &mesh, model,  mesh.T,  mesh.dT,  mesh.rhs_t, mesh.T, &particles, model.cooling_time, scaling );
        if ( model.therm_pert == 1 ) SetThermalPert( &mesh, model, scaling );
        Interp_Grid2P( particles, particles.T,    &mesh, mesh.T, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type );
        
        printf("*************************************\n");
        printf("******** Initialize pressure ********\n");
        printf("*************************************\n");
        
        // Strain rate field
        if (model.cpc==-1) CountPartCell_BEN( &particles, &mesh, model, topo, 0, scaling );
        if (model.cpc== 0) CountPartCell_Old( &particles, &mesh, model, topo, 0, scaling  );
        if (model.cpc== 1) CountPartCell    ( &particles, &mesh, model, topo, topo_ini, 0, scaling  );
        StrainRateComponents( &mesh, scaling, &model );
        
        // Free surface - subgrid density correction
        if ( model.free_surf == 1 ) {
            SurfaceDensityCorrection( &mesh, model, topo, scaling  );
        }
        else {
            ArrayEqualArray( mesh.rho_app_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );
            ArrayEqualArray( mesh.rho_app_s, mesh.rho_s, (mesh.Nx)*(mesh.Nz) );
        }
        
        // Lithostatic pressure for initial visco-plastic viscosity field
        ComputeLithostaticPressure( &mesh, &model, materials.rho[0], scaling, 0 );
        Interp_Grid2P( particles, particles.P,    &mesh, mesh.p_lith, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
        ArrayEqualArray( mesh.p_in, mesh.p_lith,  (mesh.Nx-1)*(mesh.Nz-1) );
        
        printf("*************************************\n");
        printf("******* Initialize grain size *******\n");
        printf("*************************************\n");
        
        // Grain size
        Interp_P2C ( particles, particles.d,  &mesh, mesh.d, mesh.xg_coord,  mesh.zg_coord, 1, 0 );
        ArrayEqualArray( mesh.d0, mesh.d,  (mesh.Nx-1)*(mesh.Nz-1) );
        
        printf("*************************************\n");
        printf("******** Initialize density *********\n");
        printf("*************************************\n");
        
        if ( model.eqn_state > 0 ) {
            UpdateDensity( &mesh, &particles, &materials, &model, &scaling );
        }
        ArrayEqualArray( mesh.rho_app_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );
        ArrayEqualArray( mesh.rho_app_s, mesh.rho_s, (mesh.Nx)*(mesh.Nz) );
        
        printf("*************************************\n");
        printf("******* Initialize viscosity ********\n");
        printf("*************************************\n");
        
        if ( model.iselastic == 1 ) ShearModulusGrid( &mesh, materials, model, scaling );
        
        // Compute cohesion and friction angle on the grid
        CohesionFrictionGrid( &mesh, materials, model, scaling );
        
        Interp_Grid2P( particles, particles.P,    &mesh, mesh.p_in, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
        Interp_Grid2P( particles, particles.T,    &mesh, mesh.T,    mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type );
        // Some stuff to be put on vertices
        Interp_TPdphi_centroid2vertices (&mesh, &model );
        if ( model.rheo_on_cells == 0 ) NonNewtonianViscosityGrid (     &mesh, &materials, &model, Nmodel, &scaling, 1 );
        if ( model.rheo_on_cells == 1 ) NonNewtonianViscosityCells(     &mesh, &materials, &model, Nmodel, &scaling, 1 );
        
        printf("Number of phases : %d\n", model.Nb_phases);
        MinMaxArrayTag( mesh.p_lith,     scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P litho   ", mesh.BCp.type );
        MinMaxArrayTag( mesh.p_in,       scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P initial ", mesh.BCp.type );
        MinMaxArrayTag( mesh.T,          scaling.T,   (mesh.Nx-1)*(mesh.Nz-1), "T         ", mesh.BCp.type );
        MinMaxArrayTag( mesh.d,          scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );
        MinMaxArrayTag( mesh.mu_s,       scaling.S,   (mesh.Nx)*(mesh.Nz),     "mu_s      ", mesh.BCg.type );
        MinMaxArrayTag( mesh.mu_n,       scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "mu_n      ", mesh.BCp.type );
        MinMaxArrayTag( mesh.eta_s,      scaling.eta, (mesh.Nx)*(mesh.Nz),     "eta_s     ", mesh.BCg.type );
        MinMaxArrayTag( mesh.eta_n,      scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_n     ", mesh.BCp.type );
        MinMaxArrayTag( mesh.eta_phys_s, scaling.eta, (mesh.Nx)*(mesh.Nz),     "eta_phys_s", mesh.BCg.type );
        MinMaxArrayTag( mesh.eta_phys_n, scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_phys_n", mesh.BCp.type );
        MinMaxArrayTag( mesh.rho_s,      scaling.rho, (mesh.Nx)*(mesh.Nz),     "rho_s     ", mesh.BCg.type );
        MinMaxArrayTag( mesh.rho_n,      scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
        
        printf("*************************************\n");
        printf("******** Initialize timestep ********\n");
        printf("*************************************\n");
        
        DefineInitialTimestep( &model, &mesh, particles, materials, scaling );
        
        if ( model.rec_T_P_x_z == 1 ) {
            ArrayEqualArray( particles.T0,   particles.T, particles.Nb_part );
            ArrayEqualArray( particles.P0,   particles.P, particles.Nb_part );
            ArrayEqualArray( particles.x0,   particles.x, particles.Nb_part );
            ArrayEqualArray( particles.z0,   particles.z, particles.Nb_part );
            ArrayEqualArray( particles.Tmax, particles.T, particles.Nb_part );
            ArrayEqualArray( particles.Pmax, particles.P, particles.Nb_part );
        }
        
        printf("*************************************\n");
        printf("*** Write initial file or restart ***\n");
        printf("*************************************\n");
        
        // Write initial output
#ifndef _VG_
        if ( writer == 1 ) {
            WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, model, "Output",  materials, scaling );
            if ( model.write_markers == 1 ) WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, model, "Particles",  materials, scaling );
        }
#endif
        // Set initial stresses and pressure to zero
        //        Initialise1DArrayDouble( mesh.p_in,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.sxxd,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.szzd,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.sxz,   (mesh.Nx)  *(mesh.Nz)  , 0.0 );
        // Generate deformation maps
        if ( model.def_maps == 1 ) GenerateDeformationMaps( &mesh, &materials, &model, Nmodel, &scaling, 1 );
    }
    else {
        // Which step do we restart from (BreakpointXXXX.dat)
        model.step = istep;
        printf("Restarting from step number %05d...\n", model.step);
        LoadBreakpointParticles( &particles, &mesh, &topo_chain, &topo_chain_ini, &model, &topo, &topo_ini, scaling  );
        SetGridCoordinates( &mesh, &model, model.Nx, model.Nz ); // Overwrite previous grid
    }
    
    // Evaluate rhs functions
    EvaluateRHS( &mesh, model, scaling, materials.rho[0] );
    
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    //                                             TIME LOOP : en avant les Doud'series !
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    
    model.step += 1;
    
    for (; model.step<=model.Nt; model.step++) {
        
        printf("*****************************************************\n");
        printf("****************** Time step %05d ******************\n", model.step);
        printf("*****************************************************\n");


        //------------------------------------------------------------------------------------------------------------------------------//
       
        t_omp_step = (double)omp_get_wtime();

        // Define new time step
        EvaluateCourantCriterion_BEN( mesh.u_in, mesh.v_in, &model, scaling, &mesh, 0 );

        // Save initial dt
        model.dt0 = model.dt;

        Initialise1DArrayInt( particles.generation,   particles.Nb_part  , 0 );

        // Initial solution fields (Fine mesh)
        if ( (model.ispureshear_ale == 1 || model.isperiodic_x == 1) && model.HomoFields == 1) InitialiseSolutionFields( &mesh, &model );

        //------------------------------------------------------------------------------------------------------------------------------//

        // Remove particles that would be above the surface
        if ( model.free_surf == 1 ) {
            CleanUpSurfaceParticles( &particles, &mesh, topo, scaling );
            CellFlagging( &mesh, model, topo, scaling );
        }

        // Interpolate material properties from particles to nodes
        t_omp = (double)omp_get_wtime();

        // Energy - interpolate thermal parameters and advected energy
        if ( model.isthermal == 1 ) {

            // Get energy and related material parameters from particles
            Interp_P2C ( particles, materials.Cv,   &mesh, mesh.Cv,   mesh.xg_coord, mesh.zg_coord,  0, 0 );
            Interp_P2C ( particles, materials.Qr,   &mesh, mesh.Qr,   mesh.xg_coord, mesh.zg_coord,  0, 0 );

            Interp_P2U ( particles, materials.k_eff,    &mesh, mesh.kz,   mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz,   0, mesh.BCv.type , &model);
            Interp_P2U ( particles, materials.k_eff,    &mesh, mesh.kx,   mesh.xg_coord, mesh.zvx_coord,  mesh.Nx, mesh.Nz+1,   0, mesh.BCu.type , &model);

            // Get T from previous step from particles
            Interp_P2C ( particles, particles.T, &mesh, mesh.T, mesh.xg_coord, mesh.zg_coord,  1, 0 );
        }

        // Get physical properties that are constant throughout each timestep
        if ( model.eqn_state  > 0 ) {
            UpdateDensity( &mesh, &particles, &materials, &model, &scaling );
        }
        else {
            Interp_P2N ( particles, materials.rho,  &mesh, mesh.rho_s, mesh.xg_coord,  mesh.zg_coord, 0, 0, &model );
            Interp_P2C ( particles, materials.rho,  &mesh, mesh.rho_n, mesh.xg_coord,  mesh.zg_coord, 0, 0 );
        }

        Interp_P2C ( particles, materials.alp, &mesh, mesh.alp,      mesh.xg_coord, mesh.zg_coord, 0, 0 );
        Interp_P2C ( particles, materials.bet, &mesh, mesh.bet,      mesh.xg_coord, mesh.zg_coord, 0, 0 );

        // Free surface - subgrid density correction
        if ( model.free_surf == 1 ) {
            SurfaceDensityCorrection( &mesh, model, topo, scaling  );
        }
        else {
            ArrayEqualArray( mesh.rho_app_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );
            ArrayEqualArray( mesh.rho_app_s, mesh.rho_s, (mesh.Nx  )*(mesh.Nz  ) );
        }

        // Lithostatic pressure
        ComputeLithostaticPressure( &mesh, &model, materials.rho[0], scaling, 1 );

        // Elasticity - interpolate advected/rotated stresses
        if  (model.iselastic == 1 ) {

            // Get old stresses from particles
            Interp_P2C ( particles, particles.sxxd, &mesh, mesh.sxxd0,   mesh.xg_coord, mesh.zg_coord, 1, 0 );
            Interp_P2N ( particles, particles.sxxd, &mesh, mesh.sxxd0_s, mesh.xg_coord, mesh.zg_coord, 1, 0, &model );
            Interp_P2C ( particles, particles.szzd, &mesh, mesh.szzd0,   mesh.xg_coord, mesh.zg_coord, 1, 0 );
            Interp_P2N ( particles, particles.szzd, &mesh, mesh.szzd0_s, mesh.xg_coord, mesh.zg_coord, 1, 0, &model );

            Interp_P2N ( particles, particles.sxz,  &mesh, mesh.sxz0,    mesh.xg_coord, mesh.zg_coord, 1, 0, &model );
            Interp_P2C ( particles, particles.sxz,  &mesh, mesh.sxz0_n,  mesh.xg_coord, mesh.zg_coord, 1, 0 );

            // Interpolate elastic energy from previous step
            Interp_P2C ( particles, particles.rhoUe0, &mesh, mesh.rhoUe0,   mesh.xg_coord, mesh.zg_coord, 1, 0 );

            // Interpolate shear modulus
            ShearModulusGrid( &mesh, materials, model, scaling );
        }

        // Interpolate Grain size
        Interp_P2C ( particles,   particles.d, &mesh, mesh.d0,   mesh.xg_coord, mesh.zg_coord, 1, 0 );
        ArrayEqualArray(  mesh.d,  mesh.d0, Ncx*Ncz );

        // Interpolate Melt fraction
        Interp_P2C ( particles, particles.phi, &mesh, mesh.phi,  mesh.xg_coord, mesh.zg_coord, 1, 0 );
        //        Interp_P2G ( particles, particles.phi,   &mesh, mesh.phi,  mesh.xg_coord,   mesh.zg_coord, Ncx, Ncz, xmin_c, zmin_c, 1, 0, &model, mesh.BCp.type );
        
        // Interpolate pressure
        Interp_P2C ( particles, particles.P, &mesh, mesh.p_in,   mesh.xg_coord, mesh.zg_coord, 1, 0 );

        if (model.isPl_soft == 1) {
            Interp_P2C ( particles, particles.strain_pl, &mesh, mesh.strain_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
            Interp_P2N ( particles, particles.strain_pl, &mesh, mesh.strain_s, mesh.xg_coord, mesh.zg_coord, 1, 0, &model );
        }

        // Compute cohesion and friction angle on the grid
        CohesionFrictionGrid( &mesh, materials, model, scaling );
        
        // Detect compressible cells
        if (model.compressible == 1) DetectCompressibleCells ( &mesh, &model );
        
        // Some stuff to be put on vertices
        Interp_TPdphi_centroid2vertices (&mesh, &model );

        // Min/Max interpolated fields
        MinMaxArrayTag( mesh.rho_s,    scaling.rho, (mesh.Nx)*(mesh.Nz),     "rho_s   ", mesh.BCg.type );
        MinMaxArrayTag( mesh.rho_n,    scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n   ", mesh.BCp.type );
        MinMaxArrayTag( mesh.sxz0,     scaling.S,   (mesh.Nx)*(mesh.Nz),     "sxz0    ", mesh.BCg.type );
        MinMaxArrayTag( mesh.sxxd0,    scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "sxx0    ", mesh.BCp.type );
        MinMaxArrayTag( mesh.szzd0,    scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "szz0    ", mesh.BCp.type );
        MinMaxArrayTag( mesh.mu_s,     scaling.S,   (mesh.Nx)*(mesh.Nz),     "mu_s    ", mesh.BCg.type );
        MinMaxArrayTag( mesh.mu_n,     scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "mu_n    ", mesh.BCp.type );
        MinMaxArrayTag( mesh.C_s,      scaling.S,   (mesh.Nx)*(mesh.Nz),     "C_s     ", mesh.BCg.type );
        MinMaxArrayTag( mesh.C_n,      scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "C_n     ", mesh.BCp.type );
        MinMaxArrayTag( mesh.phi_s,    180.0/M_PI,  (mesh.Nx)*(mesh.Nz),     "phi_s   ", mesh.BCg.type );
        MinMaxArrayTag( mesh.phi_n,    180.0/M_PI,  (mesh.Nx-1)*(mesh.Nz-1), "phi_n   ", mesh.BCp.type );
        MinMaxArrayTag( mesh.strain_s,   1.0,       (mesh.Nx)*(mesh.Nz),     "strain_s", mesh.BCg.type );
        MinMaxArrayTag( mesh.strain_n,   1.0,       (mesh.Nx-1)*(mesh.Nz-1), "strain_n", mesh.BCp.type );
        MinMaxArrayTag( mesh.T,      scaling.T,     (mesh.Nx-1)*(mesh.Nz-1), "T       ", mesh.BCt.type );
        MinMaxArray(particles.T, scaling.T, particles.Nb_part, "T part  ");
        printf("** Time for particles interpolations I = %lf sec\n",  (double)((double)omp_get_wtime() - t_omp) );

        // Allocate and initialise solution and RHS vectors
#ifdef _NEW_INPUT_
        SetBCs_new( &mesh, &model, scaling , &particles, &materials );
#else
        SetBCs( &mesh, &model, scaling , &particles, &materials );
#endif
        EvalNumberOfEquations( &mesh, &Stokes );
        if ( model.Newton == 1 ) EvalNumberOfEquations( &mesh, &Jacob  );
        SAlloc( &Stokes,  Stokes.neq );
        if ( model.Newton == 1 ) SAlloc(  &Jacob,   Jacob.neq );
        if ( model.decoupled_solve == 1 ) {
            SAlloc( &StokesA, Stokes.neq_mom );
            SAlloc( &StokesB, Stokes.neq_mom );
            SAlloc( &StokesC, Stokes.neq_cont);
            SAlloc( &StokesD, Stokes.neq_cont );
        }
        if ( model.Newton == 1 ) {
            SAlloc( &JacobA, Stokes.neq_mom );
            SAlloc( &JacobB, Stokes.neq_mom );
            SAlloc( &JacobC, Stokes.neq_cont);
            SAlloc( &JacobD, Stokes.neq_cont );
        }
        printf( "Linear systems allocated\n");
        printf( "neq_tot = %d, neq_mom = %d, neq_cont = %d\n", Stokes.neq, Stokes.neq_mom, Stokes.neq_cont );

        // Set vector x = [u;p]
        InitialiseSolutionVector( &mesh, &Stokes, &model );

        //------------------------------------------------------------------------------------------------------------------------------//

        // Picard iteration cycle
        Nmodel.nit       = 0;
        Nmodel.stagnated = 0;
        nstag            = 0;
        Nmax_picard      = Nmodel.nit_max;
        //        dummy            = model.free_surf_stab;

        ArrayEqualArray( mesh.p_start,    mesh.p_in,      (mesh.Nx-1)*(mesh.Nz-1) );
        ArrayEqualArray( mesh.u_start,    mesh.u_in,      (mesh.Nx)  *(mesh.Nz+1) );
        ArrayEqualArray( mesh.v_start,    mesh.v_in,      (mesh.Nx+1)*(mesh.Nz)   );
        ArrayEqualArray( mesh.rho_app_n0, mesh.rho_app_n, (mesh.Nx-1)*(mesh.Nz-1) );
        ArrayEqualArray( mesh.rho_app_s0, mesh.rho_app_s, (mesh.Nx)  *(mesh.Nz)   );

        // Set up solver context
        if ( model.decoupled_solve == 1 ) {
            cholmod_start( &CholmodSolver.c );
            if ( Nmodel.nit==0 ) CholmodSolver.Analyze = 1;
            printf("Run CHOLMOD analysis yes/no: %d \n", CholmodSolver.Analyze);
        }
        
        int Nmax_picard = Nmodel.nit_max;
        
        while ( Nmodel.nit <= Nmax_picard ) {
            
            printf("******************************************\n");
            printf("*** Picard it. %02d of %02d (step = %05d) ***\n", Nmodel.nit, Nmodel.nit_max, model.step);
            printf("******************************************\n");
            
            // Non-Linearity
            UpdateNonLinearity( &mesh, &particles, &topo_chain, &topo, materials, &model, &Nmodel, scaling, 0, 0.0 );
            
//            MinMaxArrayTag( mesh.exxd, scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "exx", mesh.BCp.type );
//            MinMaxArrayTag( mesh.exz_n, scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "exz_n", mesh.BCp.type );
            MinMaxArrayTag( mesh.exz, scaling.E, (mesh.Nx-0)*(mesh.Nz-0), "exz", mesh.BCg.type );
//
//            MinMaxArrayTag( mesh.exxd_s, scaling.E, (mesh.Nx-0)*(mesh.Nz-0), "exx_s", mesh.BCg.type );
//
            MinMaxArrayTag( mesh.eta_phys_n, scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_phys_n", mesh.BCp.type );
            MinMaxArrayTag( mesh.eta_phys_s, scaling.eta, (mesh.Nx-0)*(mesh.Nz-0), "eta_phys_s", mesh.BCg.type );
            MinMaxArrayTag( mesh.eta_n, scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),   "eta_eff_n", mesh.BCp.type );
            MinMaxArrayTag( mesh.eta_s, scaling.eta, (mesh.Nx)*(mesh.Nz),       "eta_eff_s", mesh.BCg.type );
            
            // Stokes solver
            if ( model.ismechanical == 1 ) {

                // Build discrete system of equations - MATRIX
                if ( model.decoupled_solve == 0 ) BuildStokesOperator           ( &mesh, model, 0, mesh.p_in, mesh.u_in, mesh.v_in, &Stokes, 1 );
                if ( model.decoupled_solve == 1 ) BuildStokesOperatorDecoupled  ( &mesh, model, 0, mesh.p_in, mesh.u_in, mesh.v_in, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, 1 );
                if ( model.Newton          == 1 ) ComputeViscosityDerivatives_FD( &mesh, &materials, &model, Nmodel, &scaling, 1 );
                if ( model.Newton          == 1 ) RheologicalOperators( &mesh, &model, &scaling, 1 );
                if ( model.Newton          == 1 ) BuildJacobianOperatorDecoupled( &mesh, model, 0, mesh.p_in, mesh.u_in, mesh.v_in,  &Jacob,  &JacobA,  &JacobB,  &JacobC,   &JacobD, 1 );
                
                // Diagonal scaling
                if ( model.diag_scaling ) {
                    if ( model.Newton          == 0 )  ExtractDiagonalScale( &StokesA, &StokesB, &StokesC, &StokesD );
                    if ( model.Newton          == 1 )  ExtractDiagonalScale( &JacobA,  &JacobB,  &JacobC,   &JacobD );
                    if ( model.Newton          == 1 )  ArrayEqualArray(StokesA.d, JacobA.d, StokesA.neq);
                    if ( model.Newton          == 1 )  ArrayEqualArray(StokesC.d, JacobC.d, StokesC.neq);
                    ScaleMatrix( &JacobA,  &JacobB,  &JacobC,  &JacobD  );
                    ScaleMatrix( &StokesA, &StokesB, &StokesC, &StokesD );
                }
                
                // If iteration > 0 ---> Evaluate non linear residual and test
                if ( model.DefectCorrectionForm == 1 || (Nmodel.nit>0 && model.DefectCorrectionForm == 0 ) ) {

                    printf("---- Non-linear residual ----\n");
                    RheologicalOperators( &mesh, &model, &scaling, 0 );
                    if ( model.decoupled_solve == 0 ) EvaluateStokesResidual( &Stokes, &Nmodel, &mesh, model, scaling, 0 );
                    if ( model.decoupled_solve == 1 ) EvaluateStokesResidualDecoupled( &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &Nmodel, &mesh, model, scaling, 0 );
                    if ( model.Newton          == 1 ) {
                        ArrayEqualArray( JacobA.F, StokesA.F,  StokesA.neq );
                        ArrayEqualArray( JacobC.F, StokesC.F,  StokesC.neq );
                    }
                    Nmodel.resx_f = Nmodel.resx;
                    Nmodel.resz_f = Nmodel.resz;
                    Nmodel.resp_f = Nmodel.resp;
#ifndef _VG_
                    if ( model.write_debug == 1 ) WriteResiduals( mesh, model, Nmodel, scaling );
#endif
                    printf("---- Non-linear residual ----\n");

                    // if pass --> clear matrix break
                    if ( (Nmodel.resx < Nmodel.tol_u) && (Nmodel.resz < Nmodel.tol_u) && (Nmodel.resp < Nmodel.tol_p) ) {

                        printf( "Non-linear solver converged to tol_u = %2.2e tol_p = %2.2e\n", Nmodel.tol_u, Nmodel.tol_p );
                        if ( model.decoupled_solve == 0 ) { FreeMat( &Stokes ); }
                        if ( model.decoupled_solve == 1 ) {
                            FreeMat( &StokesA );
                            FreeMat( &StokesB );
                            FreeMat( &StokesC );
                            FreeMat( &StokesD );
                        }
                        if ( model.Newton == 1 ) {
                            FreeMat( &JacobA );
                            FreeMat( &JacobB );
                            FreeMat( &JacobC );
                            FreeMat( &JacobD );
                        }
                        break;
                    }
                }


                // Direct solve
                t_omp = (double)omp_get_wtime();
                if ( Nmodel.nit == 0 && model.DefectCorrectionForm == 0  ) {
                    if ( model.decoupled_solve == 0 ) SolveStokes( &Stokes, &CholmodSolver );
                    if ( model.decoupled_solve == 1 ) SolveStokesDecoupled( &StokesA, &StokesB, &StokesC, &StokesD, &Stokes, &CholmodSolver, model, &mesh, scaling );
                }
                else {
                    if ( model.decoupled_solve == 0 ) SolveStokesDefect( &Stokes, &CholmodSolver, &Nmodel, &mesh, &model, &particles, &topo_chain, &topo, materials, scaling );
                    if ( model.decoupled_solve == 1 ) {
                        if ( model.Newton==0 ) SolveStokesDefectDecoupled( &StokesA, &StokesB, &StokesC, &StokesD, &Stokes, &CholmodSolver, &Nmodel, &mesh, &model, &particles, &topo_chain, &topo, materials, scaling, &StokesA, &StokesB, &StokesC );
                        if ( model.Newton==1 ) SolveStokesDefectDecoupled( &StokesA, &StokesB, &StokesC, &StokesD, &Stokes, &CholmodSolver, &Nmodel, &mesh, &model, &particles, &topo_chain, &topo, materials, scaling, &JacobA, &JacobB, &JacobC );
                        
                    }

                    if ( Nmodel.stagnated == 1 && safe == 0 ) {
                        printf( "Non-linear solver stagnated to res_u = %2.2e res_z = %2.2e\n", Nmodel.resx_f, Nmodel.resz_f );

                        if ( model.decoupled_solve == 0 ) { FreeMat( &Stokes ); }
                        if ( model.decoupled_solve == 1 ) {
                            FreeMat( &StokesA );
                            FreeMat( &StokesB );
                            FreeMat( &StokesC );
                            FreeMat( &StokesD );
                        }
                        if ( model.Newton == 1 ) {
                            FreeMat( &JacobA );
                            FreeMat( &JacobB );
                            FreeMat( &JacobC );
                            FreeMat( &JacobD );
                        }
                        break;
                    }
                    if ( Nmodel.stagnated == 1 && model.iselastic == 1 && safe == 1 ) {
                        printf( "\e[1;31mWARNING : Non-linear solver stagnated (tol_u = %2.2e tol_p = %2.2e)\e[m\n", Nmodel.tol_u, Nmodel.tol_p );
                        printf( "\e[1;31mReducing the timestep, and restart the iterations cycle...\e[m\n");
                        Nmodel.nit = -1;
                        printf( "Restart solutions\n");
                        ArrayEqualArray( mesh.p_in, mesh.p_start,  (mesh.Nx-1)*(mesh.Nz-1) );
                        ArrayEqualArray( mesh.u_in, mesh.u_start,  (mesh.Nx)  *(mesh.Nz+1) );
                        ArrayEqualArray( mesh.v_in, mesh.v_start,  (mesh.Nx+1)*(mesh.Nz)   );
                        nstag++;
                    }
                }

                if ( Nmodel.stagnated == 0 ) {
                    printf("EXTRACT SOLS\n");
                    Initialise2DArray( mesh.u_in, mesh.Nx, (mesh.Nz+1), 0.0 );
                    Initialise2DArray( mesh.v_in, (mesh.Nx+1), mesh.Nz, 0.0 );
                    ExtractSolutions( &Stokes, &mesh, &model );
                    //                    }
                    if ( Nmodel.nit == 0  ) {
                        printf("---- Direct solve residual ----\n");
                        RheologicalOperators( &mesh, &model, &scaling, 0 );
                        if ( model.decoupled_solve == 0 ) EvaluateStokesResidual( &Stokes, &Nmodel, &mesh, model, scaling, 0 );
                        if ( model.decoupled_solve == 1 ) EvaluateStokesResidualDecoupled( &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &Nmodel, &mesh, model, scaling, 0 );
                        printf("---- Direct solve residual ----\n");
                        
                        
                        
#ifndef _VG_
                        if ( model.write_debug == 1 ) WriteResiduals( mesh, model, Nmodel, scaling );
#endif
                        
                        //                        if ( Nmodel.resx>1e-8/(scaling.F/pow(scaling.L,3)) || Nmodel.resz>1e-8/(scaling.F/pow(scaling.L,3)) || Nmodel.resp>1e-10/scaling.E ) {
                        //                            printf("Direct solve residuals are too large !\n Stopping now !");
                        //                            exit(22);
                        //                        }
                    }
                }
                if ( model.decoupled_solve == 0 ) { FreeMat( &Stokes ); }
                if ( model.decoupled_solve == 1 ) {
                    FreeMat( &StokesA );
                    FreeMat( &StokesB );
                    FreeMat( &StokesC );
                    FreeMat( &StokesD );
                }
                if ( model.Newton == 1 ) {
                    FreeMat( &JacobA );
                    FreeMat( &JacobB );
                    FreeMat( &JacobC );
                    FreeMat( &JacobD );
                }
            } // end is_mechanical
            
            Nmodel.nit++;
        }
        
        // Clean solver context
        if ( model.decoupled_solve == 1 && Nmodel.nit>0 ) {
            printf("Cleaning up Cholesky factors --- nit = %02d\n",  Nmodel.nit);
            cholmod_free_factor ( &CholmodSolver.Lfact, &CholmodSolver.c);
            cholmod_finish( &CholmodSolver.c );
        }
        
        // Update rheology
        UpdateNonLinearity( &mesh, &particles, &topo_chain, &topo, materials, &model, &Nmodel, scaling, 0, 1.0 );
        
        // Min/Max velocities
        MinMaxArray( mesh.u_in,  scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
        MinMaxArray( mesh.v_in,  scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
        MinMaxArray( mesh.p_in,  scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );
        MinMaxArray( mesh.div_u, scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "  div(V)" );
        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        // Particle velocities
        VelocitiesToParticles( &mesh, &particles, particles.Vx, particles.Vz, model, scaling );
        
        // Free surface - interpolate velocity components on the free surface
        if ( model.free_surf == 1 ) {
            SurfaceVelocity( &mesh, model, &topo, &topo_chain, scaling );
            VelocitiesToParticles( &mesh, &topo_chain_ini, topo_chain_ini.Vx, topo_chain_ini.Vz, model, scaling );
            
            MinMaxArray( topo_chain.Vx,  scaling.V, topo_chain.Nb_part,   "Vx surf." );
            MinMaxArray( topo_chain.Vz,  scaling.V, topo_chain.Nb_part,   "Vz surf." );
            MinMaxArray( topo_chain_ini.Vx,  scaling.V, topo_chain_ini.Nb_part,   "Vx surf. ini." );
            MinMaxArray( topo_chain_ini.Vz,  scaling.V, topo_chain_ini.Nb_part,   "Vz surf. ini." );
        }
        
        // Update stresses on markers
        if (model.iselastic == 1 ) {
            // Compute stress changes on the grid and update stress on the particles
            UpdateParticleStress(  &mesh, &particles, &model, &materials, &scaling );
        }
        else {
            Interp_Grid2P( particles, particles.sxxd, &mesh, mesh.sxxd, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
            Interp_Grid2P( particles, particles.szzd, &mesh, mesh.szzd, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
            Interp_Grid2P( particles, particles.sxz,  &mesh, mesh.sxz , mesh.xg_coord,  mesh.zg_coord,  mesh.Nx  , mesh.Nz, mesh.BCg.type   );
        }
        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        if (model.isthermal == 1 ) {
            
            printf("*************************************\n");
            printf("*********** Thermal solver **********\n");
            printf("*************************************\n");
            
            t_omp = (double)omp_get_wtime();
            
            // Matrix assembly and direct solve
            EnergyDirectSolve( &mesh, model,  mesh.T,  mesh.dT,  mesh.rhs_t, mesh.T, &particles, model.dt, model.shear_heat, model.adiab_heat, scaling, 1 );
            
            // Update energy on particles
            UpdateParticleEnergy( &mesh, scaling, model, &particles, &materials );
            
            if ( model.iselastic == 1 ) Interp_Grid2P( particles, particles.rhoUe0, &mesh, mesh.rhoUe0, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
            
            // Calculate energies
            if (model.write_energies==1) Energies( &mesh, model, scaling );
            
            printf("** Time for Thermal solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
        }
        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        // Grain size evolution
        UpdateParticleGrainSize( &mesh, scaling, model, &particles, &materials );
        MinMaxArrayTag( mesh.d0    , scaling.L, (mesh.Nx-1)*(mesh.Nz-1), "d0", mesh.BCp.type );
        MinMaxArrayTag( mesh.d     , scaling.L, (mesh.Nx-1)*(mesh.Nz-1), "d ", mesh.BCp.type );
        MinMaxArrayPart( particles.d, scaling.L, particles.Nb_part, "d on markers", particles.phase ) ;
        
        // Update pressure on the particles
        Interp_Grid2P( particles, particles.P, &mesh, mesh.p_in, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
        //UpdateParticlePressure( &mesh, scaling, model, &particles, &materials );
        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        // Update maximum pressure and temperature on markers
        if ( model.rec_T_P_x_z == 1 )  UpdateMaxPT( scaling, model, &particles );
        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        if ( model.advection == 1 ) {

            printf("*************************************\n");
            printf("************** Advection ************\n");
            printf("*************************************\n");

            t_omp = (double)omp_get_wtime();
            EvaluateCourantCriterion_BEN( mesh.u_in, mesh.v_in, &model, scaling, &mesh, 0 );

//            double whole_dt = model.dt;
//            int nsub, isub;
//
//            if ( model.ispureshear_ale > 0 ) nsub = 1;
//            else nsub = 1;
//            model.dt = whole_dt/nsub;
//
//            // Loop on substeps
//            for (isub=0;isub<nsub;isub++) {
//
//                printf("Advection step %03d of %03d: dtsub = %2.2e\n", isub, nsub, model.dt*scaling.t );

                // Advect domain boundaries
                if ( model.ispureshear_ale > 0 ) {
                    PureShearALE( &model, &mesh, &topo_chain, scaling );
                }

                // Advect free surface
                if ( model.free_surf == 1 ) {
//                    RogerGuntherII( &topo_chain,     model, mesh, 1, scaling );
//                    RogerGuntherII( &topo_chain_ini, model, mesh, 1, scaling );
                    AdvectFreeSurf( &topo_chain,     model, scaling );
                    AdvectFreeSurf( &topo_chain_ini, model, scaling );
                    MinMaxArray( topo_chain.z,      scaling.L, topo_chain.Nb_part,       "z surf." );
                    MinMaxArray( topo_chain_ini.z,  scaling.L, topo_chain_ini.Nb_part,   "z surf. ini." );
                }

                // Correction for particle inflow 0
                if (model.ispureshear_ale == -1 && model.isperiodic_x == 0) ParticleInflowCheck( &particles, &mesh, model, topo, 0 );

                // Advect fluid particles
                RogerGuntherII( &particles, model, mesh, 1, scaling );

                // Correction for particle inflow 1
                if (model.ispureshear_ale == -1 && model.isperiodic_x == 0) ParticleInflowCheck( &particles, &mesh, model, topo, 1 );

                // Roration of stresses (visco-elastic flow)
                if ( model.iselastic == 1 ) RotateStresses( mesh, &particles, model, &scaling );

                // Update accumulated strain
//                AccumulatedStrain( &mesh, scaling, model, &particles );
                AccumulatedStrainII( &mesh, scaling, model, &particles,  mesh.xc_coord,  mesh.zc_coord, mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );

                // Update deformation gradient tensor components
                if ( model.fstrain == 1 ) DeformationGradient( mesh, scaling, model, &particles );

                if ( model.free_surf == 1 ) {

                    if ( model.surf_remesh == 1 ) {
                        // Get current topography
                        ProjectTopography( &topo,     &topo_chain,     model, mesh, scaling, mesh.xg_coord, 0 );
                        ProjectTopography( &topo_ini, &topo_chain_ini, model, mesh, scaling, mesh.xg_coord, 0 );
                        MarkerChainPolyFit( &topo,     &topo_chain,     model, mesh );
                        MarkerChainPolyFit( &topo_ini, &topo_chain_ini, model, mesh );
                        // Remesh free surface I
                        RemeshMarkerChain( &topo_chain,     &topo,     model, scaling, &mesh, 1 );
                        RemeshMarkerChain( &topo_chain_ini, &topo_ini, model, scaling, &mesh, 1 );
                    }

                    // Project topography on vertices
                    ProjectTopography( &topo,     &topo_chain,     model, mesh, scaling, mesh.xg_coord, 0 );
                    ProjectTopography( &topo_ini, &topo_chain_ini, model, mesh, scaling, mesh.xg_coord, 0 );
                    ArrayEqualArray( topo.height0, topo.height, mesh.Nx );
                    ArrayEqualArray( topo_ini.height0, topo_ini.height, mesh.Nx );

                    //                    // Diffuse topography
                    if ( model.surf_processes >= 1 )  DiffuseAlongTopography( &mesh, model, scaling, topo.height, mesh.Nx, 0.0, model.dt );
                    //
                    //                    // Marker chain polynomial fit
                    MarkerChainPolyFit( &topo,     &topo_chain,     model, mesh );
                    CorrectTopoIni( &particles, materials, &topo_chain_ini, &topo, model, scaling, &mesh);
                    MarkerChainPolyFit( &topo_ini, &topo_chain_ini, model, mesh );

                    // Sedimentation
                    if ( model.surf_processes >= 1 ) {
                        AddPartSed( &particles, materials, &topo_chain, &topo, model, scaling, &mesh);
                        if (model.cpc==-1) CountPartCell_BEN( &particles, &mesh, model, topo, 0, scaling );
                        if (model.cpc== 0) CountPartCell_Old( &particles, &mesh, model, topo, 0, scaling );
                        if (model.cpc== 1) CountPartCell    ( &particles, &mesh, model, topo, topo_ini, 0, scaling );
                    }

#ifdef _HDF5_
                    if ( model.write_debug == 1 ) {
                        WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, Mmodel, model, "Outputx", materials, scaling );
                        WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, Mmodel, model, "Particlesx", materials, scaling );
                    }
#endif
                    MinMaxArray( topo_chain.z,      scaling.L, topo_chain.Nb_part,       "z surf." );
                    MinMaxArray( topo_chain_ini.z,  scaling.L, topo_chain_ini.Nb_part,   "z surf. ini." );

                    // Remesh free surface II
                    RemeshMarkerChain( &topo_chain,     &topo,     model, scaling, &mesh, 2 );
                    RemeshMarkerChain( &topo_chain_ini, &topo_ini, model, scaling, &mesh, 2 );
                    CorrectTopoIni( &particles, materials, &topo_chain_ini, &topo, model, scaling, &mesh);
                    MarkerChainPolyFit( &topo_ini, &topo_chain_ini, model, mesh );


                    MinMaxArray( topo_chain.z,      scaling.L, topo_chain.Nb_part,       "z surf." );
                    MinMaxArray( topo_chain_ini.z,  scaling.L, topo_chain_ini.Nb_part,   "z surf. ini." );

                    // Remove particles that are above the surface
                    CleanUpSurfaceParticles( &particles, &mesh, topo, scaling );

                    // Call cell flagging routine for free surface calculations
                    CellFlagging( &mesh, model, topo, scaling );
                }

                printf("** Time for advection solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp) );

#ifdef _HDF5_
                if ( model.write_debug == 1 ) {WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, Mmodel, model, "Outputxx", materials, scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, Mmodel, model, "Particlesxx", materials, scaling );
                }
#endif

            // Count the number of particle per cell
            printf("Before re-seeding : number of particles = %d\n", particles.Nb_part);
            t_omp = (double)omp_get_wtime();
            
            // Count the number of particle per cell
            t_omp = (double)omp_get_wtime();
            if (model.cpc==-1) CountPartCell_BEN( &particles, &mesh, model, topo, 1, scaling );
            if (model.cpc==-1) CountPartCell_BEN( &particles, &mesh, model, topo, 0, scaling );
            if (model.cpc== 0) CountPartCell_Old( &particles, &mesh, model, topo, 1, scaling );
            if (model.cpc== 0) CountPartCell_Old( &particles, &mesh, model, topo, 0, scaling );
            if (model.cpc== 1) CountPartCell    ( &particles, &mesh, model, topo, topo_ini, 1, scaling );
            if (model.cpc== 1) CountPartCell    ( &particles, &mesh, model, topo, topo_ini, 0, scaling );
            printf("** Time for CountPartCell = %lf sec\n", (double)((double)omp_get_wtime() - t_omp) );
            
            printf("After re-seeding : number of particles = %d\n", particles.Nb_part);
            printf("** Time for CountPartCell = %lf sec\n", (double)((double)omp_get_wtime() - t_omp) );
            
            // Remove particles that would be above the surface
            if ( model.free_surf == 1 ) {
                CleanUpSurfaceParticles( &particles, &mesh, topo, scaling );
                CellFlagging( &mesh, model, topo, scaling );
            }


////            }
//            model.dt = whole_dt;
        }
        
        // Update time
        model.time += model.dt;
        
        // Free solution arrays
        SFree( &Stokes );
        if ( model.Newton == 1 ) SFree( &Jacob );
        if ( model.decoupled_solve == 1 ) {
            SFree( &StokesA );
            SFree( &StokesB );
            SFree( &StokesC );
            SFree( &StokesD );
        }
        if ( model.Newton == 1 ) {
            SFree( &JacobA );
            SFree( &JacobB );
            SFree( &JacobC );
            SFree( &JacobD );
        }
        DoodzFree( Stokes.eqn_u );
        DoodzFree( Stokes.eqn_v );
        DoodzFree( Stokes.eqn_p );
        if ( model.Newton == 1 ) {
            DoodzFree( Jacob.eqn_u );
            DoodzFree( Jacob.eqn_v );
            DoodzFree( Jacob.eqn_p );
        }

        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        // Write output data
        if ( writer == 1 && model.step % writer_step == 0 ) {
            
            printf("*************************************\n");
            printf("********* Write output files ********\n");
            printf("*************************************\n");
            
            // Breakpoint file
            t_omp = (double)omp_get_wtime();
            if ( model.rm_break == 1 ) DeletePreviousBreakpoint( model.step, writer_step  );
            MakeBreakpointParticles( &particles, &mesh, &topo_chain, &topo_chain_ini, model, &topo, &topo_ini, scaling );
            UpdateInputFile( fin_name, model.step);
            printf("** Time for Breakpoint file write = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
            
            // Visualisation file
#ifndef _VG_
            t_omp = (double)omp_get_wtime();
            WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, model, "Output", materials, scaling );
            if ( model.write_markers == 1 ) WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, model, "Particles", materials, scaling );
            printf("** Time for Output file write = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
#endif
        }
        
        printf("** Total timestep calculation time = %lf sec\n", (double)((double)omp_get_wtime() - t_omp_step) );
        printf("** Model time = %2.2e sec\n", model.time*scaling.t );
        printf("** Current dt = %2.2e sec, Old dt = %2.2e sec\n", model.dt*scaling.t, model.dt0*scaling.t );

    }
    
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    //                                           END TIME LOOP : c'est fini les Doud'series...
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    
//    // Free markers chains
    if ( model.free_surf == 1 )    FreeMarkerChain( &topo,     &topo_chain     );
    if ( model.free_surf == 1 )    FreeMarkerChain( &topo_ini, &topo_chain_ini );

    // Free Phase diagrams if activated
    if ( model.isPD == 1 ) FreePhaseDiagrams( &model );
//
    // Free particle arrays
    PartFree( &particles, &model );

    // Free arrays
    GridFree( &mesh, &model );

    // Free char*'s
    free(model.input_file);
    free(fin_name);
#ifdef _NEW_INPUT_
    free(PartFileName);
#endif

    printf("\n********************************************************\n");
    printf("************* Ending MDOODZ 6.0 simulation *************\n");
    printf("********************************************************\n");

    // Exit success signal
    return(EXIT_SUCCESS);
}
