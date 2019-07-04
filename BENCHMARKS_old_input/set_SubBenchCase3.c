// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2018  MDOODZ Developper team
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
#include "header_MDOODZ.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {
    
    int k;
    double TopoLevel = 0.0e3/scaling.L; // sets zero initial topography
    
    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        topo_chain->z[k]     = TopoLevel;
        topo_chain->phase[k] = 0;
    }
    
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {

	int np;
    double Lx = (double) (model.xmax - model.xmin);
    double Lz = (double) (model.zmax - model.zmin);
    double Tbg  = 773.15/scaling.T;                         // reference temperature
    double gsbg = 2e-3/scaling.L;                           // reference grain size
    double H    = 1.2e-2/scaling.L;                         // plate thickness

    // set polygon vertices coordinates
    double ax0 = 20.0e-2/scaling.L,        ay0 = 0e-2/scaling.L;
    double ax1 = ax0 - 4.974e-2/scaling.L, ay1 = ay0 - 3.3552e-2/scaling.L;
    double ax2 = ax1 + 0.0049/scaling.L,   ay2 = ay1 - 0.0110/scaling.L;
    double ax3 = ax0 + 0.0049/scaling.L,   ay3 = ay0 - 0.0110/scaling.L;
    
    // Upper slab limit
    double aa0 = (ay0-ay1)/(ax0-ax1);
    double bb0 = ay1  - (ax1*aa0);
    
    // Lower slab limit
    double aa2 = (ay3-ay2)/(ax3-ax2);
    double bb2 = ay2  - (ax2*aa2);
    
    // Left slab limit
    double aa1 = (ay2-ay1)/(ax2-ax1);
    double bb1 = ay2  - (ax2*aa1);
    
    // Right slab limit
    double aa3 = (ay3-ay0)/(ax3-ax0);
    double bb3 = ay0  - (ax0*aa3);

    // Loop over particles
    for( np=0; np<particles->Nb_part; np++ ) {
        
        particles->Vx[np]    = -1.0*particles->x[np]*model.EpsBG;
        particles->Vz[np]    =  particles->z[np]*model.EpsBG;
        particles->d[np]     = gsbg;                          // same grain size everywhere
        particles->phi[np]   = 0.0;                           // zero porosity everywhere
        particles->X[np]     = 0.0;                           // X set to 0
        particles->T[np]     = Tbg;                           // same temperature everywhere

        particles->phase[np] = 1;
        
        // Draw dipping slap
        if ( particles->z[np] < aa0*particles->x[np]+bb0 && particles->z[np] > aa2*particles->x[np]+bb2 && particles->z[np] > aa1*particles->x[np]+bb1 && particles->z[np] < aa3*particles->x[np]+bb3 ) {
            particles->phase[np] = 0;
        }

        // Lower plate
        if (particles->x[np] > 20.0e-1 && particles->z[np] > -H) {
            particles->phase[np] = 0;
        }
        
        }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials ) {
    
    int   k, l, c, c1;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ;
    
    int StressBC_W=0, StressBC_E=0;
    if ( model->user0 == 1 ) StressBC_E = 1;
    double TN = 273.15/scaling.T, TS = 1330/scaling.T;
    double TW = 273.15/scaling.T, TE = 1330/scaling.T;
    
    NX  = mesh->Nx;
    NZ  = mesh->Nz;
    NCX = NX-1;
    NCZ = NZ-1;
    
    X  = malloc (NX*sizeof(double));
    Z  = malloc (NZ*sizeof(double));
    XC = malloc (NCX*sizeof(double));
    ZC = malloc (NCZ*sizeof(double));
    
    for (k=0; k<NX; k++) {
        X[k] = mesh->xg_coord[k];
    }
    for (k=0; k<NCX; k++) {
        XC[k] = mesh->xc_coord[k];
    }
    for (l=0; l<NZ; l++) {
        Z[l] = mesh->zg_coord[l];
    }
    for (l=0; l<NCZ; l++) {
        ZC[l] = mesh->zc_coord[l];
    }
    
    /* --------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for Vx on all grid levels                                                                   */
    /* Type  0: Dirichlet point that matches the physical boundary (Vx: left/right, Vz: bottom/top)            */
    /* Type 11: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)       */
    /* Type  2: Neumann point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)         */
    /* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)              */
    /* Type -2: periodic in the x direction (matches the physical boundary)                                    */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
    /* --------------------------------------------------------------------------------------------------------*/
    
    //--------------------------------------------------------------------------------
    
    for (l=0; l<mesh->Nz+1; l++) {
        for (k=0; k<mesh->Nx; k++) {
            
            c = k + l*(mesh->Nx);
            
            if ( mesh->BCu.type[c] != 30 ) {
                
                // Internal points:  -1
                mesh->BCu.type[c] = -1;
                mesh->BCu.val[c]  =  0;
                
                // Matching BC nodes WEST
                if (k==0 && l>0 && l<mesh->Nz ) { //
                    
                    if ( StressBC_W==1 ) {
                        mesh->BCu.type[c] = 2;
                        c1 = k + (l-1)*(mesh->Nx-1) ;
                        mesh->BCu.val[c]  = -mesh->p_lith[c1];
                    }
                    else {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = 0.0;
                    }
                }
                
                // Matching BC nodes EAST
                if (k==mesh->Nx-1 && l>0 && l<mesh->Nz) { //
                    if ( StressBC_E==1 ) {
                        mesh->BCu.type[c] = 2;
                        c1 = k + (l-1)*(mesh->Nx-1) -1;
                        mesh->BCu.val[c]  = -mesh->p_lith[c1];
                    }
                    else {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = 0.0;
                    }
                }
                
                
                // Free slip
                if (l==0  ) {
                    mesh->BCu.type[c] = 13;
                    mesh->BCu.val[c]  =  0;
                }
                
                if ( l==mesh->Nz ) {
                    mesh->BCu.type[c] = 13;
                    mesh->BCu.val[c]  =  0;
                }
            }
            
        }
    }
    
    
    /* --------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for Vz on all grid levels                                                                   */
    /* Type  0: Dirichlet point that matches the physical boundary (Vx: left/right, Vz: bottom/top)            */
    /* Type 11: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)       */
    /* Type  2: Neumann point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)         */
    /* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)              */
    /* Type -2: periodic in the x direction (does not match the physical boundary)                             */
    /* Type-10: useless point (set to zero)                                                                    */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
    /* --------------------------------------------------------------------------------------------------------*/

    for (l=0; l<mesh->Nz; l++) {
        for (k=0; k<mesh->Nx+1; k++) {
            
            c  = k + l*(mesh->Nx+1);
            
            if ( mesh->BCv.type[c] != 30 ) {

                // Internal points:  -1
                mesh->BCv.type[c] = -1;
                mesh->BCv.val[c]  =  0;

                // Matching BC nodes SOUTH
                if (l==0  && k>0 && k<mesh->Nx) {
                    mesh->BCv.type[c] = 0;
                    mesh->BCv.val[c]  = 0.0;
                }

                // Matching BC nodes NORTH
                if (l==mesh->Nz-1  && k>0 && k<mesh->Nx) {
                    mesh->BCv.type[c] = 0;
                    mesh->BCv.val[c]  = 0.0;
                }

                // Non-matching boundary points
                if ( (k==0) ) {
                    if ( StressBC_W==1 ) {
                        mesh->BCv.type[c] =   11;
                        mesh->BCv.val[c]  =   0;
                    }
                    else {
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0;
                    }
                }

                // Non-matching boundary points
                if ( (k==mesh->Nx) ) {
                    if ( StressBC_E==1 ) {
                        mesh->BCv.type[c] =   11;
                        mesh->BCv.val[c]  =   0;
                    }
                    else{
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0;
                    }
                }

                //                    // Normal stress applied to the free surface
                //                    if (l<mesh->Nz-1 && k>0 && k<mesh->Nx) {
                //                        c1 = k + (l)*(mesh->Nx-1)-1;
                //                        if (mesh->BCp.type[0][c1]==31) mesh->BCv.val[c]  =  -0*mesh->p_lith[c1-(mesh->Nx-1)];
                //                    }
            }
        }
    }

    
    /* --------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for P on all grid levels                                                                    */
    /* Type  0: Dirichlet within the grid                                                                      */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
    /* Type 31: surface pressure (Dirichlet)                                                                   */
    /* --------------------------------------------------------------------------------------------------------*/
    
    for (l=0; l<NCZ; l++) {
        for (k=0; k<NCX; k++) {
            
            c  = k + l*(NCX);
            
            if (mesh->BCt.type[c] != 30) {
                
                // Internal points:  -1
                mesh->BCp.type[c] = -1;
                mesh->BCp.val[c]  =  0;
            }
        }
    }
    
    /* -------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for T on all grid levels                                                                   */
    /* Type  1: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)      */
    /* Type  0: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)             */
    /* Type -2: periodic in the x direction (matches the physical boundary)                                   */
    /* Type -1: not a BC point (tag for inner points)                                                         */
    /* Type 30: not calculated (part of the "air")                                                            */
    /* -------------------------------------------------------------------------------------------------------*/
    
    for (l=0; l<mesh->Nz-1; l++) {
        for (k=0; k<mesh->Nx-1; k++) {
            
            c = k + l*(NCX);
            
            if ( mesh->BCt.type[c] != 30 ) {
                
                // WEST
                if ( k==0 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.typW[l] = 0;
                    mesh->BCt.valW[l] = TW;
                }
                
                // EAST
                if ( k==NCX-1 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.typE[l] = 0;
                    mesh->BCt.valE[l] = TE;
                }
                
                // SOUTH
                if ( l==0 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.typS[k] = 1;
                    mesh->BCt.valS[k] = TS;
                }
                
                // NORTH
                if ( l==NCZ-1 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.typN[k] = 1;
                    mesh->BCt.valN[k] = TN;
                }
                
                // FREE SURFACE
                else {
                    if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c+NCX] == 30) {
                        mesh->BCt.type[c] = 1;
                        mesh->BCt.val[c]  = TN;
                    }
                }
                
                
            }
            
        }
    }
    
    free(X);
    free(Z);
    free(XC);
    free(ZC);
    printf("Velocity and pressure were initialised\n");
    printf("Boundary conditions were set up\n");
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
