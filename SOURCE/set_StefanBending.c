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

// SIMPLE
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {
    
    int k;
    double TopoLevel = 2.e-1/scaling.L; // sets zero initial topography
    double a = 100e3/scaling.L, b = 10e3/scaling.L, Rad=6370e3/scaling.L;
    double maxAngle = 18.0/2.0*M_PI/180.0;
    double maxX     = model.xmax, X, tet;
    double amp      = Rad*sin(M_PI/2.0) - Rad*sin(maxAngle+M_PI/2);
    
    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        
        if (model.polar==0) {
            topo_chain->z[k]     = Rad;//b / (1 + pow(topo_chain->x[k]/a,2.0));
        }
        if (model.polar==1) {
            // see PolarCoordinatesStuff.py
            topo_chain->z[k]     =  sqrt((Rad - topo_chain->x[k])*(Rad + topo_chain->x[k]));
        }
        printf("topo = %2.2e tet = %2.2e\n", topo_chain->z[k]*scaling.L/1e3, tet*180/M_PI);
        topo_chain->phase[k] = 0;
    }
    
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {
    int np;
    FILE *read;
    int s1, s2;
    
    // Define dimensions;
    double Lx = (double) (model.xmax - model.xmin) ;
    double Lz = (double) (model.zmax - model.zmin) ;
    double T_init = (model.user2 + zeroC)/scaling.T;
    int    cyl = (int)model.user1;
    double X, Z, xc = 0.0, zc = model.user0/scaling.L;
    double Tgrad = model.user3/(scaling.T/scaling.L);
    double Tsurf = 293.0/(scaling.T);
    double a = 100e3/scaling.L, b = 10e3/scaling.L, Rad=6370e3/scaling.L;
    double maxAngle = 18.0/2.0*M_PI/180.0;
    double maxX     = model.xmax;
    double amp      = Rad*sin(M_PI/2.0) - Rad*sin(maxAngle+M_PI/2);
    
    printf("AMP = %2.2e\n", amp*scaling.L/1000);

    // Loop on particles
    for( np=0; np<particles->Nb_part; np++ ) {
        
        // Standard initialisation of particles
        particles->Vx[np]    = -1.0*particles->x[np]*model.EpsBG;               // set initial particle velocity (unused)
        particles->Vz[np]    =  particles->z[np]*model.EpsBG;                   // set initial particle velocity (unused)
        particles->phase[np] = 0;                                               // same phase number everywhere
        particles->d[np]     = 0;                                               // same grain size everywhere
        particles->phi[np]   = 0.0;                                             // zero porosity everywhere
        particles->rho[np]   = 0;
//        particles->T[np]     = T_init;
        particles->T[np]     = Tsurf + Tgrad*particles->z[np];
        X = particles->x[np]-xc;
        Z = particles->z[np]-zc;
    
        // ------------------------- //
        
//        // Mantle
//        if ( cyl==0 ) {
//            if (particles->z[np] < zc + b / (1 + pow(particles->x[np]/a,2.0)) ) {
//                particles->phase[np] = 1;
//            }
//        }
//        if ( cyl==1 ) {
//            X = particles->x[np] * maxAngle/maxX + M_PI/2.0;
//            Z = zc + Rad*sin(X) - Rad + 0*amp/2.0;
//            if (particles->z[np] < Z )  {
//                particles->phase[np] = 1;
//            }
//        }
        
        // SANITY CHECK
        if (particles->phase[np] > model.Nb_phases) {
            printf("Lazy bastard! Fix your particle phase ID! \n");
            exit(144);
        }

        //--------------------------//
        // DENSITY
        if ( model.eqn_state > 0 ) {
            particles->rho[np] = materials->rho[particles->phase[np]] * (1 -  materials->alp[particles->phase[np]] * (T_init - materials->T0[particles->phase[np]]) );
        }
        else {
            particles->rho[np] = materials->rho[particles->phase[np]];
        }
        //--------------------------//
    }
    MinMaxArray(particles->Vx, scaling.V, particles->Nb_part, "Vxp init" );
    MinMaxArray(particles->Vz, scaling.V, particles->Nb_part, "Vzp init" );
    MinMaxArray(particles->T,  scaling.T, particles->Nb_part, "Tp init"  );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials ) {

    int   kk, k, l, c, c1;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
    double dmin, VzBC, width = 1 / scaling.L, eta = 1e4 / scaling.eta ;
    double Lx, Lz, T1, T2, rate=model->EpsBG,  z_comp=-140e3/scaling.L;
    double Vx_r, Vx_l, Vz_b, Vz_t, Vx_tot, Vz_tot;
    double Lxinit = 1400e3/scaling.L, ShortSwitchV0 = 0.40;
    double Vfix = (50.0/(1000.0*365.25*24.0*3600.0))/(scaling.L/scaling.t); // [50.0 == 5 cm/yr]
    double Inflow = 0.0, VzOutflow = 0.0, x, z, V, tet, r, Vx, Vz;
    int    cyl = (int)model->user1;
    
    // Fix temperature
    double Tgrad = model->user3/(scaling.T/scaling.L);
    double Tsurf = 293.0/(scaling.T);
    
    NX  = mesh->Nx;
    NZ  = mesh->Nz;
    NCX = NX-1;
    NCZ = NZ-1;
    NXVZ = NX+1;
    NZVX = NZ+1;
    
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
    
    // Fix temperature
    for (l=0; l<mesh->Nz-1; l++) {
        for (k=0; k<mesh->Nx-1; k++) {
            c = k + l*(NCX);
            mesh->T[c]     = Tsurf + Tgrad*mesh->zc_coord[l];
        }
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
	
		NX  = mesh->Nx;
		NZ  = mesh->Nz;
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;
        
		for (l=0; l<mesh->Nz+1; l++) {
			for (k=0; k<mesh->Nx; k++) {
				
				c = k + l*(mesh->Nx);
                
                x  = mesh->xg_coord[k];
                z  = mesh->zvx_coord[l];
                Vx = -mesh->xg_coord[k] * model->EpsBG;
                
                if (cyl==1) {
                    r           = sqrt(x*x + z*z);
                    tet         = atan(z/x);
                    V           = sqrt(Vx*Vx);
                    Vx          = -V*sin(tet);
                }
                
                if ( mesh->BCu.type[c] != 30 ) {
                    
                    // Internal points:  -1
                    mesh->BCu.type[c] = -1;
                    mesh->BCu.val[c]  =  0;
                    
                    // Matching BC nodes WEST
                    if (k==0 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = Vx;
                        Inflow           += fabs(mesh->BCu.val[c])*model->dz;
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = Vx;
                        Inflow           += fabs(mesh->BCu.val[c])*model->dz;
                    }
                    
                    // Free slip SOUTH
                    if (l==0  ) {
                        mesh->BCu.type[c] = 13;
                        mesh->BCu.val[c]  =  0;
                    }
                    
                    // Free slip NORTH
                    if ( l==mesh->Nz ) {
                        mesh->BCu.type[c] = 13;
                        mesh->BCu.val[c]  =  0;
                    }
                }
                
			}
		}
    
    VzOutflow = -Inflow/(model->xmax-model->xmin);
    printf("VzOutflow = %2.2e\n", VzOutflow*scaling.V);
		
	
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
	
		NX  = mesh->Nx;
		NZ  = mesh->Nz;
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;
		
		for (l=0; l<mesh->Nz; l++) {
			for (k=0; k<mesh->Nx+1; k++) {
				
				c  = k + l*(mesh->Nx+1);
                
                x  = mesh->xvz_coord[k];
                z  = mesh->zg_coord[l];
                Vz = 0.0;
                
                if (cyl==1) {
                    r           = sqrt(x*x + z*z);
                    tet         = atan(z/x);
                    V           = sqrt(Vx*Vx);
                    Vz          =  V*cos(tet);
//                    printf("%2.2e\n", Vz*scaling.V);
                }
                
                if ( mesh->BCv.type[c] != 30 ) {
                    
                    // Internal points:  -1
                    mesh->BCv.type[c] = -1;
                    mesh->BCv.val[c]  =  0;
                    
                    // Matching BC nodes SOUTH
                    if (l==0 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = VzOutflow;//mesh->zg_coord[l] * model->EpsBG;
                    }
                    
                    // Matching BC nodes NORTH
                    if (l==mesh->Nz-1 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = mesh->zg_coord[l] * model->EpsBG;
                    }
                    
                    // Non-matching boundary WEST
                    if ( (k==0) ) {
                        mesh->BCv.type[c] =   11;
                        mesh->BCv.val[c]  =   Vz;
                    }
                    
                    // Non-matching boundary EAST
                    if ( (k==mesh->Nx) ) {
                        mesh->BCv.type[c] =   11;
                        mesh->BCv.val[c]  =   Vz;
                    }
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
    
    
        NX  = mesh->Nx;
        NZ  = mesh->Nz;
        NCX = NX-1;
        NCZ = NZ-1;
        NXVZ = NX+1;
        NZVX = NZ+1;
        
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
    
    double Ttop = 273.15/scaling.T;
    double Tbot= (model->user2 + zeroC)/scaling.T, Tleft, Tright;
    
    printf("Ttop=%2.2e Tbot=%2.2e\n", Ttop*scaling.T, Tbot*scaling.T);
    
	
		NX  = mesh->Nx;
		NZ  = mesh->Nz;
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;
		
		for (l=0; l<mesh->Nz-1; l++) {
			for (k=0; k<mesh->Nx-1; k++) {
				
				c = k + l*(NCX);
                
                if ( mesh->BCt.type[c] != 30 ) {
                    
                    // WEST
                    if ( k==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typW[l] = 0;
                        mesh->BCt.valW[l] = 0;
                    }
                    
                    // EAST
                    if ( k==NCX-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typE[l] = 0;
                        mesh->BCt.valE[l] = 0;
                    }
                    
                    // SOUTH
                    if ( l==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typS[k] = 1;
                        mesh->BCt.valS[k] = Tbot;
                    }
                    
                    // NORTH
                    if ( l==NCZ-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typN[k] = 1;
                        mesh->BCt.valN[k] = Ttop;
                    }
                    
                    // FREE SURFACE
                    else {
                        if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c+NCX] == 30) {
                            mesh->BCt.type[c] = 1;
                            mesh->BCt.val[c]  = Ttop;
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
