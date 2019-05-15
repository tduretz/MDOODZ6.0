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
#include "header_MDOODZ.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num() 0
#endif

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ApplyBC( grid* mesh, params model ) {
    
    int nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1;
    int i, j;
    
    // Vx Neumann
    for( i=0; i<nx; i++) {
        // South
        if ( mesh->BCu.type[i] == 13 ) {
            mesh->u_in[i] = mesh->u_in[i+nx];
        }
        // North
        if ( mesh->BCu.type[i + (nzvx-1)*nx] == 13 ) {
            mesh->u_in[i + (nzvx-1)*nx] = mesh->u_in[i + (nzvx-2)*nx];
        }
    }
    
    // Vx Dirichlet
    for( i=0; i<nx; i++) {
        // South
        if ( mesh->BCu.type[i] == 11 ) {
            mesh->u_in[i] = 2.0*mesh->BCu.val[i] - mesh->u_in[i+nx];
        }
        // North
        if ( mesh->BCu.type[i + (nzvx-1)*nx] == 11 ) {
            mesh->u_in[i + (nzvx-1)*nx] = 2.0*mesh->BCu.val[i + (nzvx-1)*nx] - mesh->u_in[i + (nzvx-2)*nx];
        }
    }
    
    // Vz Neumann
    for( j=0; j<nz; j++) {
        // West
        if ( mesh->BCv.type[j*nxvz] == 13 ) {
            mesh->v_in[j*nxvz] = mesh->v_in[j*nxvz + 1];
        }
        // East
        if ( mesh->BCv.type[j*nxvz+(nxvz-1)] == 13 ) {
            mesh->v_in[j*nxvz+(nxvz-1)] = mesh->v_in[j*nxvz+(nxvz-1) -1];
        }
    }
    
    // Vz Dirichlet
    for( j=0; j<nz; j++) {
        // West
        if ( mesh->BCv.type[j*nxvz] == 11 ) {
            mesh->v_in[j*nxvz] = 2.0*mesh->BCv.val[j*nxvz] - mesh->v_in[j*nxvz + 1];
        }
        // East
        if ( mesh->BCv.type[j*nxvz+(nxvz-1)] == 11 ) {
            mesh->v_in[j*nxvz + (nxvz-1)] = 2.0*mesh->BCv.val[j*nxvz+(nxvz-1)] - mesh->v_in[j*nxvz + (nxvz-1) -1];
        }
    }
    
    // Vz Periodic
    for( j=0; j<nz; j++) {
        // West
        if ( mesh->BCv.type[j*nxvz] == -12 ) {
            mesh->v_in[j*nxvz] = mesh->v_in[j*nxvz + (nxvz-1) -1];
        }
        // East
        if ( mesh->BCv.type[j*nxvz+(nxvz-1)] == -12 ) {
            mesh->v_in[j*nxvz+(nxvz-1)] = mesh->v_in[j*nxvz+1];
        }
    }
    
    // Vx Periodic
    for( j=0; j<nzvx; j++) {
        if ( mesh->BCu.type[j*nx + nx-1] == -12 ) {
            mesh->u_in[j*nx + nx-1] = mesh->u_in[j*nx];
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateNonLinearity( grid* mesh, markers* particles, markers* topo_chain, surface *topo, mat_prop materials, params *model, Nparams *Nmodel, scale scaling, int mode, double h_contin ) {
        
    int strain_rate_formulation = 1;
    if ( model->iselastic == 1 ) strain_rate_formulation = 0;
    
    // Strain rate component evaluation
    StrainRateComponents( mesh, scaling, model );
    
    //-----------------------------------------------//
    
    if ( model->rheo_on_cells == 0 ) NonNewtonianViscosityGrid ( mesh, &materials, model, *Nmodel, &scaling, strain_rate_formulation );
    if ( model->rheo_on_cells == 1 ) NonNewtonianViscosityCells( mesh, &materials, model, *Nmodel, &scaling, strain_rate_formulation );
    //    NonNewtonianViscosityMarkers(  particles, mesh, &materials, model, Emodel, Mmodel, *Nmodel, &scaling, strain_rate_formulation );
    
    //-----------------------------------------------//
    
//    // Force upper cells to be weak
//    int Nx = mesh->Nx;
//    int Nz = mesh->Nz;
//
//
//
//    for (int i=0;i<Nx-1;i++) {
//        for (int j=0;j<Nz-1;j++) {
//
//            int ii = i + j*(Nx-1);
//
//            if (j<Nz-2) {
//                if (mesh->BCp.type[ii] == 31 && mesh->BCp.type[ii+(Nx-1)] == 30  ) {
//                    //                if (mesh->BCp.type[ii] == 0 && mesh->BCp.type[ii+(Nx-1)] == 30  ) {
////                    printf("Upper cells\n");
//                    mesh->eta_n[ii] = 32/scaling.eta;
////                    mesh->eta_phys_n[ii] = 32/scaling.eta;
//                    //                                            printf("Upper cells\n");
//                }
//            }
//
//                        if (j<Nz-3) {
//                            if (mesh->BCp.type[ii] == -1 && mesh->BCp.type[ii+2*(Nx-1)] == 30  ) {
//                                mesh->eta_n[ii] = 32/scaling.eta;
//                                mesh->eta_phys_n[ii] = 32/scaling.eta;
////                                                   printf("Upper cells\n");
//                            }
//                        }
//
//        }
//    }
//
//    for (int i=0;i<Nx;i++) {
//        for (int j=0;j<Nz;j++) {
//
//            int ii = i + j*(Nx);
//
//            if (j<Nz-1) {
//                if (mesh->BCg.type[ii] == -1 && mesh->BCg.type[ii+Nx] == 30 ) {
//                    mesh->eta_s[ii] = 32/scaling.eta;
//                    mesh->eta_phys_s[ii] = 32/scaling.eta;
//                    //                                            printf("Upper vertices\n");
//                }
//            }
//
//
//                        if (j<Nz-2) {
//                            if (mesh->BCg.type[ii] == -1 && mesh->BCg.type[ii+2*Nx] == 30 ) {
//                                mesh->eta_s[ii] = 32/scaling.eta;
//                                mesh->eta_phys_s[ii] = 32/scaling.eta;
////                                                     printf("Upper vertices\n");
//                            }
//                        }
//
//        }
//    }


    //-----------------------------------------------//
    
    // Evaluate right hand side
    EvaluateRHS( mesh, *model, scaling, materials.rho[0] );
    
    //-----------------------------------------------//
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InterpCentroidsToVerticesDouble( double* CentroidArray, double* VertexArray, grid* mesh, params *model, scale *scaling ) {
    
    int k, l, Nx, Nz, Ncx, Ncz, c0, c1;
    double *temp;
    double accu = model->accu;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    // Allocate temporary swelled centroid of size = (ncx+2) * (ncz+2)
    temp = DoodzCalloc(sizeof(DoodzFP), (Ncx+2)*(Ncz+2));
    
    // Fill interior points
    for (k=0; k<Ncx; k++) {
        for (l=0; l<Ncz; l++) {
            c0 = k + l*(Ncx);
            c1 = k + (l+1)*(Ncx+2) + 1;
            temp[c1] = CentroidArray[c0];
        }
    }
    
    // Fill sides - avoid corners - assume zero flux
    for (k=1; k<Ncx+1; k++) {
        c0 = k + (0)*(Ncx+2);       // South
        c1 = k + (1)*(Ncx+2);       // up neighbour
        temp[c0] = temp[c1];
    }
    for (k=1; k<Ncx+1; k++) {
        c0 = k + (Ncz+1)*(Ncx+2);   // North
        c1 = k + (Ncz  )*(Ncx+2);   // down neighbour
        temp[c0] = temp[c1];
    }
    for (l=1; l<Ncz+1; l++) {
        c0 = 0 + (l)*(Ncx+2);       // West
        c1 = 1 + (l)*(Ncx+2);       // right neighbour
        temp[c0] = temp[c1];
    }
    for (l=1; l<Ncz+1; l++) {
        c0 = (Ncx+1) + (l)*(Ncx+2); // East
        c1 = (Ncx  ) + (l)*(Ncx+2); // left neighbour
        temp[c0] = temp[c1];
    }
    
    // Corners - assume zero flux
    c0 = (0) + (0)*(Ncx+2);         // South-West
    c1 = (1) + (1)*(Ncx+2);         // up-right neighbour
    temp[c0] = temp[c1];
    c0 = (Ncx+1) + (0)*(Ncx+2);     // South-East
    c1 = (Ncx  ) + (1)*(Ncx+2);     // up-left neighbour
    temp[c0] = temp[c1];
    c0 = (0) + (Ncz+1)*(Ncx+2);     // North-West
    c1 = (1) + (Ncz  )*(Ncx+2);     // down-right neighbour
    temp[c0] = temp[c1];
    c0 = (Ncx+1) + (Ncz+1)*(Ncx+2); // North-West
    c1 = (Ncx  ) + (Ncz  )*(Ncx+2); // down-left neighbour
    temp[c0] = temp[c1];
    
    // interpolate from temp array to actual vertices array
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz; l++) {
            
            c1 = k + (l)*(Nx);
            c0 = k + (l+1)*(Ncx+2) + 1;
            
            // Default 0 - above free surface
            VertexArray[c1] = 0.0;
            
            // Else interpolate
            if ( mesh->BCg.type[c1] != 30 ) {
                VertexArray[c1] = 0.25*( temp[c0-1-(Ncx+2)] + temp[c0-0-(Ncx+2)] + temp[c0-1] + temp[c0] );
                if ( model->cut_noise==1 ){
                    VertexArray[c1] = round(VertexArray[c1]*accu)/accu;
                }
                
            }
        }
    }
    
    // Free temporary array
    DoodzFree(temp);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Initialise grid
void SetGridCoordinates( grid *mesh, params *model, int nx, int nz ) {
    
    // Set coordinates of mesh vertices and cell centroids, Vx and Vz nodes
	
    int i;
    double dx0, dz0;
    
    // Set spacing
    mesh->dx = (model->xmax-model->xmin) / (nx-1);
    mesh->dz = (model->zmax-model->zmin) / (nz-1);
    
    model->dx = mesh->dx;
    model->dz = mesh->dz;
	
    // Gridlines positions
    mesh->xg_coord[0] = model->xmin;
    mesh->zg_coord[0] = model->zmin;
	
    for (i = 1; i<nx; i++) {
        mesh->xg_coord[i] = mesh->xg_coord[i-1] + mesh->dx;
    }
	
    for (i = 1; i<nz; i++) {
        mesh->zg_coord[i] = mesh->zg_coord[i-1] + mesh->dz;
    }
    
    //---------------
    
    // Initial mesh
    // Set spacing
    dx0 = (model->xmax0-model->xmin0) / (nx-1);
    dz0 = (model->zmax0-model->zmin0) / (nz-1);
    
    // Gridlines positions
    mesh->xg_coord0[0] = model->xmin0;
    mesh->zg_coord0[0] = model->zmin0;
    
    for (i = 1; i<nx; i++) {
        mesh->xg_coord0[i] = mesh->xg_coord0[i-1] + dx0;
    }
    
    for (i = 1; i<nz; i++) {
        mesh->zg_coord0[i] = mesh->zg_coord0[i-1] + dz0;
    }
    
    //---------------
    
    // Cell centers positions
    mesh->xc_coord[0] = model->xmin + mesh->dx/2.0;
    mesh->zc_coord[0] = model->zmin + mesh->dz/2.0;
    
    for (i = 1; i<nx-1; i++) {
        mesh->xc_coord[i]  = mesh->xc_coord[i-1]  + mesh->dx;
    }
    for (i = 1; i<nz-1; i++) {
        mesh->zc_coord[i]  = mesh->zc_coord[i-1]  + mesh->dz;
    }
    
    // Ghosted nodes
    mesh->xvz_coord[0] = model->xmin - mesh->dx/2.0;
    mesh->zvx_coord[0] = model->zmin - mesh->dz/2.0;
    
    for (i = 1; i<nx+1; i++) {
        mesh->xvz_coord[i] = mesh->xvz_coord[i-1] + mesh->dx;
    }
    for (i = 1; i<nz+1; i++) {
        mesh->zvx_coord[i] = mesh->zvx_coord[i-1] + mesh->dz;
    }
    
    // Ghosted nodes
    mesh->xg_coord_ext[0] = model->xmin - mesh->dx;
    mesh->zg_coord_ext[0] = model->zmin - mesh->dz;
    
    for (i = 1; i<nx+2; i++) {
        mesh->xg_coord_ext[i] = mesh->xg_coord_ext[i-1] + mesh->dx;
    }
    for (i = 1; i<nz+2; i++) {
        mesh->zg_coord_ext[i] = mesh->zg_coord_ext[i-1] + mesh->dz;
    }
	
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InitialiseSolutionFields( grid *mesh, params *model ) {
    
    // Set initial velocity and pressure fields
	
    int nx,nz,nxvz,nzvx;
    int k,l,c;
	
    nx  = mesh->Nx;
    nz  = mesh->Nz;
    nxvz = nx+1;
    nzvx = nz+1;

    for( l=0; l<nzvx; l++) {
        for( k=0; k<nx; k++) {
            c = k + l*nx;
            
            // Initial velocity field (zero or pure shear)
            if (model->EpsBG == 0) mesh->u_in[c]  = 0.0;
            else mesh->u_in[c]  = -mesh->xg_coord[k]*model->EpsBG;
            
            if (model->isperiodic_x == 1) mesh->u_in[c] = 2.0*mesh->zvx_coord[l]*model->EpsBG;
  
        }
    }
	
    for( l=0; l<nz; l++) {
        for( k=0; k<nxvz; k++) {
            c = k + l*nxvz;
            
            // Initial velocity field (zero or pure shear)
            if (model->EpsBG == 0) mesh->v_in[c]  = 0.0;
            else mesh->v_in[c]  = mesh->zg_coord[l]*model->EpsBG;
            
            if (model->isperiodic_x == 1) mesh->v_in[c]  = 0.0;
	
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ComputeLithostaticPressure( grid *mesh, params *model, double RHO_REF, scale scaling, int mode ) {
	
    // Compute lithostatic pressure by cumulative sum of rho*g across model thickness
    
    int nx, nz, ncx, ncz;
    int k, l, c;
    double rho_eff;
	
    nx  = mesh->Nx;
    nz  = mesh->Nz;
    ncx = nx-1;
    ncz = nz-1;
    
    Initialise1DArrayDouble( mesh->p_lith,  (mesh->Nx-1)*(mesh->Nz-1), 0.0 );

    // Cell center arrays
    for( l=ncz-2; l>=0; l--) {
        for( k=0; k<ncx; k++) {
            
            // Initialise vertices variables
            c  = k + l*ncx;
            
            mesh->p_lith[c]  = 0.0;
            
            // density
            if ( mode == 0 ) rho_eff = RHO_REF;
            if ( mode == 1 ) rho_eff = mesh->rho_app_n[c];

            // Initialise pressure variables : Compute lithostatic pressure
            if ( mesh->BCp.type[c] != 30 && mesh->BCp.type[c] != 31 ) { // First row (surface)
                mesh->p_lith[c]  = mesh->p_lith[c+ncx] -  model->gz * mesh->dz * rho_eff;
                mesh->p[c]       = mesh->p_in[c];
            }
        }
    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void GridIndices( grid* mesh) {
    
// Compute i, j indices of flattened arrays (vertices, centroids, Vx, Vz)
    
    int k,l,kk;
    
    // Cell vertices
    for( l=0; l<mesh->Nz; l++) {
        for( k=0; k<mesh->Nx; k++) {
            kk = k + l*(mesh->Nx);
            mesh->kn[kk] = k;
            mesh->ln[kk] = l;
        }
    }
    
    // Vx nodes
    for( l=0; l<mesh->Nz+1; l++) {
        for( k=0; k<mesh->Nx; k++) {
            kk = k + l*mesh->Nx;
            mesh->kvx[kk] = k;
            mesh->lvx[kk] = l;
        }
    }
    
    // Vz nodes
    for( l=0; l<mesh->Nz; l++) {
        for( k=0; k<mesh->Nx+1; k++) {
            kk = k + l*(mesh->Nx+1);
            
            mesh->kvz[kk] = k;
            mesh->lvz[kk] = l;
        }
    }
    
    // Cell Centroids
    for( l=0; l<mesh->Nz-1; l++) {
        for( k=0; k<mesh->Nx-1; k++) {
            kk = k + l*(mesh->Nx-1);
            mesh->kp[kk] = k;
            mesh->lp[kk] = l;
        }
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/