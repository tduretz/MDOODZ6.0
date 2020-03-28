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

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ApplyBC( grid* mesh, params* model ) {
    
    int nx=model->Nx, nz=model->Nz, nzvx=nz+1, nxvz=nx+1;
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
        
    // Strain rate component evaluation
    StrainRateComponents( mesh, scaling, model );
    
    //-----------------------------------------------//
    
    NonNewtonianViscosityGrid ( mesh, &materials, model, *Nmodel, &scaling );
    
    //-----------------------------------------------//
    
    // Evaluate right hand side
    EvaluateRHS( mesh, *model, scaling, materials.rho[0] );
    
    //-----------------------------------------------//
    
    // Fill up the rheological matrices arrays
    RheologicalOperators( mesh, model, &scaling, 0 );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InterpCentroidsToVerticesDouble( double* CentroidArray, double* VertexArray, grid* mesh, params *model, scale *scaling ) {
    
    int k, l, Nx, Nz, Ncx, Ncz, c0, c1;
    double *temp;
    
    int per = model->isperiodic_x;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    // Allocate temporary swelled centroid of size = (ncx+2) * (ncz+2)
    temp = DoodzCalloc((Ncx+2)*(Ncz+2),sizeof(DoodzFP));
    
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
        if (per == 0) c1 = 1           + (l)*(Ncx+2);       // right neighbour
        if (per == 1) c1 = (Ncx+2-1-1) + (l)*(Ncx+2);       // right neighbour
        temp[c0] = temp[c1];
    }
    for (l=1; l<Ncz+1; l++) {
        c0 = (Ncx+1) + (l)*(Ncx+2); // East
        if (per==0) c1 = (Ncx  ) + (l)*(Ncx+2); // left neighbour
        if (per==1) c1 = 1       + (l)*(Ncx+2); // left neighbour
        temp[c0] = temp[c1];
    }
    
    // Corners - assume zero flux
    c0 = (0) + (0)*(Ncx+2);         // South-West
    if (per==0) c1 = (1) + (1)*(Ncx+2);         // up-right neighbour
    if (per==1) c1 = (0) + (1)*(Ncx+2);         // up       neighbour
    temp[c0] = temp[c1];
    c0 = (Ncx+1) + (0)*(Ncx+2);     // South-East
    if (per==0) c1 = (Ncx  ) + (1)*(Ncx+2);     // up-left neighbour
    if (per==1) c1 = (Ncx+1) + (1)*(Ncx+2);     // up      neighbour
    temp[c0] = temp[c1];
    c0 = (0) + (Ncz+1)*(Ncx+2);     // North-West
    if (per==0) c1 = (1) + (Ncz  )*(Ncx+2);     // down-right neighbour
    if (per==1) c1 = (0) + (Ncz  )*(Ncx+2);     // down       neighbour
    temp[c0] = temp[c1];
    c0 = (Ncx+1) + (Ncz+1)*(Ncx+2); // North-West
    if (per==0) c1 = (Ncx  ) + (Ncz  )*(Ncx+2); // down-left neighbour
    if (per==1) c1 = (Ncx+1) + (Ncz  )*(Ncx+2); // down      neighbour
    temp[c0] = temp[c1];
    
#pragma omp parallel for shared( temp, VertexArray, mesh ) private( k,l,c1, c0 )  firstprivate( Ncx,Ncz, Nx, Nz )
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
  
            }
        }
    }
    
    // Free temporary array
    DoodzFree(temp);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InterpVerticesToCentroidsDouble( double* CentroidArray, double* VertexArray, grid* mesh, params *model, scale *scaling ) {
    
    int k, l, Nx, Nz, Ncx, Ncz, c0, c1;
    double *temp;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
#pragma omp parallel for shared( CentroidArray, VertexArray, mesh ) private( k,l,c1, c0 )  firstprivate( Ncx,Ncz, Nx )
    // Fill interior points
    for (k=0; k<Ncx; k++) {
        for (l=0; l<Ncz; l++) {
            c0 = k + l*(Ncx);
            c1 = k + l*(Nx);
            
            if ( mesh->BCp.type[c1] != 30 &&  mesh->BCp.type[c1] != 31 ) {
            CentroidArray[c0] = 0.25*( VertexArray[c1] + VertexArray[c1+1] + VertexArray[c1+Nx] + VertexArray[c1+1+Nx] );
            }
        }
    }
    
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
	
    int nx, nz, nxvz, nzvx, ncx, ncz;
    int k, l, c;
    double eps = 1e-13; // perturbation to avoid zero pressure that results in Nan d(eta)dP in numerical differentiation
	
    nx  = mesh->Nx;
    nz  = mesh->Nz;
    nxvz = nx+1;
    nzvx = nz+1;
    ncx  = nx-1;
    ncz  = nz-1;

    for( l=0; l<nzvx; l++) {
        for( k=0; k<nx; k++) {
            c = k + l*nx;
            
            mesh->u_in[c]  = 0.0;
            
            if ( mesh->BCu.type[c] != 30 ) {
                // Initial velocity field (zero or pure shear)
                if (model->EpsBG == 0) mesh->u_in[c]  = 0.0;
                // Pure shear
                else mesh->u_in[c]  = -mesh->xg_coord[k]*model->EpsBG;
                // Simple shear
                if (model->isperiodic_x == 1) mesh->u_in[c] = 2.0*mesh->zvx_coord[l]*model->EpsBG;
                // Force Dirichlets
                if (mesh->BCu.type[c] == 0) mesh->u_in[c]  = mesh->BCu.val[c];
            }
        }
    }
	
    for( l=0; l<nz; l++) {
        for( k=0; k<nxvz; k++) {
            c = k + l*nxvz;
            
            mesh->v_in[c]  = 0.0;
            
            if ( mesh->BCv.type[c] != 30 ) {
                // Initial velocity field (zero or pure shear)
                if (model->EpsBG == 0) mesh->v_in[c]  = 0.0;
                else mesh->v_in[c]  = mesh->zg_coord[l]*model->EpsBG;
                if (model->isperiodic_x == 1) mesh->v_in[c]  = 0.0;
                // Force Dirichlets
                if (mesh->BCv.type[c] == 0) mesh->v_in[c]  = mesh->BCv.val[c];
            }
        }
    }
    
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {
            c = k + l*ncx;
            
            // Initial pressure field
            if (model->num_deriv==0) mesh->p_in[c]  = 0.0;//eps;
            if (model->num_deriv==1) mesh->p_in[c]  = eps;
        }
    }
    
    printf("Velocity field was set to background pure shear\n");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ComputeLithostaticPressure( grid *mesh, params *model, double RHO_REF, scale scaling, int mode ) {
	
    // Compute lithostatic pressure by cumulative sum of rho*g across model thickness
    
    int nx, nz, ncx, ncz;
    int k, l, c;
    double rho_eff;
    double eps = 1e-13; // perturbation to avoid zero pressure that results in Nan d(eta)dP in numerical differentiation

    
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
  
            }
        }
    }
    
    // + eps
    // Cell center arrays
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {            
            c  = k + l*ncx;
            mesh->p_lith[c] += eps;
            mesh->p[c]       = mesh->p_in[c];
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

//void Interp_TPdphi_centroid2vertices ( grid* mesh, params *model ) {
//    
//    int k, l, Nx, Nz, Ncx, Ncz, k1, c0, c1;
//    
//    Nx = mesh->Nx;
//    Nz = mesh->Nz;
//    Ncx = Nx-1;
//    Ncz = Nz-1;
//    
//    
//#pragma omp parallel for shared( mesh ) private( k, l, k1, c1, c0 )  firstprivate(  model, Ncx, Ncz, Nx, Nz )
//    for ( k1=0; k1<Nx*Nz; k1++ ) {
//        
//        k  = mesh->kn[k1];
//        l  = mesh->ln[k1];
//        
//        c0 = k + l*(Ncx);
//        c1 = k + l*Nx;
//        
//        mesh->T_s[c1]   = 0.0;
//        mesh->P_s[c1]   = 0.0;
//        mesh->d0_s[c1]  = 0.0;
//        mesh->phi1_s[c1] = 0.0;
//        
//        // INNER average T and P
//        if (k>0 && k<Ncx && l>0 && l<Ncz) {
//             mesh->T_s[c1]  = 0.25*( mesh->T[c0] +  mesh->T[c0-Ncx] +  mesh->T[c0-Ncx-1] +  mesh->T[c0-1] );
//             mesh->P_s[c1]  = 0.25*( mesh->p_in[c0] +  mesh->p_in[c0-Ncx] +  mesh->p_in[c0-Ncx-1] +  mesh->p_in[c0-1]);
//             mesh->d0_s[c1]  = 0.25*( mesh->d0[c0] +  mesh->d0[c0-Ncx] +  mesh->d0[c0-Ncx-1] +  mesh->d0[c0-1]);
//             mesh->phi1_s[c1] = 0.25*( mesh->phi[c0] +  mesh->phi[c0-Ncx] +  mesh->phi[c0-Ncx-1] +  mesh->phi[c0-1]);
//        }
//
//        // WEST
//        if (k==0 && (l>0 && l<Ncz) && model->isperiodic_x==0) {
//             mesh->T_s[c1]  = 0.5*( mesh->T[c0] + mesh->T[c0-Ncx] );
//             mesh->P_s[c1]  = 0.5*( mesh->p_in[c0] + mesh->p_in[c0-Ncx] );
//             mesh->d0_s[c1]  = 0.5*( mesh->d0[c0] + mesh->d0[c0-Ncx] );
//             mesh->phi1_s[c1] = 0.5*( mesh->phi[c0] + mesh->phi[c0-Ncx] );
//        }
//
//        // WEST - periodic
//        if (k==0 && (l>0 && l<Ncz) && model->isperiodic_x==1) {
//             mesh->T_s[c1]  = 0.25*( mesh->T[c0] + mesh->T[c0-Ncx] + mesh->T[c0+Ncx-1] + mesh->T[c0-1] );
//             mesh->P_s[c1]  = 0.25*( mesh->p_in[c0] + mesh->p_in[c0-Ncx] + mesh->p_in[c0+Ncx-1] + mesh->p_in[c0-1] );
//             mesh->d0_s[c1]  = 0.25*( mesh->d0[c0] + mesh->d0[c0-Ncx] + mesh->d0[c0+Ncx-1] + mesh->d0[c0-1] );
//             mesh->phi1_s[c1] = 0.25*( mesh->phi[c0] + mesh->phi[c0-Ncx] + mesh->phi[c0+Ncx-1] + mesh->phi[c0-1] );
//        }
//
//        // EAST
//        if (k==Ncx && (l>0 && l<Ncz) && model->isperiodic_x==0) {
//             mesh->T_s[c1]  = 0.5*( mesh->T[c0-1] + mesh->T[c0-Ncx-1]);
//             mesh->P_s[c1]  = 0.5*( mesh->p_in[c0-1] + mesh->p_in[c0-Ncx-1] );
//             mesh->d0_s[c1]  = 0.5*( mesh->d0[c0-1] + mesh->d0[c0-Ncx-1]);
//             mesh->phi1_s[c1] = 0.5*( mesh->phi[c0-1] + mesh->phi[c0-Ncx-1] );
//        }
//
//        // EAST - periodic
//        if (k==Ncx && (l>0 && l<Ncz) && model->isperiodic_x==1) {
//             mesh->T_s[c1]  = 0.25*( mesh->T[c0-1] + mesh->T[c0-Ncx-1] + mesh->T[c0-Ncx-Ncx] + mesh->T[c0-Ncx]);
//             mesh->P_s[c1]  = 0.25*( mesh->p_in[c0-1] + mesh->p_in[c0-Ncx-1] + mesh->p_in[c0-Ncx-Ncx] + mesh->p_in[c0-Ncx] );
//             mesh->d0_s[c1]  = 0.25*( mesh->d0[c0-1] + mesh->d0[c0-Ncx-1] + mesh->d0[c0-Ncx-Ncx] + mesh->d0[c0-Ncx]);
//             mesh->phi1_s[c1] = 0.25*( mesh->phi[c0-1] + mesh->phi[c0-Ncx-1] + mesh->phi[c0-Ncx-Ncx] + mesh->phi[c0-Ncx] );
//        }
//
//
//        // SOUTH
//        if (l==0 && (k>0 && k<Ncx)) {
//             mesh->T_s[c1]  = 0.5*( mesh->T[c0]+ mesh->T[c0-1] );
//             mesh->P_s[c1]  = 0.5*( mesh->p_in[c0] + mesh->p_in[c0-1] );
//             mesh->d0_s[c1]  = 0.5*( mesh->d0[c0]+ mesh->d0[c0-1] );
//             mesh->phi1_s[c1] = 0.5*( mesh->phi[c0] + mesh->phi[c0-1] );
//        }
//
//        // NORTH
//        if (l==Ncz && (k>0 && k<Ncx)) {
//             mesh->T_s[c1]  = 0.5*( mesh->T[c0-Ncx] + mesh->T[c0-Ncx-1] );
//             mesh->P_s[c1]  = 0.5*( mesh->p_in[c0-Ncx] + mesh->p_in[c0-Ncx-1] );
//             mesh->d0_s[c1]  = 0.5*( mesh->d0[c0-Ncx] + mesh->d0[c0-Ncx-1] );
//             mesh->phi1_s[c1] = 0.5*( mesh->phi[c0-Ncx] + mesh->phi[c0-Ncx-1] );
//        }
//
//        // SOUTH-WEST
//        if (l==0 && k==0 && model->isperiodic_x==0) {
//             mesh->T_s[c1]  = mesh->T[c0];
//             mesh->P_s[c1]  = mesh->p_in[c0];
//             mesh->d0_s[c1]  = mesh->d0[c0];
//             mesh->phi1_s[c1] = mesh->phi[c0];
//        }
//
//        // SOUTH-WEST - periodic
//        if (l==0 && k==0 && model->isperiodic_x==1) {
//             mesh->T_s[c1]  = 0.5*(mesh->T[c0] + mesh->T[c0+Ncx-1]);
//             mesh->P_s[c1]  = 0.5*(mesh->p_in[c0] + mesh->p_in[c0+Ncx-1]);
//             mesh->d0_s[c1]  = 0.5*(mesh->d0[c0] + mesh->d0[c0+Ncx-1]);
//             mesh->phi1_s[c1] = 0.5*(mesh->phi[c0] + mesh->phi[c0+Ncx-1]);
//        }
//
//        // NORTH-WEST
//        if (l==Ncz && k==0 && model->isperiodic_x==0) {
//             mesh->T_s[c1]  = mesh->T[c0-Ncx];
//             mesh->P_s[c1]  = mesh->p_in[c0-Ncx];
//             mesh->d0_s[c1]  = mesh->d0[c0-Ncx];
//             mesh->phi1_s[c1] = mesh->phi[c0-Ncx];
//        }
//
//        // NORTH-WEST - periodic
//        if (l==Ncz && k==0 && model->isperiodic_x==1) {
//             mesh->T_s[c1]  = 0.5*(mesh->T[c0-Ncx] + mesh->T[c0-1]);
//             mesh->P_s[c1]  = 0.5*(mesh->p_in[c0-Ncx] + mesh->p_in[c0-1]);
//             mesh->d0_s[c1]  = 0.5*(mesh->d0[c0-Ncx] + mesh->d0[c0-1]);
//             mesh->phi1_s[c1] = 0.5*(mesh->phi[c0-Ncx] + mesh->phi[c0-1]);
//        }
//
//        // SOUTH-EAST
//        if (l==0 && k==Ncx && model->isperiodic_x==0) {
//             mesh->T_s[c1]  = mesh->T[c0-1];
//             mesh->P_s[c1]  = mesh->p_in[c0-1];
//             mesh->d0_s[c1]  = mesh->d0[c0-1];
//             mesh->phi1_s[c1] = mesh->phi[c0-1];
//        }
//
//        // SOUTH-EAST - periodic
//        if (l==0 && k==Ncx && model->isperiodic_x==1) {
//             mesh->T_s[c1]  = 0.5*(mesh->T[c0-1] + mesh->T[c0-Ncx]);
//             mesh->P_s[c1]  = 0.5*(mesh->p_in[c0-1]+ mesh->p_in[c0-Ncx]);
//             mesh->d0_s[c1]  = 0.5*(mesh->d0[c0-1] + mesh->d0[c0-Ncx]);
//             mesh->phi1_s[c1] = 0.5*(mesh->phi[c0-1]+ mesh->phi[c0-Ncx]);
//        }
//
//        // NORTH-EAST
//        if (l==Ncz && k==Ncx && model->isperiodic_x==0) {
//             mesh->T_s[c1]  = mesh->T[c0-Ncx-1];
//             mesh->P_s[c1]  = mesh->p_in[c0-Ncx-1];
//             mesh->d0_s[c1]  = mesh->d0[c0-Ncx-1];
//             mesh->phi1_s[c1] = mesh->phi[c0-Ncx-1];
//        }
//
//        // NORTH-EAST - periodic
//        if (l==Ncz && k==Ncx && model->isperiodic_x==1) {
//             mesh->T_s[c1]  = 0.5*(mesh->T[c0-Ncx-1] + mesh->T[c0-Ncx-Ncx]);
//             mesh->P_s[c1]  = 0.5*(mesh->p_in[c0-Ncx-1] + mesh->p_in[c0-Ncx-Ncx]);
//             mesh->d0_s[c1]  = 0.5*(mesh->d0[c0-Ncx-1] + mesh->d0[c0-Ncx-Ncx]);
//             mesh->phi1_s[c1] = 0.5*(mesh->phi[c0-Ncx-1] + mesh->phi[c0-Ncx-Ncx]);
//        }
//    }
//}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SetUpModel_NoMarkers ( grid* mesh, params *model, scale *scaling ) {
    
    int k, l, Nx, Nz, Ncx, Ncz, k1, c0, c1;
    double x, z;
    double radius = model->user1/scaling->L;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    printf("Setting up mode without using markers --- DEBUG !!!!\n");
    
    // Vertices
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {
        
        k  = mesh->kp[k1];
        l  = mesh->lp[k1];
        c0 = k + l*(Ncx);
        
        x = mesh->xc_coord[k];
        z = mesh->zc_coord[l];
        
        mesh->T[k1] = 0.05;
        
        mesh->phase_perc_n[0][k1] = 1.0;
        mesh->phase_perc_n[1][k1] = 0.0;
        
        if (x*x + z*z < radius*radius) {
            mesh->phase_perc_n[0][k1] = 0.0;
            mesh->phase_perc_n[1][k1] = 1.0;
        }
    }
    
    
    
    // Vertices
    for ( k1=0; k1<Nx*Nz; k1++ ) {
        
        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;
        
        x = mesh->xg_coord[k];
        z = mesh->zg_coord[l];
        
        mesh->phase_perc_s[0][k1] = 1.0;
        mesh->phase_perc_s[1][k1] = 0.0;
        
        if (x*x + z*z < radius*radius) {
            mesh->phase_perc_s[0][k1] = 0.0;
            mesh->phase_perc_s[1][k1] = 1.0;
        }
    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
