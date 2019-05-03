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
#include "header_MDOODZ.h"
#include "umfpack.h"
#include "amd.h"
#include "cs.h"
#include "cholmod.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

#define error printf

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void ExtractSolutions( SparseMat *Stokes, grid* mesh, params model ) {
    
    int cc, nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1, kk;
    
    // -------- Get U -------- //
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( nzvx, nx )
    for( cc=0; cc<nzvx*nx; cc++) {
        
        mesh->u_in[cc] = 0.0;
        
        if ( mesh->BCu.type[cc] != 30 ) {
            if ( mesh->BCu.type[cc] == 0 ) {
                mesh->u_in[cc] = mesh->BCu.val[cc];
            }
            
            if ( mesh->BCu.type[cc] < 30 && mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13  && mesh->BCu.type[cc] != -12 ) {
                mesh->u_in[cc] = Stokes->x[Stokes->eqn_u[cc]];
            }
        }
    }
    
    // Periodic
    for( cc=0; cc<nzvx; cc++) {
        kk = cc*nx + nx-1;
        if ( mesh->BCu.type[kk] ==-12) {
            mesh->u_in[kk] = mesh->u_in[kk-mesh->Nx+1];
        }
    }
    
    // -------- Get V -------- //
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( nz, nxvz )
    for( cc=0; cc<nz*nxvz; cc++) {
        
        mesh->v_in[cc] = 0.0;
        
        if ( mesh->BCv.type[cc] != 30 ) {
            if ( mesh->BCv.type[cc] == 0 ) {
                mesh->v_in[cc] = mesh->BCv.val[cc];
            }
            
            if ( mesh->BCv.type[cc] < 30 && mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12 ) {
                mesh->v_in[cc] = Stokes->x[Stokes->eqn_v[cc]];
                
            }
        }
    }
    
    // -------- Get P -------- //
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( nz, nxvz ) //reduction(+:num_p)
    
    for( cc=0; cc<ncz*ncx; cc++) {
        
        mesh->p_in[cc] = 0.0;
        if ( mesh->BCp.type[cc] != 30 ) {
            if ( mesh->BCp.type[cc] == 0 || mesh->BCp.type[cc] == 31 ) {
                mesh->p_in[cc] = mesh->BCp.val[cc];
            }
            else {
                mesh->p_in[cc] = Stokes->x[Stokes->eqn_p[cc]];
            }
        }
    }
    
    //    mean_p /=num_p;
    //
    //    for( cc=0; cc<ncz*ncx; cc++) {
    //
    //        mesh->p_in[cc] = 0.0;
    //        if ( mesh->BCp.type[cc] != 30 ) {
    //            if ( mesh->BCp.type[cc] == 0 || mesh->BCp.type[cc] == 31 ) {
    //                mesh->p_in[cc] = mesh->BCp.val[cc];
    //            }
    //            else {
    //                mesh->p_in[cc]  = 0*Stokes->x[Stokes->eqn_p[cc]] - mean_p;
    //            }
    //        }
    //    }
    
    
    
    // Apply Bc to Vx and Vz
    ApplyBC( mesh, model );
    
    //            int c1, i, j;
    //
    //    double symVx[nx*(nz+1)];
    //    double symVz[nz*(nx+1)];
    //
    //    printf("%lf\n", ceil(nx/2));
    //    for( j=0; j<nz+1; j++) {
    //        for( i=0; i<(int)ceil(nx/2); i++) {
    //            cc = i + j*(nx);
    //            c1 =-i + j*(nx) + nx-1;
    //            symVx[cc] = fabs(mesh->u_in[cc] + mesh->u_in[c1]);
    //            printf("%.2e ", symVx[cc]);
    //        }
    //        printf("\n");
    //    }
    
    
    //    printf("Vx\n");
    //
    //    for( j=0; j<nz+1; j++) {
    //        for( i=0; i<nx; i++) {
    //            cc = i + j*(nx);
    //            c1 = i + j*(nx-0);
    //            printf("%.2e ", mesh->u_in[cc]);
    //        }
    //        printf("\n");
    //    }
    
    
    
    //   printf("%d\n", (int)ceil((nx+1)/2)+0);
    //    //    int c1, i, j;
    //    for( j=0; j<nz; j++) {
    //        for( i=0; i<(int)ceil((nx+1)/2)+0; i++) {
    //            cc = i + j*(nx+1);
    //            c1 =-i + j*(nx+1) + nx;
    //            symVz[cc] = fabs(mesh->v_in[cc] - mesh->v_in[c1]);
    //            printf("%.2e ", symVz[cc]);
    //        }
    //        printf("\n");
    //    }
    
    
    //    printf("Vz\n");
    ////    int c1, i, j;
    //    for( j=0; j<nz; j++) {
    //        for( i=0; i<nx+1; i++) {
    //            cc = i + j*(nx+1);
    //            c1 = i + j*(nx-0);
    //            printf("%.6e ", mesh->v_in[cc]);
    //        }
    //        printf("\n");
    //    }
    
    //    printf("p\n");
    //    //    int c1, i, j;
    //    for( j=0; j<nz-1; j++) {
    //        for( i=0; i<nx-1; i++) {
    //            cc = i + j*(nx-1);
    //            c1 = i + j*(nx-0);
    //            printf("%.6e ", mesh->p_in[cc]);
    //        }
    //        printf("\n");
    //    }
    
    
    
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double LineSearch( SparseMat *Stokes, double *dx, grid* mesh, params *model, Nparams* Nmodel, markers* particles, markers *topo_chain, surface *topo, mat_prop materials, scale scaling ) {
    
    int ix, kk, k, cc;
    int nx= model->Nx, nz=model->Nz, ncx=nx-1, ncz=nz-1, nzvx=nz+1, nxvz=nx+1;
    double *rx, *rz, *rp, alpha, *alphav, dalpha, minx;
    double *u, *v, *p;
    Nparams residuals;
    clock_t  t_omp;
    int success = 0, niter=0, nitermax=1;
    double minalpha[nitermax], maxalpha[nitermax], frac = 1.0;
    int ntry[nitermax];
    
    Nmodel->stagnated = 0;
    
    ntry[0]     = 6;
    minalpha[0] = -2.5;
    maxalpha[0] = -0.0;
    
    //    ntry[0]     = 11;
    //    minalpha[0] = -2.25;
    //    maxalpha[0] = -0.25;
    
    //    ntry[0]     = 1;
    //    minalpha[0] = -1.0;
    //    maxalpha[0] = -1.0;
    
    t_omp = (double)omp_get_wtime();
    
    // Allocate
    u          = DoodzMalloc( nx*nzvx*sizeof(double) );
    v          = DoodzMalloc( nxvz*nz*sizeof(double) );
    p          = DoodzMalloc( ncx*ncz*sizeof(double) );
    
    // Save initial fields (before updates)
    ArrayEqualArray( u, mesh->u_in, nx*nzvx );
    ArrayEqualArray( v, mesh->v_in, nxvz*nz );
    ArrayEqualArray( p, mesh->p_in, ncx*ncz );
    
    // Do line search iterations
    while ( success == 0 && niter < nitermax ) {
        
        // allocate array
        alphav = DoodzMalloc( ntry[niter]*sizeof(double) );
        rx     = DoodzMalloc( ntry[niter]*sizeof(double) );
        rz     = DoodzMalloc( ntry[niter]*sizeof(double) );
        rp     = DoodzMalloc( ntry[niter]*sizeof(double) );
        
        alpha    = maxalpha[niter];
        dalpha   = fabs(maxalpha[niter]-minalpha[niter])/(ntry[niter]-1);
        
        // Search for optimal relaxation parameters
        for( kk=0; kk<ntry[niter]; kk++ ) {
            
            // Update alpha
            if (kk>0) alpha -= dalpha;
            alphav[kk] = alpha;
            
            // Start from initial solutions
            ArrayEqualArray( mesh->u_in, u, nx*nzvx );
            ArrayEqualArray( mesh->v_in, v, nxvz*nz );
            ArrayEqualArray( mesh->p_in, p, ncx*ncz );
            
            // Test relaxed u-v-p solutions
#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, nzvx, nx )
            for( cc=0; cc<nzvx*nx; cc++) {
                if ( mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12  ) {
                    mesh->u_in[cc] += alpha*dx[Stokes->eqn_u[cc]];
                }
            }
            
#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, nz, nxvz )
            for( cc=0; cc<nz*nxvz; cc++) {
                if ( mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12 ) {
                    mesh->v_in[cc] += alpha*dx[Stokes->eqn_v[cc]];
                }
            }
            
#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, ncx, ncz )
            for( cc=0; cc<ncz*ncx; cc++) {
                
                if ( mesh->BCp.type[cc] != 30 && mesh->BCp.type[cc] != 0  && mesh->BCp.type[cc] != 31 ) {
                    mesh->p_in[cc] += alpha*dx[Stokes->eqn_p[cc]];
                }
            }
            
            // Apply Bc to Vx and Vz
            ApplyBC( mesh, *model );
            
            //            Print2DArrayDouble( u,  nx, nzvx, scaling.V );
            //            Print2DArrayDouble( v,  nxvz, nz, scaling.V );
            //            Print2DArrayDouble( p,  ncx, ncz, scaling.S );
            
            //------------------------------------------------------------------------------------------------------
            // Update non-linearity
            UpdateNonLinearity( mesh, particles, topo_chain, topo, materials, model, Nmodel, scaling, 0, 1.0 );
            
            //------------------------------------------------------------------------------------------------------
            
            // Calculate residual
//            model->free_surf_stab = 0;
            EvaluateStokesResidual( Stokes, &residuals, mesh, *model, scaling, 1 );
//            model->free_surf_stab = dummy;
            
            rx[kk] = residuals.resx;
            rz[kk] = residuals.resz;
            rp[kk] = residuals.resp;
            printf("\e[1;34mAlpha\e[m = %lf --> rx = %2.4e rz = %2.4e rp = %2.2e\n", alpha, rx[kk]* (scaling.F/pow(scaling.L,3)), rz[kk]* (scaling.F/pow(scaling.L,3)), rp[kk]*scaling.E );
            
        }
        
        // Look for the minimum predicted residuals
        minx  = rz[0];
        ix = 0;
        for( k=1; k<ntry[niter]; k++ ) {
            if(rz[k]<minx) {
                minx = rz[k];
                ix = k;
            }
        }
        alpha = alphav[ix];
        
        // if the minmimun residuals are lower than starting ones, then success
        if ( rx[ix] < frac*Nmodel->resx_f || rz[ix]<frac*Nmodel->resz_f  ) { //|| rp[ix]<frac*Nmodel->resp
            success = 1;
            printf("\e[1;34mPredicted Residuals\e[m : alpha  = %lf --> rx = %2.4e rz = %2.4e rp = %2.4e\n", alphav[ix], rx[ix]* (scaling.F/pow(scaling.L,3)), rz[ix]* (scaling.F/pow(scaling.L,3)), rp[ix]*scaling.E);
        }
        
        DoodzFree(rx);
        DoodzFree(rz);
        DoodzFree(rp);
        DoodzFree(alphav);
        niter++;
    }
    
    if ( fabs(alpha)<1e-13 || success == 0 ) {
        printf( "Found minimum of the function -- iteration cycle stagnates\n...\n" );
        Nmodel->stagnated = 1;
        alpha = 1.0;
    }
    
    // Reset solutions
    ArrayEqualArray( mesh->u_in, u, nx*nzvx );
    ArrayEqualArray( mesh->v_in, v, nxvz*nz );
    ArrayEqualArray( mesh->p_in, p, ncx*ncz );
    
    // Free
    DoodzFree( u );
    DoodzFree( v );
    DoodzFree( p );
    printf("** Line search took = %f sec\n", (double)((double)omp_get_wtime() - t_omp) );
    
    return alpha;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double LineSearchDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, double *dx, grid* mesh, params *model, Nparams* Nmodel, markers* particles, markers *topo_chain, surface *topo, mat_prop materials, scale scaling ) {
    
    int ix, kk, k, cc;
    int nx= model->Nx, nz=model->Nz, ncx=nx-1, ncz=nz-1, nzvx=nz+1, nxvz=nx+1;
    double *rx, *rz, *rp, alpha, *alphav, dalpha, minx;
    double *u, *v, *p;
    Nparams residuals;
    clock_t  t_omp;
    int success = 0, niter=0, nitermax=1;
    double minalpha[nitermax], maxalpha[nitermax], frac = 1.0;
    int ntry[nitermax];
    
    Nmodel->stagnated = 0;
    
    ntry[0]     = 6;
    minalpha[0] = -2.5;
    maxalpha[0] = -0.0;
    
    if (model->Newton==1) minalpha[0] = -1.0;
    
    //    ntry[0]     = 11;
    //    minalpha[0] = -2.25;
    //    maxalpha[0] = -0.25;
    
    //    ntry[0]     = 1;
    //    minalpha[0] = -1.0;
    //    maxalpha[0] = -1.0;
    
    //    printf("dalpha %lf %lf %d \n", maxalpha[0], minalpha[0], ntry[0]-1);
    
    t_omp = (double)omp_get_wtime();
    
    // Allocate
    u          = DoodzMalloc( nx*nzvx*sizeof(double) );
    v          = DoodzMalloc( nxvz*nz*sizeof(double) );
    p          = DoodzMalloc( ncx*ncz*sizeof(double) );
    
    // Save initial fields (before updates)
    ArrayEqualArray( u, mesh->u_in, nx*nzvx );
    ArrayEqualArray( v, mesh->v_in, nxvz*nz );
    ArrayEqualArray( p, mesh->p_in, ncx*ncz );
    
    // Do line search iterations
    while ( success == 0 && niter < nitermax ) {
        
        // allocate array
        alphav = DoodzMalloc( ntry[niter]*sizeof(double) );
        rx     = DoodzMalloc( ntry[niter]*sizeof(double) );
        rz     = DoodzMalloc( ntry[niter]*sizeof(double) );
        rp     = DoodzMalloc( ntry[niter]*sizeof(double) );
        
        alpha    = maxalpha[niter];
        dalpha   = fabs(maxalpha[niter]-minalpha[niter])/(ntry[niter]-1);
        
        // Search for optimal relaxation parameters
        for( kk=0; kk<ntry[niter]; kk++ ) {
            
            // Update alpha
            if (kk>0) alpha -= dalpha;
            alphav[kk] = alpha;
            
            // Start from initial solutions
            ArrayEqualArray( mesh->u_in, u, nx*nzvx );
            ArrayEqualArray( mesh->v_in, v, nxvz*nz );
            ArrayEqualArray( mesh->p_in, p, ncx*ncz );
            
            // Test relaxed u-v-p solutions
#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, nzvx, nx )
            for( cc=0; cc<nzvx*nx; cc++) {
                if ( mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12 ) {
                    mesh->u_in[cc] += alpha*dx[Stokes->eqn_u[cc]];
                }
            }
            
#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, nz, nxvz )
            for( cc=0; cc<nz*nxvz; cc++) {
                if ( mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12 ) {
                    mesh->v_in[cc] += alpha*dx[Stokes->eqn_v[cc]];
                }
            }
            
#pragma omp parallel for shared( mesh, dx, Stokes ) private( cc ) firstprivate( alpha, ncx, ncz )
            for( cc=0; cc<ncz*ncx; cc++) {
                
                if ( mesh->BCp.type[cc] != 30 && mesh->BCp.type[cc] != 0  && mesh->BCp.type[cc] != 31 ) {
                    mesh->p_in[cc] += alpha*dx[Stokes->eqn_p[cc]];
                }
            }
            
            // Apply Bc to Vx and Vz
            ApplyBC( mesh, *model );
            
            //            Print2DArrayDouble( u,  nx, nzvx, scaling.V );
            //            Print2DArrayDouble( v,  nxvz, nz, scaling.V );
            //            Print2DArrayDouble( p,  ncx, ncz, scaling.S );
            
            //------------------------------------------------------------------------------------------------------
            // Update non-linearity
            UpdateNonLinearity( mesh, particles, topo_chain, topo, materials, model, Nmodel, scaling, 0, 1.0 );
            
            //------------------------------------------------------------------------------------------------------
            
            // Calculate residual
            //            model->free_surf_stab = 0;
            EvaluateStokesResidualDecoupled( Stokes, StokesA, StokesB, StokesC, StokesD, &residuals, mesh, *model, scaling, 1 );
            //            model->free_surf_stab = dummy;
            
            rx[kk] = residuals.resx;
            rz[kk] = residuals.resz;
            rp[kk] = residuals.resp;
            //            printf("\e[1;34mAlpha\e[m = %lf --> rx = %2.6e rz = %2.6e rp = %2.6e\n", alpha, rx[kk]* (scaling.F/pow(scaling.L,3)), rz[kk]* (scaling.F/pow(scaling.L,3)), rp[kk]*scaling.E );
            printf("\e[1;34mAlpha\e[m = %lf --> rx = %2.6e rz = %2.6e rp = %2.6e\n", alpha, rx[kk], rz[kk], rp[kk]);
            
            
        }
        
        //        // Look for the minimum predicted residuals
        //        minx  = rx[0];
        //        ix = 0;
        //        for( k=1; k<ntry[niter]; k++ ) {
        //            if(rx[k]<minx) {
        //                minx = rx[k];
        //                ix = k;
        //            }
        //        }
        //        alpha = alphav[ix];
        
        // Look for the minimum predicted residuals
        double r;
        minx  = sqrt( pow( rx[0],2 ) + pow( rz[0],2 ) );
        ix = 0;
        for( k=1; k<ntry[niter]; k++ ) {
            r = sqrt( pow( rx[k],2 ) + pow( rz[k],2 ) );
            if( r < minx ) {
                minx = r;
                ix = k;
            }
        }
        alpha = alphav[ix];
        
        // if the minmimun residuals are lower than starting ones, then success
        if ( rx[ix] < frac*Nmodel->resx_f || rz[ix]<frac*Nmodel->resz_f  ) { //|| rp[ix]<frac*Nmodel->resp
            success = 1;
            //            printf("\e[1;34mPredicted Residuals\e[m : alpha  = %lf --> rx = %2.4e rz = %2.4e rp = %2.4e\n", alphav[ix], rx[ix]*(scaling.F/pow(scaling.L,3)), rz[ix]* (scaling.F/pow(scaling.L,3)), rp[ix]*scaling.E);
            printf("\e[1;34mPredicted Residuals\e[m : alpha  = %lf --> rx = %2.4e rz = %2.4e rp = %2.4e\n", alphav[ix], rx[ix], rz[ix], rp[ix]);
            
        }
        
        //        if (success == 0 ) {
        //            // Look for the minimum predicted residuals
        //            minx  = rz[0];
        //            ix = 0;
        //            for( k=1; k<ntry[niter]; k++ ) {
        //                if(rz[k]<minx) {
        //                    minx = rz[k];
        //                    ix = k;
        //                }
        //            }
        //            alpha = alphav[ix];
        //        }
        
        DoodzFree(rx);
        DoodzFree(rz);
        DoodzFree(rp);
        DoodzFree(alphav);
        niter++;
    }
    
    if ( fabs(alpha)<1e-13 || success == 0 ) {
        printf( "Found minimum of the function -- cannot iterate further down\n" );
        Nmodel->stagnated = 1;
        alpha = 0.0;
    }
    
    // Reset solutions
    ArrayEqualArray( mesh->u_in, u, nx*nzvx );
    ArrayEqualArray( mesh->v_in, v, nxvz*nz );
    ArrayEqualArray( mesh->p_in, p, ncx*ncz );
    
    // Free
    DoodzFree( u );
    DoodzFree( v );
    DoodzFree( p );
    printf("** Line search took = %f sec\n", (double)((double)omp_get_wtime() - t_omp) );
    
    return alpha;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SolveStokes( SparseMat *Stokes, DirectSolver* PardisoStokes ) {
    
    clock_t t_omp;
    
    // Call direct solver
    t_omp = (double)omp_get_wtime();
    DirectStokes( Stokes, PardisoStokes, Stokes->b, Stokes->x );
    printf("** Time for direct Stokes solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SolveStokesDefect( SparseMat *Stokes, DirectSolver *PardisoStokes, Nparams* Nmodel, grid* mesh, params* model, markers* particles, markers* topo_chain, surface *topo, mat_prop materials, scale scaling ) {
    
    double alpha = -1.0;
    double *dx;
    clock_t t_omp;
    
    dx = DoodzCalloc( Stokes->neq, sizeof(double));
    
    // Call direct solver
    t_omp = (double)omp_get_wtime();
    DirectStokes( Stokes, PardisoStokes, Stokes->F, dx );
    printf("** Time for direct Stokes solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    // Line search
    if ( model->line_search == 1 ) {
        alpha = LineSearch( Stokes, dx, mesh, model, Nmodel, particles, topo_chain, topo, materials, scaling  );
    }
    
    // Return updated solution vector
    if ( Nmodel->stagnated == 0) {
        printf( "Correct x\n");
        ArrayPlusScalarArray( Stokes->x, alpha, dx, Stokes->neq );
    }
    
    DoodzFree(dx);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SolveStokesDecoupled( SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, SparseMat *Stokes, DirectSolver* PardisoStokes, params model, grid *mesh, scale scaling ) {
    
    clock_t t_omp;
    
    // Call direct solver
    t_omp = (double)omp_get_wtime();
    if (model.lsolver==0) DirectStokesDecoupled( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->b, StokesC->b, Stokes->x, model, mesh, scaling, Stokes );
    if (model.lsolver==1) KSPStokesDecoupled   ( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->b, StokesC->b, Stokes->x, model, mesh, scaling, Stokes, Stokes, StokesA, StokesB, StokesC );
    printf("** Time for direct decoupled Stokes solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SolveStokesDefectDecoupled( SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, SparseMat *Stokes, DirectSolver *PardisoStokes, Nparams* Nmodel, grid* mesh, params* model, markers* particles, markers* topo_chain, surface *topo, mat_prop materials, scale scaling, SparseMat *JacobA, SparseMat *JacobB, SparseMat *JacobC ) {
    
    double alpha = -1.0;
    double *dx, *fu, *fp;
    clock_t t_omp;
    
    dx = DoodzCalloc( Stokes->neq,      sizeof(double));
    fu = DoodzCalloc( Stokes->neq_mom,  sizeof(double));
    fp = DoodzCalloc( Stokes->neq_cont, sizeof(double));
    
    // Call direct solver
    t_omp = (double)omp_get_wtime();
    if ( model->lsolver == 0 ) DirectStokesDecoupled( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->F, StokesC->F, dx, *model, mesh, scaling, Stokes );
    if ( model->lsolver == 1 ) KSPStokesDecoupled   ( StokesA, StokesB, StokesC, StokesD, PardisoStokes, StokesA->F, StokesC->F, dx, *model, mesh, scaling, Stokes, Stokes, JacobA, JacobB, JacobC );
    printf("** Time for direct Stokes solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    // Line search
    if ( model->line_search == 1 ) {
        alpha = LineSearchDecoupled( Stokes, StokesA, StokesB, StokesC, StokesD, dx, mesh, model, Nmodel, particles, topo_chain, topo, materials, scaling  );
    }
    
    // Return updated solution vector
    if ( Nmodel->stagnated == 0) ArrayPlusScalarArray( Stokes->x, alpha, dx, Stokes->neq );
    
    //    Nparams residuals;
    //    ExtractSolutions( Stokes, mesh, *model );
    //
    //        UpdateNonLinearity( mesh, particles, topo_chain, topo, materials, model, Mmodel, Emodel, Nmodel, scaling, 0, 1.0 );
    //    EvaluateStokesResidual( Stokes, &residuals, mesh, *model, scaling, 0 );
    
    DoodzFree(dx);
    DoodzFree(fu);
    DoodzFree(fp);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void EvaluateStokesResidual( SparseMat *Stokes, Nparams *Nmodel, grid *mesh, params model, scale scaling, int quiet ) {
    
    int cc, nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1;
    double celvol, Area=0.0;
    int ndofx=0, ndofz=0, ndofp=0;
    
    celvol = model.dx*model.dz;
    
    // Residuals
    double resx = 0.0;
    double resz = 0.0;
    double resp = 0.0;
    
    // Function evaluation
    BuildStokesOperator( mesh, model, 0, mesh->p_in,  mesh->u_in,  mesh->v_in, Stokes, 0 );
    
    // Integrate residuals
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nx, nzvx ) reduction(+:resx,ndofx)
    for( cc=0; cc<nzvx*nx; cc++) {
        if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12) {
            ndofx++;
            resx += Stokes->F[Stokes->eqn_u[cc]]*Stokes->F[Stokes->eqn_u[cc]];//*celvol;
            mesh->ru[cc] = Stokes->F[Stokes->eqn_u[cc]];
        }
    }
    Nmodel->resx = resx;
    
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nz, nxvz ) reduction(+:resz,ndofz)
    for( cc=0; cc<nz*nxvz; cc++) {
        if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12) {
            ndofz++;
            resz += Stokes->F[Stokes->eqn_v[cc]]*Stokes->F[Stokes->eqn_v[cc]];//*celvol;
            mesh->rv[cc] = Stokes->F[Stokes->eqn_v[cc]];
        }
    }
    Nmodel->resz = resz;
    
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, ncz, ncx ) reduction(+:resp,ndofp,Area)
    for( cc=0; cc<ncz*ncx; cc++) {
        if ( mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30  && mesh->BCp.type[cc] != 31) {
            ndofp++;
            Area += celvol;
            resp += Stokes->F[Stokes->eqn_p[cc]]*Stokes->F[Stokes->eqn_p[cc]];//*celvol;
            mesh->rp[cc] = Stokes->F[Stokes->eqn_p[cc]];
        }
    }
    Nmodel->resp = resp;
    
    // Sqrt
    //    Nmodel->resx =  sqrt(Nmodel->resx)/Area/(ndofx);
    //    Nmodel->resz =  sqrt(Nmodel->resz)/Area/(ndofz);
    //    Nmodel->resp =  sqrt(Nmodel->resp)/Area/(ndofp);
    Nmodel->resx =  sqrt(Nmodel->resx/ndofx);
    Nmodel->resz =  sqrt(Nmodel->resz/ndofz);
    Nmodel->resp =  sqrt(Nmodel->resp/ndofp);
    
    
    if ( quiet == 0 ) {
        //        printf("Fu = %2.6e\n", Nmodel->resx * (scaling.F/pow(scaling.L,3))); // Units of momentum
        //        printf("Fv = %2.6e\n", Nmodel->resz * (scaling.F/pow(scaling.L,3))); // Units of momentum
        //        printf("Fp = %2.6e\n", Nmodel->resp * (scaling.E)); // Units of velocity gradient
        //        printf("Unscaled residuals!\n");
        printf("Fu = %2.6e\n", Nmodel->resx ); // Units of momentum
        printf("Fv = %2.6e\n", Nmodel->resz ); // Units of momentum
        printf("Fp = %2.6e\n", Nmodel->resp ); // Units of velocity gradient
    }
    
    if ( isnan(Nmodel->resx) || isnan(Nmodel->resz) || isnan(Nmodel->resp) ) {
        printf("Solve went wrong - Nan residuals...\nExiting...\n");
        exit(122);
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void EvaluateStokesResidualDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, Nparams *Nmodel, grid *mesh, params model, scale scaling, int quiet ) {
    
    int cc, nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1;
    double celvol, Area=0.0;
    int ndofx=0, ndofz=0, ndofp=0;
    
    celvol = model.dx*model.dz;
    
    //    printf("include stab %d\n", model.free_surf_stab);
    
    
    // Residuals
    double resx = 0.0;
    double resz = 0.0;
    double resp = 0.0;
    
    // Function evaluation
    BuildStokesOperatorDecoupled( mesh, model, 0, mesh->p_in,  mesh->u_in,  mesh->v_in, Stokes, StokesA, StokesB, StokesC, StokesD, 0 );
    
    // Integrate residuals
#pragma omp parallel for shared( mesh, Stokes, StokesA ) private( cc ) firstprivate( celvol, nx, nzvx ) reduction(+:resx,ndofx)
    for( cc=0; cc<nzvx*nx; cc++) {
        if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12) {
            ndofx++;
            resx += StokesA->F[Stokes->eqn_u[cc]]*StokesA->F[Stokes->eqn_u[cc]];//*celvol;
            mesh->ru[cc] = StokesA->F[Stokes->eqn_u[cc]];
        }
    }
    Nmodel->resx = resx;
    
#pragma omp parallel for shared( mesh, Stokes, StokesA ) private( cc ) firstprivate( celvol, nz, nxvz ) reduction(+:resz,ndofz)
    for( cc=0; cc<nz*nxvz; cc++) {
        if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12) {
            ndofz++;
            resz += StokesA->F[Stokes->eqn_v[cc]]*StokesA->F[Stokes->eqn_v[cc]];//*celvol;
            mesh->rv[cc] = StokesA->F[Stokes->eqn_v[cc]];
        }
    }
    Nmodel->resz = resz;
    
#pragma omp parallel for shared( mesh, Stokes, StokesD ) private( cc ) firstprivate( celvol, ncz, ncx ) reduction(+:resp,ndofp,Area)
    for( cc=0; cc<ncz*ncx; cc++) {
        if ( mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30  && mesh->BCp.type[cc] != 31) {
            ndofp++;
            Area += celvol;
            resp += StokesC->F[Stokes->eqn_p[cc]- Stokes->neq_mom]*StokesC->F[Stokes->eqn_p[cc]- Stokes->neq_mom];//*celvol;
            mesh->rp[cc] = StokesC->F[Stokes->eqn_p[cc]- Stokes->neq_mom];
        }
    }
    Nmodel->resp = resp;
    
    // Sqrt
    //    Nmodel->resx =  sqrt(Nmodel->resx)/Area/(ndofx);
    //    Nmodel->resz =  sqrt(Nmodel->resz)/Area/(ndofz);
    //    Nmodel->resp =  sqrt(Nmodel->resp)/Area/(ndofp);
    Nmodel->resx =  sqrt(Nmodel->resx/ndofx);
    Nmodel->resz =  sqrt(Nmodel->resz/ndofz);
    Nmodel->resp =  sqrt(Nmodel->resp/ndofp);
    
    int k,l;
    
    double *fu, *fv;
    fu  = DoodzCalloc( sizeof(double),mesh->Nx*(mesh->Nz+1));
    fv  = DoodzCalloc( sizeof(double),mesh->Nz*(mesh->Nx+1));
    
    for( l=0; l<nzvx; l++) {
        for( k=0; k<nx; k++) {
            cc = k + l*nx;
            fu[cc] = 0.0;
            if (  mesh->BCu.type[cc] == -1  || mesh->BCu.type[cc] == 2) {
                fu[cc] = StokesA->F[Stokes->eqn_u[cc]];
            }
        }
    }
    
    for( l=0; l<nz; l++) {
        for( k=0; k<nxvz; k++) {
            cc = k + l*nxvz;
            fv[cc] = 0.0;
            if (  mesh->BCv.type[cc] == -1  || mesh->BCv.type[cc] == 2) {
                fv[cc] = StokesA->F[Stokes->eqn_v[cc]];
            }
        }
    }
    
    //        Print2DArrayDouble( fu, mesh->Nx, mesh->Nz+1, 1 ); //(scaling.F/pow(scaling.L,3))
    //        Print2DArrayDouble( fv, mesh->Nx+1, mesh->Nz, (scaling.F/pow(scaling.L,3)) );
    //        Print2DArrayDouble( StokesC->F, mesh->Nx-1, mesh->Nz-1, (scaling.E) );
    
    
    if ( quiet == 0 ) {
        //        printf("Fu = %2.6e\n", Nmodel->resx * (scaling.F/pow(scaling.L,3))); // Units of momentum
        //        printf("Fv = %2.6e\n", Nmodel->resz * (scaling.F/pow(scaling.L,3))); // Units of momentum
        //        printf("Fp = %2.6e\n", Nmodel->resp * (scaling.E)); // Units of velocity gradient
        //        printf("Unscaled residuals! \n");
        printf("Fu = %2.6e\n", Nmodel->resx ); // Units of momentum
        printf("Fv = %2.6e\n", Nmodel->resz ); // Units of momentum
        printf("Fp = %2.6e\n", Nmodel->resp ); // Units of velocity gradient
    }
    
    if ( isnan(Nmodel->resx) || isnan(Nmodel->resz) || isnan(Nmodel->resp) ) {
        printf("Solve went wrong - Nan residuals...\nExiting...\n");
        exit(122);
    }
    
    DoodzFree(fu);
    DoodzFree(fv);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void EvaluateRHS( grid* mesh, params model, scale scaling, double RHO_REF ) {
    
    int c, k, l, c1, c2;
    int NX, NZ, NZVX, NXVZ;
    int NCX, NCZ;
    double rhoVx, rhoVz;
    double dx, dz;
    int comp = 0;
    
    NX   = mesh->Nx;
    NXVZ = mesh->Nx+1;
    NZ   = mesh->Nz;
    NZVX = mesh->Nz+1;
    NCX=NX-1; NCZ=NZ-1;
    dx   = mesh->dx;
    dz   = mesh->dz;
    
//    double eii;
//    int    pl=0;
//
//    if (model.loc_iter==2) pl=-1;
//
//    // Evaluate plastic strain rate
//    double *exx_pl, *exz_pl;
//    exx_pl = DoodzCalloc(NCX*NCZ, sizeof(double));
//    exz_pl = DoodzCalloc(NX *NZ , sizeof(double));
//
//    double *eiin, *eiis;
//    eiin = DoodzCalloc(NCX*NCZ, sizeof(double));
//    eiis = DoodzCalloc(NX *NZ , sizeof(double));
//
//    //#pragma omp parallel for shared( mesh, exx_pl, eiin ) private( eii, c ) firstprivate( NCX, NCZ )
//    for (c=0; c<NCX*NCZ; c++) {
//        eii       = sqrt(mesh->exxd[c]*mesh->exxd[c]  + mesh->exz_n[c]*mesh->exz_n[c]);
//        exx_pl[c] = mesh->eII_pl[c]*mesh->exxd[c]/eii;
//        exx_pl[c] = mesh->exx_pl[c];
//    }
//
//    //#pragma omp parallel for shared( mesh, exz_pl, eiis ) private( eii, c ) firstprivate( NX, NZ )
//    for (c=0; c<NX*NZ; c++) {
//        eii       = sqrt(mesh->exxd_s[c]*mesh->exxd_s[c] + mesh->exz[c]*mesh->exz[c]);
//        exz_pl[c] = mesh->eII_pl_s[c]*mesh->exz[c]/eii;
//        exz_pl[c] = mesh->exz_pl[c];
//    }
//    DoodzFree( eiin );
//    DoodzFree( eiis );
    
    /* --------------------------------------------------------------------*/
    /* Here we calculate the forcing term -rho*gx on the finest grid level */
    /* --------------------------------------------------------------------*/
    
    // U POINTS
    for (l=0; l<NZVX; l++) {
        for (k=0; k<NX; k++) {
            
            c  = k + l*NX;
            c2 = k + (l-1)*NCX;
            
            mesh->roger_x[c] = 0.0;
            
            if (l>0 && l<NZVX-1) {
                
                if ( mesh->BCu.type[c] == -1 || mesh->BCu.type[c] == 2 || mesh->BCu.type[c] == -2 ) {
                    
                    // Gravity force: take apparent viscosity (free surface correction)
                    rhoVx = 0.5*(mesh->rho_app_s[c] + mesh->rho_app_s[c-NX]);
                    mesh->roger_x[c]  = - model.gx * rhoVx;
                    
                    // Elastic force
                    if ( model.iselastic == 1 ) {
                        
                        // Elasticity in periodic
                        if (   l>0 && l<NZVX-1 && k==0 && k==0 && mesh->BCu.type[c] == -2 ) {
                            mesh->roger_x[c]  -= 1.0/dx * ( mesh->VE_n[c2] * mesh->sxxd0[c2] - mesh->VE_n[c2+NCX-1]  * mesh->sxxd0[c2+NCX-1] );
                            mesh->roger_x[c]  -= 1.0/dz * ( mesh->VE_s[c]  * mesh->sxz0[c]   - mesh->VE_s[c-NX]  * mesh->sxz0[c-NX]  );
                        }
                        
                        // Inner nodes
                        if (  l>0 && l<NZVX-1 && k>0 && k<NX-1 && model.subgrid_diff != 4 ) {
                            // Elasticity
                            mesh->roger_x[c]  -= 1.0/dx * ( mesh->VE_n[c2] * mesh->sxxd0[c2] - mesh->VE_n[c2-1]  * mesh->sxxd0[c2-1] );
                            mesh->roger_x[c]  -= 1.0/dz * ( mesh->VE_s[c]  * mesh->sxz0[c]   - mesh->VE_s[c-NX]  * mesh->sxz0[c-NX]  );
//                            // Plasticity
//                            mesh->roger_x[c]  -= pl*2.0/dx * ( yE*mesh->eta_n[c2] * exx_pl[c2] - yW*mesh->eta_n[c2-1]  * exx_pl[c2-1] );
//                            mesh->roger_x[c]  -= pl*2.0/dz * ( yN*mesh->eta_s[c]  * exz_pl[c]   - yS*mesh->eta_s[c-NX]  * exz_pl[c-NX]  );
                        }
                        
                        // WEST OPEN
                        if (  l>0 && l<NZVX-1 && k==0 && mesh->BCu.type[c] == 2) {
                            // Elasticity
                            mesh->roger_x[c]  -= 1.0/dx * ( mesh->VE_n[c2] * mesh->sxxd0[c2] );
                            mesh->roger_x[c]  -= 1.0/dz * ( mesh->VE_s[c]  * mesh->sxz0[c]   - mesh->VE_s[c-NX]  * mesh->sxz0[c-NX]  );
//                            // Plasticity
//                            mesh->roger_x[c]  -= pl*2.0/dx * ( yW*mesh->eta_n[c2] * exx_pl[c2] );
//                            mesh->roger_x[c]  -= pl*2.0/dz * ( yN*mesh->eta_s[c]  * exz_pl[c]   - yS*mesh->eta_s[c-NX]  * exz_pl[c-NX]  );
                        }
                        
                        // EAST OPEN
                        if (  l>0 && l<NZVX-1 && k==NX-1 && mesh->BCu.type[c] == 2) {
                            // Elasticity
                            mesh->roger_x[c]  -= 1.0/dx * ( - mesh->VE_n[c2-1]  * mesh->sxxd0[c2-1] );
                            mesh->roger_x[c]  -= 1.0/dz * ( mesh->VE_s[c]  * mesh->sxz0[c]   - mesh->VE_s[c-NX]  * mesh->sxz0[c-NX]  );
//                            // Plasticity
//                            mesh->roger_x[c]  -= pl*2.0/dx * ( - yW*mesh->eta_n[c2-1]  * exx_pl[c2-1] );
//                            mesh->roger_x[c]  -= pl*2.0/dz * ( yN*mesh->eta_s[c]  * exz_pl[c]   - yS*mesh->eta_s[c-NX]  * exz_pl[c-NX]  );
                        }
                        
                    }
                }
            }
            mesh->roger_x[c] = -mesh->roger_x[c];
        }
    }
    
    
    /* --------------------------------------------------------------------*/
    /* Here we calculate the forcing term -rho*gz on the finest grid level */
    /* --------------------------------------------------------------------*/
    
    // V POINTS
    for (l=0; l<NZ; l++) {
        for (k=0; k<NXVZ; k++) {
            
            c  = k + l*(NXVZ);
            c1 = k + l*(NX)-1;
            c2 = k-1 + (l-1)*(NX-1);
            
            mesh->roger_z[c]  = 0.0;
            
            if (k>0 && k<NXVZ-1) {
                
                if ( mesh->BCv.type[c] == -1  ) {
                    
                    // Gravity force: use apparent density (free surface correction)
                    rhoVz = 0.5 * (mesh->rho_app_s[c1] + mesh->rho_app_s[c1+1]);  // USE THIS ALWAYS
//                    rhoVz = 0.5 * (mesh->rho_app_n[c2] + mesh->rho_app_n[c2+(NX-1)]); // DO NOT USE
                    
                    // Remove background pressure
                    mesh->roger_z[c]  = - model.gz * rhoVz;
                    
                    // Elastic force
                    if  (model.iselastic == 1 ) {
                        
                        // Backward Euler
                        if ( l>0 && l<NZ-1 && k>0 && k<NXVZ-1 && model.subgrid_diff!=4  ) {
                            // Elasticity
                            mesh->roger_z[c]  -= 1.0/dz * ( mesh->VE_n[c2+(NCX)] * -(mesh->sxxd0[c2+(NCX)]) - mesh->VE_n[c2] * -(mesh->sxxd0[c2]) );
                            mesh->roger_z[c]  -= 1.0/dx * ( mesh->VE_s[c1+1]     *   mesh->sxz0[c1+1]      - mesh->VE_s[c1] *   mesh->sxz0[c1]  );
//                            // Plasticity
//                            mesh->roger_z[c]  -= pl*2.0/dz * ( yN*mesh->eta_n[c2+(NCX)] * -(exx_pl[c2+(NCX)]) - yS*mesh->eta_n[c2] * -(exx_pl[c2]) );
//                            mesh->roger_z[c]  -= pl*2.0/dx * ( yE*mesh->eta_s[c1+1]     *   exz_pl[c1+1]      - yW*mesh->eta_s[c1] *   exz_pl[c1]  );
                        }
                        // SOUTH OPEN
                        if ( l==0 && k>0 && k<NXVZ-1 ) {
                            // Elasticity
                            mesh->roger_z[c]  -= 1.0/dz * ( mesh->VE_n[c2+(NCX)] * -(mesh->sxxd0[c2+(NCX)]) );
                            mesh->roger_z[c]  -= 1.0/dx * ( mesh->VE_s[c1+1]     *   mesh->sxz0[c1+1]      - mesh->VE_s[c1] *   mesh->sxz0[c1]  );
//                            // Plasticity
//                            mesh->roger_z[c]  -= pl*2.0/dz * (yN* mesh->eta_n[c2+(NCX)] * -(exx_pl[c2+(NCX)]) );
//                            mesh->roger_z[c]  -= pl*2.0/dx * ( yE*mesh->eta_s[c1+1]     *   exz_pl[c1+1]      - yW*mesh->eta_s[c1] *   exz_pl[c1]  );
                            
                        }
                        
                        // NORTH OPEN
                        if ( l==NZ-1 && k>0 && k<NXVZ-1 ) {
                            // Elasticity
                            mesh->roger_z[c]  -= 1.0/dz * ( - mesh->VE_n[c2] * -(mesh->sxxd0[c2]) );
                            mesh->roger_z[c]  -= 1.0/dx * ( mesh->VE_s[c1+1]     *   mesh->sxz0[c1+1]      - mesh->VE_s[c1] *   mesh->sxz0[c1]  );
//                            // Plasticity
//                            mesh->roger_z[c]  -= pl*2.0/dz * ( - yS*mesh->eta_n[c2] * -(exx_pl[c2]) );
//                            mesh->roger_z[c]  -= pl*2.0/dx * ( yE*mesh->eta_s[c1+1]     *   exz_pl[c1+1]      - yW*mesh->eta_s[c1] *   exz_pl[c1]  );
                            
                        }
                    }
                }
            }
            mesh->roger_z[c] = -mesh->roger_z[c];
        }
    }
    
    /* ------------------------------------------------------------------------*/
    /* Here we calculate the forcing term of the continuity and heat equation  */
    /* ------------------------------------------------------------------------*/
    
    // P-T RHS
    for (l=0; l<NZ-1; l++) {
        for (k=0; k<NX-1; k++) {
            
            c  = k + l*(NX-1);
            
            if ( mesh->BCp.type[c] != 30 && mesh->BCp.type[c] != 31 ) {
                
                
                // Heat equation
                mesh->rhs_t[c] = mesh->T[c];
                mesh->rhs_p[c] = 0.0;
                
                // Continuity equation
                mesh->rhs_p[c] = 0.0;
                
            }
        }
    }
    
//    DoodzFree( exx_pl );
//    DoodzFree( exz_pl );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void EvalNumberOfEquations( grid* mesh, SparseMat *Stokes ) {
    int inc = 0, incP=0;
    int k,l,kk;
    
    // Allocate
    Stokes->eqn_u = DoodzMalloc((mesh->Nz+1)*mesh->Nx * sizeof(int));
    Stokes->eqn_v = DoodzMalloc(mesh->Nz*(mesh->Nx+1) * sizeof(int));
    Stokes->eqn_p = DoodzMalloc((mesh->Nz-1)*(mesh->Nx-1) * sizeof(int));
    
    // Number U dofs
    for( l=0; l<mesh->Nz+1; l++) {
        for( k=0; k<mesh->Nx; k++) {
            kk = k + l*mesh->Nx;
            Stokes->eqn_u[kk] = -1;
            if ( mesh->BCu.type[kk] ==-1 || mesh->BCu.type[kk] ==2 || mesh->BCu.type[kk] ==-2  ) {
                Stokes->eqn_u[kk] = inc;
                inc++;
            }
            // Copy equation numbers from the West side
            if ( mesh->BCu.type[kk] ==-12) {
                Stokes->eqn_u[kk] = Stokes->eqn_u[kk-mesh->Nx+1];
            }
        }
    }
    
    // Number V dofs
    for( l=0; l<mesh->Nz; l++) {
        for( k=0; k<mesh->Nx+1; k++) {
            kk = k + l*(mesh->Nx+1);
            Stokes->eqn_v[kk] = -1;
            if ( mesh->BCv.type[kk] ==-1 || mesh->BCv.type[kk] ==2 ) {
                Stokes->eqn_v[kk] = inc;
                inc++;
            }
        }
    }
    
    Stokes->neq_mom = inc;
    
    // Number P dofs
    for( l=0; l<mesh->Nz-1; l++) {
        for( k=0; k<mesh->Nx-1; k++) {
            kk = k + l*(mesh->Nx-1);
            Stokes->eqn_p[kk] = -1;
            if ( mesh->BCp.type[kk] ==-1 ) {
                Stokes->eqn_p[kk] = inc;
                inc++;
                incP++;
            }
        }
    }
    Stokes->neq = inc;
    Stokes->neq_cont = incP;
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InitialiseSolutionVector( grid* mesh, SparseMat *Stokes, params* model ) {
    int k,l,kk;
    
    for( l=0; l<mesh->Nz+1; l++) {
        for( k=0; k<mesh->Nx; k++) {
            kk = k + l*mesh->Nx;
            if ( mesh->BCu.type[kk] ==-1 || mesh->BCu.type[kk] ==2 || mesh->BCu.type[kk] ==-2  ) {
                Stokes->x[ Stokes->eqn_u[kk] ] = mesh->u_in[kk];
            }
            
        }
    }
    
    // Number V dofs
    for( l=0; l<mesh->Nz; l++) {
        for( k=0; k<mesh->Nx+1; k++) {
            kk = k + l*(mesh->Nx+1);
            if ( mesh->BCv.type[kk] ==-1 || mesh->BCv.type[kk] ==2 ) {
                Stokes->x[ Stokes->eqn_v[kk] ] = mesh->v_in[kk];
            }
        }
    }
    
    // Number P dofs
    for( l=0; l<mesh->Nz-1; l++) {
        for( k=0; k<mesh->Nx-1; k++) {
            kk = k + l*(mesh->Nx-1);
            if ( mesh->BCp.type[kk] ==-1 ) {
                Stokes->x[ Stokes->eqn_p[kk] ] = mesh->p_in[kk];
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DirectStokes( SparseMat *mat, DirectSolver *pardi, double *rhs, double *sol ) {
    
    double Info [UMFPACK_INFO], Control [UMFPACK_CONTROL];
    int status;
    void *Symbolic = NULL;
    void *Numeric = NULL;
    int msglvl = 0;
    
    clock_t t_omp;
    
    t_omp = (double)omp_get_wtime();
    
    if (msglvl==1) printf ("\nUMFPACK V%d.%d (%s) demo: _di_ version\n",  UMFPACK_MAIN_VERSION, UMFPACK_SUB_VERSION, UMFPACK_DATE) ;
    
    /* get the default control parameters */
    umfpack_di_defaults (Control) ;
    
    /* change the default print level for this demo */
    /* (otherwise, nothing will print) */
    if (msglvl==1) Control [UMFPACK_PRL] = 6 ;
    
    /* print the license agreement */
    if (msglvl==1) umfpack_di_report_status (Control, UMFPACK_OK) ;
    if (msglvl==1) Control [UMFPACK_PRL] = 5 ;
    
    Control [UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
    //    Control [UMFPACK_PIVOT_TOLERANCE] = 1;
    Control [UMFPACK_ORDERING] = 1; // not 2
    Control [UMFPACK_AGGRESSIVE] = 1;
    Control [UMFPACK_IRSTEP]   = 0;
    Control [UMFPACK_FIXQ] = 1;
    
    /* print the control parameters */
    if (msglvl==1) umfpack_di_report_control (Control) ;
    
    int n = mat->neq;
    
    /* ---------------------------------------------------------------------- */
    /* symbolic factorization */
    /* ---------------------------------------------------------------------- */
    
    status = umfpack_di_symbolic (n, n, mat->Ic, mat->J, mat->A, &Symbolic,
                                  Control, Info) ;
    if (status < 0)
    {
        umfpack_di_report_info (Control, Info) ;
        umfpack_di_report_status (Control, status) ;
        error ("umfpack_di_symbolic failed") ;
    }
    
    /* ---------------------------------------------------------------------- */
    /* numeric factorization */
    /* ---------------------------------------------------------------------- */
    
    status = umfpack_di_numeric (mat->Ic, mat->J, mat->A, Symbolic, &Numeric,
                                 Control, Info) ;
    if (status < 0)
    {
        umfpack_di_report_info (Control, Info) ;
        umfpack_di_report_status (Control, status) ;
        error ("umfpack_di_numeric failed") ;
    }
    umfpack_di_free_symbolic(&Symbolic);
    
    /* ---------------------------------------------------------------------- */
    /* solve Ax=b */
    /* ---------------------------------------------------------------------- */
    
    status = umfpack_di_solve (UMFPACK_Aat, mat->Ic, mat->J, mat->A, sol, rhs,
                               Numeric, Control, Info) ;
    if (msglvl==1)umfpack_di_report_info (Control, Info) ;
    umfpack_di_report_status (Control, status) ;
    if (status < 0)
    {
        error ("umfpack_di_solve failed") ;
    }
    umfpack_di_free_numeric(&Numeric);
    
    printf("** Time for UMFPACK = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DirectStokesDecoupled( SparseMat *matA,  SparseMat *matB,  SparseMat *matC,  SparseMat *matD, DirectSolver *pardi, double *rhs_mom, double *rhs_cont, double *sol, params model, grid *mesh, scale scaling, SparseMat *Stokes ) {
    
    double  celvol=model.dx*model.dz ;
    cs_di  A, B, D, C, *B1, *L, *Ac, *Bc,  *Cc,  *Dc, *L1;
    DoodzFP  *u0, *p0;
    int  noisy=1;
    int nitmax=20, k;
    double mone[2] = {-1.0,0.0}, one[2] = {1.0,0.0}, zero[2] = {0.0,0.0};
    DoodzFP ru, rp;
    double maxdiv0, mindiv, maxdiv, maxdivit=0, rel_tol_div=model.rel_tol_div;
    double minru0, maxru0, minru, maxru;
    
    cholmod_common c ;
    cholmod_sparse *Lcm, *Lcml, *Acm, *Bcm, *Ccm, *Dcm, *Dcm0;
    cholmod_factor *Lfact ;
    
    //    cholmod_start( &c ) ;
    //    pardi->Analyze = 1;
    
    Lfact = pardi->Lfact;
    c     = pardi->c;
    
    // --------------- Pre-process--------------- //
    
    // ************** D ************** //
    SuiteSparse_long rsize;
    double gamma  = model.penalty, penalty, min_eta, max_eta, min_G, max_G;
    rsize = Stokes->neq_cont;
    Dcm0  = cholmod_speye (rsize, rsize, CHOLMOD_REAL, &c );
    
    MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_eta, &max_eta );
    MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_G  , &max_G   );
    penalty = -min_G*model.dt *model.dt * celvol;
    //    printf("min eta = %2.2e, max eta = %2.2e, penalty try=%2.2e, penalty =%2.2e, celvol=%2.2e...\n", min_eta, max_eta, penalty, -gamma * celvol, celvol);
    
    if (model.auto_penalty == 1) {
        printf("Trying AUTO_PENALTY...\n");
        MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_eta, &max_eta );
        MinMaxArrayVal( mesh->eta_s, mesh->Nx*mesh->Nz, &min_eta, &max_eta );
        
        //        penalty = -min_eta*model.dt * celvol;
        //        penalty = -1/(min_G*model.dt *model.dt);
        //        penalty = -min_eta/model.dt / celvol;
        //         printf("min eta = %2.2e, max eta = %2.2e, penalty try=%2.2e, penalty try=%2.2e...\n", min_eta, max_eta, penalty, -gamma * celvol);
    }
    else {
        penalty = -gamma * celvol;
        printf("Penalty factor = %2.2e\n", penalty);
    }
    
    for (k=0;k<Dcm0->nzmax;k++) ((double*)Dcm0->x)[k] *= penalty; // Should be /celvol
    
    clock_t t_omp;
    t_omp = (double)omp_get_wtime();
    
    // Build initial solution vector
    u0 = DoodzCalloc( matA->neq, sizeof(double) );
    p0 = DoodzCalloc( matC->neq, sizeof(double) );
    BuildInitialSolutions( u0, p0, mesh );
    
    // Prepare A
    A.nzmax = matA->nnz;
    A.nz    = matA->nnz;
    A.m     = matA->neq;
    A.n     = A.m;
    A.p     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.i     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.x     = DoodzCalloc( A.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( A.m, matA->Ic, A.i );
    ArrayEqualArrayI( A.p, matA->J,  A.nzmax );
    ArrayEqualArray(  A.x, matA->A,  A.nzmax );
    Ac  = cs_di_compress( &A );
    cs_droptol( Ac, 1.0e-18 );
    
    // Prepare B
    B.nzmax = matB->nnz;
    B.nz    = matB->nnz;
    B.m     = matB->neq;
    B.n     = matC->neq;
    B.p     = DoodzCalloc( B.nzmax, sizeof(int) );
    B.i     = DoodzCalloc( B.nzmax, sizeof(int) );
    B.x     = DoodzCalloc( B.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( B.m, matB->Ic, B.i );
    ArrayEqualArrayI( B.p, matB->J,  B.nzmax );
    ArrayEqualArray(  B.x, matB->A,  B.nzmax );
    Bc  = cs_di_compress( &B );
    
    // Prepare C
    C.nzmax = matC->nnz;
    C.nz    = matC->nnz;
    C.m     = matC->neq;
    C.n     = matB->neq;
    C.p     = DoodzCalloc( C.nzmax, sizeof(int) );
    C.i     = DoodzCalloc( C.nzmax, sizeof(int) );
    C.x     = DoodzCalloc( C.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( C.m, matC->Ic, C.i );
    ArrayEqualArrayI( C.p, matC->J,  C.nzmax );
    ArrayEqualArray(  C.x, matC->A,  C.nzmax );
    Cc  = cs_di_compress( &C );
    
    // Prepare D
    D.nzmax = Stokes->neq_cont;
    D.nz    = Stokes->neq_cont;
    D.m     = Stokes->neq_cont;
    D.n     = Stokes->neq_cont;
    D.p     = DoodzCalloc( Stokes->neq_cont+1, sizeof(int) );
    D.i     = DoodzCalloc( Stokes->neq_cont, sizeof(int) );
    D.x     = DoodzCalloc( Stokes->neq_cont, sizeof(double) );
    copy_cholmod_to_cs_matrix( Dcm0, &D );
    cholmod_free_sparse( &Dcm0, &c );
    Dc  = cs_di_compress( &D );
    
    // --------------- Schur complement --------------- //
    
    // Matrix multiplication: D*C
    L =  cs_di_multiply( Dc, Cc );
    
    //----- test - in case C' != B (assume B is deficient)
    B1 = cs_di_transpose( Cc, 1);
    //-----
    
    // Matrix multiplication: B*(D*C)
    L1 = cs_di_multiply( B1, L);
    cs_spfree(L);
    
    // Matrix addition: L = A - B*(D*C)
    L = cs_di_add( Ac, L1, 1, -1);
    cs_spfree(L1);
    
    // --------------- Cholesky --------------- //
    
    // LL' Cholesky
    //    cholmod_common c ;
    //    cholmod_sparse *Lcm, *Lcml, *Acm, *Bcm, *Ccm, *Dcm, *A1;
    //    cholmod_factor *Lfact ;
    //    cholmod_start( &c ) ;
    //    int *perm;
    //    c.supernodal = 2;
    Lcm = cholmod_allocate_sparse (L->m, L->n, L->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Lcm, L );
    Acm = cholmod_allocate_sparse (Ac->m, Ac->n, Ac->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Acm, Ac );
    //    Bcm = cholmod_allocate_sparse (Bc->m, Bc->n, Bc->nzmax, 0, 1, 0, 1, &c) ;
    //    copy_cs_to_cholmod_matrix( Bcm, Bc );
    Bcm = cholmod_allocate_sparse (B1->m, B1->n, B1->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Bcm, B1 );
    Ccm = cholmod_allocate_sparse (Cc->m, Cc->n, Cc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Ccm, Cc );
    Dcm = cholmod_allocate_sparse (Dc->m, Dc->n, Dc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Dcm, Dc );
    
    //    int ccc = cholmod_print_sparse( Lcm, "lcm", &c);
    //    printf("ccc=%d\n",ccc);
    
    // Keep lower part
    Lcml = cholmod_copy ( Lcm, -1, 1, &c );
    
    
    if (Lcml == NULL || Lcml->stype == 0)
    {
        printf("Unsymmetric matrix\n");
        cholmod_free_sparse (&Lcml, &c) ;
        cholmod_finish (&c) ;
    }
    
    //    c.nmethods=0;
    ////    printf("permuted\n");
    ////    c.method[0].ordering=CHOLMOD_GIVEN;
    ////    c.postorder = 1;
    //    perm     = DoodzCalloc( Lcml->nrow, sizeof(int) );
    //    cholmod_amd(Lcml, NULL, 0, perm, &c );
    //    Lfact = cholmod_analyze_p( Lcml, perm, NULL, 0, &c ) ;
    //    double beta[2] = {0,0};
    //    cholmod_factorize_p( Lcml, beta, NULL, 0, Lfact, &c) ;
    
    //    c.nmethods=1;
    //    c.method[0].ordering=CHOLMOD_AMD;
    //    c.postorder = 1;
    //
    //    Lfact = cholmod_analyze( Lcml, &c ) ;
    //    cholmod_factorize( Lcml, Lfact, &c) ;
    //
    //    if ( Lfact->is_super == 0 ) cholmod_change_factor( 1, 1, Lfact->is_super, 0, 0, Lfact, &c );
    //    cholmod_print_factor(Lfact, "L", &c);
    
    if (pardi->Analyze == 1) {
        c.nmethods           = 1;
        c.method[0].ordering = CHOLMOD_AMD;
        c.postorder          = 1;
        t_omp                = (double)omp_get_wtime();
        Lfact                = cholmod_analyze( Lcml, &c ) ;
        pardi->Analyze       = 0;
        printf("** Time for Cholesky analysis = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
        
    }
    
    t_omp = (double)omp_get_wtime();
    cholmod_factorize( Lcml, Lfact, &c);
    printf("** Time for Cholesky factorization = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    //    exit(69);
    
    cholmod_dense *u, *p, *bu, *bp, *pdum, *udum, *fu, *fp;
    pdum = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    udum = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    bu   = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    bp   = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    u    = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    p    = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    fu   = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    fp   = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    
    copy_vec_to_cholmod_dense( u, u0 );
    copy_vec_to_cholmod_dense( p, p0 );
    copy_vec_to_cholmod_dense( bu, rhs_mom );
    copy_vec_to_cholmod_dense( bp, rhs_cont );
    
    printf("Initial residual:\n");
    copy_cholmod_dense_to_cholmod_dense( fu, bu );   // fu = bu
    cholmod_sdmult ( Acm, 0, mone, one, u, fu, &c) ; // fu -= A*u
    cholmod_sdmult ( Bcm, 0, mone, one, p, fu, &c) ; // fu -= B*p
    
    
    //    SumArray(fp->x, 1.0, matC->neq, "fp");
    //    SumArray(bp->x, 1.0, matC->neq, "bp");
    copy_cholmod_dense_to_cholmod_dense( fp, bp );   // fp = bp
    //    SumArray(fp->x, 1.0, matC->neq, "fp");
    cholmod_sdmult ( Ccm, 0, mone, one, u, fp, &c) ; // fp -= C*u
    //    cholmod_sdmult ( Dcm, 0, mone, one, p, fp, &c) ;
    
    //    SumArray(fu->x, 1.0, matA->neq, "fu");
    //    SumArray(fp->x, 1.0, matC->neq, "fp");
    MinMaxArrayVal( fu->x, matA->neq, &minru0, &maxru0 );
    NormResidualCholmod( &ru, &rp, fu, fp, matA->neq, matC->neq, model, scaling, 0 );
    
    MinMaxArray( bu->x, 1, matA->neq, "u ini");
    MinMaxArray( bp->x, 1, matC->neq, "p ini");
    
    for ( k=0; k<nitmax; k++) {
        
        cholmod_sdmult ( Dcm, 0, one, zero, bp, pdum, &c) ;   // pdum <-- Dcm * bp
        copy_cholmod_dense_to_cholmod_dense( udum, bu );      // udum <-- bu
        cholmod_sdmult ( Bcm, 0, mone, one, p, udum, &c) ;    // udum <-- bu - B*p
        cholmod_sdmult ( Bcm, 0, mone, one, pdum, udum, &c) ; // udum <-- bu - B*p - B*(D*bp)
        
        cholmod_free_dense( &u, &c );
        //        t_omp = (double)omp_get_wtime;
        u = cholmod_solve (CHOLMOD_A, Lfact, udum, &c) ;
        //        printf( "substitution: %lf s\n", (double)((double)omp_get_wtime() - t_omp));
        
        copy_cholmod_dense_to_cholmod_dense( pdum, bp );      // pdum <-- bp
        cholmod_sdmult ( Ccm, 0, mone, one, u, pdum, &c);     // pdum <-- bp - C*u
        cholmod_sdmult ( Dcm, 0, one, one, pdum, p, &c) ;     // p <-- p + D*(bp - C*u)
        
        copy_cholmod_dense_to_cholmod_dense( fp, bp );
        cholmod_sdmult ( Ccm, 0, mone, one, u, fp, &c) ;
        //MinMaxArray(fp->x, scaling.E, matC->neq, "divU" );
        
        if (k>0) maxdivit = maxdiv;
        MinMaxArrayVal( fp->x, matC->neq, &mindiv, &maxdiv );
        MinMaxArrayVal( fu->x, matA->neq, &minru , &maxru  );
        
        if (k==0) maxdiv0 = maxdiv;
        
        if ( noisy == 1 ) {
            //            printf("PH iteration %01d:\n", k+1 );
            copy_cholmod_dense_to_cholmod_dense( fu, bu  );
            cholmod_sdmult ( Acm, 0, mone, one, u, fu, &c);
            cholmod_sdmult ( Bcm, 0, mone, one, p, fu, &c);
            copy_cholmod_dense_to_cholmod_dense( fp, bp );
            cholmod_sdmult ( Ccm, 0, mone, one, u, fp, &c);
            printf("PH it. %01d. rel. max. div. = %2.2e -- max. abs. div = %2.2e-- max. rel. div = %2.2e  / max. mom. = %2.2e rel. max. mom. = %2.2e\n", k, fabs(maxdiv/maxdiv0), maxdiv, maxdivit, maxru, maxru/maxru0);
            
            //cholmod_sdmult ( Dcm, 0, mone, one, p, fp, &c) ;
            //            NormResidualCholmod( &ru, &rp, fu, fp, matA->neq, matC->neq, model, scaling, 0 );
        }
        //        if (fabs(maxdiv/maxdiv0)<rel_tol_div) break;
        //        if (k>0 && fabs(maxdiv)/fabs(maxdivit)>0.75) break;
        if (k>1 && (fabs(maxdiv)<model.abs_tol_div || maxdiv/maxdiv0<rel_tol_div ) ) break;
        
        //        if ( ru<1e-11/(scaling.F/pow(scaling.L,3)) && rp<1e-11/scaling.E) break;
    }
    
    //    if (fabs(maxdiv/maxdiv0)>rel_tol_div || fabs(maxdiv) > rel_tol_div ){
    if (fabs(maxdiv)>model.abs_tol_div && maxdiv/maxdiv0>rel_tol_div) {
        printf("The code has exited since the incompressibility constrain was not satisfied to abs. tol. = %2.2e and rel. tol. = %2.2e\n Try modifying the PENALTY factor or check MIN/MAX viscosities\n Good luck!\n", model.abs_tol_div, rel_tol_div);
        exit(1);
    }
    printf("** Cholesky = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    SumArray(u->x, 1.0, matC->neq, "u");
    SumArray(p->x, 1.0, matC->neq, "p");
    
    // --------------- Solution vector --------------- //
    BackToSolutionVector( u, p, sol, mesh );
    
    // separate residuals
    double *F, Area =0.0, resx=0.0,resz=0.0, resp=0.0;
    int ndofx=0, ndofz=0, ndofp=0, cc;
    int nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1;
    F = DoodzCalloc(matA->neq+matC->neq, sizeof(double));
    
    BackToSolutionVector( fu, fp, F, mesh );
    
    // Integrate residuals
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nx, nzvx ) reduction(+:resx)
    for( cc=0; cc<nzvx*nx; cc++) {
        if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12 ) {
            ndofx++;
            resx += F[Stokes->eqn_u[cc]]*F[Stokes->eqn_u[cc]];//*celvol;
        }
    }
    
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nz, nxvz ) reduction(+:resz)
    for( cc=0; cc<nz*nxvz; cc++) {
        if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 && mesh->BCv.type[cc] != -12 ) {
            ndofz++;
            resz += F[Stokes->eqn_v[cc]]*F[Stokes->eqn_v[cc]];//*celvol;
            //            if ( mesh->BCv.type[cc] == 2) printf("F=%2.2e\n", F[Stokes->eqn_v[cc]]);
        }
    }
    
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, ncz, ncx ) reduction(+:resp)
    for( cc=0; cc<ncz*ncx; cc++) {
        if ( mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30  && mesh->BCp.type[cc] != 31) {
            ndofp++;
            Area += celvol;
            resp += F[Stokes->eqn_p[cc]]*F[Stokes->eqn_p[cc]];//*celvol;
        }
    }
    
    // Sqrt
    resx =  sqrt(resx/ndofx);
    resz =  sqrt(resz/ndofz);
    resp =  sqrt(resp/ndofp);
    
    if ( noisy == 1 ) {
        printf("Fu = %2.6e\n", resx ); // Units of momentum
        printf("Fv = %2.6e\n", resz ); // Units of momentum
        printf("Fp = %2.6e\n", resp ); // Units of velocity gradient
    }
    
    // --------------- Free --------------- //
    
    cholmod_free_dense( &bu, &c );
    cholmod_free_dense( &bp, &c );
    cholmod_free_dense( &u, &c );
    cholmod_free_dense( &p, &c );
    cholmod_free_dense( &fu, &c );
    cholmod_free_dense( &fp, &c );
    cholmod_free_dense( &pdum, &c );
    cholmod_free_dense( &udum, &c );
    cholmod_free_sparse( &Lcm, &c );
    
    
    
    cholmod_free_sparse( &Lcml, &c );
    cholmod_free_sparse( &Lcm, &c );
    cholmod_free_sparse( &Acm, &c );
    cholmod_free_sparse( &Bcm, &c );
    cholmod_free_sparse( &Ccm, &c );
    cholmod_free_sparse( &Dcm, &c );
    //    cholmod_free_factor ( &Lfact, &c) ;
    //    cholmod_finish( &c ) ;
    pardi->Lfact = Lfact;
    pardi->c     = c;
    
    cs_spfree(L);
    
    DoodzFree(u0);
    DoodzFree(p0);
    
    DoodzFree(A.p);
    DoodzFree(A.i);
    DoodzFree(A.x);
    cs_spfree(Ac);
    DoodzFree(B.p);
    DoodzFree(B.i);
    DoodzFree(B.x);
    cs_spfree(Bc);
    cs_spfree(B1);
    DoodzFree(C.p);
    DoodzFree(C.i);
    DoodzFree(C.x);
    cs_spfree(Cc);
    DoodzFree(D.p);
    DoodzFree(D.i);
    DoodzFree(D.x);
    cs_spfree(Dc);
    
    DoodzFree(F);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void KSPStokesDecoupled( SparseMat *matA,  SparseMat *matB,  SparseMat *matC,  SparseMat *matD, DirectSolver *pardi, double *rhs_mom, double *rhs_cont, double *sol, params model, grid *mesh, scale scaling, SparseMat *Stokes, SparseMat *Jacobian, SparseMat *JmatA,  SparseMat *JmatB,  SparseMat *JmatC ) {
    
    
    //double Control [AMD_CONTROL], Info [AMD_INFO] ;
    cs_di  A, B, D, C, *B1, *L, *Ac, *Bc,  *Cc,  *Dc, *L1;
    cs_di  AJ, BJ, CJ, *AJc, *BJc, *CJc;
    //int    *P, msglvl = 0;
    DoodzFP  *u0, *p0, *F;
    int  noisy=1;
    int k, cc; //nitmax=5
    double celvol = model.dx*model.dz;
    
    cholmod_common c ;
    cholmod_sparse *Lcm, *Lcml, *Acm, *Bcm, *Ccm, *Dcm; //, *A1
    cholmod_sparse *AcmJ, *BcmJ, *CcmJ;
    cholmod_sparse *Dcm0; //, *BDC, *DC, *Lcm0, *Lcml0, *Acm0, *Bcm0, *Ccm0
    cholmod_factor *Lfact ;
    cholmod_sparse *M, *AB, *CD, *Iu, *Ip, *D_zero; // *bs, *bus, *bps, *M0
    cholmod_dense *b, *x, *f, *val, *s, *v;
    //    cholmod_start( &c ) ;
    double mone[2] = {-1.0,0.0}, one[2] = {1.0,0.0}, zero[2] = {0.0,0.0};
    
    Lfact = pardi->Lfact;
    c     = pardi->c;
    
    // --------------- Pre-process--------------- //
    
    // ************** D ************** //
    SuiteSparse_long rsize;
    double gamma  = model.penalty;//1e12;//1e10*model.Nx*model.Nz;
    rsize = Stokes->neq_cont;
    Dcm0  = cholmod_speye (rsize, rsize, CHOLMOD_REAL, &c );
    for (k=0;k<Dcm0->nzmax;k++) ((double*)Dcm0->x)[k] *= -gamma*celvol;
    //    printf("-gamma*celvol = %2.2e %2.2e %2.2e %d %d\n", -gamma*celvol, model.dx*scaling.L, model.dz*scaling.L, model.Nx, model.Nz);
    
    clock_t t_omp;
    t_omp = (double)omp_get_wtime();
    
    // Build initial solution vector
    F  = DoodzCalloc(matA->neq+matC->neq, sizeof(double));
    u0 = DoodzCalloc( matA->neq, sizeof(double) );
    p0 = DoodzCalloc( matC->neq, sizeof(double) );
    BuildInitialSolutions( u0, p0, mesh );
    
    // Prepare A
    A.nzmax = matA->nnz;
    A.nz    = matA->nnz;
    A.m     = matA->neq;
    A.n     = A.m;
    //    printf( "A.nzmax = %d\n", matA->nnz );
    A.p     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.i     = DoodzCalloc( A.nzmax, sizeof(int) );
    A.x     = DoodzCalloc( A.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( A.m, matA->Ic, A.i );
    ArrayEqualArrayI( A.p, matA->J,  A.nzmax );
    ArrayEqualArray(  A.x, matA->A,  A.nzmax );
    Ac  = cs_di_compress( &A );
    cs_droptol( Ac, 1.0e-15 );
    
    // Prepare AJ
    AJ.nzmax = JmatA->nnz;
    AJ.nz    = JmatA->nnz;
    AJ.m     = JmatA->neq;
    AJ.n     = AJ.m;
    //    printf( "A.nzmax = %d\n", matA->nnz );
    AJ.p     = DoodzCalloc( AJ.nzmax, sizeof(int) );
    AJ.i     = DoodzCalloc( AJ.nzmax, sizeof(int) );
    AJ.x     = DoodzCalloc( AJ.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( AJ.m, JmatA->Ic, AJ.i );
    ArrayEqualArrayI( AJ.p, JmatA->J,  AJ.nzmax );
    ArrayEqualArray(  AJ.x, JmatA->A,  AJ.nzmax );
    AJc  = cs_di_compress( &AJ );
    cs_droptol( AJc, 1.0e-15 );
    
    // Prepare B
    B.nzmax = matB->nnz;
    B.nz    = matB->nnz;
    B.m     = matB->neq;
    B.n     = matC->neq;
    B.p     = DoodzCalloc( B.nzmax, sizeof(int) );
    B.i     = DoodzCalloc( B.nzmax, sizeof(int) );
    B.x     = DoodzCalloc( B.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( B.m, matB->Ic, B.i );
    ArrayEqualArrayI( B.p, matB->J,  B.nzmax );
    ArrayEqualArray(  B.x, matB->A,  B.nzmax );
    Bc  = cs_di_compress( &B );
    
    // Prepare BJ
    BJ.nzmax = JmatB->nnz;
    BJ.nz    = JmatB->nnz;
    BJ.m     = JmatB->neq;
    BJ.n     = JmatC->neq;
    BJ.p     = DoodzCalloc( BJ.nzmax, sizeof(int) );
    BJ.i     = DoodzCalloc( BJ.nzmax, sizeof(int) );
    BJ.x     = DoodzCalloc( BJ.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( BJ.m, JmatB->Ic, BJ.i );
    ArrayEqualArrayI( BJ.p, JmatB->J,  BJ.nzmax );
    ArrayEqualArray(  BJ.x, JmatB->A,  BJ.nzmax );
    BJc  = cs_di_compress( &BJ );
    
    // Prepare C
    C.nzmax = matC->nnz;
    C.nz    = matC->nnz;
    C.m     = matC->neq;
    C.n     = matB->neq;
    C.p     = DoodzCalloc( C.nzmax, sizeof(int) );
    C.i     = DoodzCalloc( C.nzmax, sizeof(int) );
    C.x     = DoodzCalloc( C.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( C.m, matC->Ic, C.i );
    ArrayEqualArrayI( C.p, matC->J,  C.nzmax );
    ArrayEqualArray(  C.x, matC->A,  C.nzmax );
    Cc  = cs_di_compress( &C );
    
    // Prepare CJ
    CJ.nzmax = JmatC->nnz;
    CJ.nz    = JmatC->nnz;
    CJ.m     = JmatC->neq;
    CJ.n     = JmatB->neq;
    CJ.p     = DoodzCalloc( CJ.nzmax, sizeof(int) );
    CJ.i     = DoodzCalloc( CJ.nzmax, sizeof(int) );
    CJ.x     = DoodzCalloc( CJ.nzmax, sizeof(double) );
    DecompressCSRtoTriplets( CJ.m, JmatC->Ic, CJ.i );
    ArrayEqualArrayI( CJ.p, JmatC->J,  CJ.nzmax );
    ArrayEqualArray(  CJ.x, JmatC->A,  CJ.nzmax );
    CJc  = cs_di_compress( &CJ );
    
    // Prepare D
    D.nzmax = Stokes->neq_cont;
    D.nz    = Stokes->neq_cont;
    D.m     = matD->neq;
    D.n     = matD->neq;
    D.p     = DoodzCalloc( D.nzmax+1, sizeof(int) );
    D.i     = DoodzCalloc( D.nzmax, sizeof(int) );
    D.x     = DoodzCalloc( D.nzmax, sizeof(double) );
    copy_cholmod_to_cs_matrix( Dcm0, &D );
    cholmod_free_sparse( &Dcm0, &c );
    Dc  = cs_di_compress( &D );
    
    // --------------- Schur complement --------------- //
    
    // Matrix multiplication: D*C
    L =  cs_di_multiply( Dc, Cc );
    
    //----- test - in case C' != B (assume B is deficient)
    B1 = cs_di_transpose( Cc, 1);
    //-----
    
    // Matrix multiplication: B*(D*C)
    L1 = cs_di_multiply( B1, L);
    cs_spfree(L);
    
    // Matrix addition: L = A - B*(D*C)
    L = cs_di_add( Ac, L1, 1, -1);
    cs_spfree(L1);
    
    // --------------- Cholesky --------------- //
    
    int N = Ac->m + Cc->m, cc2; //NNZ=Ac->nzmax+Bc->nzmax+Cc->nzmax+Dc->nzmax
    //    c.supernodal = 2;
    
    // Prepare Jacobian preconditioner blocks
    Lcm = cholmod_allocate_sparse (L->m, L->n, L->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Lcm, L );
    Acm = cholmod_allocate_sparse (Ac->m, Ac->n, Ac->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Acm, Ac );
    Bcm = cholmod_allocate_sparse (B1->m, B1->n, B1->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Bcm, B1 );
    Ccm = cholmod_allocate_sparse (Cc->m, Cc->n, Cc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Ccm, Cc );
    Dcm = cholmod_allocate_sparse (Dc->m, Dc->n, Dc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( Dcm, Dc );
    
    // Prepare Jacobian blocks
    AcmJ = cholmod_allocate_sparse (AJc->m, AJc->n, AJc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( AcmJ, AJc );
    BcmJ = cholmod_allocate_sparse (BJc->m, BJc->n, BJc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( BcmJ, BJc );
    CcmJ = cholmod_allocate_sparse (CJc->m, CJc->n, CJc->nzmax, 0, 1, 0, 1, &c) ;
    copy_cs_to_cholmod_matrix( CcmJ, CJc );
    
    // Keep lower part
    Lcml  = cholmod_copy ( Lcm,  -1, 1, &c );
    
    if (Lcml == NULL || Lcml->stype == 0)
    {
        printf("Unsymmetric matrix\n");
        cholmod_free_sparse (&Lcml, &c) ;
        cholmod_finish (&c) ;
    }
    
    if (pardi->Analyze == 1) {
        c.nmethods           = 1;
        c.method[0].ordering = CHOLMOD_AMD;
        c.postorder          = 1;
        t_omp                = (double)omp_get_wtime();
        Lfact                = cholmod_analyze( Lcml, &c ) ;
        pardi->Analyze       = 0;
        printf("** Time for Cholesky analysis = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    }
    
    t_omp = (double)omp_get_wtime();
    cholmod_factorize( Lcml, Lfact, &c);
    printf("** Time for Cholesky factorization = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
    
    cholmod_dense *u, *p, *bu, *bp, *pdum, *udum,  *fu, *fp;
    pdum = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    udum = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    bu   = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    bp   = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    u    = cholmod_allocate_dense( matA->neq, 1, matA->neq, CHOLMOD_REAL, &c );
    p    = cholmod_allocate_dense( matC->neq, 1, matC->neq, CHOLMOD_REAL, &c );
    
    // separate residuals
    double Area =0., resx=0.,resz=0., resp=0.;
    int ndofx=0, ndofz=0, ndofp=0;
    int nx=model.Nx, nz=model.Nz, nzvx=nz+1, nxvz=nx+1, ncx=nx-1, ncz=nz-1;
    
    //----------------- DEVELOPMENT -----------------//
    
    // Initialise KSP
    int restart = 20;//6;
    int max_it  = 1000;
    int ncycles = 0;
    int its     = 0, i1, i2, success=0;
    double eps  = 1.0e-9, norm_r, rnorm0, fact, r_dot_v, nrm;//, epsa  = 1e-8;
    double **VV, **SS;
    b      = cholmod_zeros (N, 1, CHOLMOD_REAL, &c );
    x      = cholmod_zeros (N, 1, CHOLMOD_REAL, &c );
    f      = cholmod_zeros (N, 1, CHOLMOD_REAL, &c );
    D_zero = cholmod_spzeros ( Dcm->nrow, Dcm->ncol, 0, CHOLMOD_REAL, &c ) ;
    fu     = cholmod_zeros (   Ac->m,     1, CHOLMOD_REAL, &c );
    fp     = cholmod_zeros (   Dc->m,     1, CHOLMOD_REAL, &c );
    val    = cholmod_zeros ( restart,     1, CHOLMOD_REAL, &c );
    s      = cholmod_zeros (       N,     1, CHOLMOD_REAL, &c );
    v      = cholmod_zeros (       N,     1, CHOLMOD_REAL, &c );
    Iu     = cholmod_speye (   Ac->m, Ac->n, CHOLMOD_REAL, &c );
    Ip     = cholmod_speye (   Dc->m, Dc->n, CHOLMOD_REAL, &c );
    
    // Allocate temp vecs
    VV  = (double**) DoodzCalloc( restart, sizeof(double*));
    SS  = (double**) DoodzCalloc( restart, sizeof(double*));
    for (cc=0; cc<restart; cc++) {
        VV[cc]  = (double*) DoodzCalloc( N, sizeof(double));
        SS[cc]  = (double*) DoodzCalloc( N, sizeof(double));
    }
    
    // Initial fields
    copy_vec_to_cholmod_dense( u, u0 );
    copy_vec_to_cholmod_dense( p, p0 );
    copy_vec_to_cholmod_dense( bu, rhs_mom );
    copy_vec_to_cholmod_dense( bp, rhs_cont );
    
    // Concatenate b = [bu; bp]
    for (cc=0; cc<bu->nrow; cc++) {
        ((double*)b->x)[cc ] = ((double*)bu->x)[cc];
    }
    for (cc=0; cc<bp->nrow; cc++) {
        cc2=bu->nrow+cc;
        ((double*)b->x)[cc2] = ((double*)bp->x)[cc];
    }
    
    // 1. Concatenate M = [A,B;C,D]
    AB  = cholmod_horzcat ( AcmJ, BcmJ,   1, &c );    // AB = [A B]
    CD  = cholmod_horzcat ( CcmJ, D_zero, 1, &c );    // CD = [C D]
    
    // 2. Concatenate M = [A,B;C,D]
    M   = cholmod_vertcat ( AB,  CD,     1, &c );     // M  = [AB; CD]
    cholmod_free_sparse( &AB,  &c );
    cholmod_free_sparse( &CD,  &c );
    
    // Concatenate s = [u; p]
    for (cc=0; cc<u->nrow; cc++) {
        ((double*)x->x)[cc] =  ((double*)u->x)[cc];
    }
    for (cc=0; cc<p->nrow; cc++) {
        cc2=u->nrow+cc;
        ((double*)x->x)[cc2] = ((double*)p->x)[cc];
    }
    
    // Compute residual f = b - Mx
    copy_cholmod_dense_to_cholmod_dense( f, b );
    cholmod_sdmult ( M, 0, mone, one, x, f, &c) ;
    MinMaxArray(M->x, 1, M->nzmax, "M");
    MinMaxArray(x->x, 1, x->nzmax, "x");
    MinMaxArray(b->x, 1, b->nzmax, "b");
    norm_r = cholmod_norm_dense ( f, 2, &c );
    
    // Evaluate reference norm
    rnorm0 = norm_r;
    printf("       %1.4d KSP GCR Residual %1.12e %1.12e\n", 0, norm_r, norm_r/rnorm0);
    
    while ( success == 0 && its<max_it ) {
        
        for ( i1=0; i1<restart; i1++ ) {
            
            // ---------------- Apply preconditioner, s = Q^{-1} r ---------------- //
            
            // extract ru from r = [ru; rp]
#pragma omp parallel for shared( fu, f ) private( cc )
            for (cc=0; cc<fu->nrow; cc++) {
                ((double*)fu->x)[cc] = ((double*)f->x)[cc];
            }
#pragma omp parallel for shared( fp, f ) private( cc, cc2 ) firstprivate( fu )
            for (cc=0; cc<fp->nrow; cc++) {
                cc2=fu->nrow+cc;
                ((double*)fp->x)[cc] = ((double*)f->x)[cc2];
            }
            
            cholmod_sdmult ( Dcm, 0,  one, zero, fp, pdum, &c );
            cholmod_sdmult ( Bcm, 0, mone,  one, pdum, fu, &c );
            cholmod_free_dense( &u, &c );
            u = cholmod_solve ( CHOLMOD_A, Lfact, fu, &c );
            cholmod_sdmult ( Dcm, 0,  one, zero, fp,    p, &c );
            cholmod_sdmult ( Ccm, 0,  one, zero,     u, pdum, &c );
            cholmod_sdmult ( Dcm, 0, mone,  one,  pdum,    p, &c );
            
            // Concatenate s = [u; p]
#pragma omp parallel for shared( s, u ) private( cc )
            for (cc=0; cc<u->nrow; cc++) {
                ((double*)s->x)[cc] = ((double*)u->x)[cc];
            }
#pragma omp parallel for shared( s, p ) private( cc, cc2 ) firstprivate( u )
            for (cc=0; cc<p->nrow; cc++) {
                cc2=u->nrow+cc;
                ((double*)s->x)[cc2] =((double*)p->x)[cc];
            }
            
            // Simplest preconditionning: PC = I
            //            ArrayEqualArray( s->x, f->x, N );
            
            // ---------------- Apply preconditioner, s = Q^{-1} r ---------------- //
            
            // Action of Jacobian on s
            cholmod_sdmult ( M, 0, one, zero, s, v, &c) ;
            
            // Approximation of the Jv product
            for (i2=0; i2<i1; i2++) {
                ((double*)val->x)[i2] = DotProduct( ((double*)v->x), VV[i2], N );
            }
            
            for(i2=0; i2<i1+1; i2++) {
                fact = -((double*)val->x)[i2];
                ArrayPlusScalarArray( ((double*)v->x), fact, VV[i2], N );
                ArrayPlusScalarArray( ((double*)s->x), fact, SS[i2], N );
            }
            
            r_dot_v = DotProduct( ((double*)f->x), ((double*)v->x), N );
            nrm     = sqrt(  DotProduct( ((double*)v->x), ((double*)v->x), N )  );
            r_dot_v = r_dot_v / nrm;
            
            fact = 1.0/nrm;
            ArrayTimesScalar( ((double*)v->x), fact, N );
            ArrayTimesScalar( ((double*)s->x), fact, N );
            
            fact = r_dot_v;
            ArrayPlusScalarArray( ((double*)x->x), fact, ((double*)s->x), N );
            fact =-r_dot_v;
            ArrayPlusScalarArray( ((double*)f->x), fact, ((double*)v->x), N );
            
            // Check convergence
            norm_r = cholmod_norm_dense ( f, 2, &c );
            if (norm_r < eps * rnorm0 ) { // || norm_r < epsa
                success = 1;
                break;
            }
            //            printf("[%1.4d] %1.4d KSP GCR Residual %1.12e %1.12e\n", ncycles, its, norm_r, norm_r/rnorm0);
            
            
            // Store arrays
            ArrayEqualArray( VV[i1], ((double*)v->x), N );
            ArrayEqualArray( SS[i1], ((double*)s->x), N );
            its++;
        }
        //        copy_cholmod_dense_to_cholmod_dense( fp, bp );
        //        cholmod_sdmult ( Ccm, 0, mone, one, u, fp, &c) ;
        //
        //        if ( noisy == 1 ) {
        //            //            printf("PH iteration %01d:\n", k+1 );
        //            copy_cholmod_dense_to_cholmod_dense( fu, bu );
        //            cholmod_sdmult ( Acm, 0, mone, one, u, fu, &c) ;
        //            cholmod_sdmult ( Bcm, 0, mone, one, p, fu, &c) ;
        //            copy_cholmod_dense_to_cholmod_dense( fp, bp );
        //            cholmod_sdmult ( Ccm, 0, mone, one, u, fp, &c) ;
        //        }
        its++;
        ncycles++;
    }
    printf("[%1.4d] %1.4d KSP GCR Residual %1.12e %1.12e\n", ncycles, its, norm_r, norm_r/rnorm0);
    MinMaxArray(fp->x, scaling.E, matC->neq, "divU" );
    
    // Extract ru from r = [ru; rp]
#pragma omp parallel for shared( f, x, u, fu ) private( cc )
    for (cc=0; cc<u->nrow; cc++) {
        ((double*)fu->x)[cc] = ((double*)f->x)[cc];
        ((double*) u->x)[cc] = ((double*)x->x)[cc];
    }
#pragma omp parallel for shared( f, x, p, fp ) private( cc, cc2 ) firstprivate( u )
    for (cc=0; cc<p->nrow; cc++) {
        cc2=u->nrow+cc;
        ((double*) p->x)[cc] = ((double*)x->x)[cc2];
        ((double*)fp->x)[cc] = ((double*)f->x)[cc2];
    }
    
    // --------------- Solution vector --------------- //
    BackToSolutionVector( u, p, sol, mesh );
    
    BackToSolutionVector( fu, fp, F, mesh );
    
    // Integrate residuals
    resx = 0.0;
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nx, nzvx ) reduction(+:resx)
    for( cc=0; cc<nzvx*nx; cc++) {
        if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13) {
            ndofx++;
            resx += F[Stokes->eqn_u[cc]]*F[Stokes->eqn_u[cc]]*celvol;
        }
    }
    
    resz = 0.0;
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, nz, nxvz ) reduction(+:resz)
    for( cc=0; cc<nz*nxvz; cc++) {
        if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13) {
            ndofz++;
            resz += F[Stokes->eqn_v[cc]]*F[Stokes->eqn_v[cc]]*celvol;
            //            if ( mesh->BCv.type[cc] == 2) printf("F=%2.2e\n", F[Stokes->eqn_v[cc]]);
        }
    }
    
    resp = 0.0;
#pragma omp parallel for shared( mesh, Stokes ) private( cc ) firstprivate( celvol, ncz, ncx ) reduction(+:resp)
    for( cc=0; cc<ncz*ncx; cc++) {
        if ( mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30  && mesh->BCp.type[cc] != 31) {
            ndofp++;
            Area += celvol;
            resp += F[Stokes->eqn_p[cc]]*F[Stokes->eqn_p[cc]]*celvol;
        }
    }
    
    // Sqrt
    resx =  sqrt(resx)/Area/(ndofx);
    resz =  sqrt(resz)/Area/(ndofz);
    resp =  sqrt(resp)/Area/(ndofp);
    
    if ( noisy == 1 ) {
        printf("Fu = %2.4e\n", resx * (scaling.F/pow(scaling.L,3))); // Units of momentum
        printf("Fv = %2.4e\n", resz * (scaling.F/pow(scaling.L,3))); // Units of momentum
        printf("Fp = %2.4e\n", resp * (scaling.E)); // Units of velocity gradient
    }
    
    // Free temp vecs
    for (cc=0; cc<restart; cc++) {
        DoodzFree(VV[cc]);
        DoodzFree(SS[cc]);
    }
    DoodzFree(VV);
    DoodzFree(SS);
    
    // Free
    cholmod_free_sparse( &M,  &c );
    cholmod_free_dense ( &f,  &c );
    cholmod_free_dense ( &b,  &c );
    cholmod_free_dense ( &x,  &c );
    //    cholmod_free_dense( &resiu,  &c );
    //    cholmod_free_dense( &resip,  &c );
    cholmod_free_dense( &val,  &c );
    cholmod_free_dense( &s,  &c );
    cholmod_free_dense( &v,  &c );
    cholmod_free_sparse( &Iu,  &c );
    cholmod_free_sparse( &Ip,  &c );
    cholmod_free_sparse( &D_zero,  &c );
    
    //----------------- DEVELOPMENT -----------------//
    
    // --------------- Free --------------- //
    cholmod_free_dense( &bu, &c );
    cholmod_free_dense( &bp, &c );
    cholmod_free_dense( &u, &c );
    cholmod_free_dense( &p, &c );
    cholmod_free_dense( &fu, &c );
    cholmod_free_dense( &fp, &c );
    cholmod_free_dense( &pdum, &c );
    cholmod_free_dense( &udum, &c );
    cholmod_free_sparse( &Lcm, &c );
    
    
    //    cholmod_free_factor ( &Lfact, &c) ;
    cholmod_free_sparse( &Lcml, &c );
    cholmod_free_sparse( &Lcm, &c );
    cholmod_free_sparse( &Acm, &c );
    cholmod_free_sparse( &Bcm, &c );
    cholmod_free_sparse( &Ccm, &c );
    cholmod_free_sparse( &Dcm, &c );
    
    cholmod_free_sparse( &AcmJ, &c );
    cholmod_free_sparse( &BcmJ, &c );
    cholmod_free_sparse( &CcmJ, &c );
    //    cholmod_finish( &c ) ;
    
    //    cholmod_free_factor ( &Lfact, &c) ;
    //    cholmod_finish( &c ) ;
    pardi->Lfact = Lfact;
    pardi->c     = c;
    
    cs_spfree(L);
    
    DoodzFree(u0);
    DoodzFree(p0);
    
    DoodzFree(A.p);
    DoodzFree(A.i);
    DoodzFree(A.x);
    cs_spfree(Ac);
    DoodzFree(B.p);
    DoodzFree(B.i);
    DoodzFree(B.x);
    cs_spfree(Bc);
    cs_spfree(B1);
    DoodzFree(C.p);
    DoodzFree(C.i);
    DoodzFree(C.x);
    cs_spfree(Cc);
    DoodzFree(D.p);
    DoodzFree(D.i);
    DoodzFree(D.x);
    cs_spfree(Dc);
    
    DoodzFree(AJ.p);
    DoodzFree(AJ.i);
    DoodzFree(AJ.x);
    cs_spfree(AJc);
    DoodzFree(BJ.p);
    DoodzFree(BJ.i);
    DoodzFree(BJ.x);
    cs_spfree(BJc);
    DoodzFree(CJ.p);
    DoodzFree(CJ.i);
    DoodzFree(CJ.x);
    cs_spfree(CJc);
    
    
    DoodzFree(F);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
