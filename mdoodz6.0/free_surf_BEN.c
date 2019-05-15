#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
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

void CorrectTopoIni_BEN( markers *particles, mat_prop materials, markers *topo_chain, surface *topo, params model, scale scaling, grid *mesh) {
    
    int    k, Ncx=model.Nx-1;
    double distance, dx=model.dx;
    int in;
    double grid_topo;
    
    //    // Find to which cell each marker contribute / find number of topo. markers per FINE cell column (DX/res)
    //    for (k=0;k<topo_chain->Nb_part;k++) {
    //
    //        if (topo_chain->x[k]>model.xmax || topo_chain->x[k]<model.xmin  ) topo_chain->phase[k]=-1;
    //        else topo_chain->phase[k]=0;
    //
    //        // Index of the fine grid column
    //        distance        = topo_chain->x[k] - (model.xmin + dx/2/res);
    //        in              = ceil((distance/dx*res)+0.5) - 1;
    //        if (in<0)    in = 0;
    //        if (in>res*Ncx-1)in = res*Ncx-1;
    //
    //        // Topography seen from the grid
    //        hm   = (topo->b[in]  + topo->a[in]  * ( topo_chain->x[k] ));
    //        if (topo_chain->z[k]>hm) topo_chain->z[k] = hm;
    //    }
    
    for (k=0;k<topo_chain->Nb_part;k++) {
        // Index of the coarse grid column
        distance        = (topo_chain->x[k]-model.xmin-dx/2.0);
        in              = ceil((distance/dx)+0.5) - 1;
        if (in<0)    in = 0;
        if (in>Ncx-1)in = Ncx-1;
        grid_topo       = (topo->b[in] + topo->a[in] * ( topo_chain->x[k] ));
        if ( topo_chain->z[k] > grid_topo ) {
            topo_chain->z[k]   = grid_topo;
        }
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddPartSed_BEN( markers *particles, mat_prop materials, markers *topo_chain, surface *topo, params model, scale scaling, grid *mesh) {
    
    int sed_phase = model.surf_ised1;
    int finite_strain = model.fstrain;
    int rec_T_P_x_z = model.rec_T_P_x_z;
    int time_My = floor(model.time*scaling.t / (3600.0*365.25*24.0*1.0e6));
    if ( time_My % 2 > 0 ) sed_phase = model.surf_ised1;
    else                   sed_phase = model.surf_ised2;
    
    int    k, Ncx=model.Nx-1, res = 2, ip;
    double distance, dx=model.dx, xW, xE, xmin = model.xmin;
    int nb_part = 2, in, new_ind = particles->Nb_part;
    double xnew, znew, hm0, hm;
    
    for ( k=0; k<res*Ncx; k++ ) {
        
        // West/East cell boundaries
        xW = xmin + k*dx/res;
        xE = xmin + k*dx/res + dx/res;
        
        for ( ip=0; ip<nb_part; ip++ ) {
            
            if (ip == 0) xnew = xW + dx/res/nb_part/2.0;
            if (ip == 1) xnew = xE - dx/res/nb_part/2.0;
            
            // Index of the coarse grid column
            distance        = xnew - model.xmin - dx/2.0;
            in              = ceil((distance/dx)+0.5) - 1;
            if (in<0)    in = 0;
            if (in>Ncx-1)in = Ncx-1;
            
            
            // Topography computed from the coarse grid column
            hm0  = (topo->b0[in] + topo->a0[in] * ( xnew ));
            hm   = (topo->b[in]  + topo->a[in]  * ( xnew ));
            
            if ( hm > hm0 ) { //;# && (hm-hm0)>1e-6
                
                // Compute new altitude of the marker right in between the new and old position of the free surface
                znew                      = 0.5*(hm+hm0);
                particles->x[new_ind]     = xnew;
                particles->z[new_ind]     = znew;
                particles->phase[new_ind] = sed_phase;
                
                particles->sxxd[new_ind]          =  0.0;
                particles->sxz[new_ind]           =  0.0;
                particles->Vx[new_ind]            =  0.0;
                particles->Vz[new_ind]            =  0.0;
                particles->strain[new_ind]        =  0.0;
                particles->strain_el[new_ind]     =  0.0;
                particles->strain_pl[new_ind]     =  0.0;
                particles->strain_pwl[new_ind]    =  0.0;
                particles->strain_exp[new_ind]    =  0.0;
                particles->strain_lin[new_ind]    =  0.0;
                particles->strain_gbs[new_ind]    =  0.0;
                particles->d[new_ind]             =  materials.gs_ref[sed_phase];
                particles->T[new_ind]             =  zeroC/scaling.T;
                particles->P[new_ind]             =  0.0;
                
                particles->phi[new_ind]           =  0.0;
                particles->X[new_ind]             =  0.0;
                particles->rho[new_ind]           =  0.0;
                particles->sxxd[new_ind]          =  0.0;
                particles->sxz[new_ind]           =  0.0;
                
                if (finite_strain==1) {
                    particles->Fxx[new_ind]           = 1.0;
                    particles->Fxz[new_ind]           = 0.0;
                    particles->Fzx[new_ind]           = 1.0;
                    particles->Fzz[new_ind]           = 0.0;
                }
                
                if (rec_T_P_x_z==1) {
                    particles->T0[new_ind]           = zeroC/scaling.T;
                    particles->P0[new_ind]           = 0.0;
                    particles->x0[new_ind]           = particles->x[new_ind];
                    particles->z0[new_ind]           = particles->z[new_ind];
                    particles->Tmax[new_ind]         = zeroC/scaling.T;
                    particles->Pmax[new_ind]         = 0.0;
                }
                
                new_ind++;
            }
        }
    }
    particles->Nb_part = new_ind;
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AdvectFreeSurf_BEN( markers *topo_chain, params model, scale scaling ) {
    
    int k;
    
    for (k=0; k<topo_chain->Nb_part; k++) {
        topo_chain->x[k] += model.dt*topo_chain->Vx[k];
        topo_chain->z[k] += model.dt*topo_chain->Vz[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// COMPLEX
void SetTopoChainHorizontalCoords_BEN( surface *topo, markers *topo_chain, params model, grid Mmesh, scale scaling ) {

    int k, fact=24, counter;
    double dx = model.dx/fact;
    topo_chain->Nb_part = model.Nx*fact-(fact-1) - 2;
    
    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        
        topo_chain->x[k]     = model.xmin + k*dx + dx;
        topo_chain->z[k]     = 0.0/scaling.L;
//        topo_chain->x0[k]    = topo_chain->x[k];
//        topo_chain->z0[k]    = topo_chain->z[k];
        topo_chain->phase[k] = 0;
        
//        // Save indices of some selected particles
//        if ( k%fact == 0 ) {
//            topo->VertInd[counter] = k-1+fact;
//            counter++;
//        }
//        if ( k == 0 )  topo->VertInd[0] = 0;
//        if ( k == topo_chain->Nb_part-1 )  topo->VertInd[model.Nx-1] = topo_chain->Nb_part-1;
    }
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

//// COMPLEX
//void SetTopoChainHorizontalCoords( surface *topo, markers *topo_chain, params model, grid Mmesh, scale scaling ) {
//    
//    int k, fact=24, counter;
//    topo_chain->Nb_part = (model.Nx-1)*fact;
//    double dx = (Mmesh.xc_coord[0][model.Nx-2] - Mmesh.xc_coord[0][0]) / (topo_chain->Nb_part-1);
//    
//    for ( k=0; k<topo_chain->Nb_part; k++ ) {
//        
//        topo_chain->x[k]     = Mmesh.xc_coord[0][0] + k*dx;
//        topo_chain->z[k]     = 0.0/scaling.L;
//        topo_chain->x0[k]    = topo_chain->x[k];
//        topo_chain->z0[k]    = topo_chain->z[k];
//        topo_chain->phase[k] = 0;
//        
//        // Save indices of some selected particles
//        if ( k%fact == 0 ) {
//            topo->VertInd[counter] = k-1+fact;
//            counter++;
//        }
//        if ( k == 0 )  topo->VertInd[0] = 0;
//        if ( k == topo_chain->Nb_part-1 )  topo->VertInd[model.Nx-1] = topo_chain->Nb_part-1;
//    }
//    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
//}

//void MarkerChainPolyFit( surface *topo, markers *topo_chain, params model, grid mesh ) {
//    
//    int k, ic, jc, im;
//    int ncx = model.Nx-1, ncz = model.Nz-1;
//    int *iCounts, **MarkerIndices;
//    double dx=mesh.dx[0], dz=mesh.dz[0];
//    double xmin = mesh.xg_coord[0][0] + dx/2;
//    double zmin = mesh.zg_coord[0][0] + dz/2;
//    double distance;
//    
//    iCounts       = DoodzCalloc ( ncx,sizeof(int) );
//    MarkerIndices = DoodzCalloc ( ncx,sizeof(int*) );
//    
//    // count number of surface marker per cells
//    for ( k=0; k<topo_chain->Nb_part; k++ ) {
//        
//        if (topo_chain->x[k]>model.xmin && topo_chain->x[k]<model.xmax ) {
//            
//            // Get the column:
//            distance=fabs(topo_chain->x[k] - xmin);
//            ic = ceil((distance/dx)+0.5) - 1;
//            iCounts[ic]++;
//        }
//    }
//    
//    // Allocate number of marker per cell
//    for ( ic=0; ic<ncx; ic++ ) {
//        MarkerIndices[ic] = DoodzCalloc ( iCounts[ic], sizeof(int) );
//        iCounts[ic]       = 0;
//    }
//    
//    // Add marker indices corresponding to cells
//    for ( k=0; k<topo_chain->Nb_part; k++ ) {
//        
//        if (topo_chain->x[k]>model.xmin && topo_chain->x[k]<model.xmax ) {
//            
//            // Get the column:
//            distance=fabs(topo_chain->x[k] - xmin);
//            ic = ceil((distance/dx)+0.5) - 1;
//            MarkerIndices[ic][iCounts[ic]] = k;
//            iCounts[ic]++;
//        }
//    }
//    
//    double xzbar, xbar, zbar, x2bar;
//    
//    // Build interpolant
//    for ( ic=0; ic<ncx; ic++ ) {
//        
//        xzbar = 0;
//        xbar  = 0;
//        zbar  = 0;
//        x2bar = 0;
//        
//        for (k=0; k<iCounts[ic]; k++ ) {
//            
//            im = MarkerIndices[ic][k];
//            
//            xzbar += topo_chain->x[im]*topo_chain->z[im];
//            xbar  += topo_chain->x[im];
//            zbar  += topo_chain->z[im];
//            x2bar += topo_chain->x[im]*topo_chain->x[im];
//        }
//        
//        xzbar /= iCounts[ic];
//        xbar  /= iCounts[ic];
//        zbar  /= iCounts[ic];
//        x2bar /= iCounts[ic];
//        
//        topo->a[ic] = (xzbar - xbar*zbar) / (x2bar - xbar*xbar);
//        topo->b[ic] =  zbar - topo->a[ic]*xbar;
//        topo->height[ic] = topo->b[ic] + topo->a[ic] * mesh.xc_coord[0][ic];
//    }
//    
//    double X, ZLeft, ZRight, dH;
//    
//    // Connect interpolants
//    for ( ic=0; ic<ncx-1; ic++ ) {
//        // common node
//        X      = mesh.xg_coord[0][ic+1];
//        
//        // Altitude computed from left
//        ZLeft  = topo->b[ic] + topo->a[ic] * X;
//        
//        // Altitude computed from right
//        ZRight = topo->b[ic+1] + topo->a[ic+1] * X;
//        
//        // Difference
//        dH = ZRight - ZLeft;
//        
//        // Correct
//        topo->b[ic+1] -= dH;
//    }
//    
//    // Free
//    for ( k=0; k<ncx; k++ ) {
//        DoodzFree( MarkerIndices[k] );
//    }
//    DoodzFree( MarkerIndices );
//    DoodzFree( iCounts );
//    
//}

void ProjectTopography_BEN( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling, double* X_vect, int itp_type ) {
    
    int k, in, jn, Nx=mesh.Nx;
    double dx=mesh.dx, distance, dxm, mark_val, *Xc_virtual, *Wm, *BmWm;
    
    // compute Xc_virtual
    Xc_virtual = DoodzMalloc ((Nx+1)*sizeof(double));    // allocate storage for the array
    Wm         = DoodzCalloc ( Nx, sizeof(double));
    BmWm       = DoodzCalloc ( Nx, sizeof(double));
    
    Xc_virtual[0]= X_vect[0]-0.5*dx;
    for (k=0;k<Nx-1;k++) {
        Xc_virtual[k+1]= 0.5*(X_vect[k+1]+X_vect[k]);
    }
    Xc_virtual[Nx]= X_vect[Nx-1]+0.5*dx;
    
//    for (k=0;k<Nx+1;k++) printf("%lf\n", Xc_virtual[k]*scaling.L);

    
    for (k=0;k<topo_chain->Nb_part;k++) {
        
        distance=fabs(topo_chain->x[k]-X_vect[0]);
        in=ceil((distance/dx)+0.5) - 1;
        
        dxm=fabs(0.5*(Xc_virtual[in]+Xc_virtual[in+1])-topo_chain->x[k]);
        
        mark_val = topo_chain->z[k];
        
        if (itp_type==1) {
            mark_val =  1/mark_val;
        }
        if (itp_type==2) {
            mark_val =  log(mark_val);
        }
        
        Wm[in]   += (1-(dxm/dx));
        BmWm[in] += mark_val*(1-(dxm/dx));
        
    }
    
    for (k=0;k<Nx;k++) {
//        if (WM[i]<1e-30 || mesh->BCg.type[i]==30) {
//            NodeField[i] = 0.0;
//            //printf("WARNING: Need to seed more particles for interpolation (zero weights): P2N\n");
//        }
//        else {
            topo->height[k] = BmWm[k]/Wm[k];
            if (itp_type==1) {
                topo->height[k] =  1 / topo->height[k];
            }
            if (itp_type==2) {
                topo->height[k] =  exp(topo->height[k]);
            }
//        printf("%lf\n", topo->height[k]*scaling.L);
//        }
    }
    
    DoodzFree(Xc_virtual);
    DoodzFree(Wm);
    DoodzFree(BmWm);
}

// SIMPLE
void MarkerChainPolyFit_BEN( surface *topo, markers *topo_chain, params model, grid mesh ) {
    
    int k, ic, jc, im;
    int ncx = model.Nx-1, ncz = model.Nz-1;
    int *iCounts, **MarkerIndices;
    double dx=mesh.dx, dz=mesh.dz;
    double xmin = mesh.xg_coord[0] + dx/2;
    double zmin = mesh.zg_coord[0] + dz/2;
    double distance;
    
    for ( ic=0; ic<ncx; ic++ ) {
        
        topo->a[ic] =  (topo->height[ic+1] - topo->height[ic])/(mesh.xg_coord[ic+1] - mesh.xg_coord[ic]);
        topo->b[ic] =  topo->height[ic] - (mesh.xg_coord[ic])*topo->a[ic];
    }
}

void RemeshMarkerChain_BEN( markers *topo_chain, surface *topo, params model, scale scaling, grid *mesh, int flag ) {
    
    double xmin = model.xmin + model.dx/2;
    int k, fact=24, ic, nx=model.Nx;
    double dx = model.dx/fact, distance;
    topo_chain->Nb_part = model.Nx*fact - (fact-1) - 2;
    
    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        
        // Reset x coordinate
        topo_chain->x[k]     = model.xmin + k*dx + dx;
        distance=fabs(topo_chain->x[k] - xmin);
        ic = ceil((distance/model.dx)+0.5) - 1;
        
        if ( ic<0)          ic = 0;
        if ( ic>model.Nx-1) ic =model.Nx-1;
        topo_chain->z[k]     = (topo->b[ic] + topo->a[ic] * ( topo_chain->x[k] ));
    }
    
    topo_chain->x[topo_chain->Nb_part] = topo_chain->x[topo_chain->Nb_part-1];
    topo_chain->z[topo_chain->Nb_part] = topo_chain->z[topo_chain->Nb_part-1];
    
//    // Calculate smoothed slope using markers on cells vertices coordinates
//    for ( ic=0; ic<model.Nx-1; ic++ ) {
//        //        printf("%d %d\n", topo->VertInd[ic],topo->VertInd[ic+1] );
//        //        printf("%lf %lf\n", topo_chain->z[topo->VertInd[ic]],topo_chain->z[topo->VertInd[ic+1]] );
//        //        printf("%lf ",topo->a[ic]);
//        topo->a[ic] =  (topo_chain->z[topo->VertInd[ic+1]] - topo_chain->z[topo->VertInd[ic]])/(topo_chain->x[topo->VertInd[ic+1]] - topo_chain->x[topo->VertInd[ic]]);
//        //        printf("%lf\n",topo->a[ic]);
//        
//        topo->b[ic] =  topo_chain->z[topo->VertInd[ic]] - (topo_chain->x[topo->VertInd[ic]])*topo->a[ic];
//    }
//    
//    for ( k=0; k<topo_chain->Nb_part; k++ ) {
//        distance=fabs(topo_chain->x[k] - xmin);
//        ic = ceil((distance/model.dx)+0.5) - 1;
//        topo_chain->z[k]     = topo->b[ic] + topo->a[ic]* topo_chain->x[k];
////        printf("%lf %lf %lf\n", topo_chain->z[k]*scaling.L, topo->a[ic], topo->b[ic]);
//        //        topo_chain->x0[k]    = topo_chain->x[k];
//        //        topo_chain->z0[k]    = topo_chain->z[k];
//        //        topo_chain->phase[k] = 0;
//    }
    
//    for ( k=1; k<topo_chain->Nb_part-1; k+=2 ) {
//        topo_chain->z[k]     = 0.5*(topo_chain->z[k-1] + topo_chain->z[k+1]);
//        
//    }
    
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

//// SIMPLE
//void SetTopoChainHorizontalCoords( surface *topo, markers *topo_chain, params model, grid Mmesh, scale scaling ) {
//    
//    int k, fact=1;
//    double dx = model.dx/fact;
//    topo_chain->Nb_part = model.Nx*fact;
//    
//    for ( k=0; k<topo_chain->Nb_part; k++ ) {
//        
//        topo_chain->x[k]     = model.xmin + k*dx;
//        topo_chain->z[k]     = 0.0/scaling.L;
//        
//        topo_chain->x0[k]    = topo_chain->x[k];
//        topo_chain->z0[k]    = topo_chain->z[k];
//        topo_chain->phase[k] = 0;
//    }
//    
//    topo_chain->x[0] += dx/2;
//    topo_chain->x[topo_chain->Nb_part-1] -= dx/2;
//    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
//}
//
//// SIMPLE
//void MarkerChainPolyFit( surface *topo, markers *topo_chain, params model, grid mesh ) {
//
//    int k, ic, jc, im;
//    int ncx = model.Nx-1, ncz = model.Nz-1;
//    int *iCounts, **MarkerIndices;
//    double dx=mesh.dx[0], dz=mesh.dz[0];
//    double xmin = mesh.xg_coord[0][0] + dx/2;
//    double zmin = mesh.zg_coord[0][0] + dz/2;
//    double distance;
//
//    for ( ic=0; ic<ncx; ic++ ) {
//
//        topo->a[ic] =  (topo_chain->z[ic+1] - topo_chain->z[ic])/(topo_chain->x[ic+1] - topo_chain->x[ic]);
//        topo->b[ic] =  topo_chain->z[ic] - (topo_chain->x[ic])*topo->a[ic];
//    }
//}
//
//// SIMPLE
//void RemeshMarkerChain( markers *topo_chain, surface *topo, params model, scale scaling, grid *mesh  ) {
//
//    double xmin = model.xmin + model.dx/2;
//    int k, fact=1, ic, nx = model.Nx;
//    double dx = model.dx/fact, distance;
//    topo_chain->Nb_part = model.Nx*fact;
//    double topo_s[nx];
//    
//    for ( k=0; k<nx; k++ ) {
//        
//        if (k==0)    topo_s[k]  = topo->b[k] + topo->a[k]*mesh->xc_coord[0][k];
//        if (k==nx-1) topo_s[k]  = topo->b[nx-2] + topo->a[nx-2]*mesh->xc_coord[0][nx-2];
//        if (k>0 && k<nx-1) {
//            topo_s[k]  = 0.5*(topo->b[k] + topo->a[k]*mesh->xc_coord[0][k]);
//            topo_s[k] += 0.5*(topo->b[k-1] + topo->a[k-1]*mesh->xc_coord[0][k-1]);
//        }
////        if (k==0)    topo_s[k]  = topo->b[k] + topo->a[k]*topo_chain->x[k];
////        if (k==nx-1) topo_s[k]  = topo->b[nx-2] + topo->a[nx-2]*topo_chain->x[k];
////        if (k>0 && k<nx-1) {
////            topo_s[k]  = 0.5*(topo->b[k] + topo->a[k]*topo_chain->x[k]);
////            topo_s[k] += 0.5*(topo->b[k-1] + topo->a[k-1]*topo_chain->x[k]);
////        }
////
//        
////        topo_chain->z[k] =  topo_s[k];
////        printf("%lf\n", topo_s[k]*scaling.L);
//    }
//    
////    printf("\n-----------\n");
////    
////    for ( k=0; k<nx-1; k++ ) {
////        printf("%lf\n", (topo->b[k] + topo->a[k]* mesh->xc_coord[0][k])*scaling.L);
////
////    }
//    
//    for ( k=0; k<topo_chain->Nb_part; k++ ) {
//        
//        topo_chain->x[k]     = model.xmin + k*dx;
//        
//        if ( k==0 ) topo_chain->x[k] += dx/2;
//        if ( k==topo_chain->Nb_part-1 ) topo_chain->x[k] -= dx/2;
//        
//        distance=fabs(topo_chain->x[k] - xmin);
//        ic = ceil((distance/model.dx)+0.5) - 1;
//        topo_chain->z[k]     = topo->b[ic] + topo->a[ic]* topo_chain->x[k];
//        topo_chain->z[k]     = 0.5*(topo->b[ic] + topo->a[ic] * ( topo_chain->x[k] ));
//        topo_chain->z[k]    += 0.5*0.5*(topo_s[ic] + topo_s[ic+1]);
////        topo_chain->z[k]    = 0.5*(topo_s[ic] + topo_s[ic+1]);
//        
////        if (ic==0) {
////            topo_chain->z[k]     = 0.5*(topo->b[ic] + topo->a[ic]* topo_chain->x[k]);
////            topo_chain->z[k]    += 0.5*(topo->b[ic+1] + topo->a[ic+1]* topo_chain->x[k]);
////        }
////        if (ic==nx-2) {
////            topo_chain->z[k]     = 0.5*(topo->b[ic] + topo->a[ic]* topo_chain->x[k]);
////            topo_chain->z[k]    += 0.5*(topo->b[ic-1] + topo->a[ic-1]* topo_chain->x[k]);
////        }
////        if (ic==0 && ic==nx-2) {
////            topo_chain->z[k]     = 0.333333333*(topo->b[ic] + topo->a[ic]* topo_chain->x[k]);
////            topo_chain->z[k]    += 0.333333333*(topo->b[ic-1] + topo->a[ic-1]* topo_chain->x[k]);
////            topo_chain->z[k]    += 0.333333333*(topo->b[ic+1] + topo->a[ic+1]* topo_chain->x[k]);
////        }
//        
//        topo_chain->x0[k]    = topo_chain->x[k];
//        topo_chain->z0[k]    = topo_chain->z[k];
//        topo_chain->phase[k] = 0;
//    }
//    
//    for ( k=1; k<topo_chain->Nb_part-1; k+=2 ) {
//        topo_chain->z[k]     = 0.5*(topo_chain->z[k-1] + topo_chain->z[k+1]);
//
//    }
//    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
//    printf("xminc=%lf  xmaxc=%lf\n", topo_chain->x[0]*scaling.L, topo_chain->x[topo_chain->Nb_part-1]*scaling.L );
//}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AllocateMarkerChain_BEN( surface *topo, markers* topo_chain, params model ) {
    
    topo_chain->Nb_part_max = 50*model.Nx;
    topo_chain->x0          = DoodzCalloc( topo_chain->Nb_part_max,sizeof(DoodzFP) );
    topo_chain->z0          = DoodzCalloc( topo_chain->Nb_part_max,sizeof(DoodzFP) );
    topo_chain->x           = DoodzCalloc( topo_chain->Nb_part_max,sizeof(DoodzFP) );
    topo_chain->z           = DoodzCalloc( topo_chain->Nb_part_max,sizeof(DoodzFP) );
    topo_chain->Vx          = DoodzCalloc( topo_chain->Nb_part_max, sizeof(DoodzFP) );
    topo_chain->Vz          = DoodzCalloc( topo_chain->Nb_part_max, sizeof(DoodzFP) );
    topo_chain->phase       = DoodzCalloc( topo_chain->Nb_part_max, sizeof(int) );
    
    topo->height0           = DoodzCalloc( (model.Nx),sizeof(DoodzFP) );
    topo->height            = DoodzCalloc( (model.Nx),sizeof(DoodzFP) );
    topo->vx                = DoodzCalloc( (model.Nx),sizeof(DoodzFP) );
    topo->vz                = DoodzCalloc( (model.Nx+1),sizeof(DoodzFP) );
    topo->a0                = DoodzCalloc( (model.Nx-1),sizeof(DoodzFP) );
    topo->b0                = DoodzCalloc( (model.Nx-1),sizeof(DoodzFP) );
    topo->a                 = DoodzCalloc( (model.Nx-1),sizeof(DoodzFP) );
    topo->b                 = DoodzCalloc( (model.Nx-1),sizeof(DoodzFP) );
    topo->VertInd           = DoodzCalloc( model.Nx,sizeof(DoodzFP) );
    printf( "Marker chain for topography was allocated\n" );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FreeMarkerChain_BEN( surface *topo, markers* topo_chain ) {
    
    DoodzFree( topo_chain->x0 );
    DoodzFree( topo_chain->z0 );
    DoodzFree( topo_chain->x );
    DoodzFree( topo_chain->z );
    DoodzFree( topo_chain->Vx );
    DoodzFree( topo_chain->Vz );
    DoodzFree( topo_chain->phase );
    DoodzFree( topo->height0 );
    DoodzFree( topo->height  );
    DoodzFree( topo->a0 );
    DoodzFree( topo->b0 );
    DoodzFree( topo->a  );
    DoodzFree( topo->b  );
    DoodzFree( topo->vx );
    DoodzFree( topo->vz );
    DoodzFree( topo->VertInd );

    printf( "Marker chain for topography was freed\n" );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double TopoFun0_BEN( double X, double H, double L, scale scaling ) {
    double Y = -10e3/scaling.L;
    Y = 0 + 2000/scaling.L*cos(X*M_PI/L*6);
    Y = -5000/scaling.L + 4500/scaling.L*cos(X*M_PI/L*1);
    return Y;
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double TopoFun_BEN( double X, int ic, surface topo, scale scaling ) {
    double Y;
    Y = topo.b[ic] + topo.a[ic] * X;
    return Y;
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void CellFlagging_BEN( grid *mesh, params model, surface topo, scale scaling ) {
    
    // Define dimensions;
    double Lx = (double) (model.xmax - model.xmin) ;
    double Lz = (double) (model.zmax - model.zmin) ;
    int nx   = mesh->Nx;
	int nz   = mesh->Nz;
	int ncx  = nx-1;
	int ncz  = nz-1;
	int nxvz = nx+1;
	int nzvx = nz+1;
    int i, j, c1, c2, c3;
    double h0, h, dz = fabs(mesh->zg_coord[1]-mesh->zg_coord[0]);
    int PVtag0[(ncx+2)*(ncz+2)], PVtag[(ncx+2)*(ncz+2)];
    
    
    //---------------------------- METHOD 1 
    
    //------------- FLAG EXTENDED PRESSURE CELLS -------------//
    for( i=0; i<ncx+2; i++ ) {
        
        // Topographic function        
        if (i==0) {
            h = TopoFun_BEN( mesh->xvz_coord[i], 0, topo, scaling );
        }
        
        if (i==ncx+1) {
            h = TopoFun_BEN( mesh->xvz_coord[i], ncx-1, topo, scaling );
        }
        
        if (i>0 && i<ncx+1) {
            h = TopoFun_BEN( mesh->xvz_coord[i], i-1, topo, scaling );
        }
        
        //-------------------------------------------------------------//

        for( j=0; j<ncz+2; j++ ) {
            
            c1 = i + j*(ncx+2);
            
            PVtag0[c1] = -1;

            
            if (mesh->zvx_coord[j]>h) {
                // Deactivate pressure cells over the surface
                PVtag0[c1] = 30;
            }

        }
    }
    
    //---------------------------- METHOD 1
    
//    //---------------------------- METHOD 2 
//    
//    //------------- FLAG EXTENDED PRESSURE CELLS -------------//
//    for( i=0; i<ncx+2; i++ ) {
//        
//        // Topographic function
//        if (i==0) {
//            h = TopoFun_BEN( mesh->xvz_coord[0][i], 0, topo, scaling );
//        }
//        
//        if (i==ncx+1) {
//            h = TopoFun_BEN( mesh->xvz_coord[0][i], ncx-1, topo, scaling );
//        }
//        
//        if (i>0 && i<ncx+1) {
//            h = TopoFun_BEN( mesh->xvz_coord[0][i], i-1, topo, scaling );
//        }
//        
//        //-------------------------------------------------------------//
//        
//        for( j=0; j<ncz+2; j++ ) {
//            
//            c1 = i + j*(ncx+2);
//            
//            PVtag0[c1] = -1;
//            
//            if (fabs(mesh->zvx_coord[0][j]-h)<dz/2) {
//                // Deactivate pressure cells over the surface
//                PVtag0[c1] = 60;
//            }
//        }
//    }
//    
//    int found;
//    for( i=0; i<ncx+2; i++ ) {
//        found = 0;
//        for( j=0; j<ncz+2; j++ ) {
//             c1 = i + j*(ncx+2);
//            
//            if ( found == 1 ) {
//                PVtag0[c1] = 30;
//            }
//            
//            if ( PVtag0[c1] == 60 ) {
//                found = 1;
//            }
//             
//        }
//    }
//    
//    //---------------------------- METHOD 2 
    
    
    for( i=0; i<ncx+2; i++ ) {
        for( j=0; j<ncz+2; j++ ) {
            
            c1 = i + j*(ncx+2);
        
            PVtag[c1] = PVtag0[c1];

            if ( (j>0 && j<ncz+1) && (PVtag0[c1] != PVtag0[c1-(ncx+2)]) ) {
                PVtag[c1] = 60;
            }
        }
    }
    
    for( i=1; i<ncx+1; i++ ) {
        for( j=0; j<ncz+2; j++ ) {
            
            c1 = i + j*(ncx+2);
                        
            if ( PVtag[c1]==-1 && (PVtag[c1-1]==60 && PVtag[c1+1]==60) ) {
                PVtag[c1] = 60; 
            }
        }
    }
    
    
    // Sides
    for( j=0; j<ncz+2; j++ ) {
        
        c1 = 0 + j*(ncx+2);
        c1 = ncx+1 + j*(ncx+2);
        
        PVtag[c1] = PVtag[c1-1];
    }
    
   

    
    //------------- FLAG PRESSURE/TEMPERATURE CELLS -------------//
    
//    printf("SURFACE TEMP.\n");
    for( j=0; j<ncz; j++ ) {
        for( i=0; i<ncx; i++ ) {
            
            c1 = i + j*(ncx);
            c2 = i + 1 + (j+1)*(ncx+2);
            
            mesh->BCp.type[c1] = PVtag[c2];
            mesh->BCt.type[c1]    = PVtag[c2];
            
            if ( PVtag[c2] == 30 ) {
                mesh->BCt.type[c1] = 30;
                mesh->BCt.val[c1]  = 273.15/scaling.T;
                mesh->T[c1]     = 273.15/scaling.T;
            }
            
            if ( PVtag[c2] == 60 ) {
                mesh->BCp.type[c1] = 0;
                mesh->BCp.val[c1]  = 0.0;
                
//                if ( i==0 ) mesh->BCp.type[0][c1] = 0;
                
                // TEMPERATURE ALONG THE SURFACE
                mesh->BCt.type[c1] = 30;
                mesh->BCt.val[c1]  = 0.0;

            }
        }
    }
    
    //------------- FLAG Vx CELLS -------------//
    for ( i=0; i<nx; i++ ) {
        for( j=0; j<ncz+2; j++ ) {
            c2 = i + (j)*(nx);
            mesh->BCu.type[c2] = 30;
            mesh->BCu.val[c2]  = 0.0;
        }
    }
    
    for ( i=0; i<ncx+2; i++ ) {
        for( j=0; j<ncz+2; j++ ) {
            
            c1 = i + j*(ncx+2);
            c2 = i + (j)*(nx);
            
            if ( i<ncx+1 && PVtag[c1]==-1 ) {
                mesh->BCu.type[c2] = -1;
                mesh->BCu.val[c2]  = 0;
            }
            
            if ( i>0 && PVtag[c1]==-1) {
                mesh->BCu.type[c2-1] = -1;
                mesh->BCu.val[c2-1]  = 0;
            }
        }
    }
    
    
        
//            for ( i=0; i<nx; i++ ) {
//                for( j=0; j<nzvx; j++ ) {
//            
//            c2 = i + (j)*(nx);
//            
//            if ( j<nzvx && mesh->BCu.type[0][c2] == 30 && mesh->BCu.type[0][c2-nx] == -1 ) {
//                mesh->BCu.type[0][c2] = 4;
//                mesh->BCu.val[0][c2]  = 0;
//            }
//        }
//    }
    
    for (j=0; j<nzvx; j++) {
        
        c1 = 0 + j*(nx);
        c2 = nx-1 + j*(nx);
        
        mesh->BCu.type[c1] = mesh->BCu.type[c1+1];
        mesh->BCu.type[c2] = mesh->BCu.type[c2-1];
    }
    
    //------------- FLAG Vz CELLS -------------//
    for (i=0; i<ncx+2; i++) {
        for (j=0; j<nz; j++) {
            
            c2 = i + j*(nxvz);
            mesh->BCv.type[c2] = 30;
            mesh->BCv.val[c2] = 0.0;
        }
    }
    
    
    for (i=0; i<ncx+2; i++) {
        for (j=0; j<ncz+2; j++) {
            
            c1 = i + j*(ncx+2);
            c2 = i + j*(nxvz);
            
            if ( j<ncz+1 && PVtag[c1]==-1 ) {
                mesh->BCv.type[c2] = -1;
                mesh->BCv.val[c2]  = 0;
            }
            
            if ( j>0  && PVtag[c1]==-1 ) {
                mesh->BCv.type[c2-nxvz] = -1;
                mesh->BCv.val[c2-nxvz]  = 0;
            }
            
        }
    }
    
    
    for (j=0; j<nz; j++) {
        
        c1 = 0 + j*(nxvz);
        c2 = nxvz-1 + j*(nxvz);
        
        mesh->BCv.type[c1] = mesh->BCv.type[c1+1];
        mesh->BCv.type[c2] = mesh->BCv.type[c2-1];
    }
    
    //------------- FLAG VERTICES -------------//
    for( j=0; j<nz; j++ ) {
        for ( i=0; i<nx; i++ ) {
            
            c1 = i + j*nx;
            c2 = i + j*nxvz;
            
            if (i==0)    h  = topo.b[i] + topo.a[i]*mesh->xc_coord[i];
            if (i==nx-1) h  = topo.b[nx-2] + topo.a[nx-2]*mesh->xc_coord[nx-2];
            if (i>0 && i<nx-1) {
                h  = 0.5*(topo.b[i] + topo.a[i]*mesh->xc_coord[i]);
                h += 0.5*(topo.b[i-1] + topo.a[i-1]*mesh->xc_coord[i-1]);
            }
            
            
            mesh->BCg.type[c1] = -1;
            
            if (mesh->zg_coord[j] > h) {
                mesh->BCg.type[c1] = 30;
            }
//
//
//            if ( j>0 ) {
//                if ( mesh->BCu.type[0][c1] == 30 || mesh->BCu.type[0][c1-nx] == 30 ) { /////////// ???????? +nx
//                    mesh->BCg.type[c1] = 30;
//                }
//            }
//            
//            if ( i<nxvz-1 ) {
//                if ( mesh->BCv.type[0][c2] == 30 || mesh->BCv.type[0][c2+1] == 30 ) {
//                    mesh->BCg.type[c1] = 30;
//                }
//            }
            
        }
    }

    
    // set surface vetrices to zero vel
    for ( i=0; i<nx; i++ ) {
        for( j=0; j<nz; j++ ) {
            
            c1 = i + j*nx;
            
            //            if ( (i==0 || i==nx-1) && j<nz-1 ) {
            //
            //                if ( mesh->BCg.type[c1] == -1 && mesh->BCg.type[c1+nx] == 30 ) mesh->BCg.type[c1] = 30;
            //            }
            
            if ( j<nz-1 ) {
                
                if ( mesh->BCg.type[c1] == -1 && mesh->BCg.type[c1+nx] == 30 ) mesh->BCg.type[c1] = 30;
            }
            //            if ( j<nz-1 ) {
            //
            //                if ( mesh->BCg.type[c1] == -1 && mesh->BCg.type[c1+nx] == 30 && mesh->BCu.type[0][c1+nx] == 30) mesh->BCg.type[c1] = 30;
            //            }
            
        }
    }
    
    
//    //------------- FLAG VERTICES -------------//
//    for( j=0; j<nz; j++ ) {
//        for ( i=0; i<nx; i++ ) {
//            
//            c1 = i + j*nx;
//            c2 = i + j*nxvz;
//            
////            mesh->BCg.type[c1] = 30;
//            
//            if (i==0) {
//                if ( mesh->BCv.type[0][c2+1] == -1 ) mesh->BCg.type[c1] = -1;
//            }
//            if (i==nx-1) {
//                if ( mesh->BCv.type[0][c2] == -1 )   mesh->BCg.type[c1] = -1;
//            }
//            
//            if (i>0 && i<nx-1) {
//                if ((mesh->BCv.type[0][c2]==-1 && mesh->BCv.type[0][c2+1]==-1) && (mesh->BCu.type[0][c1]==-1 && mesh->BCu.type[0][c1+nx]==-1)) {
//                    mesh->BCg.type[c1] = -1;
//                }
//            }
//            
//        }
//    }
    
//    //------------- FLAG VERTICES -------------//
//    for( j=0; j<nz; j++ ) {
//        for ( i=0; i<nx; i++ ) {
//            
//            c1 = i + j*nx;
//            c2 = i + j*nxvz;
//            
//            mesh->BCg.type[c1] = -1;
//            
//            
//            if ( j>0 ) {
//                if ( mesh->BCu.type[0][c1] == 30 || mesh->BCu.type[0][c1-nx] == 30 ) { /////////// ???????? +nx
//                    mesh->BCg.type[c1] = 30;
//                }
//            }
//            
//            if ( i<nxvz-1 ) {
//                if ( mesh->BCv.type[0][c2] == 30 || mesh->BCv.type[0][c2+1] == 30 ) {
//                    mesh->BCg.type[c1] = 30;
//                }
//            }
//            
//        }
//    }
    

    
    // METHOD 2 Vertices
    for( j=0; j<nz; j++ ) {
        for ( i=0; i<nx; i++ ) {

            c1 = i + j*nx;
            c2 = i + j*ncx;
            if ( i>0 && i<nx-1 && j>0 && j<nz-1 ) {
                if ( mesh->BCg.type[c1] == 30) {
                    
                    
                    if ( mesh->BCp.type[c2-1] == -1  && mesh->BCp.type[c2-1-ncx] == -1 && mesh->BCp.type[c2] == -1 && mesh->BCp.type[c2-ncx] == -1 ) {
                        //                printf("ALLO!\n");
                        mesh->BCg.type[c1] = -1;
                }

            }
            }
        }
    }
   
//    //------------- Clean viscosity and density field -------------//
    
//    // Cells
//    for( j=0; j<ncz; j++ ) {
//        for ( i=0; i<ncx; i++ ) {
//            
//            c1 = i + j*ncx;
//            c2 = i + 1 + (j+1)*(ncx+2);
//            
//            if ( mesh->BCt.type[c1] == 30 ) {
//                mesh->eta_n[0][c1] = 0.0;
//                mesh->rho_n[0][c1] = 1.0;
//            }
//        }
//    }
//    
//    // Vertices
//    for( j=0; j<nz; j++ ) {
//        for ( i=0; i<nx; i++ ) {
//            
//            c1 = i + j*nx;
//            
//            if ( mesh->BCg.type[c1] == 30 ) {
//                mesh->eta_s[0][c1] = 0.0;
//                mesh->rho_s[0][c1] = 1.0;
//            }
//        }
//    }
    
    
    
//    int c,k,l;

//    printf("PV tags:\n");
//    for (l=0; l<ncz+2; l++) {
//        for (k=0; k<ncx+2; k++) {
//            c = k + l*(ncx+2);
//            printf("%02d ", PVtag[c]);
//            
//        }
//        printf("\n");
//    }
//    
//    printf("P tags:\n");
//    for (l=0; l<ncz; l++) {
//        for (k=0; k<ncx; k++) {
//            c = k + l*(ncx);
//            printf("%02d ", mesh->BCp.type[0][c]);
//            
//        }
//        printf("\n");
//    }
//
//    
//    printf("T tags:\n");
//    for (l=0; l<ncz; l++) {
//        for (k=0; k<ncx; k++) {
//            c = k + l*(ncx);
//            printf("%02d ", mesh->BCt.type[c]);
//            
//        }
//        printf("\n");
//    }
//    
//       
//    printf("Vx tags:\n");
//    for (l=0; l<nzvx; l++) {
//        for (k=0; k<nx; k++) {
//            c = k + l*(nx);
//            printf("%02d ", mesh->BCu.type[0][c]);
//            
//        }
//        printf("\n");
//    }
//    
//    printf("Vz tags:\n");
//    for (l=0; l<nz; l++) {
//        for (k=0; k<nxvz; k++) {
//            c = k + l*(nxvz);
//            printf("%02d ", mesh->BCv.type[0][c]);
//            
//        }
//        printf("\n");
//    }
//
//    printf("Vertex tags:\n");
//    for (l=0; l<nz; l++) {
//        for (k=0; k<nx; k++) {
//            c = k + l*(nx);
//            printf("%02d ", mesh->BCg.type[c]);
//            
//        }
//        printf("\n");
//    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void CleanUpSurfaceParticles_BEN( markers* particles, grid *mesh, surface topo, scale scaling ) {
    
    int k, ic, ncx=mesh->Nx-1;
    double h;
    double dx=mesh->dx, dz=mesh->dz;
    double xmin = mesh->xg_coord[0] + dx/2;
    double distance;
    
    int count = 0;
    
    //#pragma omp parallel for shared ( particles, topo ) private ( k, h, ic, distance, dx, xmin )
    for ( k=0; k<particles->Nb_part; k++ ) {
        
        if ( particles->phase[k] != -1 ) {
            
            // Get the column:
            distance=fabs(particles->x[k] - xmin);
            ic = ceil((distance/dx)+0.5) - 1;
            
            if (ic<0) ic = 0;
            if (ic>ncx-1) ic = ncx-1;
            
            h = topo.b[ic] + topo.a[ic]*particles->x[k];
            
            if ( particles->z[k]>h ) {
                particles->phase[k] = -1;
                count++;
            }
        }
    }
    printf("%d particle above surface, Nb_part: %d\n", count, particles->Nb_part);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void SurfaceDensityCorrection_BEN( grid *mesh, params model, surface topo, scale scaling ) {
    
    // Define dimensions;
    double Lx = (double) (model.xmax - model.xmin) ;
    double Lz = (double) (model.zmax - model.zmin) ;
    int nx   = mesh->Nx;
	int nz   = mesh->Nz;
	int ncx  = nx-1;
	int ncz  = nz-1;
	int nxvz = nx+1;
	int nzvx = nz+1;
    int i, j, c1, c2;
    double h0, h, dz = fabs(mesh->zg_coord[1]-mesh->zg_coord[0]);
    
    
    // Density on cell centers
    for( j=0; j<ncz; j++ ) {
    for( i=0; i<ncx; i++ ) {
            c1 = i + j*(ncx);
            
//            printf("%02d ", mesh->BCp.type[0][c1]);
        
            if (mesh->BCp.type[c1] == -1 && mesh->BCp.type[c1+ncx] == 0 ) {
                h  = topo.b[i] + topo.a[i]*mesh->xc_coord[i];
                h0 = fabs(h - mesh->zc_coord[j]);
                mesh->rho_n[c1] *= h0/dz;
//                printf("h0=%lf rho=%2.2e\n", h0*scaling.L, mesh->rho_n[0][c1]*scaling.rho);
            }
//        printf("%2.2e ", mesh->rho_n[0][c1]*scaling.rho);

        }
//        printf("\n");
    }
    
    // Density on cell vertices   
    for( j=0; j<nz; j++ ) {
        for( i=0; i<nx; i++ ) {
            
            c1 = i + j*(nx);
            
            if (mesh->BCg.type[c1] == -1 && mesh->BCg.type[c1+nx] == 30) {
                if (i==0)    h  = topo.b[i] + topo.a[i]*mesh->xc_coord[i];
                if (i==nx-1) h  = topo.b[nx-2] + topo.a[nx-2]*mesh->xc_coord[nx-2];
                if (i>0 && i<nx-1) {
                    h  = 0.5*(topo.b[i] + topo.a[i]*mesh->xc_coord[i]);
                    h += 0.5*(topo.b[i-1] + topo.a[i-1]*mesh->xc_coord[i-1]);
                }
                h0 = fabs(h - mesh->zg_coord[j]);
                mesh->rho_s[c1] = h0/dz*mesh->rho_s[c1];
            }
//            printf("%2.2e ", mesh->rho_s[0][c1]*scaling.rho);

        }
//        printf("\n");

    }
    
  
    
//    // Density on Vz points
//    for( j=0; j<nz; j++ ) {
//    for( i=1; i<nxvz-1; i++ ) {
//            
//            c1 = i + j*(nxvz);
//    
//            if (mesh->BCv.type[0][c1] == -1 && mesh->BCv.type[0][c1+nxvz] == 30) {
//                h  = topo.b[i-1] + topo.a[i-1]*mesh->xc_coord[0][i-1];
//                h0 = fabs(h - mesh->zg_coord[0][j]);
//                mesh->rhoVz[c1] *= h0/dz;
//            }
////        printf("%2.2e ", mesh->rhoVz[c1]*scaling.rho);
//
//        }
////        printf("\n");
//    }
//    
//    // Density on Vx points
//    for( j=1; j<nzvx-1; j++ ) {
//        for( i=0; i<nx; i++ ) {
//            
//            c1 = i + j*(nx);
//                        
//            if (mesh->BCu.type[0][c1] == -1 && mesh->BCu.type[0][c1+nx] == 30) {
//                if (i==0)    h  = topo.b[i] + topo.a[i]*mesh->xc_coord[0][i];
//                if (i==nx-1) h  = topo.b[nx-2] + topo.a[nx-2]*mesh->xc_coord[0][nx-2];
//                if (i>0 && i<nx-1) {
//                    h  = 0.5*(topo.b[i] + topo.a[i]*mesh->xc_coord[0][i]);
//                    h += 0.5*(topo.b[i-1] + topo.a[i-1]*mesh->xc_coord[0][i-1]);
//                }
//
//                h0 = fabs(h - mesh->zvx_coord[0][j]);
//                mesh->rhoVx[c1] *= h0/dz;
//            }
//        }
//    }
    
  }

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SurfaceVelocity_BEN( grid *mesh, params model, surface *topo, markers* topo_chain, scale scaling ) {
    
    // Define dimensions;
    double Lx = (double) (model.xmax - model.xmin) ;
    double Lz = (double) (model.zmax - model.zmin) ;
    int nx   = mesh->Nx;
	int nz   = mesh->Nz;
	int ncx  = nx-1;
	int ncz  = nz-1;
	int nxvz = nx+1;
	int nzvx = nz+1;
    int i, j, c1, c2, c3, k;
    double h0, h, dz = fabs(mesh->zg_coord[1]-mesh->zg_coord[0]), dx = fabs(mesh->xg_coord[1]-mesh->xg_coord[0]);
    
//    // Interpolate velocities to cell centers
//    double Vxc[ncx*ncz], Vzc[ncx*ncz];
//    
//    for( j=0; j<ncz; j++ ) {
//        for( i=0; i<ncx; i++ ) {
//            
//            c1 = i + j*(nx);
//            c2 = i + j*(ncx);
//            c3 = i + j*(nxvz);
//            
//            Vxc[c2] = 0.5 *( mesh->u_in[c1] + mesh->u_in[c1+1] );
//            Vzc[c2] = 0.5 *( mesh->v_in[c3] + mesh->v_in[c3+nxvz] );
//        }
//    }
//    Interp_Grid2P( *topo_chain, topo_chain->Vx,    mesh, Vxc, mesh->xc_coord[0],  mesh->zc_coord[0],  mesh->Nx[0]-1, mesh->Nz[0]-1, mesh->BCp.type[0] );
//    Interp_Grid2P( *topo_chain, topo_chain->Vz,    mesh, Vzc, mesh->xc_coord[0],  mesh->zc_coord[0],  mesh->Nx[0]-1, mesh->Nz[0]-1, mesh->BCp.type[0] );
    
    
//    // Interpolate velocities to cell centers
//    double Vx2[nx*nzvx], Vz2[nxvz*nz];
//    
//    for( j=0; j<nzvx; j++ ) {
//        for( i=0; i<nx; i++ ) {
//            
//            c1 = i + j*(nx);
//            Vx2[c1] = mesh->u_in[c1];
//            
//            if (j>0) {
//                if (mesh->BCu.type[0][c1] == 30 && mesh->BCu.type[0][c1-nx] == -1) {
//                    Vx2[c1] = Vx2[c1-nx];
//
//                }
//            }
//            
//            if (j>1) {
//                if (mesh->BCu.type[0][c1] == 30 && mesh->BCu.type[0][c1-2*nx] == -1) {
//                    Vx2[c1] = Vx2[c1-2*nx];
//                    
//                }
//            }
//        }
//    }
//    
//    for( j=0; j<nz; j++ ) {
//        for( i=0; i<nxvz; i++ ) {
//            
//            c2 = i + j*(nxvz);
//            
//            Vz2[c2] = mesh->v_in[c2];
//            
//            if (i==0) {Vz2[c2] = Vz2[c2+1];printf("ALLO 1?\n");}
//            if (i==nxvz-1) {Vz2[c2] = Vz2[c2-1];printf("ALLO 1?\n");}
//            
//            if (j>0) {
//                if (mesh->BCv.type[0][c2] == 30 && mesh->BCv.type[0][c2-nxvz] == -1) {
//                    Vz2[c2] = Vz2[c2-nxvz];
//                }
//            }
////            if (j>1) {
////                if (mesh->BCv.type[0][c2] == 30 && mesh->BCv.type[0][c2-2*nxvz] == -1) {
////                    Vz2[c2] = Vz2[c2-2*nxvz];
////                }
////            }
//        }
//    }
//    
//    Interp_Grid2P( *topo_chain, topo_chain->Vx,    mesh, Vx2, mesh->xg_coord[0],  mesh->zvx_coord[0],  mesh->Nx[0], mesh->Nz[0]+1, mesh->BCu.type[0] );
//    Interp_Grid2P( *topo_chain, topo_chain->Vz,    mesh, Vz2, mesh->xvz_coord[0],  mesh->zg_coord[0],  mesh->Nx[0]+1, mesh->Nz[0], mesh->BCv.type[0] );
    
    // Build surface velocity vectors from the mesh
    for( j=0; j<nzvx; j++ ) {
        for( i=0; i<nx; i++ ) {
            
            c1 = i + j*(nx);
            
            if (j>0) {
                if (mesh->BCu.type[c1] == 30 && mesh->BCu.type[c1-nx] == -1) {
                    topo->vx[i] = mesh->u_in[c1-nx];
                }
            }
        }
    }
//    topo->vx[1]    = 0.5*(topo->vx[0]+topo->vx[3]);
//    topo->vx[nx-2] = 0.5*(topo->vx[nx]+topo->vx[nx-3]);
    
    for( j=0; j<nz; j++ ) {
        for( i=0; i<nxvz; i++ ) {
            
            c2 = i + j*(nxvz);
            
            if (j>0) {
                if (mesh->BCv.type[c2] == 30 && mesh->BCv.type[c2-nxvz] == -1) {
                    topo->vz[i] = mesh->v_in[c2-nxvz];
                }
            }
        }
    }
//    topo->vz[1]    = 0.5*(topo->vz[0]+topo->vz[3]);
//    topo->vz[nxvz-2] = 0.5*(topo->vz[nxvz]+topo->vz[nxvz-3]);
    
    
//    DiffuseAlongTopography( mesh, model, scaling, topo->vx, nx, 1, 2*model.dt );
//    DiffuseAlongTopography( mesh, model, scaling, topo->vz, nx+1, 1, 2*model.dt );

    //---------------------------
    double distance, dxm;
    int in;
    
    for( k=0; k<topo_chain->Nb_part; k++ ) {
        
        topo_chain->Vx[k] = 0;
        
        distance=fabs(topo_chain->x[k]-mesh->xg_coord[0]);
        in=ceil((distance/dx)) - 1;
        if (in<0) {
            in = 0;
        }
        if (in>nx-2) {
            in = nx-2;
        }
        dxm = topo_chain->x[k] - mesh->xg_coord[in];
        topo_chain->Vx[k] =  (1-dxm/dx)* topo->vx[in];
        topo_chain->Vx[k] +=  (dxm/dx)* topo->vx[in+1];
        
        //---------------------------
        
        topo_chain->Vz[k] = 0;
        
        distance=fabs(topo_chain->x[k]-mesh->xvz_coord[0]);
        in=ceil((distance/dx)) - 1;
        if (in<0) {
            in = 0;
        }
        if (in>nx-2) {
            in = (nx+1)-2;
        }
        dxm = topo_chain->x[k] - mesh->xvz_coord[in];
        topo_chain->Vz[k] =  (1-dxm/dx)* topo->vz[in];
        topo_chain->Vz[k] +=  (dxm/dx)* topo->vz[in+1];

    }
    

    
    //-----------
    // Interpolate velocities to cell centers
//    double Vxc[ncx*ncz], Vzc[ncx*ncz];
//
//    for( j=0; j<ncz; j++ ) {
//        for( i=0; i<ncx; i++ ) {
//
//            c1 = i + j*(nx);
//            c2 = i + j*(ncx);
//            c3 = i + j*(nxvz);
//
//            Vxc[c2] = 0.5 *( Vx2[c1] + Vx2[c1+1] );
//            Vzc[c2] = 0.5 *( Vz2[c3+1] + Vz2[c3+nxvz+1] );
//        }
//    }
//    Interp_Grid2P( *topo_chain, topo_chain->Vx,    mesh, Vxc, mesh->xc_coord[0],  mesh->zc_coord[0],  mesh->Nx[0]-1, mesh->Nz[0]-1, mesh->BCp.type[0] );
//    Interp_Grid2P( *topo_chain, topo_chain->Vz,    mesh, Vzc, mesh->xc_coord[0],  mesh->zc_coord[0],  mesh->Nx[0]-1, mesh->Nz[0]-1, mesh->BCp.type[0] );


}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DiffuseAlongTopography_BEN( grid *mesh, params model, scale scaling, double *array, int size, double diff, double diff_time ) {
    
    // Explicit diffusion solver;
    int i, j, c1, c2, c3, k, it;
    double dx = fabs(mesh->xg_coord[1]-mesh->xg_coord[0]);
    double dt=0.4*dx*dx/diff;
    double correct[size];
    int nstep = diff_time/dt + 1;
    printf("Do %d topographic diffusion steps\n", nstep);
    
    for (it=0; it<nstep; it++) {
        for (i=1; i<size-1; i++) {
            correct[i] = 0.5*dt/dx/dx*diff*(array[i-1]+array[i+1]-2*array[i]);
        }
        for (i=size-2; i==1; i--) {
            correct[i] += 0.5*dt/dx/dx*diff*(array[i-1]+array[i+1]-2*array[i]);
        }
        for (i=1; i<size-1; i++) {
            array[i] += correct[i];
        }
    }
    
}

// SIMPLE
void BuildInitialTopography_BEN( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {
    
    int k, fact=1;
    double dx = model.dx/fact;
    topo_chain->Nb_part = model.Nx*fact;
    double TopoLevel = 660e3/scaling.L;
    double L = model.xmax - model.xmin;
    
    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        
        topo_chain->x[k]     = model.xmin + k*dx;
        topo_chain->z[k]     = 0.0/scaling.L;
        
//        topo_chain->x0[k]    = topo_chain->x[k];
//        topo_chain->z0[k]    = topo_chain->z[k];
        topo_chain->phase[k] = 0;
    }
}
