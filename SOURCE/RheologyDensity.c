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
#include "time.h"
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

void AccumulatedStrainII( grid* mesh, scale scaling, params model, markers* particles, double* X_vect, double* Z_vect, int Nx, int Nz, char *tag ) {
    
    double *strain_inc;
    int k, l, c1;
    DoodzFP *strain_inc_el, *strain_inc_pl, *strain_inc_pwl, *strain_inc_exp, *strain_inc_lin, *strain_inc_gbs;
    
    printf("Accumulating strain\n");
    
    // Allocate
    strain_inc      = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_el   = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_pl   = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_pwl  = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_exp  = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_lin  = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_gbs  = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));

    // Interpolate exz to cell centers
#pragma omp parallel for shared( mesh, model, strain_inc, strain_inc_el, strain_inc_pl, strain_inc_pwl, strain_inc_exp, strain_inc_lin, strain_inc_gbs  ) private( c1 )
    for (c1=0; c1<(mesh->Nx-1)*(mesh->Nz-1); c1++) {

        strain_inc[c1]     = model.dt*(mesh->eII_pl[c1]+mesh->eII_pwl[c1]+mesh->eII_exp[c1]+mesh->eII_lin[c1]+mesh->eII_gbs[c1]+mesh->eII_el[c1]+mesh->eII_cst[c1]);
        if (strain_inc[c1]<0) {
            printf("negative strain increment\n");
            printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n" ,mesh->eII_pl[c1],mesh->eII_pwl[c1],mesh->eII_exp[c1],mesh->eII_lin[c1],mesh->eII_gbs[c1],mesh->eII_el[c1]);
            exit(0);
        }

        strain_inc_el[c1]  = model.dt*mesh->eII_el[c1];
        strain_inc_pl[c1]  = model.dt*mesh->eII_pl[c1];
        strain_inc_pwl[c1] = model.dt*mesh->eII_pwl[c1];
        strain_inc_exp[c1] = model.dt*mesh->eII_exp[c1];
        strain_inc_lin[c1] = model.dt*mesh->eII_lin[c1];
        strain_inc_gbs[c1] = model.dt*mesh->eII_gbs[c1];
    }

    //---------------------------------------------------------------//

    double dE_tot, dE_el, dE_pl, dE_pwl, dE_exp, dE_lin, dE_gbs;
    double dx, dz, dxm, dzm, dst, sumW;
    int    i_part, j_part, iSW, iNW, iSE, iNE;
    dx=mesh->dx;
    dz=mesh->dz;

//#pragma omp parallel for shared( particles , strain_inc, strain_inc_el, strain_inc_pl, strain_inc_pwl, strain_inc_exp, strain_inc_lin, strain_inc_gbs, tag  ) private( dE_tot, dE_el, dE_pl, dE_pwl, dE_exp, dE_lin, dE_gbs, sumW, dst, dxm, dzm, iSW, iSE, iNW, iNE, i_part, j_part  ) firstprivate (dx, dz, Nx, Nz)
    for (k=0;k<particles->Nb_part;k++) {

        dE_tot = 0.0;
        dE_el  = 0.0;
        dE_pl  = 0.0;
        dE_pwl = 0.0;
        dE_exp = 0.0;
        dE_lin = 0.0;
        dE_gbs = 0.0;


        if (particles->phase[k] != -1) {

            dst = (particles->x[k]-X_vect[0]);
            j_part   = ceil((dst/dx)) - 1;
            if (j_part<0) {
                j_part = 0;
            }
            if (j_part>Nx-2) {
                j_part = Nx-2;
            }

            // Get the line:
            dst = (particles->z[k]-Z_vect[0]);
            i_part   = ceil((dst/dz)) - 1;
            if (i_part<0) {
                i_part = 0;
            }
            if (i_part>Nz-2) {
                i_part = Nz-2;
            }

            dxm = (particles->x[k] - X_vect[j_part]);
            dzm = (particles->z[k] - Z_vect[i_part]);

            iSW = j_part+i_part*Nx;
            iSE = j_part+i_part*Nx+1;
            iNW = j_part+(i_part+1)*Nx;
            iNE = j_part+(i_part+1)*Nx+1;

            dE_tot = 0.0;
            dE_el  = 0.0;
            dE_pl  = 0.0;
            dE_pwl = 0.0;
            dE_exp = 0.0;
            dE_lin = 0.0;
            dE_gbs = 0.0;

            sumW         = 0.0;

//            if (j_part>0 && j_part<Nx-2 && i_part>0 && i_part<Nz-2) {

                if (tag[iSW]!=30 && tag[iSW]!=31) {
                    dE_tot +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc[iSW];
                    dE_el  +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_el[iSW];
                    dE_pl  +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_pl[iSW];
                    dE_pwl +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_pwl[iSW];
                    dE_exp +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_exp[iSW];
                    dE_lin +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_lin[iSW];
                    dE_gbs +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_gbs[iSW];
                    sumW   +=  (1.0-dxm/dx) * (1.0-dzm/dz);
                }
                if (tag[iSE]!=30 && tag[iSE]!=31) {
                    dE_tot +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc[iSE];
                    dE_el  +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_el[iSE];
                    dE_pl  +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_pl[iSE];
                    dE_pwl +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_pwl[iSE];
                    dE_exp +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_exp[iSE];
                    dE_lin +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_lin[iSE];
                    dE_gbs +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_gbs[iSE];
                    sumW   +=  (dxm/dx) * (1.0-dzm/dz);
                }
                if (tag[iNW]!=30 && tag[iNW]!=31) {
                    dE_tot +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc[iNW];
                    dE_el  +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_el[iNW];
                    dE_pl  +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_pl[iNW];
                    dE_pwl +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_pwl[iNW];
                    dE_exp +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_exp[iNW];
                    dE_lin +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_lin[iNW];
                    dE_gbs +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_gbs[iNW];
                    sumW   +=  (1.0-dxm/dx) * (dzm/dz);
                }
                if (tag[iNE]!=30 && tag[iNE]!=31) {
                    dE_tot +=  (dxm/dx) * (dzm/dz) * strain_inc[iNE];
                    dE_el  +=  (dxm/dx) * (dzm/dz) * strain_inc_el[iNE];
                    dE_pl  +=  (dxm/dx) * (dzm/dz) * strain_inc_pl[iNE];
                    dE_pwl +=  (dxm/dx) * (dzm/dz) * strain_inc_pwl[iNE];
                    dE_exp +=  (dxm/dx) * (dzm/dz) * strain_inc_exp[iNE];
                    dE_lin +=  (dxm/dx) * (dzm/dz) * strain_inc_lin[iNE];
                    dE_gbs +=  (dxm/dx) * (dzm/dz) * strain_inc_gbs[iNE];
                    sumW += (dxm/dx)* (dzm/dz);
                }
        
//            printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n", sumW, dE_tot, strain_inc[iSW], strain_inc[iSE], strain_inc[iNW], strain_inc[iNE]);
//
//             if (dE_tot<0.0)exit(1);
//            }

//            if (j_part==0    && i_part==0   ) {
//                dE_tot = strain_inc[iNE];
//                dE_el  = strain_inc[iNE];
//                dE_pl  = strain_inc[iNE];
//                dE_pwl = strain_inc[iNE];
//                dE_exp = strain_inc[iNE];
//                dE_lin = strain_inc[iNE];
//                dE_gbs = strain_inc[iNE];
//            }
//            if (j_part==Nx-2 && i_part==0   ) {
//                dE_tot = strain_inc[iNW];
//                dE_el  = strain_inc[iNW];
//                dE_pl  = strain_inc[iNW];
//                dE_pwl = strain_inc[iNW];
//                dE_exp = strain_inc[iNW];
//                dE_lin = strain_inc[iNW];
//                dE_gbs = strain_inc[iNW];
//            }
//            if (j_part==0    && i_part==Nz-2) {
//                dE_tot = strain_inc[iSE];
//                dE_el  = strain_inc[iSE];
//                dE_pl  = strain_inc[iSE];
//                dE_pwl = strain_inc[iSE];
//                dE_exp = strain_inc[iSE];
//                dE_lin = strain_inc[iSE];
//                dE_gbs = strain_inc[iSE];
//            }
//            if (j_part==Nx-2 && i_part==Nz-2) {
//                dE_tot = strain_inc[iSW];
//                dE_el  = strain_inc[iSW];
//                dE_pl  = strain_inc[iSW];
//                dE_pwl = strain_inc[iSW];
//                dE_exp = strain_inc[iSW];
//                dE_lin = strain_inc[iSW];
//                dE_gbs = strain_inc[iSW];
//            }

//            if (j_part==0 && i_part>0 && i_part<Nz-2) {
//
//                if (tag[iSE]!=30 && tag[iSE]!=31) {
//                    dE_tot += (1.0-dzm/dz)   * strain_inc[iSE];
//                    dE_el  += (1.0-dzm/dz)   * strain_inc_el[iSE];
//                    dE_pl  += (1.0-dzm/dz)   * strain_inc_pl[iSE];
//                    dE_pwl += (1.0-dzm/dz)   * strain_inc_pwl[iSE];
//                    dE_exp += (1.0-dzm/dz)   * strain_inc_exp[iSE];
//                    dE_lin += (1.0-dzm/dz)   * strain_inc_lin[iSE];
//                    dE_gbs += (1.0-dzm/dz)   * strain_inc_gbs[iSE];
//                    sumW +=  (1.0-dzm/dz);
//                }
//                if (tag[iNE]!=30 && tag[iNE]!=31) {
//                    dE_tot += (0.0-dzm/dz)   * strain_inc[iNE];
//                    dE_el  += (0.0-dzm/dz)   * strain_inc_el[iNE];
//                    dE_pl  += (0.0-dzm/dz)   * strain_inc_pl[iNE];
//                    dE_pwl += (0.0-dzm/dz)   * strain_inc_pwl[iNE];
//                    dE_exp += (0.0-dzm/dz)   * strain_inc_exp[iNE];
//                    dE_lin += (0.0-dzm/dz)   * strain_inc_lin[iNE];
//                    dE_gbs += (0.0-dzm/dz)   * strain_inc_gbs[iNE];
//                    sumW +=  (dzm/dz);
//                }
//                if(sumW>1e-13) {
//                    dE_tot /= sumW;
//                    dE_el  /= sumW;
//                    dE_pl  /= sumW;
//                    dE_pwl /= sumW;
//                    dE_exp /= sumW;
//                    dE_lin /= sumW;
//                    dE_gbs /= sumW;
//                }
//            }
//
//            if (j_part==Nx-2 && i_part>0 && i_part<Nz-2) {
//
//                if (tag[iSW]!=30 && tag[iSW]!=31) {
//                    dE_tot += (1.0-dzm/dz)   * strain_inc[iSW];
//                    dE_el  += (1.0-dzm/dz)   * strain_inc_el[iSW];
//                    dE_pl  += (1.0-dzm/dz)   * strain_inc_pl[iSW];
//                    dE_pwl += (1.0-dzm/dz)   * strain_inc_pwl[iSW];
//                    dE_exp += (1.0-dzm/dz)   * strain_inc_exp[iSW];
//                    dE_lin += (1.0-dzm/dz)   * strain_inc_lin[iSW];
//                    dE_gbs += (1.0-dzm/dz)   * strain_inc_gbs[iSW];
//                    sumW +=  (1.0-dzm/dz);
//                }
//                if (tag[iNW]!=30 && tag[iNW]!=31) {
//                    dE_tot += (0.0-dzm/dz)   * strain_inc[iNW];
//                    dE_el  += (0.0-dzm/dz)   * strain_inc_el[iNW];
//                    dE_pl  += (0.0-dzm/dz)   * strain_inc_pl[iNW];
//                    dE_pwl += (0.0-dzm/dz)   * strain_inc_pwl[iNW];
//                    dE_exp += (0.0-dzm/dz)   * strain_inc_exp[iNW];
//                    dE_lin += (0.0-dzm/dz)   * strain_inc_lin[iNW];
//                    dE_gbs += (0.0-dzm/dz)   * strain_inc_gbs[iNW];
//                    sumW +=  (dzm/dz);
//                }
//                if(sumW>1e-13) {
//                    dE_tot /= sumW;
//                    dE_el  /= sumW;
//                    dE_pl  /= sumW;
//                    dE_pwl /= sumW;
//                    dE_exp /= sumW;
//                    dE_lin /= sumW;
//                    dE_gbs /= sumW;
//                }
//            }

//            if (i_part==0 && j_part>0 && j_part<Nx-2) {
//                if (tag[iNW]!=30 && tag[iNW]!=31) {
//                    dE_tot += (1.0-dxm/dx)   * strain_inc[iNW];
//                    dE_el  += (1.0-dxm/dx)   * strain_inc_el[iNW];
//                    dE_pl  += (1.0-dxm/dx)   * strain_inc_pl[iNW];
//                    dE_pwl += (1.0-dxm/dx)   * strain_inc_pwl[iNW];
//                    dE_exp += (1.0-dxm/dx)   * strain_inc_exp[iNW];
//                    dE_lin += (1.0-dxm/dx)   * strain_inc_lin[iNW];
//                    dE_gbs += (1.0-dxm/dx)   * strain_inc_gbs[iNW];
//                    sumW += (1.0-dxm/dx);
//                }
//                if (tag[iNE]!=30 && tag[iNE]!=31) {
//                    dE_tot += (dxm/dx)   * strain_inc[iNE];
//                    dE_el  += (dxm/dx)   * strain_inc_el[iNE];
//                    dE_pl  += (dxm/dx)   * strain_inc_pl[iNE];
//                    dE_pwl += (dxm/dx)   * strain_inc_pwl[iNE];
//                    dE_exp += (dxm/dx)   * strain_inc_exp[iNE];
//                    dE_lin += (dxm/dx)   * strain_inc_lin[iNE];
//                    dE_gbs += (dxm/dx)   * strain_inc_gbs[iNE];
//                    sumW += (dxm/dx);
//                }
//                if(sumW>1e-13) {
//                    dE_tot /= sumW;
//                    dE_el  /= sumW;
//                    dE_pl  /= sumW;
//                    dE_pwl /= sumW;
//                    dE_exp /= sumW;
//                    dE_lin /= sumW;
//                    dE_gbs /= sumW;
//                }
//            }
//            if (i_part==Nz-2 && j_part>0 && j_part<Nx-2) {
//                if (tag[iSW]!=30 && tag[iSW]!=31) {
//                    dE_tot += (1.0-dxm/dx)   * strain_inc[iSW];
//                    dE_el  += (1.0-dxm/dx)   * strain_inc_el[iSW];
//                    dE_pl  += (1.0-dxm/dx)   * strain_inc_pl[iSW];
//                    dE_pwl += (1.0-dxm/dx)   * strain_inc_pwl[iSW];
//                    dE_exp += (1.0-dxm/dx)   * strain_inc_exp[iSW];
//                    dE_lin += (1.0-dxm/dx)   * strain_inc_lin[iSW];
//                    dE_gbs += (1.0-dxm/dx)   * strain_inc_gbs[iSW];
//                    sumW += (1.0-dxm/dx);
//                }
//                if ( tag[iSE] != 30 && tag[iSE] != 31 ) {
//                    dE_tot += (dxm/dx)   * strain_inc[iSE];
//                    dE_el  += (dxm/dx)   * strain_inc_el[iSE];
//                    dE_pl  += (dxm/dx)   * strain_inc_pl[iSE];
//                    dE_pwl += (dxm/dx)   * strain_inc_pwl[iSE];
//                    dE_exp += (dxm/dx)   * strain_inc_exp[iSE];
//                    dE_lin += (dxm/dx)   * strain_inc_lin[iSE];
//                    dE_gbs += (dxm/dx)   * strain_inc_gbs[iSE];
//                    sumW += (dxm/dx);
//                }
//                if(sumW>1e-13) {
//                    dE_tot /= sumW;
//                    dE_el  /= sumW;
//                    dE_pl  /= sumW;
//                    dE_pwl /= sumW;
//                    dE_exp /= sumW;
//                    dE_lin /= sumW;
//                    dE_gbs /= sumW;
//                }
//            }

            if(sumW>1e-13) {
                dE_tot /= sumW;
                dE_el  /= sumW;
                dE_pl  /= sumW;
                dE_pwl /= sumW;
                dE_exp /= sumW;
                dE_lin /= sumW;
                dE_gbs /= sumW;
            }


        }

        particles->strain[k]      += dE_tot;
        particles->strain_el[k]   += dE_el;
        particles->strain_pl[k]   += dE_pl;
        particles->strain_pwl[k]  += dE_pwl;
        particles->strain_exp[k]  += dE_exp;
        particles->strain_lin[k]  += dE_lin;
        particles->strain_gbs[k]  += dE_gbs;
    }

    //---------------------------------------------------------------//

    // Free
    DoodzFree(strain_inc);
    DoodzFree(strain_inc_el);
    DoodzFree(strain_inc_pl);
    DoodzFree(strain_inc_pwl);
    DoodzFree(strain_inc_exp);
    DoodzFree(strain_inc_lin);
    DoodzFree(strain_inc_gbs);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DeformationGradient ( grid mesh, scale scaling, params model, markers *particles ) {
    
    int k, l, cp, cu, cv, Nx, Nz;
    double dx, dz;
    double fxx, fxz, fzx, fzz, fxx_o, fxz_o, fzx_o, fzz_o;
    
    Nx = mesh.Nx;
    Nz = mesh.Nz;
    dx = mesh.dx;
    dz = mesh.dz;
    
    double *dudx, *dudz, *dvdx, *dvdz;
    DoodzFP *pdudx, *pdudz, *pdvdx, *pdvdz;
    
    dudx   = DoodzMalloc ((Nx-1)*(Nz-1)*sizeof(double));
    dvdz   = DoodzMalloc ((Nx-1)*(Nz-1)*sizeof(double));
    dudz   = DoodzMalloc ((Nx)*(Nz)*sizeof(double));
    dvdx   = DoodzMalloc ((Nx)*(Nz)*sizeof(double));
    pdudx  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
    pdudz  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
    pdvdx  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
    pdvdz  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
    
    // Compute dudx and dvdz (cell centers)
    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            cp = k  + l*(Nx-1);
            cu = k  + l*(Nx) + Nx;
            cv = k  + l*(Nx+1) + 1;
            dudx[cp] = 1.0/dx * ( mesh.u_in[cu+1]      - mesh.u_in[cu] );
            dvdz[cp] = 1.0/dz * ( mesh.v_in[cv+(Nx+1)] - mesh.v_in[cv] );
        }
    }
    
    // Compute dudx and dvdz (cell vertices)
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz; l++) {
            cp = k  + l*(Nx-1);
            cu = k  + l*(Nx);
            cv = k  + l*(Nx+1);
            dudz[cu] = 1.0/dz * ( mesh.u_in[cu+Nx] - mesh.u_in[cu] );
            dvdx[cu] = 1.0/dx * ( mesh.v_in[cv+1]  - mesh.v_in[cv] );
        }
    }
    
    // Interpolate from grid to particles
    Interp_Grid2P( *(particles), pdudx, &mesh, dudx, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type  );
    Interp_Grid2P( *(particles), pdvdz, &mesh, dvdz, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
    Interp_Grid2P( *(particles), pdvdx, &mesh, dvdx, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
    Interp_Grid2P( *(particles), pdudz, &mesh, dudz, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
    
#pragma omp parallel for shared ( particles ) private ( k, fxx, fxz, fzx, fzz, fxx_o, fxz_o, fzx_o, fzz_o ) firstprivate( model ) schedule( static )
    for(k=0; k<particles->Nb_part; k++) {
        
        fxx = 1.0 + pdudx[k]*model.dt;
        fxz = pdudz[k]*model.dt;
        fzz = 1.0 + pdvdz[k]*model.dt;
        fzx = pdvdx[k]*model.dt;
        
        if (model.step==1) {
            fxx_o = 1.0;
            fxz_o = 0.0;
            fzx_o = 0.0;
            fzz_o = 1.0;
            
        }
        else {
            fxx_o = particles->Fxx[k];
            fxz_o = particles->Fxz[k];
            fzx_o = particles->Fzx[k];
            fzz_o = particles->Fzz[k];
        }
        particles->Fxx[k] = fxx*fxx_o + fxz*fzx_o;
        particles->Fxz[k] = fxx*fxz_o + fxz*fzz_o;
        particles->Fzx[k] = fzx*fxx_o + fzz*fzx_o;
        particles->Fzz[k] = fzx*fxz_o + fzz*fzz_o;
    }
    
    printf("Updated deformation gradient tensor\n");
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AccumulatedStrain( grid* mesh, scale scaling, params model, markers* particles ) {
    
    double *strain_inc;
    int Nx, Nz, k, l, c1;
    DoodzFP *strain_inc_mark, *strain_inc_el, *strain_inc_pl, *strain_inc_pwl, *strain_inc_exp, *strain_inc_lin, *strain_inc_gbs;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    
    printf("Accumulating strain\n");
    
    // Allocate
    //    exz_n           = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc      = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_el   = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_pl   = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_pwl  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_exp  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_lin  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_gbs  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_mark = DoodzCalloc(sizeof(DoodzFP),particles->Nb_part);
    
    // Interpolate exz to cell centers
    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            c1 = k  + l*(Nx-1);
            
            strain_inc[c1]     = model.dt*(mesh->eII_pl[c1]+mesh->eII_pwl[c1]+mesh->eII_exp[c1]+mesh->eII_lin[c1]+mesh->eII_gbs[c1]+mesh->eII_el[c1]+mesh->eII_cst[c1]);
            if (strain_inc[c1]<0) {
                printf("negative strain increment\n");
                printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n" ,mesh->eII_pl[c1],mesh->eII_pwl[c1],mesh->eII_exp[c1],mesh->eII_lin[c1],mesh->eII_gbs[c1],mesh->eII_el[c1]);
                exit(0);
            }
            
            strain_inc_el[c1]  = model.dt*mesh->eII_el[c1];
            strain_inc_pl[c1]  = model.dt*mesh->eII_pl[c1];
            strain_inc_pwl[c1] = model.dt*mesh->eII_pwl[c1];
            strain_inc_exp[c1] = model.dt*mesh->eII_exp[c1];
            strain_inc_lin[c1] = model.dt*mesh->eII_lin[c1];
            strain_inc_gbs[c1] = model.dt*mesh->eII_gbs[c1];
        }
    }
    
    // Interp increments to particles
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_el, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_el, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_pl, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_pl, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_pwl, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_pwl, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_exp, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_exp, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_lin, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_lin, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_gbs, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_gbs, strain_inc_mark, particles->Nb_part );
    
    // Free
    //    DoodzFree(exz_n);
    DoodzFree(strain_inc);
    DoodzFree(strain_inc_el);
    DoodzFree(strain_inc_pl);
    DoodzFree(strain_inc_pwl);
    DoodzFree(strain_inc_exp);
    DoodzFree(strain_inc_lin);
    DoodzFree(strain_inc_gbs);
    DoodzFree(strain_inc_mark);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void RotateStresses( grid mesh, markers* particles, params model, scale *scaling ) {
    
    int k, l, cp, cu, cv, Nx, Nz,type=0;
    double dx, dz;
    double angle, sxxr, sxzr;
    
    Nx = mesh.Nx;
    Nz = mesh.Nz;
    dx = mesh.dx;
    dz = mesh.dz;
    
    if (type==0) {
        // JAUMANN RATE
        
        // Stress correction
#pragma omp parallel for shared ( particles ) private ( angle, k, sxxr, sxzr ) firstprivate( model ) schedule( static )
        for(k=0; k<particles->Nb_part; k++) {
            
            // Filter out particles that are inactive (out of the box)
            if (particles->phase[k] != -1) {
                
                // Angle
                angle = model.dt*particles->om_p[k];
                
                // Re-belotte:
                sxxr =   ( cos(angle)*particles->sxxd[k] - sin(angle)*particles->sxz[k] ) * cos(angle) + ( particles->sxz[k]*cos(angle) +  sin(angle)*particles->sxxd[k])*sin(angle) ;
                sxzr = - ( cos(angle)*particles->sxxd[k] - sin(angle)*particles->sxz[k] ) * sin(angle) + ( particles->sxz[k]*cos(angle) +  sin(angle)*particles->sxxd[k])*cos(angle) ;
                
                particles->sxxd[k] = sxxr;
                particles->sxz[k]  = sxzr;
            }
        }
    }
    else {
        
        // UPPER CONVECTED
        
        double *dudx, *dudz, *dvdx, *dvdz;
        DoodzFP *pdudx, *pdudz, *pdvdx, *pdvdz, *VEm;
        
        dudx   = DoodzMalloc ((Nx-1)*(Nz-1)*sizeof(double));
        dvdz   = DoodzMalloc ((Nx-1)*(Nz-1)*sizeof(double));
        dudz   = DoodzMalloc ((Nx)*(Nz)*sizeof(double));
        dvdx   = DoodzMalloc ((Nx)*(Nz)*sizeof(double));
        pdudx  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
        pdudz  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
        pdvdx  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
        pdvdz  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
        VEm    = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
        
        // Compute dudx and dvdz (cell centers)
        for (k=0; k<Nx-1; k++) {
            for (l=0; l<Nz-1; l++) {
                cp = k  + l*(Nx-1);
                cu = k  + l*(Nx) + Nx;
                cv = k  + l*(Nx+1) + 1;
                dudx[cp] = 1.0/dx * ( mesh.u_in[cu+1]      - mesh.u_in[cu] );
                dvdz[cp] = 1.0/dz * ( mesh.v_in[cv+(Nx+1)] - mesh.v_in[cv] );
            }
        }
        
        // Compute dudx and dvdz (cell vertices)
        for (k=0; k<Nx; k++) {
            for (l=0; l<Nz; l++) {
                cp = k  + l*(Nx-1);
                cu = k  + l*(Nx);
                cv = k  + l*(Nx+1);
                dudz[cu] = 1.0/dz * ( mesh.u_in[cu+Nx] - mesh.u_in[cu] );
                dvdx[cu] = 1.0/dx * ( mesh.v_in[cv+1]  - mesh.v_in[cv] );
            }
        }
        
        // Interpolate from grid to particles
        Interp_Grid2P( *(particles), pdudx, &mesh, dudx, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type  );
        Interp_Grid2P( *(particles), pdvdz, &mesh, dvdz, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
        Interp_Grid2P( *(particles), pdvdx, &mesh, dvdx, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
        Interp_Grid2P( *(particles), pdudz, &mesh, dudz, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
        Interp_Grid2P( *(particles), VEm, &mesh, mesh.VE_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type  );
        
        
        // Stress correction
#pragma omp parallel for shared ( particles, model, VEm, pdudx, pdudz, pdvdx, pdvdz, sxxr, sxzr ) private ( k ) schedule( static )
        for(k=0; k<particles->Nb_part; k++) {
            sxxr = model.dt * VEm[k] * ( -2.0*particles->sxxd[k]*pdudx[k] - 2.0*particles->sxz[k]*pdudz[k]);
            sxzr = model.dt * VEm[k] * (    particles->sxxd[k]*pdudz[k] -   particles->sxxd[k]*pdvdx[k]);
            
            particles->sxxd[k] -= sxxr;
            particles->sxz[k]  -= sxzr;
            
            // Only if explicity
            if (model.subgrid_diff==4) particles->sxxd[k] *= VEm[k];
            if (model.subgrid_diff==4) particles->sxz[k]  *= VEm[k];
        }
        MinMaxArray( pdudx,  scaling->E, particles->Nb_part, "dudx p." );
        MinMaxArray( pdvdx,  scaling->E, particles->Nb_part, "dvdx p." );
        MinMaxArray( pdudz,  scaling->E, particles->Nb_part, "dudz p." );
        
        
        // clean memory
        DoodzFree(dudx);
        DoodzFree(dudz);
        DoodzFree(dvdx);
        DoodzFree(dvdz);
        DoodzFree(pdudx);
        DoodzFree(pdudz);
        DoodzFree(pdvdx);
        DoodzFree(pdvdz);
        DoodzFree(VEm);
    }
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateMaxPT ( scale scaling, params model, markers *particles ) {
    
    // Update maximum pressure and temperature on markers
    int k;
    
#pragma omp parallel for shared ( particles ) private ( k ) firstprivate( model ) schedule( static )
    for(k=0; k<particles->Nb_part; k++) {
        
        if (particles->T[k] > particles->Tmax[k] ) particles->Tmax[k] = particles->T[k];
        if (particles->P[k] > particles->Pmax[k] ) particles->Pmax[k] = particles->P[k];
        
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleGrainSize( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
    DoodzFP *d_inc_mark, *d_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    d_inc_mark = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
    d_inc_grid = DoodzCalloc(sizeof(DoodzFP), Ncx*Ncz);
    
    for (k=0;k<Ncx*Ncz;k++) {
        d_inc_grid[k] =0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) d_inc_grid[k] = mesh->d[k] - mesh->d0[k];
    }
    
    Interp_Grid2P( *particles, particles->d, mesh, mesh->d, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCt.type  );
    
#pragma omp parallel for shared ( particles )
    for (k=0;k<particles->Nb_part;k++) {
        if (particles->d[k]<1.0e-16) particles->d[k] = 1.0e-16;
    }
    
    DoodzFree(d_inc_mark);
    DoodzFree(d_inc_grid);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleEnergy( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
    DoodzFP *T_inc_mark, *Tm0, dtm, *dTms, *dTgr, *dTmr;
    double *Tg0, *dTgs, dx=model.dx, dz=model.dz, d=1.0;
    int Nx, Nz, Ncx, Ncz, k, c0, p;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    if ( model.subgrid_diff >= 1 ) {
        
        printf("Subgrid diffusion for temperature update\n");
        Tg0  = DoodzCalloc(sizeof(DoodzFP), Ncx*Ncz);
        dTgs = DoodzCalloc(sizeof(DoodzFP), Ncx*Ncz);
        dTgr = DoodzCalloc(sizeof(DoodzFP), Ncx*Ncz);
        Tm0  = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
        dTms = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
        dTmr = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
        
        // Old temperature grid
#pragma omp parallel for shared(mesh, Tg0) private(c0) firstprivate(Ncx,Ncz)
        for ( c0=0; c0<Ncx*Ncz; c0++ ) {
            if (mesh->BCt.type[c0] != 30) Tg0[c0] = mesh->T[c0] - mesh->dT[c0];
        }
        
        // Old temperature grid --> markers
        Interp_Grid2P( *particles, Tm0, mesh, Tg0, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCt.type  );
        
        // Compute subgrid temperature increments on markers
#pragma omp parallel for shared(particles,Tm0,dTms) private(k,p,dtm) firstprivate(materials,dx,dz,model,d)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) {
                p = particles->phase[k];
                dtm     = materials->Cv[p]* particles->rho[k]/ (materials->k[p] * (1.0/dx/dx + 1.0/dz/dz));
                dTms[k] = -( particles->T[k] - Tm0[k] ) * (1.0-exp(-d*model.dt/dtm));
            }
        }
        
        // Subgrid temperature increments markers --> grid
        Interp_P2C ( *particles, dTms, mesh, dTgs, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        
        // Remaining temperature increments on the grid
#pragma omp parallel for shared(mesh, dTgs, dTgr) private(c0) firstprivate(Ncx,Ncz)
        for ( c0=0; c0<Ncx*Ncz; c0++ ) {
            if (mesh->BCt.type[c0] != 30) dTgr[c0] = mesh->dT[c0] - dTgs[c0];
        }
        
        // Remaining temperature increments grid --> markers
        Interp_Grid2P( *particles, dTmr, mesh, dTgr, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCt.type  );
        
        // Final temperature update on markers
#pragma omp parallel for shared(particles,dTms,dTmr) private(k)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) particles->T[k] += dTms[k] + dTmr[k];
        }
        DoodzFree(Tg0);
        DoodzFree(Tm0);
        DoodzFree(dTms);
        DoodzFree(dTmr);
        DoodzFree(dTgs);
        DoodzFree(dTgr);
    }
    else {
        
        T_inc_mark = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
        
        // Interp increments to particles
        Interp_Grid2P( *particles, T_inc_mark, mesh, mesh->dT, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCt.type  );
        
        // Increment temperature on particles
        ArrayPlusArray( particles->T, T_inc_mark, particles->Nb_part );
        
        DoodzFree(T_inc_mark);
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticlePressure( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
    DoodzFP *P_inc_mark;
    int Nx, Nz, Ncx, Ncz, k, c0, p;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    P_inc_mark = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
    
    // Interp increments to particles
    Interp_Grid2P( *particles, P_inc_mark, mesh, mesh->dp, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type  );
    
    // Increment pressure on particles
    ArrayPlusArray( particles->P, P_inc_mark, particles->Nb_part );
    
    DoodzFree(P_inc_mark);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleStress( grid* mesh, markers* particles, params* model, mat_prop* materials, scale *scaling ) {
    
    int k, l, c0, c1, Nx, Nz, Ncx, Ncz, p;
    DoodzFP *mdsxxd, *mdsxz,*dsxxd, *dsxz, d=1.0, dtaum;
    DoodzFP *dtxxgs, *dtxzgs, *dtxxgr, *dtxzgr, *txxm0, *txzm0, *dtxxms, *dtxzms, *dtxxmr, *dtxzmr, *etam;
    
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    if ( model->subgrid_diff > -1 ) {
        
        // Alloc
        dtxxgs = DoodzCalloc(sizeof(DoodzFP), Ncx*Ncz           );
        dtxzgs = DoodzCalloc(sizeof(DoodzFP), Nx*Nz             );
        dtxxgr = DoodzCalloc(sizeof(DoodzFP), Ncx*Ncz           );
        dtxzgr = DoodzCalloc(sizeof(DoodzFP), Nx*Nz             );
        txxm0  = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
        txzm0  = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
        dtxxms = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
        dtxzms = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
        dtxxmr = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
        dtxzmr = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
        etam   = DoodzCalloc(sizeof(DoodzFP), particles->Nb_part);
        
        // Old stresses grid --> markers
        Interp_Grid2P( *particles, txxm0, mesh, mesh->sxxd0, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
        Interp_Grid2P( *particles, txzm0, mesh, mesh->sxz0 , mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type    );
        
        if ( model->subgrid_diff == 2 ) {
            
            printf("Subgrid diffusion for stress tensor component update\n");
            
            // Compute subgrid stress increments on markers
#pragma omp parallel for shared(particles,txxm0,txzm0,dtxxms,dtxzms,etam) private(k,p,dtaum) firstprivate(materials,model,d)
            for ( k=0; k<particles->Nb_part; k++ ) {
                if (particles->phase[k] != -1) {
                    p         = particles->phase[k];
                    dtaum     = etam[k]/ (materials->mu[p]);
                    dtxxms[k] = -( particles->sxxd[k] - txxm0[k]) * (1.0 - exp(-d*model->dt/dtaum));
                    dtxzms[k] = -( particles->sxz[k]  - txzm0[k]) * (1.0 - exp(-d*model->dt/dtaum));
                }
            }
            
            // Subgrid stress increments markers --> grid
            Interp_P2C ( *particles, dtxxms, mesh, dtxxgs, mesh->xg_coord, mesh->zg_coord, 1, 0 );
            Interp_P2N ( *particles, dtxzms, mesh, dtxzgs, mesh->xg_coord, mesh->zg_coord, 1, 0, model );
            
            // Remaining stress increments on the grid
#pragma omp parallel for shared(mesh,dtxxgs,dtxxgr) private(c0) firstprivate(Ncx,Ncz)
            for ( c0=0; c0<Ncx*Ncz; c0++ ) {
                if (mesh->BCp.type[c0]!=30) dtxxgr[c0] = (mesh->sxxd[c0]-mesh->sxxd0[c0]) - dtxxgs[c0];
            }
#pragma omp parallel for shared(mesh,dtxzgs,dtxzgr) private(c0) firstprivate(Nx,Nz)
            for ( c0=0; c0<Nx*Nz; c0++ ) {
                if (mesh->BCg.type[c0]   !=30) dtxzgr[c0] = (mesh->sxz[c0]-mesh->sxz0[c0]) - dtxzgs[c0];
            }
            
            // Remaining stress increments grid --> markers
            Interp_Grid2P( *particles, dtxxmr, mesh, dtxxgr, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type  );
            Interp_Grid2P( *particles, dtxzmr, mesh, dtxzgr, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type    );
            
            // Final stresses update on markers
#pragma omp parallel for shared(particles,dtxxms,dtxzms,dtxxmr,dtxzmr) private(k)
            for ( k=0; k<particles->Nb_part; k++ ) {
                if (particles->phase[k] != -1) particles->sxxd[k] += dtxxms[k] + dtxxmr[k];
                if (particles->phase[k] != -1) particles->sxz[k]  += dtxzms[k] + dtxzmr[k];
            }
        }
        
        // Free
        DoodzFree( dtxxgs );
        DoodzFree( dtxzgs );
        DoodzFree( dtxxgr );
        DoodzFree( dtxzgr );
        DoodzFree( txxm0  );
        DoodzFree( txzm0  );
        DoodzFree( dtxxms );
        DoodzFree( dtxzms );
        DoodzFree( dtxxmr );
        DoodzFree( dtxzmr );
        DoodzFree( etam   );
        
    }
    if (model->subgrid_diff==0 || model->subgrid_diff==1 || model->subgrid_diff==4){
        
        printf("No subgrid diffusion for stress tensor component update\n");
        
        // Alloc
        dsxxd  = DoodzCalloc((Nx-1)*(Nz-1),sizeof(double));
        dsxz   = DoodzCalloc((Nx)*(Nz),sizeof(double));
        mdsxxd = DoodzCalloc(particles->Nb_part,sizeof(DoodzFP));
        mdsxz  = DoodzCalloc(particles->Nb_part,sizeof(DoodzFP));
        
        // Cell: normal stress change
        for (k=0; k<Nx-1; k++) {
            for (l=0; l<Nz-1; l++) {
                c0 = k  + l*(Nx-1);
                if (mesh->BCp.type[c0] !=30 && mesh->BCp.type[c0] !=31) dsxxd[c0] = mesh->sxxd[c0] - mesh->sxxd0[c0];
            }
        }
        
        // Vertex: shear stress change
        for (k=0; k<Nx; k++) {
            for (l=0; l<Nz; l++) {
                c1 = k  + l*(Nx);
                if (mesh->BCg.type[c1] !=30 ) dsxz[c1] =  mesh->sxz[c1] - mesh->sxz0[c1];
            }
        }
        
        // Interpolate stress changes to markers
        Interp_Grid2P( *particles, mdsxxd, mesh, dsxxd, mesh->xc_coord,  mesh->zc_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type  );
        Interp_Grid2P( *particles, mdsxz,  mesh, dsxz,  mesh->xg_coord,  mesh->zg_coord,  mesh->Nx,   mesh->Nz, mesh->BCg.type    );
        
        // Update marker stresses
        ArrayPlusArray( particles->sxxd, mdsxxd, particles->Nb_part );
        ArrayPlusArray( particles->sxz,  mdsxz,  particles->Nb_part );
        
        // Free
        DoodzFree(dsxxd);
        DoodzFree(dsxz);
        DoodzFree(mdsxxd);
        DoodzFree(mdsxz);
    }
}



/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double Viscosity( int phase, double G, double T, double P, double d, double phi, double Eii_in, double tII0_in, double exx, double exz, double txx0, double txz0, mat_prop* materials, params *model, scale *scaling, int strain_rate_formulation, double *txxn, double *txzn, double* etaVE, double* VEcoeff, double* Eii_el, double* Eii_pl, double* Eii_pwl, double* Eii_exp , double* Eii_lin, double* Eii_gbs, double* Eii_cst, double* Exx_pwl, double* Exz_pwl, double* Exx_el, double* Exz_el, double* Exx_diss, double* Exz_diss, double* Exx_pl, double* Exz_pl, double *d1, double* A2_pwl, double strain_acc, double Phi, double C ) {
    
    //    C=materials->C[phase], Phi=materials->phi[phase],  G=materials->mu[phase];
    
    // General paramaters
    double eta=0.0, R=materials->R, dt=model->dt;
    double minEta=model->mineta, maxEta=model->maxeta;
    double TmaxPeierls = (1200.0+zeroC)/scaling->T;      // max. T for Peierls
    
    // Parameters for deformation map calculations
    int    local_iter = model->loc_iter, it, nitmax = 100, noisy=0;
    int    constant=0, dislocation=0, peierls=0, diffusion=0, gbs=0, elastic = model->iselastic;
    double tol = 1.0e-11, res=0.0, dfdeta=0.0, txx1=0.0, txz1=0.0, tII1=0.0, ieta_sum=0.0, tII0 = sqrt(txx0*txx0 + txz0*txz0);
    double eta_up=0.0, eta_lo=0.0, eta_0=0.0, eta_p=0.0, r_eta_pl=0.0, r_eta_0=0.0, r_eta_p=0.0;
    double eta_pwl=0.0, eta_exp=0.0, eta_pl=0.0, eta_lin=0.0, eta_el=0.0, eta_gbs=0.0, eta_cst=0.0, eta_step=0.0;
    double Exx=0.0, Exz=0.0, Eii_vis=0.0, Eii= Eii_in;
    
    // Flow law parameters from input file
    double MC_yield=0.0, yield=0.0;
    double Eapwl = materials->Qpwl[phase], Vapwl = materials->Vpwl[phase], npwl  = materials->npwl[phase], mpwl  = materials->mpwl[phase], rpwl  = materials->rpwl[phase], Apwl  = materials->Apwl[phase], fpwl = materials->fpwl[phase], apwl = materials->apwl[phase], Fpwl  = materials->Fpwl[phase], pre_factor = materials->pref_pwl[phase];
    double Ealin = materials->Qlin[phase], Valin = materials->Vlin[phase], nlin  = materials->nlin[phase], mlin  = materials->mlin[phase], rlin  = materials->rlin[phase], Alin  = materials->Alin[phase], flin = materials->flin[phase], alin = materials->alin[phase], Flin  = materials->Flin[phase];
    double Eagbs = materials->Qgbs[phase], Vagbs = materials->Vgbs[phase], ngbs  = materials->ngbs[phase], mgbs  = materials->mgbs[phase], rgbs  = materials->rlin[phase], Agbs  = materials->Agbs[phase], fgbs = materials->fgbs[phase], agbs = materials->agbs[phase], Fgbs  = materials->Fgbs[phase];
    double A1pwl=0.0, A1lin=0.0, A1exp=0.0, A1gbs=0.0, B_pwl=0.0,  B_lin=0.0, B_exp=0.0, B_gbs=0.0;
    double Eaexp = materials->Qexp[phase], Sexp  = materials->Sexp[phase], Eexp  = materials->Eexp[phase] , texp  = materials->texp[phase], Fexp=0.0;
    double gamma=materials->Gexp[phase], ST, q=materials->qexp[phase],  nexp=materials->nexp[phase];
    double Exx_lin=0.0, Exz_lin=0.0, Exx_exp=0.0, Exz_exp=0.0, Exx_gbs=0.0, Exz_gbs=0.0, Exx_cst=0.0, Exz_cst=0.0;
    int gs = materials->gs[phase];
    double pg = materials->ppzm[phase], Kg = materials->Kpzm[phase], Qg = materials->Qpzm[phase], gam = materials->Gpzm[phase], cg = materials->cpzm[phase], lambda = materials->Lpzm[phase];
    
    double etaVP = materials->eta_VP;
    
    double accu = model->accu;
    
    //------------------------------------------------------------------------//
    
    // Initialise strain rate invariants to 0
    *Eii_exp = 0.0; *Eii_lin = 0.0; *Eii_pl = 0.0; *Eii_pwl = 0.0; *Eii_el = 0.0, *Eii_gbs=0, *Eii_cst=0.0;
    *txxn=0.0; *txzn=0.0; *etaVE=0.0; *VEcoeff=0.0, *A2_pwl=0.0; *Exx_pwl=0.0; *Exz_pwl=0.0, *Exx_el=0.0, *Exz_el=0.0, *Exx_diss=0.0, *Exz_diss=0.0, *Exx_pl=0.0, *Exz_pl=0.0, *d1=0.0;
    
    // Activate deformation mechanisms
    if ( materials->cstv[phase] !=0                  ) constant    = 1;
    if ( materials->pwlv[phase] !=0                  ) dislocation = 1;
    if ( materials->expv[phase] !=0 && T<TmaxPeierls ) peierls     = 1;
    if ( materials->linv[phase] !=0                  ) diffusion   = 1;
    if ( materials->gbsv[phase] !=0                  ) gbs         = 1;
    if ( materials->gs[phase]   !=0                  ) gs          = 1;
    
    // Turn of elasticity for the initialisation step  (viscous flow stress)
    if ( model->step    == 0                         ) elastic     = 0;
    
    // Constant grain size
    if ( gs == 1 ) *d1 = materials->gs_ref[phase];
    else           *d1 = materials->gs_ref[phase];
    
    // Tensional cut-off
    if ( model->gz>0.0 && P<0.0     ) P = 0.0;
    
    // Visco-plastic limit
    if ( elastic==0                 ) G = 10.0;
    
    // Zero C limit
    if ( T< zeroC/scaling->T        ) T = zeroC/scaling->T;
    
    //------------------------------------------------------------------------//
    
    // Precomputations
    if ( dislocation == 1 ) {
        A1pwl = pre_factor * Fpwl * pow(Apwl,-1.0/npwl) * exp( (Eapwl + P*Vapwl)/R/npwl/T ) * pow(d, mpwl/npwl) * pow(fpwl, -rpwl/npwl) * exp(-apwl*phi/npwl);
        B_pwl = pow(2.0*A1pwl, -npwl);
        
    }
    if ( diffusion == 1 ) {
        if (mlin>0.0 && d<1e-13/scaling->L){
            printf("Cannot run with grain size dependent viscosity if grain size is set to 0 --> d = %2.2e!!!\n", d*scaling->L);
            exit(1);
        };
        A1lin = Flin * pow(Alin,-1.0/nlin) * exp( (Ealin + P*Valin)/R/nlin/T ) * pow(flin, -rlin/nlin) * exp(-alin*phi/nlin); // * pow(d, mlin/nlin) !!!!!!!!!!!!!!!!!!!!!!!!
        B_lin = pow(2.0*A1lin, -nlin);
    }
    if ( gbs == 1 ) {
        A1gbs = Fgbs * pow(Agbs,-1.0/ngbs) * exp( (Eagbs + P*Vagbs)/R/ngbs/T ) * pow(d, mgbs/ngbs) * pow(fgbs, -rgbs/ngbs) * exp(-agbs*phi/ngbs);
        B_gbs = pow(2.0*A1gbs, -ngbs);
    }
    if ( peierls   == 1 ) {
        ST                    = Eaexp/R/T * pow((1.0-gamma),(q-1.0)) * q*gamma;
        if ( texp == 0) Fexp  = 1.0;
        if ( texp == 1) Fexp  = 1.0/6.0*pow(2.0,1.0/(ST+nexp)) * pow(3.0,(ST+nexp-1.0)/2.0/(ST+nexp));
        if ( texp == 2) Fexp  = 1.0/4.0*pow(2,1.0/(ST+nexp));
        A1exp                   = Fexp * pow(Eexp*exp(-Eaexp/R/T*pow(1.0-gamma,2.0)), -1.0/(ST+nexp)) * pow(gamma*Sexp, ST/(ST+nexp));
        B_exp                   = pow(2.0*A1exp, -(ST+nexp));
    }
    
    MC_yield                    = C*cos(Phi) +  (P+model->PrBG)*sin(Phi);
    yield                       = MC_yield;
    yield                       = MINV(MC_yield,materials->Slim[phase]);
    
    //------------------------------------------------------------------------//
    
    // Isolated viscosities
    Exx  = exx;
    Exz  = exz;
    if ( elastic == 0 || local_iter == 0 || local_iter == 2 ) {
        Exx  = exx;
        Exz  = exz;
    }
    else {
        Exx  = exx + txx0/(2.0*G*dt);
        Exz  = exz + txz0/(2.0*G*dt);
    }
    Eii             = sqrt(Exx*Exx + Exz*Exz);
    double eII      = sqrt(exx*exx + exz*exz);
    if (Eii*scaling->E<1e-30) Eii=1e-30/scaling->E;
    
    //------------------------------------------------------------------------//
    
    // Isolated viscosities
    eta_pl                           = yield / (2.0*Eii);
    eta_el                           = G*dt;
    
    if ( constant    == 1 ) eta_cst  = materials->eta0[phase];
    if ( dislocation == 1 ) eta_pwl  = A1pwl * pow( Eii, 1.0/npwl - 1.0 );
    if ( diffusion   == 1 ) eta_lin  = A1lin * pow( Eii, 1.0/nlin - 1.0 ) * pow(d, mlin/nlin); // !!! gs - dependence !!!
    if ( gbs         == 1 ) eta_gbs  = A1gbs * pow( Eii, 1.0/ngbs - 1.0 );
    if ( peierls     == 1 ) eta_exp  = A1exp * pow( Eii, 1.0/(ST+nexp) - 1.0 );
    
    // Effective viscosity
    eta                          = 0.0;
    if ( constant    == 1 ) eta += 1.0/eta_cst;
    if ( dislocation == 1 ) eta += 1.0/eta_pwl;
    if ( diffusion   == 1 ) eta += 1.0/eta_lin;
    if ( gbs         == 1 ) eta += 1.0/eta_gbs;
    if ( peierls     == 1 ) eta += 1.0/eta_exp;
    eta                      = 1.0/eta;
    *etaVE                   = 1.0/(1.0/eta_cst + 1.0/eta_el);
    
    //------------------------------------------------------------------------//
    
    if ( local_iter == 2 ) {
        
        // Effective viscosity
        eta                          = 0.0;
        if ( constant    == 1 ) eta += 1.0/eta_cst;
        if ( dislocation == 1 ) eta += 1.0/eta_pwl;
        if ( diffusion   == 1 ) eta += 1.0/eta_lin;
        if ( gbs         == 1 ) eta += 1.0/eta_gbs;
        if ( peierls     == 1 ) eta += 1.0/eta_exp;
        eta                      = 1.0/eta;
        *etaVE                   = 1.0/(1.0/eta + 1.0/eta_el);
        *VEcoeff                 = 1.0 / (1.0 + G*model->dt/eta);
        if ( elastic== 0 ) {
            *etaVE               = eta;
            *VEcoeff             = 0.0;
        }
        
        // Stress
        txx1 = 2.0 * (*etaVE) * Exx + (*VEcoeff)*txx0;
        txz1 = 2.0 * (*etaVE) * Exz + (*VEcoeff)*txz0;
        tII1 = sqrt(txx1*txx1 + txz1*txz1);
        
        // Strain rates
        *Eii_pl                          = 0.0;
        if ( dislocation == 1 ) *Eii_pwl = tII1/2.0/eta_pwl;
        if ( diffusion   == 1 ) *Eii_lin = tII1/2.0/eta_lin * pow(d,-mlin); // !!! gs - dependence !!!
        if ( gbs         == 1 ) *Eii_gbs = tII1/2.0/eta_gbs;
        if ( peierls     == 1 ) *Eii_exp = tII1/2.0/eta_exp;
        
        if (tII1>yield) {
            *Eii_pl = ( 1.0/(1.0 + etaVP/(*etaVE)) ) * ( eII - (yield - tII0)/2.0/eta_el - yield/2.0/eta );
            *Exx_pl = *Eii_pl  * (exx/eII);
            *Exz_pl = *Eii_pl  * (exz/eII);
            txx1 = 2.0 * (*etaVE) * (Exx-(*Exx_pl)) + (*VEcoeff)*txx0;
            txz1 = 2.0 * (*etaVE) * (Exz-(*Exz_pl)) + (*VEcoeff)*txz0;
            //printf("tII1 = %2.4e Ty = %2.4e F = %2.4e Eii = %2.4e eii = %2.4e eii_pl = %2.4e  yield = %2.4e num2 = %2.4e num2 = %2.4e\n", tII1*scaling->S, yield*scaling->S, (tII1 - yield)*scaling->S, Eii*scaling->E, eII*scaling->E, (*Eii_pl)*scaling->E, yield*scaling->S, (yield - tII0)/2.0/eta_el, yield/2.0/eta);
        }
        
        *Eii_el  = eII - (*Eii_pl) - (*Eii_pwl) - (*Eii_lin) - (*Eii_gbs)  - (*Eii_exp);
        
        // Grain size update rule (Paleowattmeter)
        if ( gs          == 1 ) *d1      = exp(log( Kg*exp(-Qg/R/T) *gam/(lambda*(1.0/cg)*tII1*(*Eii_pwl + *Eii_exp + *Eii_gbs + *Eii_pl)*pg))/(1.0+pg));
        
        // Viscosity limiter
        if( eta > maxEta ) {
            eta = maxEta;
        }
        
        if( eta < minEta ) {
            eta = minEta;
        }
    }
    
    //------------------------------------------------------------------------//
    
    if ( local_iter == 0 ) {
        
        // Effective viscosity
        eta                          = 0.0;
        if ( constant    == 1 ) eta += 1.0/eta_cst;
        if ( dislocation == 1 ) eta += 1.0/eta_pwl;
        if ( diffusion   == 1 ) eta += 1.0/eta_lin;
        if ( gbs         == 1 ) eta += 1.0/eta_gbs;
        if ( peierls     == 1 ) eta += 1.0/eta_exp;
        eta                      = 1.0/eta;
        *etaVE                   = 1.0/(1.0/eta + 1.0/eta_el);
        *VEcoeff                 = 1.0 / (1.0 + G*model->dt/eta);
        if ( elastic== 0 ) {
            *etaVE               = eta;
            *VEcoeff             = 0.0;
        }
        
        // Stress
        txx1 = 2.0 * (*etaVE) * Exx + (*VEcoeff)*txx0;
        txz1 = 2.0 * (*etaVE) * Exz + (*VEcoeff)*txz0;
        tII1 = sqrt(txx1*txx1 + txz1*txz1);
        
        // Strain rates
        *Eii_pl                          = 0.0;
        if ( elastic     == 1 ) *Eii_el  = (tII1-tII0)/2.0/G/dt;
        if ( dislocation == 1 ) *Eii_pwl = tII1/2.0/eta_pwl;
        if ( diffusion   == 1 ) *Eii_lin = tII1/2.0/eta_lin * pow(d,-mlin); // !!! gs - dependence !!!
        if ( gbs         == 1 ) *Eii_gbs = tII1/2.0/eta_gbs;
        if ( peierls     == 1 ) *Eii_exp = tII1/2.0/eta_exp;
        
        //  Drucker-Prager plasticity
        if ( tII1 > yield ) {
            eta      = fabs(yield) / ( 2.0 * Eii );
            *etaVE   = eta;//1.0 / (1.0/eta + 1.0/G*model->dt);
            *VEcoeff = 0.0;//1.0 / (1.0 + G*model->dt/eta);
            txx1     = 2.0 * (*etaVE) * Exx + (*VEcoeff)*txx0;
            txz1     = 2.0 * (*etaVE) * Exz + (*VEcoeff)*txz0;
            *Eii_el  = 0.0;
            *Eii_pl  = Eii;
            *Eii_pwl = 0.0;
            *Eii_exp = 0.0;
            *Eii_lin = 0.0;
        }
        
        // Grain size update rule (Paleowattmeter)
        if ( gs          == 1 ) *d1      = exp(log( Kg*exp(-Qg/R/T) *gam/(lambda*(1.0/cg)*tII1*(*Eii_pwl + *Eii_exp + *Eii_gbs + *Eii_pl)*pg))/(1.0+pg));
        //        if ( gs          == 1 ) *d1      = pzPreComp * pow(tII1*(*Eii_pwl + *Eii_exp + *Eii_gbs), pzPow);
        //        if ( gs          == 1 ) *d1      = pow(  (Kg*exp(-Qg/RT)*cg*gam) / (pg*lambda*tII1*(*Eii_pwl + *Eii_exp + *Eii_gbs))  , 1.0/(1.0+pg));
        
        // Viscosity limiter
        if( eta > maxEta ) {
            eta = maxEta;
        }
        
        if( eta < minEta ) {
            eta = minEta;
        }
    }
    
    //------------------------------------------------------------------------//
    
    if ( local_iter == 1 ) {
        
        // Function at plastic viscosity
        tII1 = 2.0 * eta_pl * Eii;
        if ( constant    == 1 ) *Eii_cst = tII1/2.0/eta_cst;
        if ( dislocation == 1 ) *Eii_pwl = B_pwl * pow(tII1, npwl    );
        if ( gbs         == 1 ) *Eii_gbs = B_gbs * pow(tII1, ngbs    );
        if ( peierls     == 1 ) *Eii_exp = B_exp * pow(tII1, ST+nexp ); // Peierls - power law
        if ( gs          == 1 ) *d1      = exp(log( Kg*exp(-Qg/R/T) *gam/(lambda*(1.0/cg)*tII1*(*Eii_pwl + *Eii_exp + *Eii_gbs + *Eii_pl)*pg))/(1.0+pg));
        if ( diffusion   == 1 ) *Eii_lin = B_lin * pow(tII1, nlin) * pow(*d1,-mlin); // !!!!!!!!!!!!!!!
        Eii_vis                          = *Eii_pwl + *Eii_exp + *Eii_lin + *Eii_gbs + *Eii_cst;
        r_eta_pl                         = (Eii - elastic*tII1/(2.0*G*dt) - Eii_vis);
        //        r_eta_pl                         = ( 1.0/(1.0 + etaVP/(*etaVE)) ) * ( eII - (yield - tII0)/2.0/eta_el - yield/2.0/eta );
        
        if (r_eta_pl >= 0.0) {
            
            *Eii_pl    = r_eta_pl;                              // plastic strain rate
            eta_pl     = (yield + etaVP*r_eta_pl) / ( 2.0 * Eii );
            *etaVE     = eta_pl;                                // elasto-plastic viscosity
            eta_0      = eta_pl;                                // elasto-plastic viscosity
            eta        = yield / (2.0*sqrt(exx*exx + exz*exz)); //        plastic viscosity
            
        }
        
        else {
            
            // Viscoelasticity
            *Eii_pl = 0.0;
            
            // Define viscosity bounds
            eta_up   = 1.0e100/scaling->S;
            ieta_sum = 0.0;
            
            if ( constant == 1 ) {
                eta_up   = MINV(eta_up, eta_cst);
                ieta_sum += 1.0/eta_cst;
            }
            
            if ( dislocation == 1 ) {
                eta_up   = MINV(eta_up, eta_pwl);
                ieta_sum += 1.0/eta_pwl;
            }
            
            if ( elastic == 1 ) {
                eta_up   = MINV(eta_up, eta_el);
                ieta_sum += 1.0/eta_el;
            }
            
            if ( peierls == 1 ) {
                eta_up = MINV(eta_up, eta_exp);
                ieta_sum += 1.0/eta_exp;
            }
            if ( diffusion == 1 ) {
                eta_up = MINV(eta_up, eta_lin);
                ieta_sum += 1.0/eta_lin;
            }
            if ( gbs       == 1 ) {
                eta_up = MINV(eta_up, eta_gbs);
                ieta_sum += 1.0/eta_gbs;
            }
            eta_lo = 1.0/(ieta_sum);
            
            
            if ( isnan(eta_lo)==1  || isnan(eta_up)==1) printf("%d %d %d %d %d %d\n", constant, dislocation, diffusion, gbs, peierls, elastic);
            if ( isnan(eta_lo)==1  || isnan(eta_up)==1) printf("d = %2.2e eII = %2.2e  T = %lf\n", d*scaling->L, Eii*scaling->E, T*scaling->T);
            if ( isnan(eta_lo)==1  || isnan(eta_up)==1) printf("eta_el=%2.2e eta_pl=%2.2e eta_pwl=%2.2e eta_exp=%2.2e eta_lin=%2.2e eta_gbs=%2.2e eta_cst=%2.2e\n", eta_el*scaling->eta,eta_pl*scaling->eta,eta_pwl*scaling->eta,eta_exp*scaling->eta,eta_lin*scaling->eta, eta_gbs*scaling->eta, eta_cst*scaling->eta);
            if ( isnan(eta_lo)==1  || isnan(eta_up)==1) printf("eta_lo=%2.2e eta_up=%2.2e\n",eta_lo*scaling->eta, eta_up*scaling->eta);
            if ( isnan(eta_lo)==1  || isnan(eta_up)==1) printf("G=%2.2e dt=%2.2e exx=%2.2e exz=%2.2e Exx=%2.2e Exz=%2.2e txx0=%2.2e txz0=%2.2e\n", G*scaling->S, dt*scaling->t, exx*scaling->E, exz*scaling->E, Exx*scaling->E, Exz*scaling->E, txx0*scaling->S, txz0*scaling->S);
            if ( isnan(eta_lo)==1  || isnan(eta_up)==1) exit(123);
            
            // Local iterations
            for (it=0; it<nitmax; it++) {
                
                // Secant method
                if (it == 0) {
                    // Start at midpoint (initial guess)
                    eta_0     = 0.5*(eta_up+eta_lo);
                }
                
                // Function evaluation at current effective viscosity
                tII1 = 2.0 * eta_0 * Eii;
                if ( constant    == 1 ) *Eii_cst = tII1/2.0/eta_cst;
                if ( dislocation == 1 ) *Eii_pwl = B_pwl * pow(tII1, npwl    );
                if ( gbs         == 1 ) *Eii_gbs = B_gbs * pow(tII1, ngbs    );
                if ( peierls     == 1 ) *Eii_exp = B_exp * pow(tII1, ST+nexp ); // Peierls - power law
                if ( gs          == 1 ) *d1      = exp(log( Kg*exp(-Qg/R/T) *gam/(lambda*(1.0/cg)*tII1*(*Eii_pwl + *Eii_exp + *Eii_gbs + *Eii_pl)*pg))/(1.0+pg));
                if ( diffusion   == 1 ) *Eii_lin = B_lin * pow(tII1, nlin) * pow(*d1,-mlin); // !!! gs - dependence !!!
                Eii_vis                          = *Eii_pwl + *Eii_exp + *Eii_lin + *Eii_gbs + *Eii_cst;
                r_eta_0                          = Eii - elastic*tII1/(2.0*eta_el) - Eii_vis;
                
                // Residual check
                res = fabs(r_eta_0/Eii);
                if (noisy==1) printf("It. %02d, r = %2.2e\n", it, res);
                if (res < tol) break;
                
                // Perturbation
                eta_step  = tol*eta_0;
                eta_p     = eta_0 + eta_step;
                
                // Function evaluation at pertubated effective viscosity
                tII1 = 2.0 * eta_p * Eii;
                if ( constant    == 1 ) *Eii_cst = tII1/2.0/eta_cst;
                if ( dislocation == 1 ) *Eii_pwl = B_pwl * pow(tII1, npwl    );
                if ( gbs         == 1 ) *Eii_gbs = B_gbs * pow(tII1, ngbs    );
                if ( peierls     == 1 ) *Eii_exp = B_exp * pow(tII1, ST+nexp ); // Peierls - power law
                if ( gs          == 1 ) *d1      = exp(log( Kg*exp(-Qg/R/T) *gam/(lambda*(1.0/cg)*tII1*(*Eii_pwl + *Eii_exp + *Eii_gbs + *Eii_pl)*pg))/(1.0+pg));
                if ( diffusion   == 1 ) *Eii_lin = B_lin * pow(tII1, nlin) * pow(*d1,-mlin); // !!! gs - dependence !!!
                Eii_vis                          = *Eii_pwl + *Eii_exp + *Eii_lin + *Eii_gbs + *Eii_cst;
                r_eta_p                          = Eii - elastic*tII1/(2.0*eta_el) - Eii_vis;
                
                // Finite difference derivative of function
                dfdeta = (r_eta_p - r_eta_0) / (eta_p - eta_0);
                
                // Update viscosity
                eta_0  = eta_0 - r_eta_0 / dfdeta;
            }
            
            // Check for success
            if (res>tol) {
                printf("Local iterations broke down -> res = %2.2e, eII = %2.2e\n", r_eta_0*scaling->E, Eii*scaling->E );
                printf("%d eII = %2.1e  T = %lf d= %2.2e eta_pwl = %2.4e eta_lin = %2.4e eta_gbs = %2.4e eta_exp = %2.4e eta_pl = %2.4e eta_cst = %2.4e\n", phase, Eii*scaling->E, (T)*scaling->T-zeroC,  d*scaling->L, eta_pwl*scaling->eta, eta_lin*scaling->eta, eta_gbs*scaling->eta, eta_exp*scaling->eta, eta_pl*scaling->eta, eta_cst*scaling->eta);
                printf("e_pwl = %2.4e e_lin = %2.4e e_gbs = %2.4e e_exp = %2.4e e_pl = %2.4e e_cst = %2.4e sum=%2.2e\n", *Eii_pwl*scaling->E, *Eii_lin*scaling->E, *Eii_gbs*scaling->E, *Eii_exp*scaling->E, *Eii_pl*scaling->E, *Eii_cst*scaling->E,(*Eii_pwl+*Eii_lin+*Eii_gbs+*Eii_exp+*Eii_pl+*Eii_cst-eII)*scaling->E);
                exit(12);
            }
            
            // pre-calculation for analytical Jacobian
            *A2_pwl   =  A1pwl * pow( (*Eii_pwl)*(*Eii_pwl), (0.5*(1.0/npwl-1.0) - 1.0 ) )*(1.0/npwl-1.0);
            
            // Viscosity for dissipative processes (no elasticity)
            eta        = tII1/2.0/Eii_vis;
            //            if (phase==1) printf("TII=%2.2e eta_0=%2.2e  eta=%2.2e\n", tII1*scaling->S, eta_0*scaling->eta, eta*scaling->eta);
            if (isnan(eta)) {
                printf("eta_el=%2.2e eta_pl=%2.2e eta_pwl=%2.2e eta_exp=%2.2e eta_lin=%2.2e eta_gbs=%2.2e eta_cst=%2.2e\n", eta_el*scaling->eta,eta_pl*scaling->eta,eta_pwl*scaling->eta,eta_exp*scaling->eta,eta_lin*scaling->eta, eta_gbs*scaling->eta, eta_cst*scaling->eta);
                printf("eii_el=%2.2e eii_pl=%2.2e eii_pwl=%2.2e eii_exp=%2.2e eii_lin=%2.2e eii_gbs=%2.2e eii_cst=%2.2e\n", *Eii_el*scaling->E,*Eii_pl*scaling->E,*Eii_pwl*scaling->E,*Eii_exp*scaling->E,*Eii_lin*scaling->E, *Eii_gbs*scaling->E, *Eii_cst*scaling->E);
                printf("d = %2.2e eII = %2.2e  T = %lf\n", d*scaling->L, Eii*scaling->E, T*scaling->T);
                printf("eta_el=%2.2e eta_pl=%2.2e eta_pwl=%2.2e eta_exp=%2.2e eta_lin=%2.2e \n", eta_el*scaling->eta,eta_pl*scaling->eta,eta_pwl*scaling->eta,eta_exp*scaling->eta,eta_lin*scaling->eta);
                printf("G=%2.2e dt=%2.2e exx=%2.2e exz=%2.2e Exx=%2.2e Exz=%2.2e txx0=%2.2e txz0=%2.2e\n", G*scaling->S, dt*scaling->t, exx*scaling->E, exz*scaling->E, Exx*scaling->E, Exz*scaling->E, txx0*scaling->S, txz0*scaling->S);
                printf("Arrhenius lin. = %2.2e Arrhenius exp. = %2.2e\n", exp( (Ealin + P*Valin)/R/nlin/T ), pow(Eexp*exp(-Eaexp/R/T*pow(1.0-gamma,2.0)), -1.0/(ST+nexp)) );
                printf("Arrhenius lin. = %2.2e Arrhenius exp. = %2.2e\n", exp( (Ealin)/R/nlin/T ), pow(Eexp*exp(-Eaexp/R/T*pow(1.0-gamma,2.0)), -1.0/(ST+nexp)) );
                printf("\n");
                
            }
            
            if (isinf(eta)) {
                printf("eta_el=%2.2e eta_pl=%2.2e eta_pwl=%2.2e eta_exp=%2.2e eta_lin=%2.2e eta_gbs=%2.2e eta_cst=%2.2e\n", eta_el*scaling->eta,eta_pl*scaling->eta,eta_pwl*scaling->eta,eta_exp*scaling->eta,eta_lin*scaling->eta, eta_gbs*scaling->eta, eta_cst*scaling->eta);
                
                printf("eii_el=%2.2e eii_pl=%2.2e eii_pwl=%2.2e eii_exp=%2.2e eii_lin=%2.2e eii_gbs=%2.2e eii_cst=%2.2e\n", *Eii_el*scaling->E,*Eii_pl*scaling->E,*Eii_pwl*scaling->E,*Eii_exp*scaling->E,*Eii_lin*scaling->E, *Eii_gbs*scaling->E, *Eii_cst*scaling->E);
                printf("d = %2.2e eII = %2.2e  T = %lf\n", d*scaling->L, Eii*scaling->E, T*scaling->T);
                printf("eta_el=%2.2e eta_pl=%2.2e eta_pwl=%2.2e eta_exp=%2.2e eta_lin=%2.2e \n", eta_el*scaling->eta,eta_pl*scaling->eta,eta_pwl*scaling->eta,eta_exp*scaling->eta,eta_lin*scaling->eta);
                printf("G=%2.2e dt=%2.2e exx=%2.2e exz=%2.2e Exx=%2.2e Exz=%2.2e txx0=%2.2e txz0=%2.2e\n", G*scaling->S, dt*scaling->t, exx*scaling->E, exz*scaling->E, Exx*scaling->E, Exz*scaling->E, txx0*scaling->S, txz0*scaling->S);
                printf("Arrhenius lin. = %2.2e Arrhenius exp. = %2.2e\n", exp( (Ealin + P*Valin)/R/nlin/T ), pow(Eexp*exp(-Eaexp/R/T*pow(1.0-gamma,2.0)), -1.0/(ST+nexp)) );
                printf("Arrhenius lin. = %2.2e Arrhenius exp. = %2.2e\n", exp( (Ealin)/R/nlin/T ), pow(Eexp*exp(-Eaexp/R/T*pow(1.0-gamma,2.0)), -1.0/(ST+nexp)) );
                printf("\n");
                
            }
            
        }
        (*etaVE)   = eta_0;
        *VEcoeff = 1.0 / (1.0 + G*model->dt/eta);
        if (elastic==0) *VEcoeff = 0.0;
        
        // Recalculate stress components
        txx1                             = 2.0*eta_0*Exx;
        txz1                             = 2.0*eta_0*Exz;
        tII1                             = 2.0*eta_0*Eii;
        //          if (r_eta_pl >= 0.0 && model->step>0) {
        //              printf("F = %2.2e Pa tII = %2.2e, Sy = %2.2e\n", tII1 - (yield), tII1, yield);
        //              printf("F = %2.2e Pa tII = %2.2e, Sy = %2.2e\n", tII1 - (yield + etaVP*r_eta_pl), tII1, yield);
        //              exit(1);
        //          }
        
        // Post-process elastic effective strain rate
        *Exx_el = elastic*(txx1 - txx0) / (2.0*eta_el);
        *Exz_el = elastic*(txz1 - txz0) / (2.0*eta_el);
        *Eii_el = sqrt( (*Exx_el)*(*Exx_el) + (*Exz_el)*(*Exz_el) );
        *Exx_pl = *Eii_pl * (txx1/tII1);  *Exz_pl = *Eii_pl * (txz1/tII1);
        
        // Recompute strain rates and check if sum of components is equal to total strain rate components
        Exx_lin = *Eii_lin * (txx1/tII1);  Exz_lin = *Eii_lin * (txz1/tII1);
        *Exx_pwl = *Eii_pwl * (txx1/tII1); *Exz_pwl = *Eii_pwl * (txz1/tII1);
        Exx_exp = *Eii_exp * (txx1/tII1);  Exz_exp = *Eii_exp * (txz1/tII1);
        Exx_gbs = *Eii_gbs * (txx1/tII1);  Exz_gbs = *Eii_gbs * (txz1/tII1);
        Exx_cst = *Eii_cst * (txx1/tII1);  Exz_cst = *Eii_cst * (txz1/tII1);
        
        //        double chkxx, chkxz;
        //        chkxx = exx - (*Exx_el + Exx_lin + *Exx_pwl + Exx_exp + Exx_gbs + Exx_cst);
        //        chkxz = exz - (*Exz_el + Exz_lin + *Exz_pwl + Exz_exp + Exz_gbs + Exz_cst);
        //
        //        if ( model->step    > 0 ) {
        //            chkxx = exx - (*Exx_el + Exx_lin + *Exx_pwl + Exx_exp + Exx_gbs + Exx_cst + *Exx_pl);
        //            chkxz = exz - (*Exz_el + Exz_lin + *Exz_pwl + Exz_exp + Exz_gbs + Exz_cst + *Exz_pl);
        //        }
        
    }
    
    /*----------------------------------------------------*/
    /*----------------------------------------------------*/
    /*----------------------------------------------------*/
    
    // Compute dissipative strain rate components
    *Exx_diss = *Exx_pl + Exx_lin + *Exx_pwl + Exx_exp + Exx_gbs + Exx_cst;
    *Exz_diss = *Exz_pl + Exz_lin + *Exz_pwl + Exz_exp + Exz_gbs + Exz_cst;
    
    // Override viscosty at step 0 (100% visco-plastic)
    if ( model->step == 0 ) *etaVE = eta;
    *txxn = txx1;
    *txzn = txz1;
    
    // Viscosity limiter
    if( eta > maxEta ) {
        eta = maxEta;
    }
    
    if( eta < minEta ) {
        eta = minEta;
    }
    
    // Viscosity limiter
    if( *etaVE > maxEta ) {
        *etaVE = maxEta;
    }
    
    if( *etaVE < minEta ) {
        *etaVE = minEta;
    }
    
    if ( model->cut_noise==1 ){
        eta      = round(eta*accu)/accu;
        *etaVE   = round(*etaVE*accu)/accu;
        *txxn    = round(*txxn*accu)  /accu;
        *txzn    = round(*txzn*accu)/accu;
        *VEcoeff = round(*VEcoeff*accu)/accu;
        *d1      = round(*d1*accu)/accu;
    }
    
    return eta;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
void NonNewtonianViscosityCells( grid* mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling, int flag ) {
    
    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1, k1, cond;
    double eta, txx1, txz1, Pn, Tn, etaVE, VEcoeff=0.0, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1;
    double exx_pwl, exz_pwl, A2_pwl, exx_el, exz_el, exx_diss, exz_diss, exx_pl, exz_pl;
    int average = model->eta_avg;
    double accu = model->accu;
    
    printf("In NonNewtonianViscosityGrid\n");
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    // Evaluate cell center viscosities
#pragma omp parallel for shared( mesh  ) private( cond, k, l, k1, p, eta, c1, c0, Pn, Tn, txx1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1, exx_pwl, exz_pwl, A2_pwl, exx_el, exz_el, exx_diss, exz_diss, exx_pl, exz_pl ) firstprivate( materials, scaling, flag, average, model, Ncx, Ncz )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {
        
        //    for ( l=0; l<Ncz; l++ ) {
        //        for ( k=0; k<Ncx; k++ ) {
        
        k      = mesh->kp[k1];
        l      = mesh->lp[k1];
        c0 = k  + l*(Ncx);
        
        mesh->eta_n[c0]   = 0.0;
        mesh->eta_phys_n[c0] = 0.0;
        mesh->VE_n[c0]       = 0.0;
        mesh->sxxd[c0]       = 0.0;
        mesh->eII_el[c0]     = 0.0;
        mesh->eII_pl[c0]     = 0.0;
        mesh->eII_pwl[c0]    = 0.0;
        mesh->eII_exp[c0]    = 0.0;
        mesh->eII_lin[c0]    = 0.0;
        mesh->eII_gbs[c0]    = 0.0;
        mesh->eII_cst[c0]    = 0.0;
        mesh->d[c0]          = 0.0;
        mesh->exx_pwl_n[c0]  = 0.0;
        mesh->exz_pwl_n[c0]  = 0.0;
        mesh->A2_pwl_n[c0]   = 0.0;
        mesh->exx_el[c0]     = 0.0;
        mesh->exx_diss[c0]   = 0.0;
        mesh->exx_pl[c0]     = 0.0;
        
        mesh->sxz_n[c0]      = 0.0;
        mesh->exz_n_el[c0]   = 0.0;
        mesh->exz_n_pl[c0]   = 0.0;
        mesh->exz_n_diss[c0] = 0.0;
        
        // Loop on grid nodes
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
            
            // If background pressure is removed in the RHS
            //                    Pn = (mesh->p_in[c0]+mesh->p_lith[c0]);
            Pn = mesh->p_in[c0];
            Tn = mesh->T[c0];
            
            // Loop on phases
            for ( p=0; p<model->Nb_phases; p++) {
                
                cond =  fabs(mesh->phase_perc_n[p][c0])>1.0e-13;
                
                if ( cond == 1 ) {
                    eta =  Viscosity( p, mesh->mu_n[c0], Tn, Pn, mesh->d0[c0], mesh->phi[c0],  mesh->eii_n[c0],  mesh->tii0_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0], materials    , model, scaling, flag, &txx1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &exz_el, &exx_diss, &exz_diss, &exx_pl, &exz_pl, &d1, &A2_pwl, mesh->strain_n[c0], mesh->phi_n[c0], mesh->C_n[c0]);
                }
                
                // ARITHMETIC AVERAGE
                if (average == 0) {
                    if ( cond == 1 ) mesh->eta_n[c0]   += mesh->phase_perc_n[p][c0] * etaVE;
                    if ( cond == 1 ) mesh->eta_phys_n[c0] += mesh->phase_perc_n[p][c0] * eta;
                }
                if ( cond == 1 ) mesh->VE_n[c0]       += mesh->phase_perc_n[p][c0] * VEcoeff;
                if (average == 0 || average == 2) {
                    if ( cond == 1 ) mesh->sxxd[c0]       += mesh->phase_perc_n[p][c0] * txx1;
                    if ( cond == 1 ) mesh->sxz_n[c0]      += mesh->phase_perc_n[p][c0] * txz1;
                }
                if ( cond == 1 ) mesh->eII_el[c0]     += mesh->phase_perc_n[p][c0] * eII_el;
                if ( cond == 1 ) mesh->eII_pl[c0]     += mesh->phase_perc_n[p][c0] * eII_pl;
                if ( cond == 1 ) mesh->eII_pwl[c0]    += mesh->phase_perc_n[p][c0] * eII_pwl;
                if ( cond == 1 ) mesh->eII_exp[c0]    += mesh->phase_perc_n[p][c0] * eII_exp;
                if ( cond == 1 ) mesh->eII_lin[c0]    += mesh->phase_perc_n[p][c0] * eII_lin;
                if ( cond == 1 ) mesh->eII_gbs[c0]    += mesh->phase_perc_n[p][c0] * eII_gbs;
                if ( cond == 1 ) mesh->eII_cst[c0]    += mesh->phase_perc_n[p][c0] * eII_cst;
                if ( cond == 1 ) mesh->d[c0]          += mesh->phase_perc_n[p][c0] * 1.0/d1;
                
                if ( cond == 1 ) mesh->exx_pwl_n[c0]  += mesh->phase_perc_n[p][c0] * exx_pwl;
                if ( cond == 1 ) mesh->exz_pwl_n[c0]  += mesh->phase_perc_n[p][c0] * exz_pwl;
                if ( cond == 1 ) mesh->A2_pwl_n[c0]   += mesh->phase_perc_n[p][c0] * A2_pwl;
                if ( cond == 1 ) mesh->exx_el[c0]     += mesh->phase_perc_n[p][c0] * exx_el;
                if ( cond == 1 ) mesh->exx_diss[c0]   += mesh->phase_perc_n[p][c0] * exx_diss;
                if ( cond == 1 ) mesh->exx_pl[c0]     += mesh->phase_perc_n[p][c0] * exx_pl;
                
                if ( cond == 1 ) mesh->exz_n_pl[c0]     += mesh->phase_perc_n[p][c0] * exz_pl;
                if ( cond == 1 ) mesh->exz_n_el[c0]     += mesh->phase_perc_n[p][c0] * exz_el;
                if ( cond == 1 ) mesh->exz_n_diss[c0]   += mesh->phase_perc_n[p][c0] * exz_diss;
                //                    sum += mesh->phase_perc_n[p][c0];
                
                // HARMONIC AVERAGE
                if (average == 1) {
                    if ( cond == 1 ) mesh->sxxd[c0]       += mesh->phase_perc_n[p][c0] * 1.0/txx1;
                    if ( cond == 1 ) mesh->eta_n[c0]   += mesh->phase_perc_n[p][c0] * 1.0/etaVE;
                    if ( cond == 1 ) mesh->eta_phys_n[c0] += mesh->phase_perc_n[p][c0] * 1.0/eta;
                }
                
                // GEOMETRIC AVERAGE
                if (average == 2) {
                    if ( cond == 1 ) mesh->eta_n[c0]   += mesh->phase_perc_n[p][c0] * log(etaVE);
                    if ( cond == 1 ) mesh->eta_phys_n[c0] += mesh->phase_perc_n[p][c0] * log(eta);
                    
                    
                }
                
            }
            
            mesh->d[c0]          = 1.0/mesh->d[c0];
            
            // HARMONIC AVERAGE
            if (average == 1) {
                mesh->sxxd[c0]       = 1.0/mesh->sxxd[c0];
                mesh->eta_n[c0]   = 1.0/mesh->eta_n[c0];
                mesh->eta_phys_n[c0] = 1.0/mesh->eta_phys_n[c0];
                if (isinf (mesh->eta_phys_n[c0]) ) {
                    for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_n[p][c0]);
                    printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n", mesh->mu_n[c0], Tn, Pn, mesh->d0[c0], mesh->phi[c0],  mesh->eii_n[c0],  mesh->tii0_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0]);
                    printf("flag %d nb part cell = %d cell index = %d\n", mesh->BCp.type[c0],mesh->nb_part_cell[c0], c0);
                    printf("x=%2.2e z=%2.2e\n", mesh->xc_coord[k]*scaling->L/1000, mesh->zc_coord[l]*scaling->L/1000);
                    exit(1);
                }
            }
            
            // GEOMETRIC AVERAGE
            if (average == 2) {
                mesh->eta_n[c0]   = exp(mesh->eta_n[c0]);
                mesh->eta_phys_n[c0] = exp(mesh->eta_phys_n[c0]);
            }
            
            if ( model->cut_noise == 1 ) {
                mesh->eta_n[c0]   = round(mesh->eta_n[c0]*accu)/accu;
                mesh->eta_phys_n[c0] = round(mesh->eta_phys_n[c0]*accu)/accu;
                mesh->sxxd[c0]       = round(mesh->sxxd[c0]*accu)/accu;
                mesh->d[c0]          = round(mesh->d[c0]*accu)/accu;
                mesh->exz_n_pl[c0]   = round(mesh->exz_n_pl[c0]*accu)/accu;
                mesh->exz_n_el[c0]   = round(mesh->exz_n_el[c0]*accu)/accu;
                mesh->exz_n_diss[c0] = round(mesh->exz_n_diss[c0]*accu)/accu;
                
                mesh->eII_el[c0]    = round(mesh->eII_el[c0]*accu)/accu;
                mesh->eII_pl[c0]    = round(mesh->eII_pl[c0]*accu)/accu;
                mesh->eII_pwl[c0]   = round(mesh->eII_pwl[c0]*accu)/accu;
                mesh->eII_exp[c0]   = round(mesh->eII_exp[c0]*accu)/accu;
                mesh->eII_lin[c0]   = round(mesh->eII_lin[c0]*accu)/accu;
                mesh->eII_gbs[c0]   = round(mesh->eII_gbs[c0]*accu)/accu;
                mesh->eII_cst[c0]   = round(mesh->eII_cst[c0]*accu)/accu;
                
                mesh->exx_pwl_n[c0] = round(mesh->exx_pwl_n[c0]*accu)/accu;
                mesh->exz_pwl_n[c0] = round(mesh->exz_pwl_n[c0]*accu)/accu;
                mesh->A2_pwl_n[c0]  = round(mesh->A2_pwl_n[c0]*accu)/accu;
                mesh->exx_el[c0]    = round(mesh->exx_el[c0]*accu)/accu;
                mesh->exx_diss[c0]  = round(mesh->exx_diss[c0]*accu)/accu;
                mesh->exx_pl[c0]    = round(mesh->exx_pl[c0]*accu)/accu;
                
                mesh->VE_n[c0]  = round(mesh->VE_n[c0]*accu)/accu;
                mesh->sxz_n[c0]    = round(mesh->sxz_n[c0]*accu)/accu;
            }
        }
        
    }
    
    InterpCentroidsToVerticesDouble( mesh->sxz_n, mesh->sxz, mesh, model, scaling );
    InterpCentroidsToVerticesDouble( mesh->eta_n, mesh->eta_s, mesh, model, scaling );
    InterpCentroidsToVerticesDouble( mesh->eta_phys_n, mesh->eta_phys_s, mesh, model, scaling );
    InterpCentroidsToVerticesDouble( mesh->VE_n, mesh->VE_s, mesh, model, scaling );
    
    double eta_per, eta_VE_per, VE_per;
    for ( l=0; l<Nz; l++ ) {
        
        c0 = 0 + l*(Nx);
        c1 = (Nx-1) + l*(Nx);
        
        eta_per    = 0.5*( mesh->eta_phys_s[c0] + mesh->eta_phys_s[c1] );
        eta_VE_per = 0.5*( mesh->eta_s[c0] + mesh->eta_s[c1] );
        VE_per     = 0.5*( mesh->VE_s[c0] + mesh->VE_s[c1] );
        
        if ( mesh->BCg.type[c1] != 30 ) {
            
            mesh->eta_phys_s[c0] = eta_per;
            mesh->eta_phys_s[c1] = eta_per;
            mesh->eta_s[c0]   = eta_VE_per;
            mesh->eta_s[c1]   = eta_VE_per;
            mesh->VE_s[c0]       = VE_per;
            mesh->VE_s[c1]       = VE_per;
            
        }
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void NonNewtonianViscosityGrid( grid* mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling, int flag ) {
    
    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1, k1, cond;
    double eta, txx1, txz1, Pn, Tn, etaVE, VEcoeff=0.0, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1;
    double exx_pwl, exz_pwl, A2_pwl, exx_el, exz_el, exx_diss, exz_diss, exx_pl, exz_pl;
    int average = model->eta_avg;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    // Evaluate cell center viscosities
#pragma omp parallel for shared( mesh  ) private( cond, k, l, k1, p, eta, c1, c0, Pn, Tn, txx1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1, exx_pwl, exz_pwl, A2_pwl, exx_el, exz_el, exx_diss, exz_diss, exx_pl, exz_pl ) firstprivate( materials, scaling, flag, average, model, Ncx, Ncz )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {
        
        //    for ( l=0; l<Ncz; l++ ) {
        //        for ( k=0; k<Ncx; k++ ) {
        
        k      = mesh->kp[k1];
        l      = mesh->lp[k1];
        c0 = k  + l*(Ncx);
        
        mesh->eta_n[c0]   = 0.0;
        mesh->eta_phys_n[c0] = 0.0;
        mesh->VE_n[c0]       = 0.0;
        mesh->sxxd[c0]       = 0.0;
        mesh->eII_el[c0]     = 0.0;
        mesh->eII_pl[c0]     = 0.0;
        mesh->eII_pwl[c0]    = 0.0;
        mesh->eII_exp[c0]    = 0.0;
        mesh->eII_lin[c0]    = 0.0;
        mesh->eII_gbs[c0]    = 0.0;
        mesh->eII_cst[c0]    = 0.0;
        mesh->d[c0]          = 0.0;
        mesh->exx_pwl_n[c0]  = 0.0;
        mesh->exz_pwl_n[c0]  = 0.0;
        mesh->A2_pwl_n[c0]   = 0.0;
        mesh->exx_el[c0]     = 0.0;
        mesh->exx_diss[c0]   = 0.0;
        mesh->exx_pl[c0]     = 0.0;
        
        // Loop on grid nodes
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
            
            // If background pressure is removed in the RHS
            Pn = mesh->p_in[c0];
            Tn = mesh->T[c0];
            
            // Loop on phases
            for ( p=0; p<model->Nb_phases; p++) {
                
                cond =  fabs(mesh->phase_perc_n[p][c0])>1.0e-13;
                
                
                
                if ( cond == 1 ) {
                    eta =  Viscosity( p, mesh->mu_n[c0], Tn, Pn, mesh->d0[c0], mesh->phi[c0],  mesh->eii_n[c0],  mesh->tii0_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0], materials    , model, scaling, flag, &txx1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &exz_el, &exx_diss, &exz_diss, &exx_pl, &exz_pl, &d1, &A2_pwl, mesh->strain_n[c0], mesh->phi_n[c0], mesh->C_n[c0]);
                }
                
                // ARITHMETIC AVERAGE
                if (average == 0) {
                    if ( cond == 1 ) mesh->eta_n[c0]   += mesh->phase_perc_n[p][c0] * etaVE;
                    if ( cond == 1 ) mesh->eta_phys_n[c0] += mesh->phase_perc_n[p][c0] * eta;
                }
                if ( cond == 1 ) mesh->VE_n[c0]       += mesh->phase_perc_n[p][c0] * VEcoeff;
                if (average == 0 || average == 2) {
                    if ( cond == 1 ) mesh->sxxd[c0]   += mesh->phase_perc_n[p][c0] * txx1;
                }
                if ( cond == 1 ) mesh->eII_el[c0]     += mesh->phase_perc_n[p][c0] * eII_el;
                if ( cond == 1 ) mesh->eII_pl[c0]     += mesh->phase_perc_n[p][c0] * eII_pl;
                if ( cond == 1 ) mesh->eII_pwl[c0]    += mesh->phase_perc_n[p][c0] * eII_pwl;
                if ( cond == 1 ) mesh->eII_exp[c0]    += mesh->phase_perc_n[p][c0] * eII_exp;
                if ( cond == 1 ) mesh->eII_lin[c0]    += mesh->phase_perc_n[p][c0] * eII_lin;
                if ( cond == 1 ) mesh->eII_gbs[c0]    += mesh->phase_perc_n[p][c0] * eII_gbs;
                if ( cond == 1 ) mesh->eII_cst[c0]    += mesh->phase_perc_n[p][c0] * eII_cst;
                if ( cond == 1 ) mesh->d[c0]          += mesh->phase_perc_n[p][c0] * 1.0/d1;
                
                if ( cond == 1 ) mesh->exx_pwl_n[c0]  += mesh->phase_perc_n[p][c0] * exx_pwl;
                if ( cond == 1 ) mesh->exz_pwl_n[c0]  += mesh->phase_perc_n[p][c0] * exz_pwl;
                if ( cond == 1 ) mesh->A2_pwl_n[c0]   += mesh->phase_perc_n[p][c0] * A2_pwl;
                if ( cond == 1 ) mesh->exx_el[c0]     += mesh->phase_perc_n[p][c0] * exx_el;
                if ( cond == 1 ) mesh->exx_diss[c0]   += mesh->phase_perc_n[p][c0] * exx_diss;
                if ( cond == 1 ) mesh->exx_pl[c0]     += mesh->phase_perc_n[p][c0] * exx_pl;
                
                // HARMONIC AVERAGE
                if (average == 1) {
                    if ( cond == 1 ) mesh->sxxd[c0]       += mesh->phase_perc_n[p][c0] * 1.0/txx1;
                    if ( cond == 1 ) mesh->eta_n[c0]   += mesh->phase_perc_n[p][c0] * 1.0/etaVE;
                    if ( cond == 1 ) mesh->eta_phys_n[c0] += mesh->phase_perc_n[p][c0] * 1.0/eta;
                }
                
                // GEOMETRIC AVERAGE
                if (average == 2) {
                    if ( cond == 1 ) mesh->eta_n[c0]   += mesh->phase_perc_n[p][c0] * log(etaVE);
                    if ( cond == 1 ) mesh->eta_phys_n[c0] += mesh->phase_perc_n[p][c0] * log(eta);
                }
            }
            
            mesh->d[c0]          = 1.0/mesh->d[c0];
            
            // HARMONIC AVERAGE
            if (average == 1) {
                mesh->sxxd[c0]       = 1.0/mesh->sxxd[c0];
                mesh->eta_n[c0]   = 1.0/mesh->eta_n[c0];
                mesh->eta_phys_n[c0] = 1.0/mesh->eta_phys_n[c0];
                if (isinf (mesh->eta_phys_n[c0]) ) {
                    for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_n[p][c0]);
                    printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n", mesh->mu_n[c0], Tn, Pn, mesh->d0[c0], mesh->phi[c0],  mesh->eii_n[c0],  mesh->tii0_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0]);
                    printf("flag %d nb part cell = %d cell index = %d\n", mesh->BCp.type[c0],mesh->nb_part_cell[c0], c0);
                    printf("x=%2.2e z=%2.2e\n", mesh->xc_coord[k]*scaling->L/1000.0, mesh->zc_coord[l]*scaling->L/1000.0);
                    exit(1);
                }
                if (isnan (mesh->eta_phys_n[c0]) ) {
                    for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_n[p][c0]);
                    printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n", mesh->mu_n[c0], Tn, Pn, mesh->d0[c0], mesh->phi[c0],  mesh->eii_n[c0],  mesh->tii0_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0]);
                    printf("flag %d nb part cell = %d cell index = %d\n", mesh->BCp.type[c0],mesh->nb_part_cell[c0], c0);
                    printf("x=%2.2e z=%2.2e\n", mesh->xc_coord[k]*scaling->L/1000.0, mesh->zc_coord[l]*scaling->L/1000.0);
                    exit(1);
                }
            }
            // GEOMETRIC AVERAGE
            if (average == 2) {
                mesh->eta_n[c0]   = exp(mesh->eta_n[c0]);
                mesh->eta_phys_n[c0] = exp(mesh->eta_phys_n[c0]);
            }
            
        }
        
    }
    
    //    MinMaxArray( mesh->eII_el,  scaling->E, Ncx*Ncz, "eII_el" );
    //    MinMaxArray( mesh->eII_pl,  scaling->E, Ncx*Ncz, "eII_pl" );
    //    MinMaxArray( mesh->eII_pwl, scaling->E, Ncx*Ncz, "eII_pwl" );
    //    MinMaxArray( mesh->eII_exp, scaling->E, Ncx*Ncz, "eII_exp" );
    
    // Calculate vertices viscosity
    double T_s, P_s, d0s, d1s, phis;
    
#pragma omp parallel for shared( mesh, model ) private( cond, k, l, k1, p, eta, c1, c0, P_s, T_s, txx1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d0s, d1s, phis, exx_pwl, exz_pwl, A2_pwl, exx_el, exz_el, exx_diss, exz_diss, exx_pl, exz_pl ) firstprivate( materials, scaling, flag, average, Ncx, Ncz )
    for ( k1=0; k1<Nx*Nz; k1++ ) {
        
        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c0 = k + l*(Ncx);
        c1 = k + l*Nx;
        
        mesh->VE_s[c1]       = 0.0;
        mesh->sxz[c1]        = 0.0;
        mesh->eta_phys_s[c1] = 0.0;
        mesh->eta_s[c1]   = 0.0;
        mesh->eII_pwl_s[c1]  = 0.0;
        mesh->exx_pwl_s[c1]  = 0.0;
        mesh->exz_pwl_s[c1]  = 0.0;
        mesh->A2_pwl_s[c1]   = 0.0;
        mesh->exz_el[c1]     = 0.0;
        mesh->exz_diss[c1]   = 0.0;
        mesh->eII_pl_s[c1]   = 0.0;
        mesh->exz_pl[c1]     = 0.0;
        
        
        if ( mesh->BCg.type[c1] != 30 ) {
            
            // INNER average T and P
            if (k>0 && k<Ncx && l>0 && l<Ncz) {
                T_s  = 0.25*( mesh->T[c0] +  mesh->T[c0-Ncx] +  mesh->T[c0-Ncx-1] +  mesh->T[c0-1] );
                P_s  = 0.25*( mesh->p_in[c0] +  mesh->p_in[c0-Ncx] +  mesh->p_in[c0-Ncx-1] +  mesh->p_in[c0-1]);
                d0s  = 0.25*( mesh->d0[c0] +  mesh->d0[c0-Ncx] +  mesh->d0[c0-Ncx-1] +  mesh->d0[c0-1]);
                phis = 0.25*( mesh->phi[c0] +  mesh->phi[c0-Ncx] +  mesh->phi[c0-Ncx-1] +  mesh->phi[c0-1]);
            }
            
            // WEST
            if (k==0 && (l>0 && l<Ncz) && model->isperiodic_x==0) {
                T_s  = 0.5*( mesh->T[c0] + mesh->T[c0-Ncx] );
                P_s  = 0.5*( mesh->p_in[c0] + mesh->p_in[c0-Ncx] );
                d0s  = 0.5*( mesh->d0[c0] + mesh->d0[c0-Ncx] );
                phis = 0.5*( mesh->phi[c0] + mesh->phi[c0-Ncx] );
            }
            
            // WEST - periodic
            if (k==0 && (l>0 && l<Ncz) && model->isperiodic_x==1) {
                T_s  = 0.25*( mesh->T[c0] + mesh->T[c0-Ncx] + mesh->T[c0+Ncx-1] + mesh->T[c0-1] );
                P_s  = 0.25*( mesh->p_in[c0] + mesh->p_in[c0-Ncx] + mesh->p_in[c0+Ncx-1] + mesh->p_in[c0-1] );
                d0s  = 0.25*( mesh->d0[c0] + mesh->d0[c0-Ncx] + mesh->d0[c0+Ncx-1] + mesh->d0[c0-1] );
                phis = 0.25*( mesh->phi[c0] + mesh->phi[c0-Ncx] + mesh->phi[c0+Ncx-1] + mesh->phi[c0-1] );
            }
            
            // EAST
            if (k==Ncx && (l>0 && l<Ncz) && model->isperiodic_x==0) {
                T_s  = 0.5*( mesh->T[c0-1] + mesh->T[c0-Ncx-1]);
                P_s  = 0.5*( mesh->p_in[c0-1] + mesh->p_in[c0-Ncx-1] );
                d0s  = 0.5*( mesh->d0[c0-1] + mesh->d0[c0-Ncx-1]);
                phis = 0.5*( mesh->phi[c0-1] + mesh->phi[c0-Ncx-1] );
            }
            
            // EAST - periodic
            if (k==Ncx && (l>0 && l<Ncz) && model->isperiodic_x==1) {
                T_s  = 0.25*( mesh->T[c0-1] + mesh->T[c0-Ncx-1] + mesh->T[c0-Ncx-Ncx] + mesh->T[c0-Ncx]);
                P_s  = 0.25*( mesh->p_in[c0-1] + mesh->p_in[c0-Ncx-1] + mesh->p_in[c0-Ncx-Ncx] + mesh->p_in[c0-Ncx] );
                d0s  = 0.25*( mesh->d0[c0-1] + mesh->d0[c0-Ncx-1] + mesh->d0[c0-Ncx-Ncx] + mesh->d0[c0-Ncx]);
                phis = 0.25*( mesh->phi[c0-1] + mesh->phi[c0-Ncx-1] + mesh->phi[c0-Ncx-Ncx] + mesh->phi[c0-Ncx] );
            }
            
            
            // SOUTH
            if (l==0 && (k>0 && k<Ncx)) {
                T_s  = 0.5*( mesh->T[c0]+ mesh->T[c0-1] );
                P_s  = 0.5*( mesh->p_in[c0] + mesh->p_in[c0-1] );
                d0s  = 0.5*( mesh->d0[c0]+ mesh->d0[c0-1] );
                phis = 0.5*( mesh->phi[c0] + mesh->phi[c0-1] );
            }
            
            // NORTH
            if (l==Ncz && (k>0 && k<Ncx)) {
                T_s  = 0.5*( mesh->T[c0-Ncx] + mesh->T[c0-Ncx-1] );
                P_s  = 0.5*( mesh->p_in[c0-Ncx] + mesh->p_in[c0-Ncx-1] );
                d0s  = 0.5*( mesh->d0[c0-Ncx] + mesh->d0[c0-Ncx-1] );
                phis = 0.5*( mesh->phi[c0-Ncx] + mesh->phi[c0-Ncx-1] );
            }
            
            // SOUTH-WEST
            if (l==0 && k==0 && model->isperiodic_x==0) {
                T_s  = mesh->T[c0];
                P_s  = mesh->p_in[c0];
                d0s  = mesh->d0[c0];
                phis = mesh->phi[c0];
            }
            
            // SOUTH-WEST - periodic
            if (l==0 && k==0 && model->isperiodic_x==1) {
                T_s  = 0.5*(mesh->T[c0] + mesh->T[c0+Ncx-1]);
                P_s  = 0.5*(mesh->p_in[c0] + mesh->p_in[c0+Ncx-1]);
                d0s  = 0.5*(mesh->d0[c0] + mesh->d0[c0+Ncx-1]);
                phis = 0.5*(mesh->phi[c0] + mesh->phi[c0+Ncx-1]);
            }
            
            // NORTH-WEST
            if (l==Ncz && k==0 && model->isperiodic_x==0) {
                T_s  = mesh->T[c0-Ncx];
                P_s  = mesh->p_in[c0-Ncx];
                d0s  = mesh->d0[c0-Ncx];
                phis = mesh->phi[c0-Ncx];
            }
            
            // NORTH-WEST - periodic
            if (l==Ncz && k==0 && model->isperiodic_x==1) {
                T_s  = 0.5*(mesh->T[c0-Ncx] + mesh->T[c0-1]);
                P_s  = 0.5*(mesh->p_in[c0-Ncx] + mesh->p_in[c0-1]);
                d0s  = 0.5*(mesh->d0[c0-Ncx] + mesh->d0[c0-1]);
                phis = 0.5*(mesh->phi[c0-Ncx] + mesh->phi[c0-1]);
            }
            
            // SOUTH-EAST
            if (l==0 && k==Ncx && model->isperiodic_x==0) {
                T_s  = mesh->T[c0-1];
                P_s  = mesh->p_in[c0-1];
                d0s  = mesh->d0[c0-1];
                phis = mesh->phi[c0-1];
            }
            
            // SOUTH-EAST - periodic
            if (l==0 && k==Ncx && model->isperiodic_x==1) {
                T_s  = 0.5*(mesh->T[c0-1] + mesh->T[c0-Ncx]);
                P_s  = 0.5*(mesh->p_in[c0-1]+ mesh->p_in[c0-Ncx]);
                d0s  = 0.5*(mesh->d0[c0-1] + mesh->d0[c0-Ncx]);
                phis = 0.5*(mesh->phi[c0-1]+ mesh->phi[c0-Ncx]);
            }
            
            // NORTH-EAST
            if (l==Ncz && k==Ncx && model->isperiodic_x==0) {
                T_s  = mesh->T[c0-Ncx-1];
                P_s  = mesh->p_in[c0-Ncx-1];
                d0s  = mesh->d0[c0-Ncx-1];
                phis = mesh->phi[c0-Ncx-1];
            }
            
            // NORTH-EAST - periodic
            if (l==Ncz && k==Ncx && model->isperiodic_x==1) {
                T_s  = 0.5*(mesh->T[c0-Ncx-1] + mesh->T[c0-Ncx-Ncx]);
                P_s  = 0.5*(mesh->p_in[c0-Ncx-1] + mesh->p_in[c0-Ncx-Ncx]);
                d0s  = 0.5*(mesh->d0[c0-Ncx-1] + mesh->d0[c0-Ncx-Ncx]);
                phis = 0.5*(mesh->phi[c0-Ncx-1] + mesh->phi[c0-Ncx-Ncx]);
            }
            
            for ( p=0; p<model->Nb_phases; p++) {
                
                cond = fabs(mesh->phase_perc_s[p][c1])>1.0e-13;
                
                if ( cond == 1 ) {
                    eta =  Viscosity( p, mesh->mu_s[c1], T_s, P_s, d0s, phis,  mesh->eii_s[c1],  mesh->tii0_s[c1], mesh->exxd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->sxz0[c1], materials, model, scaling, flag, &txx1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &exz_el, &exx_diss, &exz_diss, &exx_pl, &exz_pl, &d1s, &A2_pwl, mesh->strain_s[c1], mesh->phi_s[c1], mesh->C_s[c1]);
                }
                
                if (average ==0) {
                    if ( cond == 1 ) mesh->eta_s[c1]   += mesh->phase_perc_s[p][c1] * etaVE;
                    if ( cond == 1 ) mesh->eta_phys_s[c1] += mesh->phase_perc_s[p][c1] * eta;
                }
                if ( cond == 1 ) mesh->VE_s[c1]       += mesh->phase_perc_s[p][c1] * VEcoeff;
                if (average ==0 || average==2 ) {
                    if ( cond == 1 ) mesh->sxz[c1]    += mesh->phase_perc_s[p][c1] * txz1;
                }
                if ( cond == 1 ) mesh->eII_pl_s[c1]   += mesh->phase_perc_s[p][c1] * eII_pl;
                if ( cond == 1 ) mesh->eII_pwl_s[c1]  += mesh->phase_perc_s[p][c1] * eII_pwl;
                if ( cond == 1 ) mesh->exx_pwl_s[c1]  += mesh->phase_perc_s[p][c1] * exx_pwl;
                if ( cond == 1 ) mesh->exz_pwl_s[c1]  += mesh->phase_perc_s[p][c1] * exz_pwl;
                if ( cond == 1 ) mesh->A2_pwl_s[c1]   += mesh->phase_perc_s[p][c1] * A2_pwl;
                if ( cond == 1 ) mesh->exz_el[c1]     += mesh->phase_perc_s[p][c1] * exz_el;
                if ( cond == 1 ) mesh->exz_diss[c1]   += mesh->phase_perc_s[p][c1] * exz_diss;
                if ( cond == 1 ) mesh->exz_pl[c1]     += mesh->phase_perc_s[p][c1] * exz_pl;
                
                
                if (average == 1) {
                    
                    if ( cond == 1  ) mesh->sxz[c1]        += mesh->phase_perc_s[p][c1] *  1.0/txz1;
                    if ( cond == 1  ) mesh->eta_s[c1]   += mesh->phase_perc_s[p][c1] *  1.0/etaVE;
                    if ( cond == 1  ) mesh->eta_phys_s[c1] += mesh->phase_perc_s[p][c1] *  1.0/eta;
                }
                if (average == 2) {
                    if ( cond == 1  ) mesh->eta_s[c1  ] += mesh->phase_perc_s[p][c1] *  log(etaVE);
                    if ( cond == 1  ) mesh->eta_phys_s[c1] += mesh->phase_perc_s[p][c1] *  log(eta);
                }
            }
            // HARMONIC AVERAGE
            if (average == 1) {
                mesh->sxz[c1]        = 1.0/mesh->sxz[c1];
                mesh->eta_s[c1]   = 1.0/mesh->eta_s[c1];
                mesh->eta_phys_s[c1] = 1.0/mesh->eta_phys_s[c1];
                if (isinf (mesh->eta_phys_s[c1]) ) {
                    for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_s[p][c1]);
                    printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n", mesh->mu_s[c1], T_s, P_s, d0s, phis,  mesh->eii_s[c1],  mesh->tii0_s[c1], mesh->exxd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->sxz0[c1]);
                    exit(1);
                }
                if (isnan (mesh->eta_phys_s[c1]) ) {
                    for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_s[p][c1]);
                    printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n", mesh->mu_s[c1], T_s, P_s, d0s, phis,  mesh->eii_s[c1],  mesh->tii0_s[c1], mesh->exxd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->sxz0[c1]);
                    exit(1);
                }
            }
            // GEOMETRIC AVERAGE
            if (average == 2) {
                mesh->eta_s[c1]   = exp(mesh->eta_s[c1]);
                mesh->eta_phys_s[c1] = exp(mesh->eta_phys_s[c1]);
            }
        }
    }
    
    //    printf("** Rheology cell centers/vertices ---> %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/



void CohesionFrictionGrid( grid* mesh, mat_prop materials, params model, scale scaling ) {
    
    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1;
    int average = 0;
    double phi, C, strain_acc;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    // Calculate cell centers cohesion and friction
    for ( l=0; l<Ncz; l++ ) {
        for ( k=0; k<Ncx; k++ ) {
            
            // Cell center index
            c0 = k  + l*(Ncx);
            
            // First - initialize to 0
            mesh->phi_n[c0] = 0.0;
            mesh->C_n[c0]   = 0.0;
            
            // Compute only if below free surface
            if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
                
                // Retrieve accumulated plastic strain
                strain_acc = mesh->strain_n[c0];
                
                // Loop on phases
                for ( p=0; p<model.Nb_phases; p++) {
                    
                    phi = materials.phi[p];
                    C   = materials.C[p];
                    
                    // Apply strain softening
                    if ( model.isPl_soft == 1 ) {
                        
                        // If we are below the lower strain limit
                        if (strain_acc < materials.pls_start[p]) {
                            //                            printf("Allo 1 - %d!\n", p);
                            
                            phi = materials.phi[p];
                            C   = materials.C[p];
                        }
                        // If we are above the upper strain limit
                        if (strain_acc >= materials.pls_end[p]) {
                            //                            printf("Allo 2 - %d!\n", p);
                            phi = materials.phi_end[p];
                            C   = materials.C_end[p];
                        }
                        // If we are in the softening strain range
                        if (strain_acc >= materials.pls_start[p] && strain_acc < materials.pls_end[p] ) {
                            //                            printf("Allo 3 - %d!\n", p);
                            phi = materials.phi[p] - (materials.phi[p] - materials.phi_end[p]) * ( strain_acc / (materials.pls_end[p] - materials.pls_start[p]) );
                            C   = materials.C[p]   + (  materials.C_end[p] -   materials.C[p]) * MINV( 1.0, strain_acc /materials.pls_end[p] );
                            
                        }
                        
                    }
                    
                    // Arithmetic
                    if (average ==0) {
                        mesh->phi_n[c0] += mesh->phase_perc_n[p][c0] * phi;
                        mesh->C_n[c0]   += mesh->phase_perc_n[p][c0] * C;
                    }
                    // Harmonic
                    if (average == 1) {
                        mesh->phi_n[c0] += mesh->phase_perc_n[p][c0] *  1.0/phi;
                        mesh->C_n[c0]   += mesh->phase_perc_n[p][c0] *  1.0/C;
                    }
                    // Geometric
                    if (average == 2) {
                        mesh->phi_n[c0] += mesh->phase_perc_n[p][c0] *  log(phi);
                        mesh->C_n[c0]   += mesh->phase_perc_n[p][c0] *  log(C);
                        
                    }
                }
                // Post-process for geometric/harmonic averages
                if ( average==1 ) mesh->phi_n[c0] = 1.0/mesh->phi_n[c0];
                if ( average==2 ) mesh->phi_n[c0] = exp(mesh->phi_n[c0]);
                if ( average==1 ) mesh->C_n[c0]   = 1.0/mesh->C_n[c0];
                if ( average==2 ) mesh->C_n[c0]   = exp(mesh->C_n[c0]);
                
            }
        }
    }
    
    
    // Calculate vertices cohesion and friction
    for ( l=0; l<Nz; l++ ) {
        for ( k=0; k<Nx; k++ ) {
            
            // Vertex index
            c1 = k + l*Nx;
            
            // First - initialize to 0
            mesh->phi_s[c1] = 0.0;
            mesh->C_s[c1]   = 0.0;
            
            // Compute only if below free surface
            if ( mesh->BCg.type[c1] != 30 ) {
                
                // Retrieve accumulated plastic strain
                strain_acc = mesh->strain_s[c1];
                
                // Loop on phases
                for ( p=0; p<model.Nb_phases; p++) {
                    
                    phi        = materials.phi[p];
                    C          = materials.C[p];
                    
                    // Apply strain softening
                    if ( model.isPl_soft == 1 ) {
                        
                        // If we are below the lower strain limit
                        if (strain_acc < materials.pls_start[p]) {
                            phi = materials.phi[p];
                            C   = materials.C[p];
                        }
                        // If we are above the upper strain limit
                        if (strain_acc >= materials.pls_end[p]) {
                            phi = materials.phi_end[p];
                            C   = materials.C_end[p];
                        }
                        // If we are in the softening strain range
                        if (strain_acc >= materials.pls_start[p] && strain_acc < materials.pls_end[p] ) {
                            phi = materials.phi[p] - (materials.phi[p] - materials.phi_end[p]) * ( strain_acc / (materials.pls_end[p] - materials.pls_start[p]) );
                            C   = materials.C[p]   - (  materials.C[p] -   materials.C_end[p]) * ( strain_acc / (materials.pls_end[p] - materials.pls_start[p]) );
                        }
                    }
                    
                    // Arithmetic
                    if (average ==0) {
                        mesh->phi_s[c1] += mesh->phase_perc_s[p][c1] * phi;
                        mesh->C_s[c1]   += mesh->phase_perc_s[p][c1] * C;
                    }
                    // Harmonic
                    if (average == 1) {
                        mesh->phi_s[c1] += mesh->phase_perc_s[p][c1] *  1.0/phi;
                        mesh->C_s[c1]   += mesh->phase_perc_s[p][c1] *  1.0/C;
                    }
                    // Geometric
                    if (average == 2) {
                        mesh->phi_s[c1] += mesh->phase_perc_s[p][c1] *  log(phi);
                        mesh->C_s[c1]   += mesh->phase_perc_s[p][c1] *  log(C);
                    }
                }
                // Post-process for geometric/harmonic averages
                if ( average==1 ) mesh->phi_s[c1] = 1.0/mesh->phi_s[c1];
                if ( average==2 ) mesh->phi_s[c1] = exp(mesh->phi_s[c1]);
                if ( average==1 ) mesh->C_s[c1]   = 1.0/mesh->C_s[c1];
                if ( average==2 ) mesh->C_s[c1]   = exp(mesh->C_s[c1]);
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ShearModulusGrid( grid* mesh, mat_prop materials, params model, scale scaling ) {
    
    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1;
    int average = model.eta_avg;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    // Calculate cell centers shear modulus
    for ( l=0; l<Ncz; l++ ) {
        for ( k=0; k<Ncx; k++ ) {
            
            // Cell center index
            c0 = k  + l*(Ncx);
            
            // First - initialize to 0
            mesh->mu_n[c0] = 0.0;
            
            // Compute only if below free surface
            if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
                
                // Loop on phases
                for ( p=0; p<model.Nb_phases; p++) {
                    
                    // Arithmetic
                    if (average ==0) {
                        mesh->mu_n[c0] += mesh->phase_perc_n[p][c0] * materials.mu[p];
                    }
                    // Harmonic
                    if (average == 1) {
                        mesh->mu_n[c0] += mesh->phase_perc_n[p][c0] *  1.0/materials.mu[p];
                    }
                    // Geometric
                    if (average == 2) {
                        mesh->mu_n[c0] += mesh->phase_perc_n[p][c0] *  log(materials.mu[p]);
                    }
                }
                // Post-process for geometric/harmonic averages
                if ( average==1 ) mesh->mu_n[c0] = 1.0/mesh->mu_n[c0];
                if ( average==2 ) mesh->mu_n[c0] = exp(mesh->mu_n[c0]);
                
            }
        }
    }
    
    
    // Calculate vertices shear modulus
    for ( l=0; l<Nz; l++ ) {
        for ( k=0; k<Nx; k++ ) {
            
            // Vertex index
            c1 = k + l*Nx;
            
            // First - initialize to 0
            mesh->mu_s[c1] = 0.0;
            
            // Compute only if below free surface
            if ( mesh->BCg.type[c1] != 30 ) {
                
                // Loop on phases
                for ( p=0; p<model.Nb_phases; p++) {
                    
                    // Arithmetic
                    if (average ==0) {
                        mesh->mu_s[c1] += mesh->phase_perc_s[p][c1] * materials.mu[p];
                    }
                    // Harmonic
                    if (average == 1) {
                        mesh->mu_s[c1] += mesh->phase_perc_s[p][c1] *  1.0/materials.mu[p];
                    }
                    // Geometric
                    if (average == 2) {
                        mesh->mu_s[c1] += mesh->phase_perc_s[p][c1] *  log(materials.mu[p]);
                    }
                }
                // Post-process for geometric/harmonic averages
                if ( average==1 ) mesh->mu_s[c1] = 1.0/mesh->mu_s[c1];
                if ( average==2 ) mesh->mu_s[c1] = exp(mesh->mu_s[c1]);
                
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateDensity( grid* mesh, markers* particles, mat_prop *materials, params *model, scale *scaling ) {
    
    int k, p, c0, Nx, Nz, Ncx, Ncz;
    int    phase_diag;
    double rho0, T0, alpha, P0, beta, drho, Tgrid, Pgrid, dT, dP;
    double percT, percP;
    double rhonew, rhop, epsi = 1e-13;
    
    int iT, iP, NT, NP, iSW, iSE, iNW, iNE;
    double dstT, dstP;
    double PW, TW;
    
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    // Get X on the cell centers
    Interp_P2C ( *particles, particles->X, mesh, mesh->X, mesh->xg_coord, mesh->zg_coord, 1, 0 );
    
#pragma omp parallel for shared(mesh,materials) private( NT, NP, iT, iP, TW, PW, iSW, iSE, iNW, iNE, phase_diag, dT, dP, Tgrid, Pgrid, c0, p, rhonew, rho0, alpha, beta, T0, P0, rhop, drho, percT, percP) firstprivate(Ncx, Ncz, model, epsi)
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {
        rhonew = 0.0;
        // Loop on phases
        for ( p=0; p<model->Nb_phases; p++) {
            
            if ( fabs(mesh->phase_perc_n[p][c0])>epsi) {
                
                // Constant density
                if ( materials->density_model[p] == 0 ) {
                    rhop = materials->rho[p];
                }
                
                // T and P dependent density based on EOS
                if ( materials->density_model[p] == 1 ) {
                    
                    rho0   = materials->rho [p];
                    drho   = materials->drho[p];
                    T0     = materials->T0 [p];
                    alpha  = materials->alp[p];
                    P0     = materials->P0 [p];
                    beta   = materials->bet[p];
                    rhop   = (1.0 -  alpha * (mesh->T[c0] - T0) ) * (1.0 +  beta * (mesh->p_in[c0] - P0) ); // EOS general
                    rhop   = ((1.0-mesh->X[c0])*rho0 + mesh->X[c0]*(rho0+drho))*rhop; // Average density based on X
                }
                
                // T and P dependent density based on phase diagrams
                if ( materials->density_model[p] == 2 ) {
                    // Identify current phase diagram
                    phase_diag = materials->phase_diagram[p];
                    // Determine T and P increments in the current phase diagram
                    NT         = model->PDMnT[phase_diag];
                    NP         = model->PDMnP[phase_diag];
                    dT         = (model->PDMTmax[phase_diag] - model->PDMTmin[phase_diag])/(NT-1.0);
                    dP         = (model->PDMPmax[phase_diag] - model->PDMPmin[phase_diag])/(NP-1.0);
                    // Pressure and temperature + correct to remain within the database bounds
                    Tgrid      = mesh->T[c0];
                    Pgrid      = mesh->p_in[c0];
                    if (Tgrid<model->PDMTmin[phase_diag]) Tgrid = model->PDMTmin[phase_diag] + 0.01*dT;
                    if (Tgrid>model->PDMTmax[phase_diag]) Tgrid = model->PDMTmax[phase_diag] - 0.01*dT;
                    if (Pgrid<model->PDMPmin[phase_diag]) Pgrid = model->PDMPmin[phase_diag] + 0.01*dP;
                    if (Pgrid>model->PDMPmax[phase_diag]) Pgrid = model->PDMPmax[phase_diag] - 0.01*dP;
                    // Find index of minimum/west temperature node
                    dstT = ( Tgrid - model->PDMTmin[phase_diag] );
                    iT   = ceil( dstT/dT - 0.0 ) - 1;
                    // Find index of minimum/west pressure node
                    dstP = ( Pgrid - model->PDMPmin[phase_diag] );
                    iP   = ceil( dstP/dP  - 0.0 ) - 1;
                    // Calculate Weighting coeeficients fo bilinear interpolant
                    TW    = (model->PDMTmin[phase_diag] + iT*dT);
                    PW    = (model->PDMPmin[phase_diag] + iP*dP);
                    percT = 1.0 - (Tgrid - TW )/dT;
                    percP = 1.0 - (Pgrid - PW )/dP;
                    // Indices of neigbours
                    iSW   = iT + iP*NT;
                    iSE   = iT + iP*NT+1;
                    iNW   = iT + (iP+1)*NT;
                    iNE   = iT + (iP+1)*NT+1;
                    // Interpolate from 4 neighbours
                    rhop  = 0.0;
                    rhop +=  (1.0-percT)* (1.0-percP) * model->PDMrho[phase_diag][iSW];
                    rhop +=  (    percT)* (1.0-percP) * model->PDMrho[phase_diag][iSE];
                    rhop +=  (1.0-percT)* (    percP) * model->PDMrho[phase_diag][iNW];
                    rhop +=  (    percT)* (    percP) * model->PDMrho[phase_diag][iNE];
                    
                }
                // Average density base on phase density and phase volume fraction
                if ( mesh->BCp.type[c0] != 30 ) rhonew += mesh->phase_perc_n[p][c0] * rhop;
            }
        }
        mesh->rho_n[c0]   = rhonew;
    }
    
    // Interpolate center values to vertices
    int l, c1;
#pragma omp parallel for shared( mesh ) private( c1, k, l, c0) firstprivate( Nx, Nz, Ncx, Ncz )
    for (c1=0; c1<Nx*Nz; c1++) {
        
        k = mesh->kn[c1];
        l = mesh->ln[c1];
        c1 = k + l*Nx;
        c0 = k + l*(Ncx);
        mesh->rho_s[c1]   = 0.0;
        
        if ( mesh->BCg.type[c1] != 30 ) {
            
            // Inner grid
            if ( k>0 && l>0 && k<Ncx && l<Ncz ) {
                mesh->rho_s[c1] = 0.25*( mesh->rho_n[c0] + mesh->rho_n[c0-1] + mesh->rho_n[c0-Ncx] + mesh->rho_n[c0-Ncx-1]  );
            }
            // Sides
            if (k==0 && (l>0 && l<Ncz)) {
                mesh->rho_s[c1] = 0.5*( mesh->rho_n[c0] + mesh->rho_n[c0-Ncx] );
            }
            if (k==Ncx && (l>0 && l<Ncz)) {
                mesh->rho_s[c1] = 0.5*( mesh->rho_n[c0-1] + mesh->rho_n[c0-Ncx-1] );
            }
            if (l==0 && (k>0 && k<Ncx)) {
                mesh->rho_s[c1] = 0.5*( mesh->rho_n[c0] + mesh->rho_n[c0-1] );
            }
            if (l==Ncz && (k>0 && k<Ncx)) {
                mesh->rho_s[c1] = 0.5*( mesh->rho_n[c0-Ncx] + mesh->rho_n[c0-Ncx-1] );
            }
            // Corners
            if (l==0 && k==0) {
                mesh->rho_s[c1] = mesh->rho_n[c0];
            }
            if (l==Ncz && k==0) {
                mesh->rho_s[c1] = mesh->rho_n[c0-Ncx];
            }
            if (l==0 && k==Ncx) {
                mesh->rho_s[c1] = mesh->rho_n[c0-1];
            }
            if (l==Ncz && k==Ncx) {
                mesh->rho_s[c1] = mesh->rho_n[c0-Ncx-1];
            }
        }
    }
    
    printf("Updated density fields:\n");
    MinMaxArrayTag(        mesh->X,            1, Ncx*Ncz,     "X", mesh->BCp.type );
    MinMaxArrayTag( mesh->rho_n, scaling->rho, Ncx*Ncz, "rho_n", mesh->BCp.type );
    MinMaxArrayTag( mesh->rho_s, scaling->rho, Nx*Nz,   "rho_s", mesh->BCg.type    );
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Strain rate
void StrainRateComponents( grid* mesh, scale scaling, params* model ) {
    
    int k, l, c0, c1, c2, Nx, Nz, Ncx, Ncz, k1;
    double dx, dz;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    dx = mesh->dx;
    dz = mesh->dz;
    
#pragma omp parallel for shared( mesh ) private( k, k1, l, c0, c1, c2  ) firstprivate( dx, dz, Nx, Ncx, Ncz )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {
        k  = mesh->kp[k1];
        l  = mesh->lp[k1];
        c0 = k  + l*(Nx-1);
        c1 = k  + l*(Nx);
        c2 = k  + l*(Nx+1);
        
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
            
            // Velocity divergence
            mesh->div_u[c0] = (mesh->u_in[c1+1+Nx] - mesh->u_in[c1+Nx])/dx + (mesh->v_in[c2+Nx+1+1] - mesh->v_in[c2+1])/dz;
            
            // Normal strain rates
            mesh->exxd[c0]  = (mesh->u_in[c1+1+Nx]     - mesh->u_in[c1+Nx] )/dx;
        }
        else {
            mesh->div_u[c0] = 0.0;
            mesh->exxd[c0]  = 0.0;
        }
    }
    
    // Shear components : we only calculate sxz because sxz = szx (stress tensor symmetry)
#pragma omp parallel for shared( mesh ) private( k, k1, l, c1, c2  ) firstprivate( dx, dz, Nx, Nz )
    for ( k1=0; k1<Nx*Nz; k1++ ) {
        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k  + l*(Nx);
        c2 = k  + l*(Nx+1);
        mesh->exz[c1] = 0.0;
        
        if ( mesh->BCg.type[c1] != 30 ) {
            if (mesh->BCu.type[c1] != 30 && mesh->BCu.type[c1+Nx] != 30) mesh->exz[c1] +=  0.5 * (mesh->u_in[c1+Nx] - mesh->u_in[c1])/dz;
            if (mesh->BCv.type[c2] != 30 && mesh->BCv.type[c2+1]  != 30) mesh->exz[c1] +=  0.5 * (mesh->v_in[c2+1]  - mesh->v_in[c2])/dx;
        }
        else {
            mesh->exz[c1] = 0.0;
        }
    }
    
    double sum=0.0;
#pragma omp parallel for shared( mesh ) private( k, k1, l, c0, c1, sum ) firstprivate( Nx, Ncx, Ncz )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {
        
        k  = mesh->kp[k1];
        l  = mesh->lp[k1];
        c0 = k  + l*(Nx-1);
        c1 = k  + l*(Nx);
        mesh->exz_n[c0] = 0.0;
        sum = 0.0;
        
        if ( mesh->BCp.type[c0]      != 30 && mesh->BCp.type[c0] != 31) {
            if (mesh->BCg.type[c1]      != 30 ) { mesh->exz_n[c0] += mesh->exz[c1];      sum++;}
            if (mesh->BCg.type[c1+1]    != 30 ) { mesh->exz_n[c0] += mesh->exz[c1+1];    sum++;}
            if (mesh->BCg.type[c1+Nx]   != 30 ) { mesh->exz_n[c0] += mesh->exz[c1+Nx];   sum++;}
            if (mesh->BCg.type[c1+Nx+1] != 30 ) { mesh->exz_n[c0] += mesh->exz[c1+Nx+1]; sum++;}
            if (sum>0) mesh->exz_n[c0] /= sum;
        }
        else {
            mesh->exz_n[c0] = 0.0;
        }
    }
    
    // Interpolate normal strain rate on vertices
#pragma omp parallel for shared( mesh ) private( k, k1, l, c1, c0, sum  ) firstprivate( dx, dz, Nx, Nz, Ncx, Ncz, model )
    for ( k1=0; k1<Nx*Nz; k1++ ) {
        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k  + l*(Nx);
        c0 = k  + l*(Nx-1);
        
        mesh->exxd_s[c1] = 0.0;
        sum = 0.0;
        
        if ( mesh->BCg.type[c1] != 30 ) {
            
            // INNER
            if ( k > 0 && k < Ncx && l > 0 && l < Ncz ) {
                if ( mesh->BCp.type[c0]       != 30 && mesh->BCp.type[c0]       != 31) { mesh->exxd_s[c1] += mesh->exxd[c0]; sum++;}
                if ( mesh->BCp.type[c0-1]     != 30 && mesh->BCp.type[c0-1]     != 31) { mesh->exxd_s[c1] += mesh->exxd[c0-1]; sum++;}
                if ( mesh->BCp.type[c0-Ncx]   != 30 && mesh->BCp.type[c0-Ncx]   != 31) { mesh->exxd_s[c1] += mesh->exxd[c0-Ncx]; sum++;}
                if ( mesh->BCp.type[c0-Ncx-1] != 30 && mesh->BCp.type[c0-Ncx-1] != 31) { mesh->exxd_s[c1] += mesh->exxd[c0-Ncx-1]; sum++;};
            }
            
            // WEST
            if ( k == 0 && l > 0 && l < Ncz && model->isperiodic_x==0  ) {
                if ( mesh->BCp.type[c0]     != 30 && mesh->BCp.type[c0]     != 31) {mesh->exxd_s[c1] += mesh->exxd[c0]; sum++;};
                if ( mesh->BCp.type[c0-Ncx] != 30 && mesh->BCp.type[c0-Ncx] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-Ncx]; sum++;};
            }
            
            // Periodic
            if ( k == 0 && l > 0 && l < Ncz && model->isperiodic_x==1 ) {
                if ( mesh->BCp.type[c0]       != 30 && mesh->BCp.type[c0]       != 31) {mesh->exxd_s[c1] += mesh->exxd[c0]; sum++;};
                if ( mesh->BCp.type[c0-Ncx]   != 30 && mesh->BCp.type[c0-Ncx]   != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-Ncx]; sum++;};
                if ( mesh->BCp.type[c0-1]     != 30 && mesh->BCp.type[c0-1]     != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-1]; sum++;};
                if ( mesh->BCp.type[c0+Ncx-1] != 30 && mesh->BCp.type[c0+Ncx-1] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0+Ncx-1]; sum++;};
            }
            
            // EAST
            if ( k == Ncx && l > 0 && l< Ncz && model->isperiodic_x==0  ) {
                if ( mesh->BCp.type[c0-1]     != 30 && mesh->BCp.type[c0-1]     != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-1]; sum++;};
                if ( mesh->BCp.type[c0-Ncx-1] != 30 && mesh->BCp.type[c0-Ncx-1] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-Ncx-1]; sum++;};
            }
            
            // Periodic
            if ( k == Ncx && l > 0 && l< Ncz && model->isperiodic_x==1 ) {
                if ( mesh->BCp.type[c0-1]       != 30 && mesh->BCp.type[c0-1]       != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-1];       sum++;};
                if ( mesh->BCp.type[c0-Ncx-1]   != 30 && mesh->BCp.type[c0-Ncx-1]   != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-Ncx-1];   sum++;};
                if ( mesh->BCp.type[c0-Ncx-Ncx] != 30 && mesh->BCp.type[c0-Ncx-Ncx] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-Ncx-Ncx]; sum++;};
                if ( mesh->BCp.type[c0-Ncx  ]   != 30 && mesh->BCp.type[c0-Ncx  ]   != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-Ncx];     sum++;};
            }
            
            // NORTH
            if ( l == Ncz && k> 0 && k < Ncx ) {
                if ( mesh->BCp.type[c0-Ncx-1] != 30 && mesh->BCp.type[c0-Ncx-1] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-Ncx-1]; sum++;};
                if ( mesh->BCp.type[c0-Ncx]   != 30 && mesh->BCp.type[c0-Ncx]   != 31) {mesh->exxd_s[c1] +=mesh->exxd[c0-Ncx]; sum++;};
            }
            
            // SOUTH
            if ( l == 0 && k > 0 && k< Ncx ) {
                if ( mesh->BCp.type[c0]   != 30 && mesh->BCp.type[c0]   != 31) {mesh->exxd_s[c1] += mesh->exxd[c0]; sum++;};
                if ( mesh->BCp.type[c0-1] != 30 && mesh->BCp.type[c0-1] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-1]; sum++;};
            }
            
            // NORTH-EAST
            if ( k == Ncx && l == Ncz && model->isperiodic_x==0  ) {
                if ( mesh->BCp.type[c0-Ncx-1] != 30 && mesh->BCp.type[c0-Ncx-1] != 31) {mesh->exxd_s[c1] = mesh->exxd[c0-Ncx-1]; sum++;};
            }
            
            if ( k == Ncx && l == Ncz && model->isperiodic_x==1 ) {
                if ( mesh->BCp.type[c0-Ncx-1]   != 30 && mesh->BCp.type[c0-Ncx-1  ] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-Ncx-1];   sum++;};
                if ( mesh->BCp.type[c0-Ncx-Ncx] != 30 && mesh->BCp.type[c0-Ncx-Ncx] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-Ncx-Ncx]; sum++;};
            }
            
            // NORTH-WEST
            if ( k == 0 && l== Ncz && model->isperiodic_x==0 ) {
                if ( mesh->BCp.type[c0-Ncx] != 30 && mesh->BCp.type[c0-Ncx] != 31) {mesh->exxd_s[c1] = mesh->exxd[c0-Ncx]; sum++;};
            }
            
            if ( k == 0 && l== Ncz && model->isperiodic_x==1  ) {
                if ( mesh->BCp.type[c0-Ncx] != 30 && mesh->BCp.type[c0-Ncx] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-Ncx]; sum++;};
                if ( mesh->BCp.type[c0-1  ] != 30 && mesh->BCp.type[c0-1  ] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-1  ]; sum++;};
            }
            
            // SOUTH-EAST
            if ( k == Ncx && l == 0 && model->isperiodic_x==0  ) {
                if ( mesh->BCp.type[c0-1] != 30 && mesh->BCp.type[c0-1] != 31) {mesh->exxd_s[c1] = mesh->exxd[c0-1]; sum++;};
            }
            
            if ( k == Ncx && l == 0 && model->isperiodic_x==1  ) {
                if ( mesh->BCp.type[c0-1]   != 30 && mesh->BCp.type[c0-1]   != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-1  ]; sum++;};
                if ( mesh->BCp.type[c0-Ncx] != 30 && mesh->BCp.type[c0-Ncx] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0-Ncx]; sum++;};
                
            }
            
            // SOUTH-WEST
            if ( k == 0 && l== 0 && model->isperiodic_x==0  ) {
                if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {mesh->exxd_s[c1] = mesh->exxd[c0]; sum++;};
            }
            
            if ( k == 0 && l== 0 && model->isperiodic_x==1  ) {
                if ( mesh->BCp.type[c0]       != 30 && mesh->BCp.type[c0      ] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0      ]; sum++;};
                if ( mesh->BCp.type[c0+Ncx-1] != 30 && mesh->BCp.type[c0+Ncx-1] != 31) {mesh->exxd_s[c1] += mesh->exxd[c0+Ncx-1]; sum++;};
                
            }
            
            if (sum>0) mesh->exxd_s[c1] /= sum;
        }
        else {
            mesh->exxd_s[c1] = 0.0;
        }
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void GenerateDeformationMaps( grid* mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling, int flag ) {
    
    // This functions generates deformation maps and writem to disk
    
    //Definition of parameters and allocation of memory
    int nT = model->nT, nE = model->nE, nd = model->nd, ix, iy, iz;
    double stepT, stepE, stepd;
    double Tmin = model->Tmin;
    double Tmax = model->Tmax;
    double Emin = model->Emin;
    double Emax = model->Emax;
    double dmin = model->dmin;
    double dmax = model->dmax;
    
    double *T, *E, *d, *dlog, *Elog, gs_ref;
    double txx1, txz1, etaVE, VEcoeff, eII_el=0.0, eII_pl=0.0, eII_pwl=0.0, eII_exp=0.0, eII_lin=0.0, eII_gbs=0.0, eII_cst=0.0, d1;
    double exx_pwl, exz_pwl, A2_pwl, exx_el, exz_el, exx_diss, exz_diss, exx_pl, exz_pl;
    double Pn = model->Pn , eta;
    int    *dom_mech, ind, mech, loud=0, k;
    double *stress, *visco, t_omp;
    
    T    = malloc(nT*sizeof(double));
    Elog = malloc(nE*sizeof(double));
    E    = malloc(nE*sizeof(double));
    dlog = malloc(nd*sizeof(double));
    d    = malloc(nd*sizeof(double));
    
    //Initialization of arrays
    stepT = (Tmax-Tmin)/(nT-1);
    stepE = (Emax-Emin)/(nE-1);
    stepd = (dmax-dmin)/(nd-1);
    
    // Temperature vector
    for ( ix=0; ix<nT; ix++) {
        if (ix==0) T[ix] = Tmin;
        else       T[ix] = T[ix-1] + stepT;
    }
    // Strain rate vector
    for ( iy=0; iy<nE; iy++) {
        if (iy==0) Elog[iy] = Emin;
        else       Elog[iy] = Elog[iy-1] + stepE;
        E[iy] = pow(10.0,Elog[iy])/scaling->E;
    }
    // Grain size vector
    for (iz=0; iz<nd; iz++) {
        if (iz==0) dlog[iz] = dmin;
        else       dlog[iz] = dlog[iz-1] + stepd;
        d[iz] = pow(10.0,dlog[iz])/scaling->L;
    }
    
    double Flin, Fgbs, Fpwl, texp;
    
    // Loop on all phases
    for (k=0;k<model->Nb_phases; k++) {
        
        printf("Generating deformation map of phase %2d\n", k);
        
        // Save real values
        Flin   = materials->Flin[k];
        Fgbs   = materials->Fgbs[k];
        Fpwl   = materials->Fpwl[k];
        texp   = materials->texp[k];
        gs_ref = materials->gs_ref[k];
        
        // Set no correction for deformation maps
        materials->Flin[k]   = 1.0;
        materials->Fgbs[k]   = 1.0;
        materials->Fpwl[k]   = 1.0;
        materials->texp[k]   = 0.0;
        
        // Allocate maps
        dom_mech = malloc(nT*nE*nd*sizeof(int));
        stress   = malloc(nT*nE*nd*sizeof(double));
        visco    = malloc(nT*nE*nd*sizeof(double));
        
        t_omp = (double)omp_get_wtime();
        
        // Boucles spatiales 2D pour créations des carte de déformation
        
        for ( iz=0; iz<nd; iz++) {
            
            // Force grain size
            materials->gs_ref[k] = d[iz];
            
            for ( ix=0; ix<nT; ix++) {
                for ( iy=0; iy<nE; iy++) {
                    
                    // Evaluate viscosity and stress
                    eta =  Viscosity( k, 0.0, T[ix], Pn, d[iz], 0.0, E[iy], 0, E[iy], 0.0, 0.0, 0.0, materials, model, scaling, 1, &txx1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &exz_el, &exx_diss, &exz_diss, &exx_pl, &exz_pl,  &d1, &A2_pwl, 0.0, materials->phi[k], materials->C[k]);
                    
                    // Select mechanism
                    if (eII_pwl>eII_el && eII_pwl>eII_pl  && eII_pwl>eII_exp && eII_pwl>eII_lin && eII_pwl>eII_gbs && eII_pwl>eII_cst) mech = 1; // dislocation
                    if (eII_lin>eII_el && eII_lin>eII_pl  && eII_lin>eII_exp && eII_lin>eII_pwl && eII_lin>eII_gbs && eII_lin>eII_cst) mech = 2; // diffusion
                    if (eII_gbs>eII_el && eII_gbs>eII_pl  && eII_gbs>eII_lin && eII_gbs>eII_pwl && eII_gbs>eII_exp && eII_gbs>eII_cst) mech = 3; // gbs
                    if (eII_exp>eII_el && eII_exp>eII_pl  && eII_exp>eII_lin && eII_exp>eII_pwl && eII_exp>eII_gbs && eII_exp>eII_cst) mech = 4; // peierls
                    if (eII_pl >eII_el && eII_pl >eII_pwl && eII_pl >eII_exp && eII_pl >eII_lin && eII_pl >eII_gbs && eII_pl >eII_cst) mech = 5; // plastic
                    if (eII_cst>eII_pl && eII_cst>eII_pwl && eII_cst>eII_exp && eII_cst>eII_lin && eII_cst>eII_gbs && eII_cst>eII_el ) mech = 6; // constant viscosity
                    if (eII_el >eII_pl && eII_el >eII_pwl && eII_el >eII_exp && eII_el >eII_lin && eII_el >eII_gbs && eII_el >eII_cst) mech = 7; // elastic
                    
                    // Global index for 3D matrix flattened in a 1D vector
                    ind           = ix + iy*nT + iz*nE*nT;
                    dom_mech[ind] = mech;
                    stress[ind]   = txx1;
                    visco[ind]    = eta;
                }
            }
        }
        
        printf("** Deformation map %2d ---> %lf sec\n", k, (double)((double)omp_get_wtime() - t_omp));
        
        
        // Print deformation maps to screen
        if (loud ==1) {
            for ( iz=0; iz<nd; iz++) {
                printf("GS = %2.2e\n", d[iz]*scaling->L);
                for ( ix=0; ix<nT; ix++) {
                    if (ix==0) printf("           ");
                    printf("%2.2lf ", T[ix]*scaling->T);
                }
                printf("\n");
                for ( iy=0; iy<nE; iy++) {
                    printf("%2.2e   ", E[iy]*scaling->E);
                    for ( ix=0; ix<nT; ix++) {
                        ind = ix + iy*nT + iz*nd*nT;
                        
                        printf("%6d ", dom_mech[ind] );
                        
                    }
                    printf("\n");
                }
                printf("\n");
            }
        }
        
        // Print deformation maps to file

        char *filename;
        asprintf( &filename, "DefMap%02d.gzip.h5", k );
        double *CT, *Cd, *CE, *Cstress, *Cvisco;
        
        CT = DoodzMalloc( sizeof(double)*nT);
        ArrayEqualArray( CT, T, nT );
        //        DoubleToFloat( T, CT, nT );
        ScaleBackD( CT, scaling->T, nT );
        
        CE = DoodzMalloc( sizeof(double)*nE);
        ArrayEqualArray( CE, E, nE );
        //        DoubleToFloat( E, CE, nE );
        ScaleBackD( CE, scaling->E, nE );
        
        Cd = DoodzMalloc( sizeof(double)*nd);
        ArrayEqualArray( Cd, d, nd );
        //        DoubleToFloat( d, Cd, nd );
        ScaleBackD( Cd, scaling->L, nd );
        
        Cstress = DoodzMalloc( sizeof(double)*nT*nE*nd);
        ArrayEqualArray( Cstress, stress, nT*nE*nd );
        //        DoubleToFloat( stress, Cstress, nT*nE*nd );
        ScaleBackD( Cstress, scaling->S, nT*nE*nd );
        
        Cvisco= DoodzMalloc( sizeof(double)*nT*nE*nd);
        ArrayEqualArray( Cvisco, visco, nT*nE*nd );
        //        DoubleToFloat( visco, Cvisco, nT*nE*nd );
        ScaleBackD( Cvisco, scaling->eta, nT*nE*nd );
        
        // Fill in DD data structure
        double params[3];
        params[0] = nT;
        params[1] = nE;
        params[2] = nd;
        
        // Send data to file
        create_output_hdf5( filename );
        AddGroup_to_hdf5( filename, "model" );
        AddGroup_to_hdf5( filename, "arrays" );
        
        AddFieldToGroup_generic( _TRUE_, filename, "model", "params" , 'd',  3, params,  1 );
        AddFieldToGroup_generic( _TRUE_, filename, "arrays", "T"     , 'd', nT, CT,  1 );
        AddFieldToGroup_generic( _TRUE_, filename, "arrays", "E"     , 'd', nE, CE,  1 );
        AddFieldToGroup_generic( _TRUE_, filename, "arrays", "d"     , 'd', nd, Cd,  1 );
        AddFieldToGroup_generic( _TRUE_, filename, "arrays", "stress", 'd', nT*nE*nd, Cstress,  1 );
        AddFieldToGroup_generic( _TRUE_, filename, "arrays", "visco" , 'd', nT*nE*nd, Cvisco,  1 );
        AddFieldToGroup_generic( _TRUE_, filename, "arrays", "map"   , 'i', nT*nE*nd, dom_mech,  1 );
        
        
        free(filename);
        DoodzFree(CT);
        DoodzFree(CE);
        DoodzFree(Cd);
        DoodzFree(Cstress);
        DoodzFree(Cvisco);

        // Free memory 1
        free(dom_mech);
        free(stress);
        free(visco);
        
        // Set no correction for deformation maps
        materials->Flin[k]   = Flin;
        materials->Fgbs[k]   = Fgbs;
        materials->Fpwl[k]   = Fpwl;
        materials->texp[k]   = texp;
        materials->gs_ref[k] = gs_ref;
    }
    // Free memory 2
    free(T);
    free(E);
    free(d);
    free(dlog);
    free(Elog);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void RotateDirectorVector( grid mesh, markers* particles, params model, scale *scaling ) {
    
    int k;
    double angle, nx, nz, norm;
    
#pragma omp parallel for shared ( particles ) private ( angle, k, nx, nz, norm ) firstprivate( model ) schedule( static )
        for(k=0; k<particles->Nb_part; k++) {
            
            // Filter out particles that are inactive (out of the box)
            if (particles->phase[k] != -1) {
                
                // Angle
                angle = -model.dt*particles->om_p[k];
                
                nx = particles->nx[k] * cos(angle) - particles->nz[k] * sin(angle) ;
                nz = particles->nx[k] * sin(angle) + particles->nz[k] * cos(angle) ;
                norm             = sqrt( nx*nx + nz*nz );
                particles->nx[k] = nx/norm;
                particles->nz[k] = nz/norm;
                
            }
        }
    
}

    /*--------------------------------------------------------------------------------------------------------------------*/
    /*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
    /*--------------------------------------------------------------------------------------------------------------------*/