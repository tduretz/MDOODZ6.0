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
double UpdateReactionProgress ( double ttrans0, double Pold, double Pc, double treac, double Preac, double dt ,scale *scaling) {
    
    double ttrans= 0.0;
    
    //printf("Update Reaction Progress\n");
    //printf("Pc = %2.2e; Preac = %2.2e; Pold = %2.2e\n", Pc*scaling->S, Preac*scaling->S, Pold*scaling->S);
    
    double dtr = 0.0;
    
    if ((Pc >= Preac) && (Pold >= Preac)) dtr = dt;
    if ((Pc <  Preac) && (Pold <  Preac)) dtr = 0.0;
    if ((Pc >= Preac) && (Pold <  Preac)) dtr = dt*(Pc-Preac)/(Pc-Pold);
    if ((Pc <  Preac) && (Pold >= Preac)) dtr = dt*(Preac-Pc)/(Pold-Pc);
    
    //printf("dt = %2.2e; dtr = %2.2e\n", dt*scaling->t, dtr*scaling->t);
    ttrans = ttrans0 + dtr;
    
    //printf("ttrans = %2.2e; ttrans0 = %2.2e\n", ttrans*scaling->t, ttrans0*scaling->t);
    
    
    
    return ttrans;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void RheologicalOperators( grid* mesh, params* model, scale* scaling, int Jacobian ) {
    
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    double nx, nz, deta=5e21/scaling->eta;
    
    if (Jacobian == 0  && model->aniso == 0) {
        
        // Loop on cell centers
#pragma omp parallel for shared( mesh )
        for (k=0; k<Ncx*Ncz; k++) {
            
            if ( mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) {
                mesh->D11_n[k] = 2*mesh->eta_n[k];
                mesh->D22_n[k] = 2*mesh->eta_n[k];
            }
            else {
                mesh->D11_n[k] = 0.0;
                mesh->D22_n[k] = 0.0;
            }
        }
        
        // Loop on cell vertices
#pragma omp parallel for shared( mesh )
        for (k=0; k<Nx*Nz; k++) {
            
            if ( mesh->BCg.type[k] != 30 ) {
                mesh->D33_s[k] = mesh->eta_s[k];
            }
            else {
                mesh->D33_s[k] = 0.0;
            }
        }
        
    }
    
    if (Jacobian == 1  && model->aniso == 0) {
        
        double etae;
        double yn = 1.0;
        int el = model->iselastic;
        // Loop on cell centers
#pragma omp parallel for shared( mesh ) private ( etae ) firstprivate ( el, yn )
        for (k=0; k<Ncx*Ncz; k++) {
            
            if ( mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) {
                switch ( el ) {
                    case 0:
                        etae = model->dt*mesh->mu_n[k];
                        mesh->D11_n[k] = 2.0*mesh->eta_n[k] + yn*2.0*mesh->detadexx_n[k]*mesh->exxd[k];
                        mesh->D12_n[k] =                      yn*2.0*mesh->detadezz_n[k]*mesh->exxd[k];
                        mesh->D13_n[k] =                      yn*2.0*mesh->detadgxz_n[k]*mesh->exxd[k];
                        mesh->D14_n[k] =                      yn*2.0*mesh->detadp_n[k]  *mesh->exxd[k];
                        
                        mesh->D21_n[k] =                      yn*2.0*mesh->detadexx_n[k]*mesh->ezzd[k];
                        mesh->D22_n[k] = 2.0*mesh->eta_n[k] + yn*2.0*mesh->detadezz_n[k]*mesh->ezzd[k];
                        mesh->D23_n[k] =                      yn*2.0*mesh->detadgxz_n[k]*mesh->ezzd[k];
                        mesh->D24_n[k] =                      yn*2.0*mesh->detadp_n[k]  *mesh->ezzd[k];
                        
                        break;
                    case 1:
                        etae = model->dt*mesh->mu_n[k];
                        mesh->D11_n[k] = 2.0*mesh->eta_n[k] + yn*2.0*mesh->detadexx_n[k]*mesh->exxd[k] + mesh->sxxd0[k]*mesh->detadexx_n[k] / etae;
                        mesh->D12_n[k] =                      yn*2.0*mesh->detadezz_n[k]*mesh->exxd[k] + mesh->sxxd0[k]*mesh->detadezz_n[k] / etae;
                        mesh->D13_n[k] =                      yn*2.0*mesh->detadgxz_n[k]*mesh->exxd[k] + mesh->sxxd0[k]*mesh->detadgxz_n[k] / etae;
                        mesh->D14_n[k] =                      yn*2.0*mesh->detadp_n[k]  *mesh->exxd[k] + mesh->sxxd0[k]*mesh->detadp_n[k]   / etae;
                        
                        mesh->D21_n[k] =                      yn*2.0*mesh->detadexx_n[k]*mesh->ezzd[k] + mesh->szzd0[k]*mesh->detadexx_n[k] / etae;
                        mesh->D22_n[k] = 2.0*mesh->eta_n[k] + yn*2.0*mesh->detadezz_n[k]*mesh->ezzd[k] + mesh->szzd0[k]*mesh->detadezz_n[k] / etae;
                        mesh->D23_n[k] =                      yn*2.0*mesh->detadgxz_n[k]*mesh->ezzd[k] + mesh->szzd0[k]*mesh->detadgxz_n[k] / etae;
                        mesh->D24_n[k] =                      yn*2.0*mesh->detadp_n[k]  *mesh->ezzd[k] + mesh->szzd0[k]*mesh->detadp_n[k]   / etae;
                        break;
                }
            }
            else {
                mesh->D11_n[k] = 0.0;
                mesh->D12_n[k] = 0.0;
                mesh->D13_n[k] = 0.0;
                mesh->D14_n[k] = 0.0;
                
                mesh->D21_n[k] = 0.0;
                mesh->D22_n[k] = 0.0;
                mesh->D23_n[k] = 0.0;
                mesh->D24_n[k] = 0.0;
            }
        }
        
        // Loop on cell vertices
#pragma omp parallel for shared( mesh )  private ( etae ) firstprivate ( el )
        for (k=0; k<Nx*Nz; k++) {
            
            if ( mesh->BCg.type[k] != 30 ) {
                switch ( el ) {
                    case 0:
                        mesh->D31_s[k] =                  yn*mesh->detadexx_s[k]*2.0*mesh->exz[k];  // Factor 2 is important!!
                        mesh->D32_s[k] =                  yn*mesh->detadezz_s[k]*2.0*mesh->exz[k];
                        mesh->D33_s[k] = mesh->eta_s[k] + yn*mesh->detadgxz_s[k]*2.0*mesh->exz[k];
                        mesh->D34_s[k] =                  yn*mesh->detadp_s[k]  *2.0*mesh->exz[k];
                        break;
                    case 1:
                        etae = model->dt*mesh->mu_s[k];
                        mesh->D31_s[k] =                  yn*mesh->detadexx_s[k]*2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadexx_s[k] / etae;  // Factor 2 is important!!
                        mesh->D32_s[k] =                  yn*mesh->detadezz_s[k]*2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadezz_s[k] / etae;
                        mesh->D33_s[k] = mesh->eta_s[k] + yn*mesh->detadgxz_s[k]*2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadgxz_s[k] / etae;
                        mesh->D34_s[k] =                  yn*mesh->detadp_s[k]  *2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadp_s[k]   / etae;
                        break;
                }
                
            }
            else {
                mesh->D31_s[k] = 0.0;
                mesh->D32_s[k] = 0.0;
                mesh->D33_s[k] = 0.0;
                mesh->D34_s[k] = 0.0;
            }
            if (isnan(mesh->D34_s[k])) exit(1);
            if (isinf(mesh->D34_s[k])) exit(1);
            
        }
        
    }
    
    if ( Jacobian == 0 && model->aniso == 1 ) {
        
        printf("Computing anisotropic viscosity tensor\n");
        
        // Loop on cell centers
#pragma omp parallel for shared( mesh ) private ( nx, nz ) firstprivate (deta)
        for (k=0; k<Ncx*Ncz; k++) {
            
            // Director
            nx = mesh->nx_n[k];
            nz = mesh->nz_n[k];
            
            if ( mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) {
                mesh->D11_n[k] = 2.0*mesh->eta_n[k] - 4.0*deta*pow(nx, 2.0)*pow(nz, 2.0);
                mesh->D12_n[k] =                      4.0*deta*pow(nx, 2.0)*pow(nz, 2.0);
                mesh->D13_n[k] =                      1.0*deta*(pow(nx, 3.0)*nz - pow(nx, 4.0));
                mesh->D14_n[k] =                      0.0;
                
                mesh->D21_n[k] =                      4.0*deta*pow(nx, 2.0)*pow(nz, 2.0);
                mesh->D22_n[k] = 2.0*mesh->eta_n[k] - 4.0*deta*pow(nx, 2.0)*pow(nz, 2.0);
                mesh->D23_n[k] =                      1.0*deta*(-pow(nx, 3.0)*nz + pow(nx, 4.0));
                mesh->D24_n[k] =                      0.0;
            }
            else {
                mesh->D11_n[k] = 0.0;
                mesh->D12_n[k] = 0.0;
                mesh->D13_n[k] = 0.0;
                mesh->D14_n[k] = 0.0;
                
                mesh->D21_n[k] = 0.0;
                mesh->D22_n[k] = 0.0;
                mesh->D23_n[k] = 0.0;
                mesh->D24_n[k] = 0.0;
            }
        }
        
        // Loop on cell vertices
#pragma omp parallel for shared( mesh )  private ( nx, nz ) firstprivate (deta)
        for (k=0; k<Nx*Nz; k++) {
            
            // Director
            nx = mesh->nx_s[k];
            nz = mesh->nz_s[k];
            
            if ( mesh->BCg.type[k] != 30 ) {
                mesh->D31_s[k] =                  2.0*deta*(pow(nx, 3.0)*nz - pow(nx, 4.0));
                mesh->D32_s[k] =                  2.0*deta*(-pow(nx, 3.0)*nz + pow(nx, 4.0));
                mesh->D33_s[k] = mesh->eta_s[k] + deta*(2.0*pow(nx, 2.0)*pow(nz, 2.0) - 0.5);
                mesh->D34_s[k] =                  0.0;
            }
            else {
                mesh->D31_s[k] = 0.0;
                mesh->D32_s[k] = 0.0;
                mesh->D33_s[k] = 0.0;
                mesh->D34_s[k] = 0.0;
            }
        }
        
    }
    
    
    //    double D31_per, D32_per, D33_per, D34_per;
    //    for ( int l=0; l<Nz; l++ ) {
    //
    //        int c0 = 0 + l*(Nx);
    //        int c1 = (Nx-1) + l*(Nx);
    //
    //        D31_per     = 0.5*( mesh->D31_s[c0] + mesh->D31_s[c1] );
    //        D32_per     = 0.5*( mesh->D32_s[c0] + mesh->D32_s[c1] );
    //        D33_per     = 0.5*( mesh->D33_s[c0] + mesh->D33_s[c1] );
    //        D34_per     = 0.5*( mesh->D34_s[c0] + mesh->D34_s[c1] );
    //
    //        if ( mesh->BCg.type[c1] != 30 ) {
    //
    //            mesh->D31_s[c0] = D31_per;
    //            mesh->D31_s[c1] = D31_per;
    //            mesh->D32_s[c0] = D32_per;
    //            mesh->D32_s[c1] = D32_per;
    //            mesh->D33_s[c0] = D33_per;
    //            mesh->D33_s[c1] = D33_per;
    //            mesh->D34_s[c0] = D34_per;
    //            mesh->D34_s[c1] = D34_per;
    //
    //        }
    //    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AccumulatedStrainII( grid* mesh, scale scaling, params model, markers* particles, double* X_vect, double* Z_vect, int Nx, int Nz, char *tag ) {
    
    double *strain_inc;
    int k, l, c1;
    DoodzFP *strain_inc_el, *strain_inc_pl, *strain_inc_pwl, *strain_inc_exp, *strain_inc_lin, *strain_inc_gbs;
    
    //    printf("Accumulating strain\n");
    
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
   
    // Free
    DoodzFree( dudx );
    DoodzFree( dvdz );
    DoodzFree( dudz );
    DoodzFree( dvdx );
    DoodzFree( pdudx );
    DoodzFree( pdudz );
    DoodzFree( pdvdx );
    DoodzFree( pdvdz );
 
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
    strain_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    
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

//void RotateStresses( grid mesh, markers* particles, params model, scale *scaling ) {
//    
//    int k, l, cp, cu, cv, Nx, Nz,type=0;
//    double dx, dz;
//    double angle, txx, tzz, txz;
//    double sxxr, sxzr;
//    
//    Nx = mesh.Nx;
//    Nz = mesh.Nz;
//    dx = mesh.dx;
//    dz = mesh.dz;
//    
//    if (type==0) {
//        // JAUMANN RATE
//        
//        // Stress correction
//#pragma omp parallel for shared ( particles ) private ( angle, k, txx, tzz, txz ) firstprivate( model ) schedule( static )
//        for(k=0; k<particles->Nb_part; k++) {
//            
//            // Filter out particles that are inactive (out of the box)
//            if (particles->phase[k] != -1) {
//                
//                // Angle
//                angle = model.dt*particles->om_p[k];
//                
//                txx = particles->sxxd[k];
//                tzz = particles->szzd[k];
//                txz = particles->sxz[k];
//                
//                particles->sxxd[k] = (txx*cos(angle) - txz*sin(angle))*cos(angle) - (txz*cos(angle) - tzz*sin(angle))*sin(angle);
//                particles->szzd[k] = (txx*sin(angle) + txz*cos(angle))*sin(angle) + (txz*sin(angle) + tzz*cos(angle))*cos(angle);
//                particles->sxz[k]  = (txx*cos(angle) - txz*sin(angle))*sin(angle) + (txz*cos(angle) - tzz*sin(angle))*cos(angle);
//                
////                                // Re-belotte:
////                                double sxxr =   ( cos(angle)*particles->sxxd[k] - sin(angle)*particles->sxz[k] ) * cos(angle) + ( particles->sxz[k]*cos(angle) +  sin(angle)*particles->sxxd[k])*sin(angle) ;
////                                double sxzr = - ( cos(angle)*particles->sxxd[k] - sin(angle)*particles->sxz[k] ) * sin(angle) + ( particles->sxz[k]*cos(angle) +  sin(angle)*particles->sxxd[k])*cos(angle) ;
////
////                                particles->sxxd[k] = sxxr;
////                                particles->szzd[k] =-sxxr;
////                                particles->sxz[k]  = sxzr;
//            }
//        }
//    }
//        else {
//    
//            // UPPER CONVECTED
//    
//            double *dudx, *dudz, *dvdx, *dvdz;
//            DoodzFP *pdudx, *pdudz, *pdvdx, *pdvdz, *VEm;
//    
//            dudx   = DoodzMalloc ((Nx-1)*(Nz-1)*sizeof(double));
//            dvdz   = DoodzMalloc ((Nx-1)*(Nz-1)*sizeof(double));
//            dudz   = DoodzMalloc ((Nx)*(Nz)*sizeof(double));
//            dvdx   = DoodzMalloc ((Nx)*(Nz)*sizeof(double));
//            pdudx  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
//            pdudz  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
//            pdvdx  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
//            pdvdz  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
//            VEm    = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
//    
//            // Compute dudx and dvdz (cell centers)
//            for (k=0; k<Nx-1; k++) {
//                for (l=0; l<Nz-1; l++) {
//                    cp = k  + l*(Nx-1);
//                    cu = k  + l*(Nx) + Nx;
//                    cv = k  + l*(Nx+1) + 1;
//                    dudx[cp] = 1.0/dx * ( mesh.u_in[cu+1]      - mesh.u_in[cu] );
//                    dvdz[cp] = 1.0/dz * ( mesh.v_in[cv+(Nx+1)] - mesh.v_in[cv] );
//                }
//            }
//    
//            // Compute dudx and dvdz (cell vertices)
//            for (k=0; k<Nx; k++) {
//                for (l=0; l<Nz; l++) {
//                    cp = k  + l*(Nx-1);
//                    cu = k  + l*(Nx);
//                    cv = k  + l*(Nx+1);
//                    dudz[cu] = 1.0/dz * ( mesh.u_in[cu+Nx] - mesh.u_in[cu] );
//                    dvdx[cu] = 1.0/dx * ( mesh.v_in[cv+1]  - mesh.v_in[cv] );
//                }
//            }
//    
//            // Interpolate from grid to particles
//            Interp_Grid2P( *(particles), pdudx, &mesh, dudx, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type  );
//            Interp_Grid2P( *(particles), pdvdz, &mesh, dvdz, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
//            Interp_Grid2P( *(particles), pdvdx, &mesh, dvdx, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
//            Interp_Grid2P( *(particles), pdudz, &mesh, dudz, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
//            Interp_Grid2P( *(particles), VEm, &mesh, mesh.VE_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type  );
//    
//    
//            // Stress correction
//    #pragma omp parallel for shared ( particles, model, VEm, pdudx, pdudz, pdvdx, pdvdz, sxxr, sxzr ) private ( k ) schedule( static )
//            for(k=0; k<particles->Nb_part; k++) {
//                sxxr = model.dt * VEm[k] * ( -2.0*particles->sxxd[k]*pdudx[k] - 2.0*particles->sxz[k]*pdudz[k]);
//                sxzr = model.dt * VEm[k] * (      particles->sxxd[k]*pdudz[k] -     particles->sxxd[k]*pdvdx[k] - particles->sxz[k]*(pdudx[k] + pdvdz[k]) );
////
////                sxxr = model.dt * VEm[k] * ( -2.0*particles->sxxd[k]*pdudx[k] - 2.0*particles->sxz[k]*pdudz[k]);
////                sxzr = model.dt * VEm[k] * (      particles->sxxd[k]*pdudz[k] -     particles->sxxd[k]*pdvdx[k] - particles->sxz[k]*(pdudx[k] + pdvdz[k]) );
//    
//                particles->sxxd[k] -= sxxr;
//                particles->szzd[k]  =-sxxr;
//                particles->sxz[k]  -= sxzr;
//    
////                // Only if explicity
////                if (model.subgrid_diff==4) particles->sxxd[k] *= VEm[k];
////                if (model.subgrid_diff==4) particles->sxz[k]  *= VEm[k];
//            }
//            MinMaxArray( pdudx,  scaling->E, particles->Nb_part, "dudx p." );
//            MinMaxArray( pdvdx,  scaling->E, particles->Nb_part, "dvdx p." );
//            MinMaxArray( pdudz,  scaling->E, particles->Nb_part, "dudz p." );
//    
//    
//            // clean memory
//            DoodzFree(dudx);
//            DoodzFree(dudz);
//            DoodzFree(dvdx);
//            DoodzFree(dvdz);
//            DoodzFree(pdudx);
//            DoodzFree(pdudz);
//            DoodzFree(pdvdx);
//            DoodzFree(pdvdz);
//            DoodzFree(VEm);
//        }
//}


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

void UpdateParticleDensity( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
    DoodzFP *rho_inc_mark, *rho_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    rho_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    rho_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    
    for (k=0;k<Ncx*Ncz;k++) {
        rho_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) rho_inc_grid[k] = mesh->rho_n[k] - mesh->rho0_n[k];
    }
    
    // Interp increments to particles
    Interp_Grid2P( *particles, rho_inc_mark, mesh, rho_inc_grid, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCt.type  );
    
    // Increment temperature on particles
    ArrayPlusArray( particles->rho, rho_inc_mark, particles->Nb_part );
    
    DoodzFree(rho_inc_grid);
    DoodzFree(rho_inc_mark);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleGrainSize( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
    DoodzFP *d_inc_mark, *d_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    d_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    d_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    
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
        Tg0  = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dTgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dTgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        Tm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        dTms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        dTmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        
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
        
        T_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        
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
    
    P_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    
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
    DoodzFP *mdsxxd, *mdszzd, *mdsxz, *dsxxd, *dszzd, *dsxz, d=1.0, dtaum;
    DoodzFP *dtxxgs, *dtzzgs, *dtxzgs, *dtxxgr, *dtzzgr, *dtxzgr, *txxm0, *tzzm0, *txzm0, *dtxxms, *dtzzms, *dtxzms, *dtxxmr, *dtzzmr, *dtxzmr, *etam;
    
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    if ( model->subgrid_diff > -1 ) {
        
        // Alloc
        dtxxgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dtzzgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dtxzgs = DoodzCalloc(Nx*Nz,   sizeof(DoodzFP));
        dtxxgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dtzzgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dtxzgr = DoodzCalloc(Nx*Nz,   sizeof(DoodzFP));
        txxm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        tzzm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        txzm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        dtxxms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        dtzzms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        dtxzms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        dtxxmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        dtzzmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        dtxzmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        etam   = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        
        // Old stresses grid --> markers
        Interp_Grid2P( *particles, txxm0, mesh, mesh->sxxd0, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
        Interp_Grid2P( *particles, tzzm0, mesh, mesh->szzd0, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
        Interp_Grid2P( *particles, txzm0, mesh, mesh->sxz0 , mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type     );
        Interp_Grid2P( *particles, etam,  mesh, mesh->eta_phys_s, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type);      
 
        if ( model->subgrid_diff == 2 ) {
            
            printf("Subgrid diffusion for stress tensor component update\n");
            
            // Compute subgrid stress increments on markers
#pragma omp parallel for shared(particles,txxm0,tzzm0,txzm0,dtxxms,dtzzms,dtxzms,etam) private(k,p,dtaum) firstprivate(materials,model,d)
            for ( k=0; k<particles->Nb_part; k++ ) {
                if (particles->phase[k] != -1) {
                    p         = particles->phase[k];
                    dtaum     = etam[k]/ (materials->mu[p]);
                    dtxxms[k] = -( particles->sxxd[k] - txxm0[k]) * (1.0 - exp(-d*model->dt/dtaum));
                    dtzzms[k] = -( particles->szzd[k] - tzzm0[k]) * (1.0 - exp(-d*model->dt/dtaum));
                    dtxzms[k] = -( particles->sxz[k]  - txzm0[k]) * (1.0 - exp(-d*model->dt/dtaum));
                }
            }
            
            // Subgrid stress increments markers --> grid
            Interp_P2C ( *particles, dtxxms, mesh, dtxxgs, mesh->xg_coord, mesh->zg_coord, 1, 0 );
            Interp_P2C ( *particles, dtzzms, mesh, dtzzgs, mesh->xg_coord, mesh->zg_coord, 1, 0 );
            Interp_P2N ( *particles, dtxzms, mesh, dtxzgs, mesh->xg_coord, mesh->zg_coord, 1, 0, model );
            
            // Remaining stress increments on the grid
#pragma omp parallel for shared(mesh,dtxxgs,dtxxgr,dtzzgs,dtzzgr) private(c0) firstprivate(Ncx,Ncz)
            for ( c0=0; c0<Ncx*Ncz; c0++ ) {
                if (mesh->BCp.type[c0]!=30) dtxxgr[c0] = (mesh->sxxd[c0]-mesh->sxxd0[c0]) - dtxxgs[c0];
                if (mesh->BCp.type[c0]!=30) dtzzgr[c0] = (mesh->szzd[c0]-mesh->szzd0[c0]) - dtzzgs[c0];
            }
#pragma omp parallel for shared(mesh,dtxzgs,dtxzgr) private(c0) firstprivate(Nx,Nz)
            for ( c0=0; c0<Nx*Nz; c0++ ) {
                if (mesh->BCg.type[c0]   !=30) dtxzgr[c0] = (mesh->sxz[c0]-mesh->sxz0[c0]) - dtxzgs[c0];
            }
            
            // Remaining stress increments grid --> markers
            Interp_Grid2P( *particles, dtxxmr, mesh, dtxxgr, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type  );
            Interp_Grid2P( *particles, dtzzmr, mesh, dtzzgr, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type  );
            Interp_Grid2P( *particles, dtxzmr, mesh, dtxzgr, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type  );
            
            // Final stresses update on markers
#pragma omp parallel for shared(particles,dtxxms,dtzzms,dtxzms,dtxxmr,dtzzmr,dtxzmr) private(k)
            for ( k=0; k<particles->Nb_part; k++ ) {
                if (particles->phase[k] != -1) particles->sxxd[k] += dtxxms[k] + dtxxmr[k];
                if (particles->phase[k] != -1) particles->szzd[k] += dtzzms[k] + dtzzmr[k];
                if (particles->phase[k] != -1) particles->sxz[k]  += dtxzms[k] + dtxzmr[k];
            }
        }
        
        // Free
        DoodzFree( dtxxgs );
        DoodzFree( dtzzgs );
        DoodzFree( dtxzgs );
        DoodzFree( dtxxgr );
        DoodzFree( dtzzgr );
        DoodzFree( dtxzgr );
        DoodzFree( txxm0  );
        DoodzFree( tzzm0  );
        DoodzFree( txzm0  );
        DoodzFree( dtxxms );
        DoodzFree( dtzzms );
        DoodzFree( dtxzms );
        DoodzFree( dtxxmr );
        DoodzFree( dtzzmr );
        DoodzFree( dtxzmr );
        DoodzFree( etam   );
        
    }
    if (model->subgrid_diff==0 || model->subgrid_diff==1 || model->subgrid_diff==4){
        
        printf("No subgrid diffusion for stress tensor component update\n");
        
        // Alloc
        dsxxd  = DoodzCalloc((Nx-1)*(Nz-1),sizeof(double));
        dszzd  = DoodzCalloc((Nx-1)*(Nz-1),sizeof(double));
        dsxz   = DoodzCalloc((Nx)*(Nz),sizeof(double));
        mdsxxd = DoodzCalloc(particles->Nb_part,sizeof(DoodzFP));
        mdszzd = DoodzCalloc(particles->Nb_part,sizeof(DoodzFP));
        mdsxz  = DoodzCalloc(particles->Nb_part,sizeof(DoodzFP));
        
        // Cell: normal stress change
        for (k=0; k<Nx-1; k++) {
            for (l=0; l<Nz-1; l++) {
                c0 = k  + l*(Nx-1);
                if (mesh->BCp.type[c0] !=30 && mesh->BCp.type[c0] !=31) dsxxd[c0] = mesh->sxxd[c0] - mesh->sxxd0[c0];
                if (mesh->BCp.type[c0] !=30 && mesh->BCp.type[c0] !=31) dszzd[c0] = mesh->szzd[c0] - mesh->szzd0[c0];
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
        Interp_Grid2P( *particles, mdszzd, mesh, dszzd, mesh->xc_coord,  mesh->zc_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type  );
        Interp_Grid2P( *particles, mdsxz,  mesh, dsxz,  mesh->xg_coord,  mesh->zg_coord,  mesh->Nx,   mesh->Nz,   mesh->BCg.type  );
        
        // Update marker stresses
        ArrayPlusArray( particles->sxxd, mdsxxd, particles->Nb_part );
        ArrayPlusArray( particles->szzd, mdszzd, particles->Nb_part );
        ArrayPlusArray( particles->sxz,  mdsxz,  particles->Nb_part );
        
        // Free
        DoodzFree(dsxxd);
        DoodzFree(dszzd);
        DoodzFree(dsxz);
        DoodzFree(mdsxxd);
        DoodzFree(mdszzd);
        DoodzFree(mdsxz);
    }
}

double Viscosity( int phase, double G, double T, double P, double d, double phi, double exx, double ezz, double exz, double Txx0, double Tzz0, double Txz0, mat_prop* materials, params *model, scale *scaling, double *txxn, double *tzzn, double *txzn, double* etaVE, double* VEcoeff, double* Eii_el, double* Eii_pl, double* Eii_pwl, double* Eii_exp , double* Eii_lin, double* Eii_gbs, double* Eii_cst, double* Exx_pwl, double* Exz_pwl, double* Exx_el, double* Ezz_el, double* Exz_el, double* Exx_diss, double* Ezz_diss, double* Exz_diss, double* Exx_pl, double* Exz_pl, double *d1, double strain_acc, double Phi, double C, double *detadexx, double *detadezz, double *detadexz, double *detadp, double P0, double ttrans0, double *Xreac, double *ttrans  ) {

    // General paramaters
    double eta=0.0, R=materials->R, dt=model->dt;
    double minEta=model->mineta, maxEta=model->maxeta;
    double TmaxPeierls = (1200.0+zeroC)/scaling->T;      // max. T for Peierls

    // Parameters for deformation map calculations
    int    local_iter = model->loc_iter, it, nitmax = 100, noisy=0;
    int    constant=0, dislocation=0, peierls=0, diffusion=0, gbs=0, elastic = model->iselastic;
    double tol = 1.0e-13, res=0.0, res0=0.0, dfdeta=0.0, Txx=0.0, Tzz=0.0, Txz=0.0, Tii=0.0, ieta_sum=0.0, Tii0 = sqrt(Txx0*Txx0 + Txz0*Txz0);
    double eta_up=0.0, eta_lo=0.0, eta_ve=0.0, eta_p=0.0, r_eta_pl=0.0, r_eta_ve=0.0, r_eta_p=0.0;
    double eta_pwl=0.0, eta_exp=0.0, eta_vep=0.0, eta_lin=0.0, eta_el=0.0, eta_gbs=0.0, eta_cst=0.0, eta_step=0.0;
    double Exx=0.0, Ezz=0.0, Exz=0.0, Eii_vis=0.0, Eii= 0.0, eII=0.0;

    // Flow law parameters from input file
    double Tyield=0.0, F_trial = 0.0, F_corr = 0.0, gdot = 0.0, dQdtxx = 0.0, dQdtzz= 0.0, dQdtxz= 0.0;
    int    is_pl = 0;
    double Ea_pwl = materials->Qpwl[phase], Va_pwl = materials->Vpwl[phase], n_pwl  = materials->npwl[phase], m_pwl  = materials->mpwl[phase], r_pwl  = materials->rpwl[phase], A_pwl  = materials->Apwl[phase], f_pwl = materials->fpwl[phase], a_pwl = materials->apwl[phase], F_pwl  = materials->Fpwl[phase], pre_factor = materials->pref_pwl[phase];
    double Ea_lin = materials->Qlin[phase], Va_lin = materials->Vlin[phase], n_lin  = materials->nlin[phase], m_lin  = materials->mlin[phase], r_lin  = materials->rlin[phase], A_lin  = materials->Alin[phase], f_lin = materials->flin[phase], a_lin = materials->alin[phase], F_lin  = materials->Flin[phase];
    double Ea_gbs = materials->Qgbs[phase], Va_gbs = materials->Vgbs[phase], n_gbs  = materials->ngbs[phase], m_gbs  = materials->mgbs[phase], r_gbs  = materials->rgbs[phase], A_gbs  = materials->Agbs[phase], f_gbs = materials->fgbs[phase], a_gbs = materials->agbs[phase], F_gbs  = materials->Fgbs[phase];
    double B_pwl=0.0, B_lin=0.0, B_exp=0.0, B_gbs=0.0, C_pwl=0.0,  C_lin=0.0, C_exp=0.0, C_gbs=0.0;
    double Ea_exp = materials->Qexp[phase], S_exp  = materials->Sexp[phase], E_exp  = materials->Eexp[phase] , t_exp  = materials->texp[phase], F_exp=0.0;
    double gamma=materials->Gexp[phase], ST, q=materials->qexp[phase],  n_exp=materials->nexp[phase];
    double Exx_lin=0.0, Ezz_lin=0.0, Exz_lin=0.0, Exx_exp=0.0, Ezz_exp=0.0, Exz_exp=0.0, Exx_gbs=0.0, Ezz_gbs=0.0, Exz_gbs=0.0, Exx_cst=0.0, Ezz_cst=0.0, Exz_cst=0.0;
    double Ezz_pl=0.0, Ezz_pwl=0.0;
    int gs = materials->gs[phase];
    double pg = materials->ppzm[phase], Kg = materials->Kpzm[phase], Qg = materials->Qpzm[phase], gam = materials->Gpzm[phase], cg = materials->cpzm[phase], lambda = materials->Lpzm[phase];
    double eta_vp = materials->eta_vp[phase];    
    double  detadTxx=0.0, detadTzz=0.0, detadTxz=0.0, deta_ve_dExx=0.0, deta_ve_dEzz=0.0, deta_ve_dExz=0.0, deta_ve_dP=0.0;
    double  dFdExx=0.0, dFdEzz=0.0, dFdExz=0.0, dFdP=0.0, g=0.0, dlamdExx=0.0, dlamdEzz=0.0, dlamdExz=0.0, dlamdP=0.0, a=0.0, deta_vep_dExx=0.0, deta_vep_dEzz=0.0, deta_vep_dExz=0.0, deta_vep_dP=0.0;

    //------------------------------------------------------------------------//

    // Initialise strain rate invariants to 0
    *Eii_exp = 0.0; *Eii_lin = 0.0; *Eii_pl = 0.0; *Eii_pwl = 0.0; *Eii_el = 0.0, *Eii_gbs=0, *Eii_cst=0.0;
    *txxn=0.0; *txzn=0.0; *etaVE=0.0; *VEcoeff=0.0; *Exx_pwl=0.0; *Exz_pwl=0.0, *Exx_el=0.0, *Ezz_el=0.0, *Exz_el=0.0, *Exx_diss=0.0, *Ezz_diss=0.0, *Exz_diss=0.0, *Exx_pl=0.0, *Exz_pl=0.0, *d1=0.0;
    *detadexx=0.0; *detadezz=0.0; *detadexz=0.0; *detadp=0.0;
    *Xreac=0.0; *ttrans=0.0;

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
    if ( model->gz>0.0 && P<0.0     ) { P = 0.0; printf("Aie aie aie P < 0 !!!\n"); exit(122);}

    // Visco-plastic limit
    if ( elastic==0                 ) G = 10.0;

    // Zero C limit
    if ( T< zeroC/scaling->T        ) T = zeroC/scaling->T;

    //------------------------------------------------------------------------//

    // Precomputations
    if ( dislocation == 1 ) {
        B_pwl = pre_factor * F_pwl * pow(A_pwl,-1.0/n_pwl) * exp( (Ea_pwl + P*Va_pwl)/R/n_pwl/T ) * pow(d, m_pwl/n_pwl) * pow(f_pwl, -r_pwl/n_pwl) * exp(-a_pwl*phi/n_pwl);
        C_pwl = pow(2.0*B_pwl, -n_pwl);
    }
    if ( diffusion == 1 ) {
        if (m_lin>0.0 && d<1e-13/scaling->L){
            printf("Cannot run with grain size dependent viscosity if grain size is set to 0 --> d = %2.2e!!!\n", d*scaling->L);
            exit(1);
        };
        B_lin = F_lin * pow(A_lin,-1.0/n_lin) * exp( (Ea_lin + P*Va_lin)/R/n_lin/T ) * pow(f_lin, -r_lin/n_lin) * exp(-a_lin*phi/n_lin); // * pow(d, m_lin/n_lin) !!!!!!!!!!!!!!!!!!!!!!!!
        C_lin = pow(2.0*B_lin, -n_lin);
    }
    if ( gbs == 1 ) {
        B_gbs = F_gbs * pow(A_gbs,-1.0/n_gbs) * exp( (Ea_gbs + P*Va_gbs)/R/n_gbs/T ) * pow(d, m_gbs/n_gbs) * pow(f_gbs, -r_gbs/n_gbs) * exp(-a_gbs*phi/n_gbs);
        C_gbs = pow(2.0*B_gbs, -n_gbs);
    }
    if ( peierls   == 1 ) {
        ST                    = Ea_exp/R/T * pow((1.0-gamma),(q-1.0)) * q*gamma;
        if ( t_exp == 0) F_exp  = 1.0;
        if ( t_exp == 1) F_exp  = 1.0/6.0*pow(2.0,1.0/(ST+n_exp)) * pow(3.0,(ST+n_exp-1.0)/2.0/(ST+n_exp));
        if ( t_exp == 2) F_exp  = 1.0/4.0*pow(2,1.0/(ST+n_exp));
        B_exp                   = F_exp * pow(E_exp*exp(-Ea_exp/R/T*pow(1.0-gamma,2.0)), -1.0/(ST+n_exp)) * pow(gamma*S_exp, ST/(ST+n_exp));
        C_exp                   = pow(2.0*B_exp, -(ST+n_exp));
    }

    Tyield                      = C*cos(Phi) +  ( P + model->PrBG)*sin(Phi);
    Tyield                      = MINV( Tyield, materials->Slim[phase] );

    //------------------------------------------------------------------------//
    // Reaction stuff:
    
        //printf("isreac = %d\n",materials->Reac[phase]);
        if ( materials->Reac[phase] > 0 && model->step > 0) {
        
        //printf("ttrans0 = %2.2e\n",ttrans0*scaling->t);
        //printf("dt      = %2.2e\n",dt*scaling->t);
        *ttrans = UpdateReactionProgress ( ttrans0, (P0+model->PrBG), (P+model->PrBG), materials->treac[phase], materials->Preac[phase], dt , scaling);
        //printf("*ttrans = %2.2e\n",*ttrans*scaling->t);
        
        //*Xreac  = 0.0; // f(ttrans)
        //*Xreac = 0.5*(1.0 + erf((*ttrans - materials->treac[phase]/2.0) / (1./6.* materials->treac[phase]) /sqrt(2)));
        
        double time_reaction, moy, sig;
        
        time_reaction = materials->treac[phase];
        moy  = time_reaction/2.0;
        sig = 1.0/6.0*time_reaction;
        *Xreac=0.5*(1.0 + erf((*ttrans-moy)/sig/sqrt(2.0)));
        
        double npwlreac  = materials->npwl[phase];
        double Apwlreac  = materials->Apwl[phase];
        double Eapwlreac = materials->Qpwl[phase];
        double Fpwlreac  = 1.0/6.0*pow(2.0,1.0/npwlreac) * pow(3.0,(npwlreac-1.0)/2.0/npwlreac); // axial compression
        
        
        double npwlprod  = materials->npwl[materials->Reac[phase]];
        double Apwlprod  = materials->Apwl[materials->Reac[phase]];
        double Eapwlprod = materials->Qpwl[materials->Reac[phase]];
        double Fpwlprod  = 1.0/6.0*pow(2.0,1.0/npwlprod) * pow(3.0,(npwlprod-1.0)/2.0/npwlprod); // axial compression
        
        
        //        npwl_mix   = npwlreac;
        //        Apwl_mix   = Apwlreac;
        //        Eapwl_mix  = Eapwlreac;

        // Huet et al 2014 ---------------
        double npwl_mix, Apwl_mix,Eapwl_mix, Fpwl_mix;
        
        double f1 = 1.0-*Xreac;
        double f2 = *Xreac;
        
        // (1) Calcul des ai
        
        double a1 = npwlprod + 1.0;
        double a2 = npwlreac + 1.0;
        
        // (2) Calcul de n bulk:
        
        double sum_up   = f1*a1*npwlreac + f2*a2*npwlprod;
        double sum_down = f1*a1+f2*a2;
        npwl_mix = sum_up/sum_down;
        
        // (3) Calcul de Q bulk:
        
        sum_up   = f1*a1*Eapwlreac + f2*a2*Eapwlprod;
        sum_down = f1*a1+f2*a2;
        Eapwl_mix = sum_up/sum_down;
        
        // (4) Calcul de A bulk:
        
        sum_down = f1*a1 + f2*a2;
        double Prod1 = pow(Apwlreac,(f1*a1/sum_down)) * pow(Apwlprod,(f2*a2/sum_down));
        double sum_n = f1*npwlreac/(npwlreac + 1.0) + f2*npwlprod/(npwlprod + 1.0);
        //double Prod2 = (npwlreac/(npwlreac + 1.0))^(f1*a1*npwlreac/sum_down) * (npwlprod/(npwlprod+1.0))^(f2*a2*npwlprod/sum_down);
        
        double Prod2 = pow(npwlreac/(npwlreac + 1.0),f1*a1*npwlreac/sum_down) * pow(npwlprod/(npwlprod+1.0),f2*a2*npwlprod/sum_down);
        
        Apwl_mix = Prod1 * pow(sum_n,-npwl_mix) * Prod2;
        
        // -------------------------------
        
        Fpwl_mix   = 1.0/6.0*pow(2.0,1.0/npwl_mix) * pow(3.0,(npwl_mix-1.0)/2.0/npwl_mix);
        
        if ( dislocation == 1 ) {
            B_pwl  = pre_factor * Fpwl_mix * pow(Apwl_mix,-1.0/npwl_mix) * exp( (Eapwl_mix)/R/npwl_mix/T );
            C_pwl  = pow(2.0*B_pwl, -npwl_mix);
            n_pwl  = npwl_mix;
            Ea_pwl = Apwl_mix;
            A_pwl  = Apwl_mix;
        }
        
        
    }
    
    //------------------------------------------------------------------------//
    
    // Isolated viscosities
    Exx  = exx;
    Ezz  = ezz;
    Exz  = exz;
    if ( elastic == 0 || local_iter == 0 || local_iter == 2 ) {
        Exx  = exx;
        Ezz  = ezz;
        Exz  = exz;
    }
    else {
        Exx  = exx + Txx0/(2.0*G*dt);
        Ezz  = ezz + Tzz0/(2.0*G*dt);
        Exz  = exz + Txz0/(2.0*G*dt);
    }
    Eii  = sqrt(1.0/2.0*Exx*Exx + 1.0/2.0*Ezz*Ezz + Exz*Exz);
    eII  = sqrt(1.0/2.0*exx*exx + 1.0/2.0*ezz*ezz + exz*exz);
    if (Eii*scaling->E<1e-30) Eii=1e-30/scaling->E;

    //------------------------------------------------------------------------//

    // Isolated viscosities
    eta_el                           = G*dt;
    if ( constant    == 1 ) eta_cst  = materials->eta0[phase];
    if ( dislocation == 1 ) eta_pwl  = B_pwl * pow( Eii, 1.0/n_pwl - 1.0 );
    if ( diffusion   == 1 ) eta_lin  = B_lin * pow( Eii, 1.0/n_lin - 1.0 ) * pow(d, m_lin/n_lin); // !!! gs - dependence !!!
    if ( gbs         == 1 ) eta_gbs  = B_gbs * pow( Eii, 1.0/n_gbs - 1.0 );
    if ( peierls     == 1 ) eta_exp  = B_exp * pow( Eii, 1.0/(ST+n_exp) - 1.0 );
    
    //------------------------------------------------------------------------//

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

    //------------------------------------------------------------------------//

    // Initial guess
    eta_ve                  = 0.5*(eta_up+eta_lo);
    
    // Local iterations
    for (it=0; it<nitmax; it++) {
        
        // Function evaluation at current effective viscosity
        Tii = 2.0 * eta_ve * Eii;
        if ( constant    == 1 ) *Eii_cst = Tii/2.0/eta_cst;
        if ( dislocation == 1 ) *Eii_pwl = C_pwl * pow(Tii, n_pwl    );
        if ( gbs         == 1 ) *Eii_gbs = C_gbs * pow(Tii, n_gbs    );
        if ( peierls     == 1 ) *Eii_exp = C_exp * pow(Tii, ST+n_exp ); // Peierls - power law
        if ( gs          == 1 ) *d1      = exp(log( Kg*exp(-Qg/R/T) *gam/(lambda*(1.0/cg)* Tii *(*Eii_pwl + *Eii_exp + *Eii_gbs + *Eii_pl)*pg))/(1.0+pg));
        if ( diffusion   == 1 ) *Eii_lin = C_lin * pow(Tii, n_lin) * pow(*d1,-m_lin); // !!! gs - dependence !!!
        Eii_vis                          = *Eii_pwl + *Eii_exp + *Eii_lin + *Eii_gbs + *Eii_cst;
        r_eta_ve                         = Eii - elastic*Tii/(2.0*eta_el) - Eii_vis;
        
        // Residual check
        res = fabs(r_eta_ve/Eii);
        if (it==0) res0 = res;
        if (noisy==1) printf("It. %02d, r = %2.2e\n", it, res);
        if (res < tol) break;
        
        // UNDER CONSTRUCTION
        
        // Analytical derivative of function
        dfdeta = -Eii/eta_el;
        if ( dislocation == 1 ) dfdeta += -(*Eii_pwl)*n_pwl/eta_ve;
        if ( constant    == 1 ) dfdeta += -Eii/eta_cst;
        
        // Update viscosity
        eta_ve -= r_eta_ve / dfdeta;
    }

    // Recalculate stress components
    Txx                  = 2.0*eta_ve*Exx;
    Tzz                  = 2.0*eta_ve*Ezz;
    Txz                  = 2.0*eta_ve*Exz;
    Tii                  = 2.0*eta_ve*Eii;

    // Partial derivatives VE
    detadTxx     = -    C_pwl * pow(Tii, n_pwl - 3.0) * Txx*(n_pwl - 1.0) * pow(eta_ve,2.0);
    detadTzz     = -    C_pwl * pow(Tii, n_pwl - 3.0) * Tzz*(n_pwl - 1.0) * pow(eta_ve,2.0);
    detadTxz     = -2.0*C_pwl * pow(Tii, n_pwl - 3.0) * Txz*(n_pwl - 1.0) * pow(eta_ve,2.0);
    deta_ve_dExx = detadTxx * 2.0*eta_ve / (1.0 - 2.0*(detadTxx*Exx + detadTzz*Ezz + detadTxz*Exz));
    deta_ve_dEzz = detadTzz * 2.0*eta_ve / (1.0 - 2.0*(detadTxx*Exx + detadTzz*Ezz + detadTxz*Exz));
    deta_ve_dExz = detadTxz * 2.0*eta_ve / (1.0 - 2.0*(detadTxx*Exx + detadTzz*Ezz + detadTxz*Exz));
    deta_ve_dP   = 0.0;

    // Check yield stress
    F_trial = Tii - Tyield;
//    if (F>0.0){ exit(1); }

    if (F_trial > 1e-17) {
        is_pl   = 1;
//        eta_vp  = 0.0;
        gdot    = F_trial / ( eta_ve + eta_vp );
        dQdtxx  = Txx/2.0/Tii;
        dQdtzz  = Tzz/2.0/Tii;
        dQdtxz  = Txz/1.0/Tii;
        Txx     = 2.0*eta_ve*(Exx - gdot*dQdtxx    );
        Tzz     = 2.0*eta_ve*(Ezz - gdot*dQdtzz    );
        Txz     = 2.0*eta_ve*(Exz - gdot*dQdtxz/2.0);
        Tii     = sqrt( 0.5*( pow(Txx,2.0) + pow(Tzz,2.0) ) + pow(Txz,2.0)  );
        Tyield += gdot*eta_vp;
        eta_vep = Tii / (2.0*Eii);
        F_corr  = Tii - Tyield;
        *Eii_pl = gdot/2.0;
        Tii     = 2.0*eta_vep*Eii;
        F_corr  = Tii - Tyield;
//        printf("%2.2e %2.2e %2.2e %2.2e %2.2e\n", Tii, Tyield, F_trial, F_corr, eta_vp*gdot*scaling->S);
    }

    //------------------------------------------------------------------------//

    // Partial derivatives VEP
    if (is_pl==1) {
        dFdExx        =     Exx*eta_ve/Eii + 2.0*Eii*deta_ve_dExx;
        dFdEzz        =     Ezz*eta_ve/Eii + 2.0*Eii*deta_ve_dEzz;
        dFdExz        = 2.0*Exz*eta_ve/Eii + 2.0*Eii*deta_ve_dExz;
        dFdP          = -sin(Phi);
        g             = 1.0 / ( eta_ve + eta_vp );
        dlamdExx      = g * (dFdExx - gdot *deta_ve_dExx);
        dlamdEzz      = g * (dFdEzz - gdot *deta_ve_dEzz);
        dlamdExz      = g * (dFdExz - gdot *deta_ve_dExz);
        dlamdP        = g * (dFdP);
        a             =  eta_vp;
        deta_vep_dExx = -0.5*Exx*Tyield/(2.0*pow(Eii,3.0)) + (a*dlamdExx)/(2.0*Eii);
        deta_vep_dEzz = -0.5*Ezz*Tyield/(2.0*pow(Eii,3.0)) + (a*dlamdEzz)/(2.0*Eii);
        deta_vep_dExz =     -Exz*Tyield/(2.0*pow(Eii,3.0)) + (a*dlamdExz)/(2.0*Eii);
        deta_vep_dP   =                           (sin(Phi) +  a*dlamdP )/(2.0*Eii);
    }

    double inv_eta_diss = 1.0/eta_pwl;
    if (is_pl == 0) {
        (*detadexx) = deta_ve_dExx;
        (*detadezz) = deta_ve_dEzz;
        (*detadexz) = deta_ve_dExz;
        (*detadp)   = deta_ve_dP;
        (*etaVE)    = eta_ve;
    }
    else {
        (*detadexx)   = deta_vep_dExx;
        (*detadezz)   = deta_vep_dEzz;
        (*detadexz)   = deta_vep_dExz;
        (*detadp)     = deta_vep_dP;
        (*etaVE)      = eta_vep;
        inv_eta_diss += 1.0/eta_vep;
    }

    /*----------------------------------------------------*/
    /*----------------------------------------------------*/
    /*----------------------------------------------------*/

    eta_pwl  = pow(2.0*C_pwl,-1.0) * pow(Tii, 1.0-n_pwl);
    *Exx_el =  (Txx-Txx0)/2.0/eta_el;
    *Ezz_el =  (Tzz-Tzz0)/2.0/eta_el;
    *Exz_el =  (Txz-Txz0)/2.0/eta_el;
    *Exx_pl = gdot*dQdtxx;
     Ezz_pl = gdot*dQdtzz;
    *Exz_pl = gdot*dQdtxz/2.0;
    *Exx_pwl = Txx/2.0/eta_pwl;
     Ezz_pwl = Tzz/2.0/eta_pwl;
    *Exz_pwl = Txz/2.0/eta_pwl;


    // Compute dissipative strain rate components
    *Exx_diss = *Exx_pl + Exx_lin + *Exx_pwl + Exx_exp + Exx_gbs + Exx_cst;
    *Ezz_diss =  Ezz_pl + Ezz_lin +  Ezz_pwl + Exx_exp + Exx_gbs + Exx_cst;
    *Exz_diss = *Exz_pl + Exz_lin + *Exz_pwl + Exz_exp + Exz_gbs + Exz_cst;

    *Eii_el  = fabs(Tii-Tii0)/2.0/eta_el;
    *Eii_pwl =  Tii/2.0/eta_pwl;

    // Viscosity for dissipative processes (no elasticity)
    eta        = 1.0/(inv_eta_diss);//Tii/2.0/Eii_vis;
    *VEcoeff   = eta_ve/eta_el;//1.0 / (1.0 + G*model->dt/eta); //eta_ve/eta_el;//  //
    if (elastic==0) *VEcoeff = 0.0;

    // Override viscosty at step 0 (100% visco-plastic)
    if ( model->step == 0 ) *etaVE = eta;
    *txxn = Txx;
    *tzzn = Tzz;
    *txzn = Txz;

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

    return eta;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void NonNewtonianViscosityGrid( grid *mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling ) {
    
    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1, k1, cond;
    double eta, txx1, tzz1, txz1, Pn, Tn, etaVE, VEcoeff=0.0, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1;
    double exx_pwl, exz_pwl, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, exx_pl, exz_pl;
    int average = model->eta_avg;
    double detadexx, detadezz, detadexz, detadp;
    double Xreac, ttrans;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    InterpCentroidsToVerticesDouble( mesh->T,    mesh->T_s,   mesh, model, scaling );
    InterpCentroidsToVerticesDouble( mesh->p_in, mesh->P_s,   mesh, model, scaling );
    InterpCentroidsToVerticesDouble( mesh->d0,   mesh->d0_s,  mesh, model, scaling );
    InterpCentroidsToVerticesDouble( mesh->phi,  mesh->phi_s, mesh, model, scaling ); // ACHTUNG NOT FRICTION ANGLE
    
//    MinMaxArrayTag( mesh->p_in, scaling->S, Ncx*Ncz, "Pn", mesh->BCp.type );
//    MinMaxArrayTag( mesh->P_s, scaling->S, Nx*Nz,   "Ps", mesh->BCg.type    );
//    MinMaxArrayTag( mesh->fric_n, 180.0/M_PI, Ncx*Ncz, "Phin", mesh->BCp.type );
//    MinMaxArrayTag( mesh->fric_s, 180.0/M_PI, Nx*Nz,   "Phis", mesh->BCg.type    );
//    MinMaxArrayTag( mesh->C_n, scaling->S, Ncx*Ncz, "Cn", mesh->BCp.type );
//    MinMaxArrayTag( mesh->C_s, scaling->S, Nx*Nz,   "Cs", mesh->BCg.type    );

    
    // Evaluate cell center viscosities
#pragma omp parallel for shared( mesh  ) private( cond, k, l, k1, p, eta, c1, c0, Pn, Tn, txx1, tzz1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1, exx_pwl, exz_pwl, exx_el, exz_el, exx_diss, exz_diss, exx_pl, exz_pl, detadexx, detadezz, detadexz, detadp, Xreac, ttrans ) firstprivate( materials, scaling, average, model, Ncx, Ncz )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {
        
        //    for ( l=0; l<Ncz; l++ ) {
        //        for ( k=0; k<Ncx; k++ ) {
        
        k      = mesh->kp[k1];
        l      = mesh->lp[k1];
        c0     = k  + l*(Ncx);
        
        mesh->eta_n[c0]      = 0.0;
        mesh->eta_phys_n[c0] = 0.0;
        mesh->VE_n[c0]       = 0.0;
        mesh->sxxd[c0]       = 0.0;
        mesh->szzd[c0]       = 0.0;
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
        mesh->exx_el[c0]     = 0.0;
        mesh->exx_diss[c0]   = 0.0;
        mesh->ezz_el[c0]     = 0.0;
        mesh->ezz_diss[c0]   = 0.0;
        mesh->exx_pl[c0]     = 0.0;
        mesh->detadexx_n[c0] = 0.0;
        mesh->detadezz_n[c0] = 0.0;
        mesh->detadgxz_n[c0] = 0.0;
        mesh->detadp_n[c0]   = 0.0;
        
        mesh->Xreac_n[c0]    = 0.0;
        mesh->ttrans_n[c0]   = 0.0;
        
        // Loop on grid nodes
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
            
            // If background pressure is removed in the RHS
            Pn = mesh->p_in[c0];
            Tn = mesh->T[c0];
            
            // Loop on phases
            for ( p=0; p<model->Nb_phases; p++) {
                
                cond =  fabs(mesh->phase_perc_n[p][c0])>1.0e-13;
                
                if ( cond == 1 ) {
                    eta =  Viscosity( p, mesh->mu_n[c0], Tn, Pn, mesh->d0[c0], mesh->phi[c0], mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1, mesh->strain_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], mesh->ttrans0_n[c0], &Xreac, &ttrans);
                }
                
                // ARITHMETIC AVERAGE
                if (average == 0) {
                    if ( cond == 1 ) mesh->eta_n[c0]      += mesh->phase_perc_n[p][c0] * etaVE;
                    if ( cond == 1 ) mesh->eta_phys_n[c0] += mesh->phase_perc_n[p][c0] * eta;
                    if ( cond == 1 ) mesh->detadexx_n[c0] += mesh->phase_perc_n[p][c0] * detadexx;
                    if ( cond == 1 ) mesh->detadezz_n[c0] += mesh->phase_perc_n[p][c0] * detadezz;
                    if ( cond == 1 ) mesh->detadgxz_n[c0] += mesh->phase_perc_n[p][c0] * detadexz/2.0;
                    if ( cond == 1 ) mesh->detadp_n[c0]   += mesh->phase_perc_n[p][c0] * detadp;
                }
                if (average == 0 || average == 2 ) {
                    if ( cond == 1 ) mesh->sxxd[c0]   += mesh->phase_perc_n[p][c0] * txx1;
                    if ( cond == 1 ) mesh->szzd[c0]   += mesh->phase_perc_n[p][c0] * tzz1;
                }
                if ( cond == 1 ) mesh->VE_n[c0]       += mesh->phase_perc_n[p][c0] * VEcoeff;
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
                if ( cond == 1 ) mesh->exx_el[c0]     += mesh->phase_perc_n[p][c0] * exx_el;
                if ( cond == 1 ) mesh->exx_diss[c0]   += mesh->phase_perc_n[p][c0] * exx_diss;
                if ( cond == 1 ) mesh->ezz_el[c0]     += mesh->phase_perc_n[p][c0] * ezz_el;
                if ( cond == 1 ) mesh->ezz_diss[c0]   += mesh->phase_perc_n[p][c0] * ezz_diss;
                if ( cond == 1 ) mesh->exx_pl[c0]     += mesh->phase_perc_n[p][c0] * exx_pl;
                
                if ( cond == 1 ) mesh->Xreac_n[c0]     += mesh->phase_perc_n[p][c0] * Xreac;
                if ( cond == 1 ) mesh->ttrans_n[c0]    += mesh->phase_perc_n[p][c0] * ttrans;
                
                // HARMONIC AVERAGE
                if (average == 1) {
                    if ( cond == 1 ) mesh->sxxd[c0]       += mesh->phase_perc_n[p][c0] * 1.0/txx1;
                    if ( cond == 1 ) mesh->szzd[c0]       += mesh->phase_perc_n[p][c0] * 1.0/tzz1;
                    if ( cond == 1 ) mesh->eta_n[c0]      += mesh->phase_perc_n[p][c0] * 1.0/etaVE;
                    if ( cond == 1 ) mesh->eta_phys_n[c0] += mesh->phase_perc_n[p][c0] * 1.0/eta;
                    if ( cond == 1 ) mesh->detadexx_n[c0] += mesh->phase_perc_n[p][c0] * detadexx     / pow(etaVE,2.0);
                    if ( cond == 1 ) mesh->detadezz_n[c0] += mesh->phase_perc_n[p][c0] * detadezz     / pow(etaVE,2.0);
                    if ( cond == 1 ) mesh->detadgxz_n[c0] += mesh->phase_perc_n[p][c0] * detadexz/2.0 / pow(etaVE,2.0);
                    if ( cond == 1 ) mesh->detadp_n[c0]   += mesh->phase_perc_n[p][c0] * detadp       / pow(etaVE,2.0);
                    
                }
                
                // GEOMETRIC AVERAGE
                if (average == 2) {
                    if ( cond == 1 ) mesh->eta_n[c0]      += mesh->phase_perc_n[p][c0] * log(etaVE);
                    if ( cond == 1 ) mesh->eta_phys_n[c0] += mesh->phase_perc_n[p][c0] * log(eta);
                    if ( cond == 1 ) mesh->detadexx_n[c0] += mesh->phase_perc_n[p][c0] * detadexx     / etaVE;
                    if ( cond == 1 ) mesh->detadezz_n[c0] += mesh->phase_perc_n[p][c0] * detadezz     / etaVE;
                    if ( cond == 1 ) mesh->detadgxz_n[c0] += mesh->phase_perc_n[p][c0] * detadexz/2.0 / etaVE;
                    if ( cond == 1 ) mesh->detadp_n[c0]   += mesh->phase_perc_n[p][c0] * detadp       / etaVE;
//                    if ( cond == 1 ) mesh->sxxd[c0]       += mesh->phase_perc_n[p][c0] * log(txx1);
//                    if ( cond == 1 ) mesh->szzd[c0]       += mesh->phase_perc_n[p][c0] * log(tzz1);
                }
            }
            
            mesh->d[c0]          = 1.0/mesh->d[c0];
            
            // HARMONIC AVERAGE
            if (average == 1) {
                mesh->sxxd[c0]       = 1.0/mesh->sxxd[c0];
                mesh->szzd[c0]       = 1.0/mesh->szzd[c0];
                mesh->eta_n[c0]      = 1.0/mesh->eta_n[c0];
                mesh->eta_phys_n[c0] = 1.0/mesh->eta_phys_n[c0];
                mesh->detadexx_n[c0] *= pow(mesh->eta_n[c0],2.0);
                mesh->detadezz_n[c0] *= pow(mesh->eta_n[c0],2.0);
                mesh->detadgxz_n[c0] *= pow(mesh->eta_n[c0],2.0);
                mesh->detadp_n[c0]   *= pow(mesh->eta_n[c0],2.0);
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
                mesh->eta_n[c0]      = exp(mesh->eta_n[c0]);
                mesh->eta_phys_n[c0] = exp(mesh->eta_phys_n[c0]);
                mesh->detadexx_n[c0] *= mesh->eta_n[c0];
                mesh->detadezz_n[c0] *= mesh->eta_n[c0];
                mesh->detadgxz_n[c0] *= mesh->eta_n[c0];
                mesh->detadp_n[c0]   *= mesh->eta_n[c0];
            }
//            printf("eta_n = %2.2e\n", mesh->eta_n[c0]*scaling->eta);

        }
        //if (mesh->eta_n[c0]<1e-8) exit(1);
    }
    ;
    
    // Calculate vertices viscosity
    double d1s; // dummy variable that stores updated grain size on current vertice
    
#pragma omp parallel for shared( mesh, model ) private( cond, k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1s, exx_pwl, exz_pwl, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, exx_pl, exz_pl, detadexx, detadezz, detadexz, detadp , Xreac, ttrans) firstprivate( materials, scaling, average, Nx, Nz  )
    for ( k1=0; k1<Nx*Nz; k1++ ) {
        
        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;
        
        mesh->VE_s[c1]       = 0.0;
        mesh->sxz[c1]        = 0.0;
        mesh->eta_phys_s[c1] = 0.0;
        mesh->eta_s[c1]      = 0.0;
        mesh->eII_pwl_s[c1]  = 0.0;
        mesh->exx_pwl_s[c1]  = 0.0;
        mesh->exz_pwl_s[c1]  = 0.0;
        mesh->A2_pwl_s[c1]   = 0.0;
        mesh->exz_el[c1]     = 0.0;
        mesh->exz_diss[c1]   = 0.0;
        mesh->eII_pl_s[c1]   = 0.0;
        mesh->exz_pl[c1]     = 0.0;
        mesh->detadexx_s[c1] = 0.0;
        mesh->detadezz_s[c1] = 0.0;
        mesh->detadgxz_s[c1] = 0.0;
        mesh->detadp_s[c1]   = 0.0;
        
//        mesh->Xreac_s[c1]    = 0.0;
//        mesh->ttrans_s[c1]   = 0.0;
        
        if ( mesh->BCg.type[c1] != 30 ) {
            
            
            for ( p=0; p<model->Nb_phases; p++) {
                
                cond = fabs(mesh->phase_perc_s[p][c1])>1.0e-13;
                
                if ( cond == 1 ) {
                    
                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi_s[c1], mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1s, mesh->strain_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp,mesh->p0_s[c1], mesh->ttrans0_s[c1], &Xreac, &ttrans);
                    
//                    if (c1==546) printf("eta_s = %2.2e index = %d\n", etaVE*scaling->eta, c1);

                }
                
                if (average ==0) {
                    if ( cond == 1 ) mesh->eta_s[c1]      += mesh->phase_perc_s[p][c1] * etaVE;
                    if ( cond == 1 ) mesh->eta_phys_s[c1] += mesh->phase_perc_s[p][c1] * eta;
                    if ( cond == 1 ) mesh->detadexx_s[c1] += mesh->phase_perc_s[p][c1] * detadexx;
                    if ( cond == 1 ) mesh->detadezz_s[c1] += mesh->phase_perc_s[p][c1] * detadezz;
                    if ( cond == 1 ) mesh->detadgxz_s[c1] += mesh->phase_perc_s[p][c1] * detadexz/2.0;
                    if ( cond == 1 ) mesh->detadp_s[c1]   += mesh->phase_perc_s[p][c1] * detadp;
                }
                if (average ==0 || average==2) {
                    if ( cond == 1 ) mesh->sxz[c1]    += mesh->phase_perc_s[p][c1] * txz1;
                }
                if ( cond == 1 ) mesh->VE_s[c1]       += mesh->phase_perc_s[p][c1] * VEcoeff;
                if ( cond == 1 ) mesh->eII_pl_s[c1]   += mesh->phase_perc_s[p][c1] * eII_pl;
                if ( cond == 1 ) mesh->eII_pwl_s[c1]  += mesh->phase_perc_s[p][c1] * eII_pwl;
                if ( cond == 1 ) mesh->exx_pwl_s[c1]  += mesh->phase_perc_s[p][c1] * exx_pwl;
                if ( cond == 1 ) mesh->exz_pwl_s[c1]  += mesh->phase_perc_s[p][c1] * exz_pwl;
                if ( cond == 1 ) mesh->exz_el[c1]     += mesh->phase_perc_s[p][c1] * exz_el;
                if ( cond == 1 ) mesh->exz_diss[c1]   += mesh->phase_perc_s[p][c1] * exz_diss;
                if ( cond == 1 ) mesh->exz_pl[c1]     += mesh->phase_perc_s[p][c1] * exz_pl;
                
//                if ( cond == 1 ) mesh->Xreac_s[c1]    += mesh->phase_perc_s[p][c1] * Xreac;
//                if ( cond == 1 ) mesh->ttrans_s[c1]   += mesh->phase_perc_s[p][c1] * ttrans;
            
                if (average == 1) {
                    if ( cond == 1 ) mesh->sxz[c1]        += mesh->phase_perc_s[p][c1] * 1.0/txz1;
                    if ( cond == 1 ) mesh->eta_s[c1]      += mesh->phase_perc_s[p][c1] * 1.0/etaVE;
                    if ( cond == 1 ) mesh->eta_phys_s[c1] += mesh->phase_perc_s[p][c1] * 1.0/eta;
                    if ( cond == 1 ) mesh->detadexx_s[c1] += mesh->phase_perc_s[p][c1] * detadexx      / pow(etaVE,2.0);
                    if ( cond == 1 ) mesh->detadezz_s[c1] += mesh->phase_perc_s[p][c1] * detadezz      / pow(etaVE,2.0);
                    if ( cond == 1 ) mesh->detadgxz_s[c1] += mesh->phase_perc_s[p][c1] * detadexz/2.0  / pow(etaVE,2.0);
                    if ( cond == 1 ) mesh->detadp_s[c1]   += mesh->phase_perc_s[p][c1] * detadp        / pow(etaVE,2.0);
                }
                if (average == 2) {
                    if ( cond == 1 ) mesh->eta_s[c1]      += mesh->phase_perc_s[p][c1] * log(etaVE);
                    if ( cond == 1 ) mesh->eta_phys_s[c1] += mesh->phase_perc_s[p][c1] * log(eta);
                    if ( cond == 1 ) mesh->detadexx_s[c1] += mesh->phase_perc_s[p][c1] * detadexx      / etaVE;
                    if ( cond == 1 ) mesh->detadezz_s[c1] += mesh->phase_perc_s[p][c1] * detadezz      / etaVE;
                    if ( cond == 1 ) mesh->detadgxz_s[c1] += mesh->phase_perc_s[p][c1] * detadexz/2.0  / etaVE;
                    if ( cond == 1 ) mesh->detadp_s[c1]   += mesh->phase_perc_s[p][c1] * detadp        / etaVE;
//                    if ( cond == 1  ) mesh->sxz[c1] += mesh->phase_perc_s[p][c1]        *  log(txz1);
                }
            }
            // HARMONIC AVERAGE
            if (average == 1) {
                mesh->sxz[c1]        = 1.0/mesh->sxz[c1];
                mesh->eta_s[c1]      = 1.0/mesh->eta_s[c1];
                mesh->eta_phys_s[c1] = 1.0/mesh->eta_phys_s[c1];
                mesh->detadexx_s[c1] *= pow(mesh->eta_s[c1],2.0);
                mesh->detadezz_s[c1] *= pow(mesh->eta_s[c1],2.0);
                mesh->detadgxz_s[c1] *= pow(mesh->eta_s[c1],2.0);
                mesh->detadp_s[c1]   *= pow(mesh->eta_s[c1],2.0);
                if (isinf (mesh->eta_phys_s[c1]) ) {
                    for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_s[p][c1]);
                    printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e \n", mesh->mu_s[c1],  mesh->eii_s[c1],  mesh->tii0_s[c1], mesh->exxd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->sxz0[c1]);
                    printf("x=%2.2e z=%2.2e\n", mesh->xg_coord[k]*scaling->L/1000, mesh->zg_coord[l]*scaling->L/1000);
                    
                    exit(1);
                }
                if (isnan (mesh->eta_phys_s[c1]) ) {
                    for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_s[p][c1]);
                    printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e \n", mesh->mu_s[c1],  mesh->eii_s[c1],  mesh->tii0_s[c1], mesh->exxd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->sxz0[c1]);
                    printf("x=%2.2e z=%2.2e\n", mesh->xg_coord[k]*scaling->L/1000, mesh->zg_coord[l]*scaling->L/1000);
                    
                    exit(1);
                }
            }
            // GEOMETRIC AVERAGE
            if (average == 2) {
                mesh->eta_s[c1]       = exp(mesh->eta_s[c1]);
                mesh->eta_phys_s[c1]  = exp(mesh->eta_phys_s[c1]);
                mesh->detadexx_s[c1] *= mesh->eta_s[c1];
                mesh->detadezz_s[c1] *= mesh->eta_s[c1];
                mesh->detadgxz_s[c1] *= mesh->eta_s[c1];
                mesh->detadp_s[c1]   *= mesh->eta_s[c1];
            }
//            if (c1==546) printf("eta_s = %2.2e index = %d\n", etaVE*scaling->eta, c1);
        }
    }
    
//    MinMaxArrayTag( mesh->eta_n, scaling->eta, Ncx*Ncz, "eta_n", mesh->BCp.type );
//    MinMaxArrayTag( mesh->eta_s, scaling->eta, Nx*Nz,   "eta_s", mesh->BCg.type );
//    MinMaxArrayTag( mesh->detadexx_n, scaling->eta/scaling->E, Ncx*Ncz, "detadexx_n", mesh->BCp.type );
//    MinMaxArrayTag( mesh->detadexx_s, scaling->eta/scaling->E, Nx*Nz,   "detadexx_s", mesh->BCg.type );
//    MinMaxArrayTag( mesh->fric_n, 180.0/M_PI, Ncx*Ncz, "Phin", mesh->BCp.type );
//    MinMaxArrayTag( mesh->fric_s, 180.0/M_PI, Nx*Nz,   "Phis", mesh->BCg.type    );
//    MinMaxArrayTag( mesh->C_n, scaling->S, Ncx*Ncz, "Cn", mesh->BCp.type );
//    MinMaxArrayTag( mesh->C_s, scaling->S, Nx*Nz,   "Cs", mesh->BCg.type    );

    
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
            mesh->fric_n[c0] = 0.0;
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
                            //   printf("Allo 1 - %d!\n", p);
                            phi = materials.phi[p];
                            C   = materials.C[p];
                        }
                        // If we are above the upper strain limit
                        if (strain_acc >= materials.pls_end[p]) {
                            // printf("Allo 2 - %d!\n", p);
                            phi = materials.phi_end[p];
                            C   = materials.C_end[p];
                        }
                        // If we are in the softening strain range
                        if (strain_acc >= materials.pls_start[p] && strain_acc < materials.pls_end[p] ) {
                            //   printf("Allo 3 - %d!\n", p);
                            phi = materials.phi[p] - (materials.phi[p] - materials.phi_end[p]) * ( strain_acc / (materials.pls_end[p] - materials.pls_start[p]) );
                            C   = materials.C[p]   + (  materials.C_end[p] -   materials.C[p]) * MINV( 1.0, strain_acc /materials.pls_end[p] );
                        }
                        
                    }
                    
                    // Arithmetic
                    if (average ==0) {
                        mesh->fric_n[c0] += mesh->phase_perc_n[p][c0] * phi;
                        mesh->C_n[c0]   += mesh->phase_perc_n[p][c0] * C;
                    }
                    // Harmonic
                    if (average == 1) {
                        mesh->fric_n[c0] += mesh->phase_perc_n[p][c0] *  1.0/phi;
                        mesh->C_n[c0]   += mesh->phase_perc_n[p][c0] *  1.0/C;
                    }
                    // Geometric
                    if (average == 2) {
                        mesh->fric_n[c0] += mesh->phase_perc_n[p][c0] *  log(phi);
                        mesh->C_n[c0]   += mesh->phase_perc_n[p][c0] *  log(C);
                        
                    }
                }
                // Post-process for geometric/harmonic averages
                if ( average==1 ) mesh->fric_n[c0] = 1.0/mesh->fric_n[c0];
                if ( average==2 ) mesh->fric_n[c0] = exp(mesh->fric_n[c0]);
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
            mesh->fric_s[c1] = 0.0;
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
                        mesh->fric_s[c1] += mesh->phase_perc_s[p][c1] * phi;
                        mesh->C_s[c1]   += mesh->phase_perc_s[p][c1] * C;
                    }
                    // Harmonic
                    if (average == 1) {
                        mesh->fric_s[c1] += mesh->phase_perc_s[p][c1] *  1.0/phi;
                        mesh->C_s[c1]   += mesh->phase_perc_s[p][c1] *  1.0/C;
                    }
                    // Geometric
                    if (average == 2) {
                        mesh->fric_s[c1] += mesh->phase_perc_s[p][c1] *  log(phi);
                        mesh->C_s[c1]   += mesh->phase_perc_s[p][c1] *  log(C);
                    }
                }
                // Post-process for geometric/harmonic averages
                if ( average==1 ) mesh->fric_s[c1] = 1.0/mesh->fric_s[c1];
                if ( average==2 ) mesh->fric_s[c1] = exp(mesh->fric_s[c1]);
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
    int average = 1;//%model.eta_avg; // SHOULD NOT BE ALLOWED TO BE ELSE THAN 1
    
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
void UpdateParticlettrans( grid *mesh, scale *scaling, params model, markers *particles, mat_prop *materials ) {
    
    DoodzFP *ttrans_inc_mark, *ttrans_inc_grid;
    double dx=model.dx, dz=model.dz;
    int Nx, Nz, Ncx, Ncz, k, c0, p;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
//    printf("=============== UpdateParticlettrans ============\n");
//    MinMaxArrayTag( mesh->ttrans0_n,   scaling->t,    Ncx*Ncz,   "ttrans0_n",   mesh->BCp.type );
//    MinMaxArrayTag( mesh->ttrans_n,   scaling->t,    Ncx*Ncz,   "ttrans_n",   mesh->BCp.type );
    
    //printf("Pc = %2.2e; Preac = %2.2e; Pold = %2.2e\n", Pc*scaling->S, Preac*scaling->S, Pold*scaling->S);
    
    ttrans_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    ttrans_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    
    for (k=0;k<Ncx*Ncz;k++) {
        ttrans_inc_grid[k] =0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) ttrans_inc_grid[k] = mesh->ttrans_n[k] - mesh->ttrans0_n[k];
    }
    
//    MinMaxArrayTag( mesh->ttrans0_n,   scaling->t,    Ncx*Ncz,   "ttrans0_n",   mesh->BCp.type );
//    MinMaxArrayTag( mesh->ttrans_n,   scaling->t,    Ncx*Ncz,   "ttrans_n",   mesh->BCp.type );
    
    Interp_Grid2P( *particles, ttrans_inc_mark, mesh, ttrans_inc_grid, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCt.type  );
    ArrayPlusArray( particles->ttrans, ttrans_inc_mark, particles->Nb_part );
    
    DoodzFree(ttrans_inc_mark);
    DoodzFree(ttrans_inc_grid);
    printf("=============== UpdateParticlettrans END============\n");
    
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
    

#pragma omp parallel for shared(mesh,materials) private( NT, NP, iT, iP, TW, PW, iSW, iSE, iNW, iNE, phase_diag, dT, dP, Tgrid, Pgrid, c0, p, rhonew, rho0, alpha, beta, T0, P0, rhop, drho, percT, percP) firstprivate(Ncx, Ncz, model, epsi)
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {
        rhonew = 0.0;
        // Loop on phases
        for ( p=0; p<model->Nb_phases; p++) {
            
            if ( fabs(mesh->phase_perc_n[p][c0])>epsi) {
                
                // Constant density
                if ( materials->density_model[p] == 3 ) {
                    rho0   = materials->rho[p];
                    P0     = 0;//materials->P0 [p];
                    beta   = materials->bet[p];
                    rhop   = rho0 * (1.0 +  beta * (mesh->p_in[c0] - P0) );
                    rhop   = rho0*exp(beta * mesh->p_in[c0] );
//                    rhop   = rho0*(1 + beta * mesh->p_in[c0] );
//                    rhop   = rho0;
//                    rhop   = rho0*(1.0 +  beta * (mesh->p_in[c0] - P0));
//                    printf("DOING %2.2e  %2.2e %2.2e %2.2e %2.2e\n",  rhop, mesh->p_in[c0], beta, rhop, P0);
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
    
    InterpCentroidsToVerticesDouble( mesh->rho_n, mesh->rho_s, mesh, model, scaling );
    
//    // Interpolate center values to vertices
//    int l, c1;
//#pragma omp parallel for shared( mesh ) private( c1, k, l, c0) firstprivate( Nx, Nz, Ncx, Ncz )
//    for (c1=0; c1<Nx*Nz; c1++) {
//
//        k = mesh->kn[c1];
//        l = mesh->ln[c1];
//        c1 = k + l*Nx;
//        c0 = k + l*(Ncx);
//        mesh->rho_s[c1]   = 0.0;
//
//        if ( mesh->BCg.type[c1] != 30 ) {
//
//            // Inner grid
//            if ( k>0 && l>0 && k<Ncx && l<Ncz ) {
//                mesh->rho_s[c1] = 0.25*( mesh->rho_n[c0] + mesh->rho_n[c0-1] + mesh->rho_n[c0-Ncx] + mesh->rho_n[c0-Ncx-1]  );
//            }
//            // Sides
//            if (k==0 && (l>0 && l<Ncz)) {
//                mesh->rho_s[c1] = 0.5*( mesh->rho_n[c0] + mesh->rho_n[c0-Ncx] );
//            }
//            if (k==Ncx && (l>0 && l<Ncz)) {
//                mesh->rho_s[c1] = 0.5*( mesh->rho_n[c0-1] + mesh->rho_n[c0-Ncx-1] );
//            }
//            if (l==0 && (k>0 && k<Ncx)) {
//                mesh->rho_s[c1] = 0.5*( mesh->rho_n[c0] + mesh->rho_n[c0-1] );
//            }
//            if (l==Ncz && (k>0 && k<Ncx)) {
//                mesh->rho_s[c1] = 0.5*( mesh->rho_n[c0-Ncx] + mesh->rho_n[c0-Ncx-1] );
//            }
//            // Corners
//            if (l==0 && k==0) {
//                mesh->rho_s[c1] = mesh->rho_n[c0];
//            }
//            if (l==Ncz && k==0) {
//                mesh->rho_s[c1] = mesh->rho_n[c0-Ncx];
//            }
//            if (l==0 && k==Ncx) {
//                mesh->rho_s[c1] = mesh->rho_n[c0-1];
//            }
//            if (l==Ncz && k==Ncx) {
//                mesh->rho_s[c1] = mesh->rho_n[c0-Ncx-1];
//            }
//        }
//    }
    
    printf("Updated density fields:\n");
//    MinMaxArrayTag(        mesh->X,            1, Ncx*Ncz,     "X", mesh->BCp.type );
//    MinMaxArrayTag( mesh->rho_n, scaling->rho, Ncx*Ncz, "rho_n", mesh->BCp.type );
//    MinMaxArrayTag( mesh->rho_s, scaling->rho, Nx*Nz,   "rho_s", mesh->BCg.type    );
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Strain rate
void  StrainRateComponents( grid* mesh, scale scaling, params* model ) {
    
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
            mesh->exxd[c0]  = (mesh->u_in[c1+1+Nx]     - mesh->u_in[c1+Nx] )/dx - 1.0/3.0*mesh->div_u[c0];
            mesh->ezzd[c0]  = (mesh->v_in[c2+1+(Nx+1)] - mesh->v_in[c2+1]  )/dz - 1.0/3.0*mesh->div_u[c0];
        }
        else {
            mesh->div_u[c0] = 0.0;
            mesh->exxd[c0]  = 0.0;
            mesh->ezzd[c0]  = 0.0;
        }
    }
    
    // Shear components we only calculate sxz because sxz = szx (stress tensor symmetry)
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
    InterpCentroidsToVerticesDouble( mesh->exxd, mesh->exxd_s, mesh, model, &scaling );
    InterpCentroidsToVerticesDouble( mesh->ezzd, mesh->ezzd_s, mesh, model, &scaling );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void GenerateDeformationMaps( grid* mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling ) {
    
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
    double txx1, tzz1, txz1, etaVE, VEcoeff, eII_el=0.0, eII_pl=0.0, eII_pwl=0.0, eII_exp=0.0, eII_lin=0.0, eII_gbs=0.0, eII_cst=0.0, d1;
    double exx_pwl, exz_pwl, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, exx_pl, exz_pl;
    double Pn = model->Pn , eta;
    int    *dom_mech, ind, mech, loud=0, k;
    double *stress, *visco, t_omp;
    double detadexx, detadezz, detadexz, detadp;
    double Xreac, ttrans;
    
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
                    eta =  Viscosity( k, 0.0, T[ix], Pn, d[iz], 0.0, E[iy], E[iy], 0.0, 0.0, 0.0, 0.0, materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl,  &d1, 0.0, materials->phi[k], materials->C[k], &detadexx, &detadezz, &detadexz, &detadp, 0.0, 0.0, &Xreac, &ttrans);
                    
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

//void RotateDirectorVector( grid mesh, markers* particles, params model, scale *scaling ) {
//
//    int k;
//    double angle, nx, nz, norm;
//
//#pragma omp parallel for shared ( particles ) private ( angle, k, nx, nz, norm ) firstprivate( model ) schedule( static )
//    for(k=0; k<particles->Nb_part; k++) {
//
//        // Filter out particles that are inactive (out of the box)
//        if (particles->phase[k] != -1) {
//
//            // Angle
//            angle = -model.dt*particles->om_p[k];
//
//            nx = particles->nx[k] * cos(angle) - particles->nz[k] * sin(angle) ;
//            nz = particles->nx[k] * sin(angle) + particles->nz[k] * cos(angle) ;
//            norm             = sqrt( nx*nx + nz*nz );
//            particles->nx[k] = nx/norm;
//            particles->nz[k] = nz/norm;
//
//        }
//    }
//
//}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ComputeViscosityDerivatives_FD( grid* mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling ) {
    
    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1, k1, cond;
    double eta, txx1, tzz1, txz1, Pn, Tn, Ps, Ts, etaVE, VEcoeff=0.0, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1;
    double exx_pwl, exz_pwl, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, exx_pl, exz_pl;
    int average = model->eta_avg;
    double detadexx, detadezz, detadexz, detadp;
    double Xreac, ttrans;
    
    double eta_exx, eta_ezz, eta_exz, eta_p;
    double etaVE_exx, etaVE_ezz, etaVE_exz, etaVE_p;
    double eps = 1e-5, pert_xx, pert_zz, pert_xz , pert_p;
    double eps1=1e-13, eii;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;
    
    printf("---> Computing numerical derivatives\n");
    
    // Evaluate cell center viscosities
#pragma omp parallel for shared( mesh  ) private( eii, cond, k, l, k1, p, eta, c1, c0, Pn, Tn, txx1, tzz1, txz1, etaVE, etaVE_exx, etaVE_ezz, etaVE_exz, etaVE_p, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1, exx_pwl, exz_pwl, exx_el, exz_el, exx_diss, exz_diss, exx_pl, exz_pl, eta_exx, eta_ezz, eta_exz, eta_p, pert_xx, pert_zz, pert_xz , pert_p, detadexx, detadezz, detadexz, detadp, Xreac, ttrans ) firstprivate( materials, scaling, average, model, Ncx, Ncz, eps, eps1 )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {
        
        k      = mesh->kp[k1];
        l      = mesh->lp[k1];
        c0     = k  + l*(Ncx);
        
        mesh->detadexx_n[c0]      = 0.0;
        mesh->detadezz_n[c0]      = 0.0;
        mesh->detadgxz_n[c0]      = 0.0;
        mesh->detadp_n[c0]        = 0.0;
        
        eta_exx = 0.0;
        eta_ezz = 0.0;
        eta_exz = 0.0;
        eta_p   = 0.0;
        
        // Loop on grid nodes
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
            
            Pn = mesh->p_in[c0];
            Tn = mesh->T[c0];
            
            eii = sqrt(0.5*pow(mesh->exxd[c0],2)+0.5*pow(mesh->ezzd[c0],2)+pow(mesh->exz_n[c0],2));
            pert_xx = eps*eii;
            pert_zz = eps*eii;
            pert_xz = eps*eii;
            pert_p  = eps*Pn;//eps1;
            pert_p  = eps1;
            
            // Loop on phases
            for ( p=0; p<model->Nb_phases; p++) {
                
                cond =  fabs(mesh->phase_perc_n[p][c0])>1.0e-13;
                
                if ( cond == 1 ) {
                    
                    eta =  Viscosity( p, mesh->mu_n[c0], Tn, Pn, mesh->d0[c0], mesh->phi[c0], mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1, mesh->strain_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], mesh->ttrans0_n[c0], &Xreac, &ttrans );
                    
                    eta =  Viscosity( p, mesh->mu_n[c0], Tn, Pn, mesh->d0[c0], mesh->phi[c0], mesh->exxd[c0]+pert_xx, mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling, &txx1, &tzz1, &txz1, &etaVE_exx, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1, mesh->strain_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], mesh->ttrans0_n[c0], &Xreac, &ttrans);
                    
                    eta =  Viscosity( p, mesh->mu_n[c0], Tn, Pn, mesh->d0[c0], mesh->phi[c0], mesh->exxd[c0], mesh->ezzd[c0]+pert_zz, mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling,  &txx1, &tzz1, &txz1, &etaVE_ezz, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1, mesh->strain_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], mesh->ttrans0_n[c0], &Xreac, &ttrans);
                    
                    eta =  Viscosity( p, mesh->mu_n[c0], Tn, Pn, mesh->d0[c0], mesh->phi[c0], mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0]+pert_xz, mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling,  &txx1, &tzz1, &txz1, &etaVE_exz, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1, mesh->strain_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], mesh->ttrans0_n[c0], &Xreac, &ttrans);
                    
                    eta =  Viscosity( p, mesh->mu_n[c0], Tn, Pn+pert_p, mesh->d0[c0], mesh->phi[c0], mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling,  &txx1, &tzz1, &txz1, &etaVE_p  , &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1, mesh->strain_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], mesh->ttrans0_n[c0], &Xreac, &ttrans);
                    
                }
                
                // ARITHMETIC AVERAGE
                if (average == 0) {
                    if ( cond == 1 ) eta_exx   += mesh->phase_perc_n[p][c0] * etaVE_exx;
                    if ( cond == 1 ) eta_ezz   += mesh->phase_perc_n[p][c0] * etaVE_ezz;
                    if ( cond == 1 ) eta_exz   += mesh->phase_perc_n[p][c0] * etaVE_exz;
                    if ( cond == 1 ) eta_p     += mesh->phase_perc_n[p][c0] * etaVE_p;
                }
                
                // HARMONIC AVERAGE
                if (average == 1) {
                    if ( cond == 1 ) eta_exx   += mesh->phase_perc_n[p][c0] * 1.0/etaVE_exx;
                    if ( cond == 1 ) eta_ezz   += mesh->phase_perc_n[p][c0] * 1.0/etaVE_ezz;
                    if ( cond == 1 ) eta_exz   += mesh->phase_perc_n[p][c0] * 1.0/etaVE_exz;
                    if ( cond == 1 ) eta_p     += mesh->phase_perc_n[p][c0] * 1.0/etaVE_p;
                }
                
                // GEOMETRIC AVERAGE
                if (average == 2) {
                    if ( cond == 1 ) eta_exx   += mesh->phase_perc_n[p][c0] * log(etaVE_exx);
                    if ( cond == 1 ) eta_ezz   += mesh->phase_perc_n[p][c0] * log(etaVE_ezz);
                    if ( cond == 1 ) eta_exz   += mesh->phase_perc_n[p][c0] * log(etaVE_exz);
                    if ( cond == 1 ) eta_p     += mesh->phase_perc_n[p][c0] * log(etaVE_p);
                }
                
                //                // General FD
                //                mesh->detadexx_n[c0]     += mesh->phase_perc_n[p][c0] * (etaVE_exx - etaVE) / pert_xx;
                //                mesh->detadezz_n[c0]     += mesh->phase_perc_n[p][c0] * (etaVE_ezz - etaVE) / pert_zz;
                //                mesh->detadgxz_n[c0]     += mesh->phase_perc_n[p][c0] * (etaVE_exz - etaVE) / pert_xz / 2.0;
                //                mesh->detadp_n[c0]       += mesh->phase_perc_n[p][c0] * (etaVE_p   - etaVE) / pert_p;
                
            }
            
            // HARMONIC AVERAGE
            if (average == 1) {
                eta_exx      = 1.0/eta_exx;
                eta_ezz      = 1.0/eta_ezz;
                eta_exz      = 1.0/eta_exz;
                eta_p        = 1.0/eta_p;
            }
            // GEOMETRIC AVERAGE
            if (average == 2) {
                eta_exx   = exp(eta_exx);
                eta_ezz   = exp(eta_ezz);
                eta_exz   = exp(eta_exz);
                eta_p     = exp(eta_p);
            }
            
            // General FD
            mesh->detadexx_n[c0]      = (eta_exx - mesh->eta_n[c0]) / pert_xx;
            mesh->detadezz_n[c0]      = (eta_ezz - mesh->eta_n[c0]) / pert_zz;
            mesh->detadgxz_n[c0]      = (eta_exz - mesh->eta_n[c0]) / pert_xz / 2.0;
            mesh->detadp_n[c0]        = (eta_p   - mesh->eta_n[c0]) / pert_p;
            
            //        printf("pert_p = %2.2e Pn = %2.2e\n", pert_p, Pn);
            //            if (isnan(mesh->detadp_n[c0])) {
            //                printf("%2.2e %2.2e %2.2e\n", pert_p, Pn, mesh->detadp_n[c0] );
            //                exit(1);
            //            }
            
        }
        
    }
    
    
    // Calculate vertices viscosity
    double d1s; // dummy variable that stores updated grain size on current vertice
#pragma omp parallel for shared( mesh, model ) private( eii, cond, k, l, k1, p, eta, c1, c0, Ps, Ts, txx1, tzz1, txz1, etaVE, etaVE_exx, etaVE_ezz, etaVE_exz, etaVE_p, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1s, exx_pwl, exz_pwl, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, exx_pl, exz_pl, eta_exx, eta_ezz, eta_exz, eta_p, pert_xx, pert_zz, pert_xz , pert_p, detadexx, detadezz, detadexz, detadp, Xreac, ttrans ) firstprivate( materials, scaling, average, Nx, Nz, eps, eps1  )
    for ( k1=0; k1<Nx*Nz; k1++ ) {
        
        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;
        
        mesh->detadexx_s[c1]      = 0.0;
        mesh->detadezz_s[c1]      = 0.0;
        mesh->detadgxz_s[c1]      = 0.0;
        mesh->detadp_s[c1]        = 0.0;
        
        eta_exx = 0.0;
        eta_ezz = 0.0;
        eta_exz = 0.0;
        eta_p   = 0.0;
        
        if ( mesh->BCg.type[c1] != 30 ) {
            
            Ps = mesh->P_s[c1];
            Ts = mesh->T_s[c1];
            
            eii = sqrt(pow(0.5*mesh->exxd_s[c1],2)+0.5*pow(mesh->ezzd_s[c1],2)+pow(mesh->exz[c1],2));
            pert_xx = eps*eii;
            pert_zz = eps*eii;
            pert_xz = eps*eii;
            pert_p  = eps*Ps;//eps1;
            pert_p  = eps1;
            
            //            pert_xx = eps1;
            //            pert_zz = eps1;
            //            pert_xz = eps1;
            //            pert_p = eps1;
            
            for ( p=0; p<model->Nb_phases; p++) {
                
                cond = fabs(mesh->phase_perc_s[p][c1])>1.0e-13;
                
                if ( cond == 1 ) {
                    
                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi_s[c1], mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1s, mesh->strain_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], mesh->ttrans0_n[c0], &Xreac, &ttrans);
                    
                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi_s[c1], mesh->exxd_s[c1]+pert_xx, mesh->ezzd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE_exx, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1s, mesh->strain_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], mesh->ttrans0_n[c0], &Xreac, &ttrans);
                    
                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi_s[c1], mesh->exxd_s[c1], mesh->ezzd_s[c1]+pert_zz, mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE_ezz, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1s, mesh->strain_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], mesh->ttrans0_n[c0], &Xreac, &ttrans);
                    
                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi_s[c1], mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1]+pert_xz, mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE_exz, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1s, mesh->strain_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], mesh->ttrans0_n[c0], &Xreac, &ttrans);
                    
                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1]+pert_p, mesh->d0_s[c1], mesh->phi_s[c1], mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE_p, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_pwl, &exz_pwl, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &exx_pl, &exz_pl, &d1s, mesh->strain_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], mesh->ttrans0_n[c0], &Xreac, &ttrans);
                    
                }
                
                if (average ==0) {
                    if ( cond == 1 ) eta_exx   += mesh->phase_perc_s[p][c1] * etaVE_exx;
                    if ( cond == 1 ) eta_ezz   += mesh->phase_perc_s[p][c1] * etaVE_ezz;
                    if ( cond == 1 ) eta_exz   += mesh->phase_perc_s[p][c1] * etaVE_exz;
                    if ( cond == 1 ) eta_p     += mesh->phase_perc_s[p][c1] * etaVE_p;
                }
                
                // HARMONIC AVERAGE
                if (average == 1) {
                    if ( cond == 1 ) eta_exx   += mesh->phase_perc_s[p][c1] * 1.0/etaVE_exx;
                    if ( cond == 1 ) eta_ezz   += mesh->phase_perc_s[p][c1] * 1.0/etaVE_ezz;
                    if ( cond == 1 ) eta_exz   += mesh->phase_perc_s[p][c1] * 1.0/etaVE_exz;
                    if ( cond == 1 ) eta_p     += mesh->phase_perc_s[p][c1] * 1.0/etaVE_p;
                }
                
                // GEOMETRIC AVERAGE
                if (average == 2) {
                    if ( cond == 1 ) eta_exx   += mesh->phase_perc_s[p][c1] * log(etaVE_exx);
                    if ( cond == 1 ) eta_ezz   += mesh->phase_perc_s[p][c1] * log(etaVE_ezz);
                    if ( cond == 1 ) eta_exz   += mesh->phase_perc_s[p][c1] * log(etaVE_exz);
                    if ( cond == 1 ) eta_p     += mesh->phase_perc_s[p][c1] * log(etaVE_p);
                }
            }
            // HARMONIC AVERAGE
            if (average == 1) {
                eta_exx      = 1.0/eta_exx;
                eta_ezz      = 1.0/eta_ezz;
                eta_exz      = 1.0/eta_exz;
                eta_p        = 1.0/eta_p;
            }
            // GEOMETRIC AVERAGE
            if (average == 2) {
                eta_exx   = exp(eta_exx);
                eta_ezz   = exp(eta_ezz);
                eta_exz   = exp(eta_exz);
                eta_p     = exp(eta_p);
            }
            
            // General FD
            mesh->detadexx_s[c1]      = (eta_exx - mesh->eta_s[c1]) / pert_xx;
            mesh->detadezz_s[c1]      = (eta_ezz - mesh->eta_s[c1]) / pert_zz;
            mesh->detadgxz_s[c1]      = (eta_exz - mesh->eta_s[c1]) / pert_xz / 2.0;
            mesh->detadp_s[c1]        = (eta_p   - mesh->eta_s[c1]) / pert_p;
            
        }
    }
    
    MinMaxArrayTag( mesh->detadexx_n, scaling->eta/scaling->E, Ncx*Ncz, "detadexx_n", mesh->BCp.type );
    MinMaxArrayTag( mesh->detadexx_s, scaling->eta/scaling->E, Nx*Nz,   "detadexx_s", mesh->BCg.type );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
