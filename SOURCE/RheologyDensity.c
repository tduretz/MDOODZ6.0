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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
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

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InitialiseGrainSizeParticles( markers* particles, mat_prop *materials ){
    
    int phase;
#pragma omp parallel for shared( particles, materials ) private(phase)
    for( int k=0; k<particles->Nb_part; k++ ) {
        
        phase           = particles->phase[k];
        if ( phase != -1) particles->d[k] = materials->gs_ref[phase];
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Old deviatoric stress field and pressure
void  OldDeviatoricStressesPressure( grid* mesh, markers* particles, scale scaling, params* model ) {

    int k, l, c0, c1, c2, Nx, Nz, Ncx, Ncz, k1;
    double *sxx0,  *syy0, *szz0, *sxz0;
    int    cent=1, vert=0, prop=1, interp=0;

    Nx  = mesh->Nx;
    Nz  = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;

    sxx0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    syy0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    szz0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    sxz0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));

    P2Mastah( model, *particles, particles->sxxd,    mesh, sxx0,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
    P2Mastah( model, *particles, particles->szzd,    mesh, szz0,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
    P2Mastah( model, *particles, particles->syy,     mesh, syy0,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
    P2Mastah( model, *particles, particles->sxz,     mesh, mesh->sxz0,   mesh->BCg.type,  1, 0, interp, vert, model->itp_stencil);

#pragma omp parallel for shared( mesh, sxx0, syy0, szz0 ) private( k, k1, l, c0, c1, c2 ) firstprivate( Nx, Ncx, Ncz )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {
        k  = mesh->kp[k1];
        l  = mesh->lp[k1];
        c0 = k  + l*(Nx-1);
        c1 = k  + l*(Nx);
        c2 = k  + l*(Nx+1);

        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
            mesh->p0_n[c0]       = -1.0/3.0*(sxx0[c0] + syy0[c0] + szz0[c0] );
            mesh->sxxd0[c0]      =  mesh->p0_n[c0] + sxx0[c0];
            mesh->szzd0[c0]      =  mesh->p0_n[c0] + szz0[c0];
        }
    }

    InterpCentroidsToVerticesDouble( mesh->sxxd0,  mesh->sxxd0_s, mesh, model );
    InterpCentroidsToVerticesDouble( mesh->szzd0,  mesh->szzd0_s, mesh, model );
    InterpCentroidsToVerticesDouble( mesh->p0_n,   mesh->p0_s,    mesh, model );
    InterpVerticesToCentroidsDouble( mesh->sxz0_n, mesh->sxz0,    mesh, model );

    DoodzFree( sxx0 );
    DoodzFree( szz0 );
    DoodzFree( sxz0 );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Total stress field
void TotalStresses( grid* mesh, markers* particles, scale scaling, params* model ) {

    int k, l, c0, c1, c2, Nx, Nz, Ncx, Ncz, k1;
    double *dsxx, *dsyy, *dszz, *dsxz, sxx, syy, szz, tyy, tyy0, sxx0, syy0, szz0;
    double *mdsxx, *mdszz, *mdsyy, *mdsxz;

    Nx  = mesh->Nx;
    Nz  = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;

    dsxx = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    dsyy = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    dszz = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    dsxz = DoodzCalloc(Nx *Nz , sizeof(DoodzFP));

    mdsxx = DoodzCalloc(particles->Nb_part,sizeof(DoodzFP));
    mdsyy = DoodzCalloc(particles->Nb_part,sizeof(DoodzFP));
    mdszz = DoodzCalloc(particles->Nb_part,sizeof(DoodzFP));
    mdsxz = DoodzCalloc(particles->Nb_part,sizeof(DoodzFP));

#pragma omp parallel for shared( mesh, dsxx, dsyy, dszz ) private( k, k1, l, c0, c1, c2, sxx, syy, szz, tyy, tyy0, sxx0, syy0, szz0  ) firstprivate( Nx, Ncx, Ncz )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {
        k  = mesh->kp[k1];
        l  = mesh->lp[k1];
        c0 = k  + l*(Nx-1);
        c1 = k  + l*(Nx);
        c2 = k  + l*(Nx+1);

        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
            tyy      = -(mesh->sxxd[c0]  + mesh->szzd[c0]);
            tyy0     = -(mesh->sxxd0[c0] + mesh->szzd0[c0]);
            sxx0     = -mesh->p0_n[c0] + mesh->sxxd0[c0];
            szz0     = -mesh->p0_n[c0] + mesh->szzd0[c0];
            syy0     = -mesh->p0_n[c0] + tyy0;
            sxx      = -mesh->p_in[c0] + mesh->sxxd[c0];
            szz      = -mesh->p_in[c0] + mesh->szzd[c0];
            syy      = -mesh->p_in[c0] + tyy;
            dsxx[c0] = sxx - sxx0;
            dsyy[c0] = syy - syy0;
            dszz[c0] = szz - szz0;
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
    Interp_Grid2P_centroids2( *particles, mdsxx,  mesh, dsxx, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
    Interp_Grid2P_centroids2( *particles, mdszz,  mesh, dszz, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
    Interp_Grid2P_centroids2( *particles, mdsyy,  mesh, dsyy, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
    Interp_Grid2P( *particles, mdsxz,  mesh, dsxz, mesh->xg_coord,  mesh->zg_coord,  mesh->Nx,   mesh->Nz,   mesh->BCg.type  );

    // Update marker stresses
    ArrayPlusArray( particles->sxxd, mdsxx, particles->Nb_part );
    ArrayPlusArray( particles->szzd, mdszz, particles->Nb_part );
    ArrayPlusArray( particles->syy , mdsyy, particles->Nb_part );
    ArrayPlusArray( particles->sxz,  mdsxz, particles->Nb_part );

    // Freedom
    DoodzFree( dsxx );
    DoodzFree( dsyy );
    DoodzFree( dszz );
    DoodzFree( dsxz );
    DoodzFree(mdsxx);
    DoodzFree(mdszz);
    DoodzFree(mdsyy);
    DoodzFree(mdsxz);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void RheologicalOperators( grid* mesh, params* model, scale* scaling, int Jacobian ) {

    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    double nx, nz, deta, d0, d1;
    int aniso_fstrain = model->aniso_fstrain;
    double etae, K, dt = model->dt;
    int el = model->iselastic;

    if (Jacobian == 0  && model->aniso == 0) {

        // Loop on cell centers
#pragma omp parallel for shared( mesh )
        for (k=0; k<Ncx*Ncz; k++) {

            if ( mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) {
                mesh->D11_n[k] = 2.0*mesh->eta_n[k];
                mesh->D22_n[k] = 2.0*mesh->eta_n[k];
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

        // Loop on cell centers
#pragma omp parallel for shared( mesh ) private ( etae, K ) firstprivate ( el, dt )
        for (k=0; k<Ncx*Ncz; k++) {

            if ( mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) {
                switch ( el ) {
                    case 0:
                        etae = model->dt*mesh->mu_n[k];
                        mesh->D11_n[k] = 2.0*mesh->eta_n[k] + 2.0*mesh->detadexx_n[k]*mesh->exxd[k];
                        mesh->D12_n[k] =                      2.0*mesh->detadezz_n[k]*mesh->exxd[k];
                        mesh->D13_n[k] =                      2.0*mesh->detadgxz_n[k]*mesh->exxd[k];
                        mesh->D14_n[k] =                      2.0*mesh->detadp_n[k]  *mesh->exxd[k];

                        mesh->D21_n[k] =                      2.0*mesh->detadexx_n[k]*mesh->ezzd[k];
                        mesh->D22_n[k] = 2.0*mesh->eta_n[k] + 2.0*mesh->detadezz_n[k]*mesh->ezzd[k];
                        mesh->D23_n[k] =                      2.0*mesh->detadgxz_n[k]*mesh->ezzd[k];
                        mesh->D24_n[k] =                      2.0*mesh->detadp_n[k]  *mesh->ezzd[k];

                        break;
                    case 1:
                        etae = model->dt*mesh->mu_n[k];
                        K    = 1.0/mesh->bet_n[k];
                        //                        mesh->D11_n[k] = 2.0*mesh->eta_n[k] + 2.0*mesh->detadexx_n[k]*mesh->exxd[k] + mesh->sxxd0[k]*mesh->detadexx_n[k] / etae - K*dt*(2.0*mesh->ddivpdexx_n[k]+mesh->ddivpdezz_n[k]);
                        //                        mesh->D12_n[k] =                      2.0*mesh->detadezz_n[k]*mesh->exxd[k] + mesh->sxxd0[k]*mesh->detadezz_n[k] / etae - K*dt*(2.0*mesh->ddivpdezz_n[k]+mesh->ddivpdexx_n[k]);
                        mesh->D11_n[k] = 2.0*mesh->eta_n[k] + 2.0*mesh->detadexx_n[k] * ( mesh->exxd[k] + mesh->sxxd0[k]/etae/2.0 ) - K*dt*mesh->ddivpdexx_n[k];
                        mesh->D12_n[k] =                      2.0*mesh->detadezz_n[k] * ( mesh->exxd[k] + mesh->sxxd0[k]/etae/2.0 ) - K*dt*mesh->ddivpdezz_n[k];
                        mesh->D13_n[k] =                      2.0*mesh->detadgxz_n[k] * ( mesh->exxd[k] + mesh->sxxd0[k]/etae/2.0 ) - K*dt*mesh->ddivpdgxz_n[k];
                        mesh->D14_n[k] =                      2.0*mesh->detadp_n[k]   * ( mesh->exxd[k] + mesh->sxxd0[k]/etae/2.0 ) - K*dt*mesh->ddivpdp_n[k];

                        //                        mesh->D21_n[k] =                      2.0*mesh->detadexx_n[k]*mesh->ezzd[k] + mesh->szzd0[k]*mesh->detadexx_n[k] / etae - K*dt*(2.0*mesh->ddivpdexx_n[k]+mesh->ddivpdezz_n[k]);
                        //                        mesh->D22_n[k] = 2.0*mesh->eta_n[k] + 2.0*mesh->detadezz_n[k]*mesh->ezzd[k] + mesh->szzd0[k]*mesh->detadezz_n[k] / etae - K*dt*(2.0*mesh->ddivpdezz_n[k]+mesh->ddivpdexx_n[k]);
                        mesh->D21_n[k] =                      2.0*mesh->detadexx_n[k] * ( mesh->ezzd[k] + mesh->szzd0[k]/etae/2.0 ) - K*dt*mesh->ddivpdexx_n[k];
                        mesh->D22_n[k] = 2.0*mesh->eta_n[k] + 2.0*mesh->detadezz_n[k] * ( mesh->ezzd[k] + mesh->szzd0[k]/etae/2.0 ) - K*dt*mesh->ddivpdezz_n[k];
                        mesh->D23_n[k] =                      2.0*mesh->detadgxz_n[k] * ( mesh->ezzd[k] + mesh->szzd0[k]/etae/2.0 ) - K*dt*mesh->ddivpdgxz_n[k];
                        mesh->D24_n[k] =                      2.0*mesh->detadp_n[k]   * ( mesh->ezzd[k] + mesh->szzd0[k]/etae/2.0 ) - K*dt*mesh->ddivpdp_n[k];
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
                        mesh->D31_s[k] =                  mesh->detadexx_s[k]*2.0*mesh->exz[k];  // Factor 2 is important!!
                        mesh->D32_s[k] =                  mesh->detadezz_s[k]*2.0*mesh->exz[k];
                        mesh->D33_s[k] = mesh->eta_s[k] + mesh->detadgxz_s[k]*2.0*mesh->exz[k];
                        mesh->D34_s[k] =                  mesh->detadp_s[k]  *2.0*mesh->exz[k];
                        break;
                    case 1:
                        etae = model->dt*mesh->mu_s[k];
                        mesh->D31_s[k] =                  mesh->detadexx_s[k]*2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadexx_s[k] / etae;  // Factor 2 is important!!
                        mesh->D32_s[k] =                  mesh->detadezz_s[k]*2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadezz_s[k] / etae;
                        mesh->D33_s[k] = mesh->eta_s[k] + mesh->detadgxz_s[k]*2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadgxz_s[k] / etae;
                        mesh->D34_s[k] =                  mesh->detadp_s[k]  *2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadp_s[k]   / etae;
                        break;
                }

            }
            else {
                mesh->D31_s[k] = 0.0;
                mesh->D32_s[k] = 0.0;
                mesh->D33_s[k] = 0.0;
                mesh->D34_s[k] = 0.0;
            }
            if (isnan(mesh->D34_s[k])) { printf("EXIT: D34 is NAN!\n"); exit(1); }
            if (isinf(mesh->D34_s[k])) { printf("EXIT: D34 is INF!\n"); exit(1); }

        }

    }

    if ( Jacobian == 0 && model->aniso == 1 ) {

        printf("Computing anisotropic viscosity tensor\n");

        // Loop on cell centers
#pragma omp parallel for shared( mesh ) private ( nx, nz, deta, d0, d1, etae )  firstprivate ( aniso_fstrain, el )
        for (k=0; k<Ncx*Ncz; k++) {

            // Director
            nx = mesh->nx0_n[k];
            nz = mesh->nz0_n[k];

            if ( mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) {
                // See Anisotropy_v2.ipynb
                if ( aniso_fstrain  == 0 ) deta =  (mesh->eta_n[k] - mesh->eta_n[k] / mesh->aniso_factor_n[k]);
                if ( aniso_fstrain  == 1 ) deta =  (mesh->eta_n[k] - mesh->eta_n[k] / mesh->FS_AR_n[k]);
                d0   = 2.0*pow(nx, 2.0)*pow(nz, 2.0);
                d1   = nx*nz*(-pow(nx, 2.0) + pow(nz, 2.0));
                //                printf("deta = %2.2e d1 = %2.2e\n",deta, d1);
                mesh->D11_n[k] = 2.0*mesh->eta_n[k] - 2.0*deta*d0;
                mesh->D12_n[k] =                      2.0*deta*d0;
                mesh->D13_n[k] =                      2.0*deta*d1;
                mesh->D14_n[k] =                      0.0;

                mesh->D21_n[k] =                      2.0*deta*d0;
                mesh->D22_n[k] = 2.0*mesh->eta_n[k] - 2.0*deta*d0;
                mesh->D23_n[k] =                     -2.0*deta*d1;
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

            //            if ( mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) {
            //                switch ( el ) {
            //                    case 0:
            //                        etae = model->dt*mesh->mu_n[k];
            //                        mesh->D11_n[k] += 2.0*mesh->eta_n[k] + yn*2.0*mesh->detadexx_n[k]*mesh->exxd[k];
            //                        mesh->D12_n[k] +=                      yn*2.0*mesh->detadezz_n[k]*mesh->exxd[k];
            //                        mesh->D13_n[k] +=                      yn*2.0*mesh->detadgxz_n[k]*mesh->exxd[k];
            //                        mesh->D14_n[k] +=                      yn*2.0*mesh->detadp_n[k]  *mesh->exxd[k];
            //
            //                        mesh->D21_n[k] +=                      yn*2.0*mesh->detadexx_n[k]*mesh->ezzd[k];
            //                        mesh->D22_n[k] += 2.0*mesh->eta_n[k] + yn*2.0*mesh->detadezz_n[k]*mesh->ezzd[k];
            //                        mesh->D23_n[k] +=                      yn*2.0*mesh->detadgxz_n[k]*mesh->ezzd[k];
            //                        mesh->D24_n[k] +=                      yn*2.0*mesh->detadp_n[k]  *mesh->ezzd[k];
            //
            //                        break;
            //                    case 1:
            //                        etae = model->dt*mesh->mu_n[k];
            //                        mesh->D11_n[k] += 2.0*mesh->eta_n[k] + yn*2.0*mesh->detadexx_n[k]*mesh->exxd[k] + mesh->sxxd0[k]*mesh->detadexx_n[k] / etae;
            //                        mesh->D12_n[k] +=                      yn*2.0*mesh->detadezz_n[k]*mesh->exxd[k] + mesh->sxxd0[k]*mesh->detadezz_n[k] / etae;
            //                        mesh->D13_n[k] +=                      yn*2.0*mesh->detadgxz_n[k]*mesh->exxd[k] + mesh->sxxd0[k]*mesh->detadgxz_n[k] / etae;
            //                        mesh->D14_n[k] +=                      yn*2.0*mesh->detadp_n[k]  *mesh->exxd[k] + mesh->sxxd0[k]*mesh->detadp_n[k]   / etae;
            //
            //                        mesh->D21_n[k] +=                      yn*2.0*mesh->detadexx_n[k]*mesh->ezzd[k] + mesh->szzd0[k]*mesh->detadexx_n[k] / etae;
            //                        mesh->D22_n[k] += 2.0*mesh->eta_n[k] + yn*2.0*mesh->detadezz_n[k]*mesh->ezzd[k] + mesh->szzd0[k]*mesh->detadezz_n[k] / etae;
            //                        mesh->D23_n[k] +=                      yn*2.0*mesh->detadgxz_n[k]*mesh->ezzd[k] + mesh->szzd0[k]*mesh->detadgxz_n[k] / etae;
            //                        mesh->D24_n[k] +=                      yn*2.0*mesh->detadp_n[k]  *mesh->ezzd[k] + mesh->szzd0[k]*mesh->detadp_n[k]   / etae;
            //                        break;
            //                }
            //            }
        }

        //        int k, l;

        // Loop on cell vertices
#pragma omp parallel for shared( mesh )  private ( nx, nz, deta, d0, d1, etae ) firstprivate ( aniso_fstrain, el )
        for (k=0; k<Nx*Nz; k++) {

            // Director
            nx = mesh->nx0_s[k];
            nz = mesh->nz0_s[k];

            if ( mesh->BCg.type[k] != 30 ) {
                // See Anisotropy_v2.ipynb
                if ( aniso_fstrain  == 0 ) deta =  (mesh->eta_s[k] - mesh->eta_s[k] / mesh->aniso_factor_s[k]);
                if ( aniso_fstrain  == 1 ) deta =  (mesh->eta_s[k] - mesh->eta_s[k] / mesh->FS_AR_s[k]);
                d0   =  2.0*pow(nx, 2.0)*pow(nz, 2.0);
                d1   = nx*nz*(-pow(nx, 2.0) + pow(nz, 2.0));

                mesh->D31_s[k] =                  2.0*deta*d1;
                mesh->D32_s[k] =                 -2.0*deta*d1;
                mesh->D33_s[k] = mesh->eta_s[k] + 2.0*deta*(d0 - 0.5);
                mesh->D34_s[k] =                  0.0;
            }
            else {
                mesh->D31_s[k] = 0.0;
                mesh->D32_s[k] = 0.0;
                mesh->D33_s[k] = 0.0;
                mesh->D34_s[k] = 0.0;
            }

            //            if ( mesh->BCg.type[k] != 30 ) {
            //                switch ( el ) {
            //                    case 0:
            //                        mesh->D31_s[k] +=                  yn*mesh->detadexx_s[k]*2.0*mesh->exz[k];  // Factor 2 is important!!
            //                        mesh->D32_s[k] +=                  yn*mesh->detadezz_s[k]*2.0*mesh->exz[k];
            //                        mesh->D33_s[k] += mesh->eta_s[k] + yn*mesh->detadgxz_s[k]*2.0*mesh->exz[k];
            //                        mesh->D34_s[k] +=                  yn*mesh->detadp_s[k]  *2.0*mesh->exz[k];
            //                        break;
            //                    case 1:
            //                        etae = model->dt*mesh->mu_s[k];
            //                        mesh->D31_s[k] +=                  yn*mesh->detadexx_s[k]*2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadexx_s[k] / etae;  // Factor 2 is important!!
            //                        mesh->D32_s[k] +=                  yn*mesh->detadezz_s[k]*2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadezz_s[k] / etae;
            //                        mesh->D33_s[k] += mesh->eta_s[k] + yn*mesh->detadgxz_s[k]*2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadgxz_s[k] / etae;
            //                        mesh->D34_s[k] +=                  yn*mesh->detadp_s[k]  *2.0*mesh->exz[k] + mesh->sxz0[k]*mesh->detadp_s[k]   / etae;
            //                        break;
            //                }
            //
            //            }


        }
    }
    if (Jacobian == 1  && model->aniso == 1) {
        printf("\nShould be here some time !!!!!!\n");
    }

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
    Interp_Grid2P_centroids2( *(particles), pdudx, &mesh, dudx, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model  );
    Interp_Grid2P_centroids2( *(particles), pdvdz, &mesh, dvdz, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model );
    Interp_Grid2P( *(particles), pdvdx, &mesh, dvdx, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
    Interp_Grid2P( *(particles), pdudz, &mesh, dudz, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );

#pragma omp parallel for shared ( particles ) private ( k, fxx, fxz, fzx, fzz, fxx_o, fxz_o, fzx_o, fzz_o ) firstprivate( model ) schedule( static )
    for(k=0; k<particles->Nb_part; k++) {

        fxx = 1.0 + pdudx[k]*model.dt;
        fxz = pdudz[k]*model.dt;
        fzz = 1.0 + pdvdz[k]*model.dt;
        fzx = pdvdx[k]*model.dt;
        fxx_o = particles->Fxx[k];
        fxz_o = particles->Fxz[k];
        fzx_o = particles->Fzx[k];
        fzz_o = particles->Fzz[k];
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

    printf("-----> Deformation gradient tensor Updated\n");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FiniteStrainAspectRatio ( grid *mesh, scale scaling, params model, markers *particles ) {

    int k, l, cp, cu, cv;
    double *FS_AR, CGxx, CGxz, CGzx, CGzz;
    double Tr, Det, U0xx, U0xz, U0zx, U0zz, s, t, e1, e2;
    int    cent=1, vert=0, prop=1, interp=0;

    FS_AR = DoodzCalloc (particles->Nb_part, sizeof(double));

    //#pragma omp parallel for shared ( particles ) private ( k, CGxx, CGxz, CGzx, CGzz, Tr, Det, U0xx, U0xz, U0zx, U0zz, s, t, e1, e2  ) schedule( static )
    for (k=0; k<particles->Nb_part; k++) {

        // Cauchy-Green
        CGxx = particles->Fxx[k]*particles->Fxx[k] + particles->Fzx[k]*particles->Fzx[k];
        CGxz = particles->Fxx[k]*particles->Fxz[k] + particles->Fzx[k]*particles->Fzz[k];
        CGzx = particles->Fxx[k]*particles->Fxz[k] + particles->Fzx[k]*particles->Fzz[k];
        CGzz = particles->Fxz[k]*particles->Fxz[k] + particles->Fzz[k]*particles->Fzz[k];

        // U0 = Sqrt(CG)
        Tr  = CGxx + CGzz;
        Det = CGxx*CGzz - CGxz*CGzx;
        s   = sqrt(Det);
        t   = sqrt(Tr + 2.0*s);
        U0xx = 1/t*(CGxx + s);
        U0xz = 1/t*CGxz;
        U0zx = 1/t*CGzx;
        U0zz = 1/t*(CGzz + s);

        // eigs(U0)
        Tr   = U0xx + U0zz;
        Det  = U0xx*U0zz - U0xz*U0zx;
        e1   = Tr/2.0 + sqrt( Tr*Tr/4.0 - Det );
        e2   = Tr/2.0 - sqrt( Tr*Tr/4.0 - Det );

        // aspect ratio
        FS_AR[k] = e1/e2;
    }

    P2Mastah( &model, *particles, FS_AR, mesh, mesh->FS_AR_n,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
    P2Mastah( &model, *particles, FS_AR, mesh, mesh->FS_AR_s,   mesh->BCg.type,  1, 0, interp, vert, model.itp_stencil);

    DoodzFree(FS_AR);

    printf("-----> Finite strain aspect ratio updated\n");
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

    printf("-----> Accumulated strain updated\n");
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

void UpdateParticleDensity( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

    DoodzFP *rho_inc_mark, *rho_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    rho_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    rho_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));

//    Interp_Grid2P_centroids2( *particles, particles->rho, mesh, mesh->rho_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );


    for (k=0;k<Ncx*Ncz;k++) {
        rho_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) rho_inc_grid[k] = mesh->rho_n[k] - mesh->rho0_n[k];
    }

    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, rho_inc_mark, mesh, rho_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    // Increment temperature on particles
    ArrayPlusArray( particles->rho, rho_inc_mark, particles->Nb_part );

    DoodzFree(rho_inc_grid);
    DoodzFree(rho_inc_mark);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticlePhi( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

    DoodzFP *phi_inc_mark, *phi_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    phi_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    phi_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));

    for (k=0;k<Ncx*Ncz;k++) {
        phi_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) phi_inc_grid[k] = mesh->phi_n[k] - mesh->phi0_n[k];
    }

    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, phi_inc_mark, mesh, phi_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    // Increment temperature on particles
    ArrayPlusArray( particles->phi, phi_inc_mark, particles->Nb_part );

    // Bounds: numerical capote
#pragma omp parallel for shared ( particles ) private( k )
    for (k=0; k<particles->Nb_part; k++) {
        if (particles->phi[k] <= 0.0) particles->phi[k] = 0.0;
        if (particles->phi[k] >= 1.0) particles->phi[k] = 1.0;
    }

    DoodzFree(phi_inc_grid);
    DoodzFree(phi_inc_mark);

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleX( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

    DoodzFP *X_inc_mark, *X_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    // // Total update (do not use!!!)
    //    Interp_Grid2P( *particles, particles->X, mesh, mesh->X_n, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    // Incremental update
    X_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    X_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
#pragma omp parallel for shared ( X_inc_grid, mesh ) private( k ) firstprivate( Ncx, Ncz )
    for (k=0; k<Ncx*Ncz; k++) {
        X_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) X_inc_grid[k] = mesh->X_n[k] - mesh->X0_n[k];
    }

    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, X_inc_mark, mesh, X_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );

    // Increment temperature on particles
    ArrayPlusArray( particles->X, X_inc_mark, particles->Nb_part );

    // Bounds: numerical capote
#pragma omp parallel for shared ( particles ) private( k )
    for (k=0; k<particles->Nb_part; k++) {
        if (particles->X[k] <= 0.0) particles->X[k] = 0.0;
        if (particles->X[k] >= 1.0) particles->X[k] = 1.0;
    }

    DoodzFree(X_inc_grid);
    DoodzFree(X_inc_mark);
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
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) d_inc_grid[k] = mesh->d_n[k] - mesh->d0_n[k];
    }

    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, d_inc_mark, mesh, d_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    // Increment grain size on particles
    ArrayPlusArray( particles->d, d_inc_mark, particles->Nb_part );

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

    DoodzFP *T_inc_mark, *Tm0, dtm, *dTms, *dTgr, *dTmr, *rho_part;
    double *Tg0, *dTgs, dx=model.dx, dz=model.dz, d=1.0;
    int    cent=1, vert=0, prop=1, interp=0;
    int Nx, Nz, Ncx, Ncz, k, c0, p;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    // Allocations
    Tm0        = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    Tg0        = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    T_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

    // Old temperature grid
#pragma omp parallel for shared(mesh, Tg0) private(c0) firstprivate(Ncx,Ncz)
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {
        if (mesh->BCt.type[c0] != 30) Tg0[c0] = mesh->T[c0] - mesh->dT[c0];
    }
    Interp_Grid2P_centroids2( *particles, Tm0, mesh, Tg0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    // SUBGRID
    if ( model.subgrid_diff >= 1 ) { /* CASE WITH SUBGRID DIFFUSION */

        printf("Subgrid diffusion for temperature update\n");
        dTgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dTgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dTms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        dTmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        rho_part = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

        // Get density
        Interp_Grid2P_centroids2( *particles, rho_part, mesh, mesh->rho_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

        // Compute subgrid temperature increments on markers
#pragma omp parallel for shared(particles,Tm0,dTms) private(k,p,dtm) firstprivate(materials,dx,dz,model,d)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) {
                p = particles->phase[k];
                dtm     = materials->Cv[p] * rho_part[k]/ (materials->k[p] * (1.0/dx/dx + 1.0/dz/dz));
                dTms[k] = -( particles->T[k] - Tm0[k] ) * (1.0-exp(-d*model.dt/dtm));
            }
        }

        // Subgrid temperature increments markers --> grid
        P2Mastah( &model, *particles, dTms,     mesh, dTgs,   mesh->BCp.type,  1, 0, interp, cent, 1);

        // Remaining temperature increments on the grid
#pragma omp parallel for shared(mesh, dTgs, dTgr) private(c0) firstprivate(Ncx,Ncz)
        for ( c0=0; c0<Ncx*Ncz; c0++ ) {
            if (mesh->BCt.type[c0] != 30) dTgr[c0] = mesh->dT[c0] - dTgs[c0];
        }

        // Remaining temperature increments grid --> markers
        Interp_Grid2P_centroids2( *particles, dTmr, mesh, dTgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

        // Final temperature update on markers
#pragma omp parallel for shared(particles,dTms,dTmr,T_inc_mark) private(k)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) T_inc_mark[k]    = dTms[k] + dTmr[k];
            if (particles->phase[k] != -1) particles->T[k]  = particles->T[k] + T_inc_mark[k];
        }
        DoodzFree(dTms);
        DoodzFree(dTmr);
        DoodzFree(dTgs);
        DoodzFree(dTgr);
        DoodzFree(rho_part);
    }
    else {  /* CASE WITHOUT SUBGRID DIFFUSION: INTERPOLATE INCREMENT DIRECTLY */

        // Interp increments to particles
        Interp_Grid2P_centroids2( *particles, T_inc_mark, mesh, mesh->dT, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

        // Increment temperature on particles
        ArrayPlusArray( particles->T, T_inc_mark, particles->Nb_part );
    }

    // Freedom
    DoodzFree(Tg0);
    DoodzFree(Tm0);
    DoodzFree(T_inc_mark);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticlePressure( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

    DoodzFP *P_inc_mark;
    int Nx, Nz, Ncx, Ncz, k, c0, p, ptrick=model.Plith_trick;
    double d=1.0, dtm;
    int    cent=1, vert=0, prop=1, interp=0;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;

    // Compute increment
#pragma omp parallel for shared(mesh) private(c0) firstprivate( ptrick )
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {
        mesh->dp[c0] = 0.0;
        if (mesh->BCp.type[c0] != 30 ) {
            if ( ptrick == 1 ) mesh->dp[c0] = (mesh->p_in[c0] - mesh->p_lith[c0]) - (mesh->p0_n[c0] - mesh->p_lith0[c0]); // dp dynamic pressure
            if ( ptrick == 0 ) mesh->dp[c0] = (mesh->p_in[c0]-mesh->p0_n[c0]);
        }
    }


//    double *p_s = DoodzCalloc(Nx*Nz, sizeof(DoodzFP));
//    double *dp = DoodzCalloc(Nx*Nz, sizeof(DoodzFP));
//    InterpCentroidsToVerticesDouble( mesh->p_in, p_s,mesh, &model );
//
//    // Compute increment
//    #pragma omp parallel for shared(mesh) private(c0) firstprivate( ptrick )
//        for ( c0=0; c0<Nx*Nz; c0++ ) {
//            dp[c0] = 0.0;
//            if (mesh->BCg.type[c0] != 30 ) {
//                dp[c0] = (p_s[c0]-mesh->p0_s[c0]);
//            }
//        }

    if ( model.subgrid_diff >= 2 ) {

        printf("Subgrid diffusion for pressure update\n");
        double *Pg0  = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        double *dPgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        double *dPgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        double *Pm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        double *dPms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        double *dPmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        double *etam = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

        Interp_Grid2P_centroids2( *particles, etam,  mesh, mesh->eta_phys_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1  , Nz-1 , mesh->BCp.type, &model);

        /* -------------- */
        // Old Pressure grid
#pragma omp parallel for shared(mesh, Pg0) private(c0) firstprivate(Ncx,Ncz) firstprivate( ptrick )
        for ( c0=0; c0<Ncx*Ncz; c0++ ) {
            if (mesh->BCt.type[c0] != 30 && ptrick == 1 ) Pg0[c0] = mesh->p0_n[c0] - mesh->p_lith0[c0];
            if (mesh->BCt.type[c0] != 30 && ptrick == 0 ) Pg0[c0] = mesh->p0_n[c0];
        }
        Interp_Grid2P_centroids2( *particles, Pm0, mesh, Pg0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );
        /* -------------- */

//        Interp_Grid2P_centroids( *particles, Pm0, mesh, mesh->p0_n, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );



        // Old temperature grid --> markers
        //        Interp_Grid2P_centroids( *particles, Pm0, mesh, mesh->p0_n, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );

        // Compute subgrid temperature increments on markers
#pragma omp parallel for shared(particles,Pm0,dPms) private(k,p,dtm) firstprivate(materials,model,d)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) {
                p         = particles->phase[k];
                dtm       = etam[k]/ (materials->bet[p]);
                dPms[k]   = -( particles->P[k] - Pm0[k]) * (1.0 - exp(-d*model.dt/dtm));
            }
        }

        // Subgrid temperature increments markers --> grid
        P2Mastah( &model, *particles, dPms,     mesh, dPgs,   mesh->BCp.type,  1, 0, interp, cent, 1);

        // Remaining temperature increments on the grid
#pragma omp parallel for shared(mesh, dPgs, dPgr) private(c0) firstprivate(Ncx,Ncz)
        for ( c0=0; c0<Ncx*Ncz; c0++ ) {
            if (mesh->BCp.type[c0] != 30 ) dPgr[c0] = mesh->dp[c0] - dPgs[c0];
        }

        // Remaining temperature increments grid --> markers
        Interp_Grid2P_centroids2( *particles, dPmr, mesh, dPgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );

        // Final temperature update on markers
#pragma omp parallel for shared(particles,dPms,dPmr) private(k)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) particles->P[k] = particles->P[k] + dPms[k] + dPmr[k];
        }

        DoodzFree(Pg0);
        DoodzFree(Pm0);
        DoodzFree(dPms);
        DoodzFree(dPmr);
        DoodzFree(dPgs);
        DoodzFree(dPgr);
        DoodzFree(etam);
    }
    else {

        P_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

        // Interp increments to particles
        Interp_Grid2P_centroids2( *particles, P_inc_mark, mesh, mesh->dp, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );

        // Increment pressure on particles
        ArrayPlusArray( particles->P, P_inc_mark, particles->Nb_part );

        DoodzFree(P_inc_mark);

    }

}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleStress( grid* mesh, markers* particles, params* model, mat_prop* materials, scale *scaling ) {

    int k, l, c0, c1, c2, Nx, Nz, Ncx, Ncz, p, k1, c3;
    DoodzFP *mdsxxd, *mdszzd, *mdsxz, *dsxxd, *dszzd, *dsxz, d=1.0, dtaum;
    DoodzFP *dtxxgs, *dtzzgs, *dtxzgs, *dtxxgr, *dtzzgr, *dtxzgr, *txxm0, *tzzm0, *txzm0, *dtxxms, *dtzzms, *dtxzms, *dtxxmr, *dtzzmr, *dtxzmr, *etam;
    double *dudx_n, *dvdz_n, *dudz_s, *dvdx_s, *om_s, *om_n, *dudz_n, *dvdx_n, *dudx_s, *dvdz_s;
    double angle, tzz, txx, txz, dx, dz, dt;
    double *txz_n, *txx_s, *tzz_s, *dtxxg0, *dtzzg0, *dtxzg0;
    int style = model->StressUpdate;
    int    cent=1, vert=0, prop=1, interp=0;

    Nx = model->Nx;
    Nz = model->Nz;
    dx = model->dx;
    dz = model->dz;
    dt = model->dt;

    om_s   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    om_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    txz_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    txx_s   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    tzz_s   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));

    dudx_n = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    dvdz_n = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    dudz_s = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    dvdx_s = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    dudz_n = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    dvdx_n = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    dudx_s = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    dvdz_s = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));

#pragma omp parallel for shared ( mesh, om_s, dudz_s, dvdx_s ) \
private ( k, l, k1, c1, c3 )                                   \
firstprivate( model )
    for ( k1=0; k1<Nx*Nz; k1++ ) {
        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;
        c3 = k + l*(Nx+1);
        if ( mesh->BCg.type[c1] != 30 ) {
            om_s[c1]   = 0.5*( (mesh->v_in[c3+1] - mesh->v_in[c3])/dx - (mesh->u_in[c1+Nx] - mesh->u_in[c1])/dz);
            dudz_s[c1] = (mesh->u_in[c1+Nx] - mesh->u_in[c1])/dz;
            dvdx_s[c1] = (mesh->v_in[c3+1       ] - mesh->v_in[c3])/dx;
        }
    }

#pragma omp parallel for shared ( mesh, dudx_n, dvdz_n ) \
private ( k, l, k1, c0, c1, c2 )                         \
firstprivate( model )
    for ( k1=0; k1<(Nx-1)*(Nz-1); k1++ ) {
        k  = mesh->kp[k1];
        l  = mesh->lp[k1];
        c0 = k  + l*(Nx-1);
        c1 = k  + l*(Nx);
        c2 = k  + l*(Nx+1);
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
            dudx_n[c0]  = (mesh->u_in[c1+1+Nx]     - mesh->u_in[c1+Nx] )/dx;
            dvdz_n[c0]  = (mesh->v_in[c2+1+(Nx+1)] - mesh->v_in[c2+1]        )/dz;
        }
    }

    InterpCentroidsToVerticesDouble( dudx_n, dudx_s, mesh, model );
    InterpCentroidsToVerticesDouble( dvdz_n, dvdz_s, mesh, model );
    InterpVerticesToCentroidsDouble( dudz_n, dudz_s, mesh, model );
    InterpVerticesToCentroidsDouble( dvdx_n, dvdx_s, mesh, model );
    InterpCentroidsToVerticesDouble( mesh->sxxd, txx_s, mesh, model );
    InterpCentroidsToVerticesDouble( mesh->szzd, tzz_s, mesh, model );
    InterpVerticesToCentroidsDouble( txz_n, mesh->sxz, mesh, model );
    InterpVerticesToCentroidsDouble( om_n, om_s, mesh, model );

    // Rotate stresses and director
    double nx, nz, ndotx, ndotz, *dnx_n, *dnz_n, *dnx_s, *dnz_s;

    dnx_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    dnz_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    dnx_s   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    dnz_s   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));

#pragma omp parallel for shared ( mesh, dudz_n, dvdx_n, dudx_n, dvdz_n, om_n ) \
private ( k1, txx, tzz, txz, angle, nx, nz, ndotx, ndotz )                     \
firstprivate( model, dt )
    for ( k1=0; k1<(Nx-1)*(Nz-1); k1++ ) {

        txx   = mesh->sxxd[k1];
        tzz   = mesh->szzd[k1];
        txz   = txz_n[k1];
        if (model->StressRotation==1) {
            angle = dt*om_n[k1];
            mesh->sxxd[k1] = (txx*cos(angle) - txz*sin(angle))*cos(angle) - (txz*cos(angle) - tzz*sin(angle))*sin(angle);
            mesh->szzd[k1] = (txx*sin(angle) + txz*cos(angle))*sin(angle) + (txz*sin(angle) + tzz*cos(angle))*cos(angle);
        }
        if (model->StressRotation==2) {
            mesh->sxxd[k1] = mesh->sxxd[k1] - dt * mesh->VE_n[k1] * ( -2.0*txx*dudx_n[k1] - 2.0*txz*dudz_n[k1]);
            mesh->szzd[k1] = mesh->szzd[k1] - dt * mesh->VE_n[k1] * ( -2.0*tzz*dvdz_n[k1] - 2.0*txz*dvdx_n[k1]);
        }

        if ( model->aniso == 1 ) { // Director vector rotation/deformation
            nx        = mesh->nx0_n[k1];// = 0.0;
            nz        = mesh->nz0_n[k1];// = 1.0;
            ndotx     = (-(dudx_n[k1] - dvdz_n[k1])*nx*nz - dvdx_n[k1]*nz*nz + dudz_n[k1]*nx*nx)*nz;
            ndotz     = ( (dudx_n[k1] - dvdz_n[k1])*nx*nz + dvdx_n[k1]*nz*nz - dudz_n[k1]*nx*nx)*nx;
            dnx_n[k1] = ndotx*dt;
            dnz_n[k1] = ndotz*dt;
        }

    }

#pragma omp parallel for shared ( mesh, dudz_s, dvdx_s, dudx_s, dvdz_s, om_s ) \
private ( k1, txx, tzz, txz, angle, nx, nz, ndotx, ndotz )                     \
firstprivate( model, dt )
    for ( k1=0; k1<(Nx-0)*(Nz-0); k1++ ) {

        txx   = txx_s[k1];
        tzz   = tzz_s[k1];
        txz   = mesh->sxz[k1];
        if (model->StressRotation==1) {
            angle = dt*om_s[k1];
            mesh->sxz[k1] = (txx*cos(angle) - txz*sin(angle))*sin(angle) + (txz*cos(angle) - tzz*sin(angle))*cos(angle);
        }
        if (model->StressRotation==2) {
            mesh->sxz[k1] = mesh->sxz[k1] - dt * mesh->VE_s[k1] * (      txx*dudz_s[k1] -     txx*dvdx_s[k1] - txz*(dudx_s[k1]+ dvdz_s[k1]) );
        }

        if ( model->aniso == 1 ) { // Director vector rotation/deformation
            nx        = mesh->nx0_s[k1];
            nz        = mesh->nz0_s[k1];
            ndotx     = (-(dudx_s[k1] - dvdz_s[k1])*nx*nz - dvdx_s[k1]*nz*nz + dudz_s[k1]*nx*nx)*nz;
            ndotz     = ( (dudx_s[k1] - dvdz_s[k1])*nx*nz + dvdx_s[k1]*nz*nz - dudz_s[k1]*nx*nx)*nx;
            dnx_s[k1] = ndotx*dt;
            dnz_s[k1] = ndotz*dt;
        }
    }

    // interpolate increments of director vectors to particles
    if ( model->aniso == 1 ) {
        //    if ( model->aniso == 1 ) Interp_Grid2P_centroids( *particles, particles->dnx, mesh, dnx_n, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
        //    if ( model->aniso == 1 ) Interp_Grid2P_centroids( *particles, particles->dnz, mesh, dnz_n, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
        Interp_Grid2P ( *particles, particles->dnx, mesh, dnx_s, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type         );
        Interp_Grid2P ( *particles, particles->dnz, mesh, dnz_s, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type         );
        ArrayPlusArray(  particles->nx, particles->dnx, particles->Nb_part );
        ArrayPlusArray(  particles->nz, particles->dnz, particles->Nb_part );
    }

    DoodzFree(txx_s);
    DoodzFree(tzz_s);
    DoodzFree(txz_n);
    DoodzFree(om_s);
    DoodzFree(om_n);
    DoodzFree(dudx_n);
    DoodzFree(dvdz_n);
    DoodzFree(dvdx_s);
    DoodzFree(dudz_s);
    DoodzFree(dudz_n);
    DoodzFree(dvdx_n);
    DoodzFree(dvdz_s);
    DoodzFree(dudx_s);
    DoodzFree(dnx_n);
    DoodzFree(dnz_n);
    DoodzFree(dnx_s);
    DoodzFree(dnz_s);
    //

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
        Interp_Grid2P_centroids2( *particles, txxm0, mesh, mesh->sxxd0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
        Interp_Grid2P_centroids2( *particles, tzzm0, mesh, mesh->szzd0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
        Interp_Grid2P( *particles, txzm0, mesh, mesh->sxz0 , mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type     );
        Interp_Grid2P( *particles, etam,  mesh, mesh->eta_phys_s, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type);

        MinMaxArray(mesh->eta_phys_s, scaling->eta, Nx*Nz, "eta grid  ");
        MinMaxArrayTag( mesh->eta_phys_s, scaling->eta, Nx*Nz,   "eta_s", mesh->BCg.type );
        MinMaxArray(etam, scaling->eta, particles->Nb_part, "eta phys part  ");


        if ( model->subgrid_diff == 2 ) {

            printf("Subgrid diffusion for stress tensor component update\n");

            // Compute subgrid stress increments on markers
#pragma omp parallel for shared(particles,txxm0,tzzm0,txzm0,dtxxms,dtzzms,dtxzms,etam) private(k,p,dtaum) firstprivate(materials,model,d)
            for ( k=0; k<particles->Nb_part; k++ ) {
                if (particles->phase[k] != -1) {
                    p         = particles->phase[k];
                    dtaum     = etam[k] / materials->mu[p];
                    dtxxms[k] = -( particles->sxxd[k] - txxm0[k]) * (1.0 - exp(-d*dt/dtaum));
                    dtzzms[k] = -( particles->szzd[k] - tzzm0[k]) * (1.0 - exp(-d*dt/dtaum));
                    dtxzms[k] = -( particles->sxz[k]  - txzm0[k]) * (1.0 - exp(-d*dt/dtaum));
                    if (isinf(dtxxms[k])) {
                       printf("Infinite dtxxms[k]: %2.2e %2.2e %2.2e\n", particles->sxxd[k], txxm0[k], exp(-d*dt/dtaum));
                        printf("%2.2e %2.2e %2.2e %2.2e %2.2e", d, dt, dtaum, etam[k]*scaling->eta, materials->mu[p]*scaling->S );
                       exit(1);
                   }
                    if (isnan(dtxxms[k])) {
                        printf("Infinite dtxxms[k]: %2.2e %2.2e %2.2e\n", particles->sxxd[k], txxm0[k], exp(-d*dt/dtaum));
                        exit(1);
                    }
                }
            }

            // Subgrid stress increments markers --> grid
            P2Mastah( model, *particles, dtxxms,     mesh, dtxxgs,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
            P2Mastah( model, *particles, dtzzms,     mesh, dtzzgs,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
            P2Mastah( model, *particles, dtxzms,     mesh, dtxzgs,   mesh->BCg.type,  1, 0, interp, vert, model->itp_stencil);

            // Remaining stress increments on the grid
#pragma omp parallel for shared(mesh,dtxxgs,dtxxgr,dtzzgs,dtzzgr) private(c0) firstprivate(Ncx,Ncz)
            for ( c0=0; c0<Ncx*Ncz; c0++ ) {
                if (mesh->BCp.type[c0]!=30 && mesh->BCp.type[c0]!=31) dtxxgr[c0] = (mesh->sxxd[c0] - mesh->sxxd0[c0]) - dtxxgs[c0];
                if (mesh->BCp.type[c0]!=30 && mesh->BCp.type[c0]!=31) dtzzgr[c0] = (mesh->szzd[c0] - mesh->szzd0[c0]) - dtzzgs[c0];
            }
#pragma omp parallel for shared(mesh,dtxzgs,dtxzgr) private(c0) firstprivate(Nx,Nz)
            for ( c0=0; c0<Nx*Nz; c0++ ) {
                if (mesh->BCg.type[c0]!=30) dtxzgr[c0] = (mesh->sxz[c0]-mesh->sxz0[c0]) - dtxzgs[c0];
            }
            
            // Remaining stress increments grid --> markers
            Interp_Grid2P_centroids2( *particles, dtxxmr, mesh, dtxxgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
            Interp_Grid2P_centroids2( *particles, dtzzmr, mesh, dtzzgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
            Interp_Grid2P(           *particles, dtxzmr, mesh, dtxzgr, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type         );

            // Final stresses update on markers
#pragma omp parallel for shared(particles,dtxxms,dtzzms,dtxzms,dtxxmr,dtzzmr,dtxzmr) private(k) firstprivate(style)
            for ( k=0; k<particles->Nb_part; k++ ) {
                if (particles->phase[k] != -1) particles->sxxd[k]  = particles->sxxd[k] + dtxxms[k] + dtxxmr[k];
                if (particles->phase[k] != -1) particles->szzd[k]  = particles->szzd[k] + dtzzms[k] + dtzzmr[k];
                if (particles->phase[k] != -1) particles->sxz[k]   = particles->sxz[k]  + dtxzms[k] + dtxzmr[k];
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
        Interp_Grid2P_centroids2( *particles, mdsxxd, mesh, dsxxd, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
        Interp_Grid2P_centroids2( *particles, mdszzd, mesh, dszzd, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
        Interp_Grid2P( *particles, mdsxz,  mesh, dsxz,  mesh->xg_coord,  mesh->zg_coord,  mesh->Nx,   mesh->Nz,   mesh->BCg.type  );

        // Update marker stresses
        ArrayPlusArray(  particles->sxxd, mdsxxd, particles->Nb_part );
        ArrayPlusArray(  particles->szzd, mdszzd, particles->Nb_part );
        ArrayPlusArray(  particles->sxz,  mdsxz,  particles->Nb_part );

        // Free
        DoodzFree(dsxxd);
        DoodzFree(dszzd);
        DoodzFree(dsxz);
        DoodzFree(mdsxxd);
        DoodzFree(mdszzd);
        DoodzFree(mdsxz);
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/




double Viscosity( int phase, double G, double T, double P, double d, double phi, double X0, double exx, double ezz, double exz, double Txx0, double Tzz0, double Txz0, mat_prop* materials, params *model, scale *scaling, double *txxn, double *tzzn, double *txzn, double* etaVE, double* VEcoeff, double* Eii_el, double* Eii_pl, double* Eii_pwl, double* Eii_exp , double* Eii_lin, double* Eii_gbs, double* Eii_cst, double* Exx_el, double* Ezz_el, double* Exz_el, double* Exx_diss, double* Ezz_diss, double* Exz_diss, double *d1, double strain_acc, double dil, double fric, double C, double *detadexx, double *detadezz, double *detadexz, double *detadp, double P0,  double *X1, double *OverS, double *ddivpdexx, double *ddivpdezz, double *ddivpdexz, double *ddivpdp, double *Pcorr, double *drhodp, double *rho, double beta, double div, double *div_el, double *div_pl, double *div_r, double *Wtot, double *Wdiss, double* Wel ) {

    // General paramaters
    double eta=0.0, R=materials->R, dt=model->dt;
    double minEta=model->mineta, maxEta=model->maxeta;
    double TmaxPeierls = (1200.0+zeroC)/scaling->T;      // max. T for Peierls

    // Parameters for deformation map calculations
    int    local_iter = model->loc_iter, it, nitmax = 20, noisy = 0;
    int    constant=0, dislocation=0, peierls=0, diffusion=0, gbs=0, elastic = model->iselastic, comp = model->compressible, VolChangeReac = model->VolChangeReac;
    double tol = 1.0e-11, res=0.0, res0=0.0, dfdeta=0.0, Txx=0.0, Tzz=0.0, Txz=0.0, Tii=0.0, ieta_sum=0.0, Tii0 = sqrt(Txx0*Txx0 + Txz0*Txz0);
    double eta_up=0.0, eta_lo=0.0, eta_ve=0.0, eta_p=0.0, r_eta_pl=0.0, r_eta_ve=0.0, r_eta_p=0.0;
    double eta_pwl=0.0, eta_exp=0.0, eta_vep=0.0, eta_lin=0.0, eta_el=0.0, eta_gbs=0.0, eta_cst=0.0, eta_step=0.0;
    double Exx=0.0, Ezz=0.0, Exz=0.0, Eii_vis=0.0, Eii= 0.0, eII=0.0;
    double f1=0.0, f2=0.0, X = X0;

    // Flow law parameters from input file
    double Tyield=0.0, F_trial = 0.0, F_corr = 0.0, gdot = 0.0, dQdtxx = 0.0, dQdtyy = 0.0, dQdtzz= 0.0, dQdtxz= 0.0;
    int    is_pl = 0;
    double Ea_pwl = materials->Qpwl[phase], Va_pwl = materials->Vpwl[phase], n_pwl  = materials->npwl[phase], m_pwl  = materials->mpwl[phase], r_pwl  = materials->rpwl[phase], A_pwl  = materials->Apwl[phase], f_pwl = materials->fpwl[phase], a_pwl = materials->apwl[phase], F_pwl  = materials->Fpwl[phase], pre_factor = materials->pref_pwl[phase], t_pwl  = materials->tpwl[phase];
    double Ea_lin = materials->Qlin[phase], Va_lin = materials->Vlin[phase], n_lin  = materials->nlin[phase], m_lin  = materials->mlin[phase], r_lin  = materials->rlin[phase], A_lin  = materials->Alin[phase], f_lin = materials->flin[phase], a_lin = materials->alin[phase], F_lin  = materials->Flin[phase];
    double Ea_gbs = materials->Qgbs[phase], Va_gbs = materials->Vgbs[phase], n_gbs  = materials->ngbs[phase], m_gbs  = materials->mgbs[phase], r_gbs  = materials->rgbs[phase], A_gbs  = materials->Agbs[phase], f_gbs = materials->fgbs[phase], a_gbs = materials->agbs[phase], F_gbs  = materials->Fgbs[phase];
    double B_pwl=0.0, B_lin=0.0, B_exp=0.0, B_gbs=0.0, C_pwl=0.0,  C_lin=0.0, C_exp=0.0, C_gbs=0.0;
    double Ea_exp = materials->Qexp[phase], S_exp  = materials->Sexp[phase], E_exp  = materials->Eexp[phase] , t_exp  = materials->texp[phase], F_exp=0.0;
    double gamma=materials->Gexp[phase], ST, q=materials->qexp[phase],  n_exp=materials->nexp[phase];
    double Exx_lin=0.0, Ezz_lin=0.0, Exz_lin=0.0, Exx_exp=0.0, Ezz_exp=0.0, Exz_exp=0.0, Exx_gbs=0.0, Ezz_gbs=0.0, Exz_gbs=0.0, Exx_cst=0.0, Ezz_cst=0.0, Exz_cst=0.0;
    double Exx_pl=0.0, Exx_pwl=0.0, Ezz_pl=0.0, Ezz_pwl=0.0, Exz_pl=0.0, Exz_pwl=0.0;
    int gs = materials->gs[phase];
    double pg = materials->ppzm[phase], Kg = materials->Kpzm[phase], Qg = materials->Qpzm[phase], gam = materials->Gpzm[phase], cg = materials->cpzm[phase], lambda = materials->Lpzm[phase];
    double eta_vp0 = materials->eta_vp[phase], n_vp = materials->n_vp[phase], eta_vp;
    double  detadTxx=0.0, detadTzz=0.0, detadTxz=0.0, deta_ve_dExx=0.0, deta_ve_dEzz=0.0, deta_ve_dExz=0.0, deta_ve_dP=0.0;
    double  dFdExx=0.0, dFdEzz=0.0, dFdExz=0.0, dFdP=0.0, K = 1.0/beta, dQdP=0.0, g=0.0, dlamdExx=0.0, dlamdEzz=0.0, dlamdExz=0.0, dlamdP=0.0, a=0.0, deta_vep_dExx=0.0, deta_vep_dEzz=0.0, deta_vep_dExz=0.0, deta_vep_dP=0.0, deta_vp_dExx, deta_vp_dEzz, deta_vp_dExz, deta_vp_dP, deta;

    double F_trial0 = F_trial;
    double dFdgdot, divp=0.0, Pc = P, dummy;

    // Mix of constant viscosity
    int phase_two    = materials->phase_two[phase];
    int constant_mix = materials->phase_mix[phase];
    int mix_avg      = model->diffuse_avg;
    double ndis1, Adis1, Qdis1, rho1;
    double ndis2, Adis2, Qdis2, rho2;
    int ProgressiveReaction = materials->reac_soft[phase], NoReturn = model->NoReturn;
    int StaticReaction = model->diffuse_X==1;
    double tau_kin = materials->tau_kin[phase], Pr = materials->Pr[phase], dPr = materials->dPr[phase];
    double dXdP=0.0, d2XdP2=0.0, dEadP=0.0, dfdn=0.0, dndP=0.0, dfdP=0.0, dAdP=0.0, retro = 1.0;
    double ddivrdpc = 0.0, divr = 0.0, drhodX = 0.0;
    double dCgbsdP = 0.0, dClindP = 0.0, dCpwldP = 0.0, Jii, alpha, rho_ref, drho_ref_dP;
    
    rho1   = materials->rho[phase];
    alpha  = materials->alp[phase];

    if (model->diffuse_X==0) constant_mix  = 0;

    //------------------------------------------------------------------------//

    // Initialise strain rate invariants to 0
    *Eii_exp = 0.0; *Eii_lin = 0.0; *Eii_pl = 0.0; *Eii_pwl = 0.0; *Eii_el = 0.0, *Eii_gbs=0, *Eii_cst=0.0;
    *txxn=0.0; *tzzn=0.0; *txzn=0.0; *etaVE=0.0; *VEcoeff=0.0, *Exx_el=0.0, *Ezz_el=0.0, *Exz_el=0.0, *Exx_diss=0.0, *Ezz_diss=0.0, *Exz_diss=0.0, *d1=0.0;
    *detadexx=0.0;  *detadezz=0.0;  *detadexz=0.0;  *detadp=0.0;
    *ddivpdexx=0.0; *ddivpdezz=0.0; *ddivpdexz=0.0; *ddivpdp=0.0;
    *X1=0.0; *drhodp = 0.0; *rho = 0.0;
    *OverS = 0.0; *Pcorr = 0.0; *div_el = 0.0; *div_pl = 0.0; *div_r = 0.0; *Wtot = 0.0; *Wdiss = 0.0; *Wel = 0.0;

    // P corr will be corrected if plasticity feedbacks on pressure (dilation)
    *Pcorr = P;

    // Activate deformation mechanisms
    if ( materials->cstv[phase] !=0                  ) constant    = 1;
    if ( materials->pwlv[phase] !=0                  ) dislocation = 1;
    if ( materials->expv[phase] !=0 && T<TmaxPeierls ) peierls     = 1;
    if ( materials->linv[phase] !=0                  ) diffusion   = 1;
    if ( materials->gbsv[phase] !=0                  ) gbs         = 1;
    if ( materials->gs[phase]   !=0                  ) gs          = 1;

    // Turn of elasticity for the initialisation step  (viscous flow stress)
    if ( model->step    == 0                         ) elastic     = 0;

    // Constant grain size initially
    *d1 = materials->gs_ref[phase];

    // Tensional cut-off
//    if ( model->gz<0.0 && P<0.0     ) { P = 0.0; printf("Aie aie aie P < 0 !!!\n"); exit(122);}
    
    // Visco-plastic limit
    if ( elastic==0                 ) { G = 1e1; K = 1e1; dil = 0.0;};

    // Zero C limit
    if ( T< zeroC/scaling->T        ) T = zeroC/scaling->T;

    //------------------------------------------------------------------------//

    // Precomputations
    if ( dislocation == 1 ) {
        B_pwl = pre_factor * F_pwl * pow(A_pwl,-1.0/n_pwl) * exp( (Ea_pwl + P*Va_pwl)/R/n_pwl/T ) * pow(d, m_pwl/n_pwl) * pow(f_pwl, -r_pwl/n_pwl) * exp(-a_pwl*phi/n_pwl);
        C_pwl   = pow(2.0*B_pwl, -n_pwl);
        dCpwldP = -C_pwl*Va_pwl/R/T;
    }
    if ( diffusion == 1 ) {
        if (m_lin>0.0 && d<1e-13/scaling->L){
            printf("Cannot run with grain size dependent viscosity if grain size is set to 0 --> d = %2.2e!!!\n", d*scaling->L);
            exit(1);
        };
        B_lin = F_lin * pow(A_lin,-1.0/n_lin) * exp( (Ea_lin + P*Va_lin)/R/n_lin/T ) * pow(f_lin, -r_lin/n_lin) * exp(-a_lin*phi/n_lin); // * pow(d, m_lin/n_lin) !!!!!!!!!!!!!!!!!!!!!!!!
        C_lin = pow(2.0*B_lin, -n_lin);
        dClindP = -C_lin*Va_lin/R/T;
    }
    if ( gbs == 1 ) {
        B_gbs = F_gbs * pow(A_gbs,-1.0/n_gbs) * exp( (Ea_gbs + P*Va_gbs)/R/n_gbs/T ) * pow(d, m_gbs/n_gbs) * pow(f_gbs, -r_gbs/n_gbs) * exp(-a_gbs*phi/n_gbs);
        C_gbs = pow(2.0*B_gbs, -n_gbs);
        dCgbsdP = -C_gbs*Va_gbs/R/T;
    }
    if ( peierls   == 1 ) {
        ST                           = Ea_exp/R/T * pow((1.0-gamma),(q-1.0)) * q*gamma;
        if ( (int)t_exp == 0) F_exp  = 1.0;
        if ( (int)t_exp == 1) F_exp  = 1.0/6.0*pow(2.0,1.0/(ST+n_exp)) * pow(3.0,(ST+n_exp-1.0)/2.0/(ST+n_exp));
        if ( (int)t_exp == 2) F_exp  = 1.0/4.0*pow(2,1.0/(ST+n_exp));
        //        B_exp                   = F_exp * pow(E_exp*exp(-Ea_exp/R/T*pow(1.0-gamma,2.0)), -1.0/(ST+n_exp)) * pow(gamma*S_exp, ST/(ST+n_exp));
        //        C_exp                   = pow(2.0*B_exp, -(ST+n_exp));
        C_exp = E_exp *exp(-Ea_exp/R/T * pow(1.0-gamma,2.0)) * pow(gamma*S_exp,-ST);  // ajouter Fexp
        B_exp = 0.5*pow(C_exp, -1./(n_exp+ST) );
    }

    double cos_fric = cos(fric);
    double sin_fric = sin(fric);
    double sin_dil = sin(dil);
    Tyield     = C*cos_fric + P*sin_fric;

    // Von-Mises cut-off
    if (materials->Slim[phase] < Tyield) {
        sin_fric   = 0.0;
        cos_fric   = 1.0;
        sin_dil    = 0.0;
        C          = materials->Slim[phase];
        Tyield     = materials->Slim[phase];
    }
    
    // Tension cut-off
    int    is_tensile = materials->is_tensile[phase], tens = 0;
    double P_tens = -25e6/scaling->S;
    double sin_fric_tens = -C*cos_fric/P_tens; // Compute appropriate slope of the tension yield (SimpleYields.m)
    
    if (C*cos_fric + P*sin_fric_tens < Tyield && is_tensile==1 ) {
        if (noisy>0) printf("Switching to tension yield stress, P = %2.2e\n", P*scaling->S);
        sin_fric   = sin_fric_tens;
        Tyield     = C*cos_fric + P*sin_fric_tens;
        tens       = 1;
    }

    //------------------------------------------------------------------------//
    // Reaction stuff: 1. Update reaction progress
    X = X0;

    // Reaction stuff: 2. Mixture rheology (Huet et al., 2014)
    if (  ProgressiveReaction == 1 ) {

        if ( model->UnsplitDiffReac == 0 ) {
            
            X     = (X0*tau_kin - 0.5*dt*erfc((P - Pr)/dPr) + dt)/(dt + tau_kin);
            dXdP   = retro*dt*exp(-pow(P - Pr, 2.0)/pow(dPr, 2.0))/(sqrt(M_PI)*dPr*(dt + tau_kin));
            d2XdP2 = -retro*dt*(2.0*P - 2.0*Pr)*exp(-pow(P - Pr, 2.0)/pow(dPr, 2.0))/(sqrt(M_PI)*pow(dPr, 3.0)*(dt + tau_kin));

            if ( X<X0 && NoReturn == 1 ) {
                retro  = 0.0;
                X      = X0;
                dXdP   = 0.0;
                d2XdP2 = 0.0;
            }

        }

        ndis1  = materials->npwl[phase];
        Adis1  = materials->Apwl[phase];
        Qdis1  = materials->Qpwl[phase];

        ndis2  = materials->npwl[materials->reac_phase[phase]];
        Adis2  = materials->Apwl[materials->reac_phase[phase]];
        Qdis2  = materials->Qpwl[materials->reac_phase[phase]];
        rho2   = materials->rho[materials->reac_phase[phase]];

        // Huet et al 2014 ---------------
        f1 = 1.0-X;
        f2 = X;

        // (1) Calcul des ai
        double a1 = ndis2 + 1.0;
        double a2 = ndis1 + 1.0;

        // (2) Calcul de n bulk:
        double sum_up   = f1*a1*ndis1 + f2*a2*ndis2;
        double sum_down = f1*a1+f2*a2;
        n_pwl     = sum_up/sum_down;

        // (3) Calcul de Q bulk:
        sum_up    = f1*a1*Qdis1 + f2*a2*Qdis2;
        sum_down  = f1*a1+f2*a2;
        Ea_pwl    = sum_up/sum_down;

        // (4) Calcul de A bulk:
        sum_down     = f1*a1 + f2*a2;
        double Prod1 = pow(Adis1,(f1*a1/sum_down)) * pow(Adis2,(f2*a2/sum_down));
        double sum_n = f1*ndis1/(ndis1 + 1.0) + f2*ndis2/(ndis2 + 1.0);
        double Prod2 = pow(ndis1/(ndis1 + 1.0),f1*a1*ndis1/sum_down) * pow(ndis2/(ndis2+1.0),f2*a2*ndis2/sum_down);
        A_pwl        = Prod1 * pow(sum_n,-n_pwl) * Prod2;

        // Partial derivatives of mixing power law parameters w.r.t. pressure
        dndP = n_pwl*(a1*dXdP - a2*dXdP)/sum_down + (-a1*dXdP*ndis1 + a2*dXdP*ndis2)/sum_down;
        dAdP = A_pwl*(((a1*dXdP*ndis1 - a2*dXdP*ndis2)/sum_down - (a1*dXdP - a2*dXdP)*(X*a2*ndis2 + a1*ndis1*(1.0 - X))/pow(sum_down, 2))*log(X*ndis2/a1 + ndis1*(1.0 - X)/a2) - (X*a2*ndis2 + a1*ndis1*(1.0 - X))*(-dXdP*ndis1/a2 + dXdP*ndis2/a1)/(sum_down*(X*ndis2/a1 + ndis1*(1.0 - X)/a2))) + A_pwl*(-a1*dXdP/sum_down + a1*(1.0 - X)*(a1*dXdP - a2*dXdP)/pow(sum_down, 2))*log(Adis1) + A_pwl*(X*a2*(a1*dXdP - a2*dXdP)/pow(sum_down, 2) + a2*dXdP/sum_down)*log(Adis2) + A_pwl*(-a1*dXdP*ndis1/sum_down + a1*ndis1*(1.0 - X)*(a1*dXdP - a2*dXdP)/pow(sum_down, 2))*log(ndis1/a2) + A_pwl*(X*a2*ndis2*(a1*dXdP - a2*dXdP)/pow(sum_down, 2) + a2*dXdP*ndis2/sum_down)*log(ndis2/a1);
        dQdP = Ea_pwl*(a1*dXdP - a2*dXdP)/sum_down + (-Qdis1*a1*dXdP + Qdis2*a2*dXdP)/sum_down;

        // Proper choices of corrections factors
        if ( (int)t_pwl == 0 ) {
            F_pwl = 1.0;
            dfdn = 0.0;
        }
        if ( (int)t_pwl == 1 ) {
            F_pwl = 1.0/6.0*pow(2.0,1.0/n_pwl) * pow(3.0,(n_pwl-1.0)/2.0/n_pwl);
            dfdn = -1.0/6.0*pow(2, 1.0/n_pwl)*pow(3, ((1.0/2.0)*n_pwl - 1.0/2.0)/n_pwl)*M_LN2/pow(n_pwl, 2) + F_pwl*((1.0/2.0)/n_pwl + (1.0/2.0 - 1.0/2.0*n_pwl)/pow(n_pwl, 2))*log(3);
        }
        if ( (int)t_pwl == 2 ) {
            F_pwl = 1.0/4.0*pow(2,1.0/n_pwl);
            dfdn = -0.25*pow(2, 1.0/n_pwl)*M_LN2/pow(n_pwl, 2);

        }
        dfdP    = dfdn*dndP;

        // Override power-law flow law parameters
        B_pwl  = pre_factor * F_pwl * pow(A_pwl,-1.0/n_pwl) * exp( (Ea_pwl)/R/n_pwl/T );
        C_pwl  = pow(2.0*B_pwl, -n_pwl);

    }

    // set pointer value
    *X1 = X;
    //------------------------------------------------------------------------//

    // Out-of-plane components
    double Eyy, Tyy;
    double eyy  = -(  exx +  ezz ); // definition of deviatoric tensor
    double Tyy0 = -( Txx0 + Tzz0 ); // definition of deviatoric tensor

    // Isolated viscosities
    Exx  = exx;
    Ezz  = ezz;
    Exz  = exz;
    Eyy  = eyy;
    if ( elastic == 0 || local_iter == 0 || local_iter == 2 ) {
        Exx  = exx;
        Ezz  = ezz;
        Exz  = exz;
        Eyy  = eyy;
    }
    else {
        Exx  = exx + Txx0/(2.0*G*dt);
        Ezz  = ezz + Tzz0/(2.0*G*dt);
        Exz  = exz + Txz0/(2.0*G*dt);
        Eyy  = eyy + Tyy0/(2.0*G*dt);
    }
    Eii  = sqrt(1.0/2.0*(Exx*Exx + Ezz*Ezz + Eyy*Eyy) + Exz*Exz);
    eII  = sqrt(1.0/2.0*(exx*exx + ezz*ezz + eyy*eyy) + exz*exz);
    if (Eii*scaling->E<1e-30) Eii=1e-30/scaling->E;

    //------------------------------------------------------------------------//

    // Isolated viscosities
    eta_el                           = G*dt;
    if ( constant    == 1 ) eta_cst  = materials->eta0[phase];
    if ( constant_mix== 1 && mix_avg==0) eta_cst  =     X0*materials->eta0[phase] + (1-X0)*materials->eta0[phase_two];
    if ( constant_mix== 1 && mix_avg==1) eta_cst  = pow(X0/materials->eta0[phase] + (1-X0)/materials->eta0[phase_two], -1.0);
    if ( constant_mix== 1 && mix_avg==2) eta_cst  = exp(X0*log(materials->eta0[phase]) + (1-X0)*log(materials->eta0[phase_two]));
    if ( dislocation == 1 ) eta_pwl  = B_pwl * pow( Eii, 1.0/n_pwl - 1.0 );
    if ( diffusion   == 1 ) eta_lin  = B_lin * pow( Eii, 1.0/n_lin - 1.0 ) * pow(d, m_lin/n_lin); // !!! gs - dependence !!!
    if ( gbs         == 1 ) eta_gbs  = B_gbs * pow( Eii, 1.0/n_gbs - 1.0 );
    if ( peierls     == 1 ) eta_exp  = B_exp * pow( Eii, 1.0/(ST+n_exp) - 1.0 );
    
    //------------------------------------------------------------------------//

    // Viscoelasticity
    *Eii_pl = 0.0;

    // Define viscosity bounds
    eta_up   = 1.0e100/scaling->eta;
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
//        if (noisy>=1) printf("Visco-Elastic iterations It. %02d, r abs. = %2.2e r rel. = %2.2e tol = %2.2e eta_ve = %2.2e eta_lo = %2.2e eta_up = %2.2e eta_lin = %2.2e eta_pwl = %2.2e eta_exp = %2.2e %d\n", it, res, res/res0, tol, eta_ve, eta_lo, eta_up, eta_lin, eta_pwl, eta_exp, nitmax);
//        printf("Visco-Elastic iterations It. %02d, r abs. = %2.2e r rel. = %2.2e tol = %2.2e eta_ve = %2.2e eta_lo = %2.2e eta_up = %2.2e eta_lin = %2.2e eta_pwl = %2.2e eta_exp = %2.2e, exz=%2.2e\n", it, res, res/res0, tol, eta_ve, eta_lo, eta_up, eta_lin, eta_pwl, eta_exp, exz*scaling->E);
        if (res < tol/100) break;

        // Analytical derivative of function
        dfdeta  = 0.0;
        if ( elastic     == 1 ) dfdeta += -Eii/eta_el;
        if ( peierls     == 1 ) dfdeta += -(*Eii_exp)*(ST+n_exp)/eta_ve;
        if ( diffusion   == 1 ) dfdeta += -(*Eii_lin)*n_lin/eta_ve;
        if ( dislocation == 1 ) dfdeta += -(*Eii_pwl)*n_pwl/eta_ve;
        if ( constant    == 1 ) dfdeta += -Eii/eta_cst;

        // Update viscosity
        eta_ve -= r_eta_ve / dfdeta;
    }
    
    if ( it==nitmax-1 && res > tol ) { printf("Visco-Elastic iterations failed!\n"); exit(0);}

    // Recalculate stress components
    Txx                  = 2.0*eta_ve*Exx;
    Tyy                  = 2.0*eta_ve*Eyy;  // out-of-plane
    Tzz                  = 2.0*eta_ve*Ezz;
    Txz                  = 2.0*eta_ve*Exz;
    Tii                  = 2.0*eta_ve*Eii;

    /*-----------------------------------------------*/

    if ( dislocation == 1 ) eta_pwl  = pow(2.0*C_pwl,-1.0) * pow(Tii, 1.0-n_pwl);
    if ( diffusion   == 1 ) eta_lin  = B_lin * pow( *Eii_lin, 1.0/n_lin - 1.0 ) * pow(*d1, m_lin/n_lin);
    if ( peierls     == 1 ) eta_exp  = B_exp * pow( *Eii_exp, 1.0/(ST+n_exp) - 1.0 );
    if ( elastic     == 1 ) *Exx_el =  (double)elastic*(Txx-Txx0)/2.0/eta_el;
    if ( elastic     == 1 ) *Ezz_el =  (double)elastic*(Tzz-Tzz0)/2.0/eta_el;
    if ( elastic     == 1 ) *Exz_el =  (double)elastic*(Txz-Txz0)/2.0/eta_el;
    Exx_pl = 0.0;
    Ezz_pl = 0.0;
    Exz_pl = 0.0;
    if ( dislocation == 1 ) Exx_pwl = Txx/2.0/eta_pwl;
    if ( dislocation == 1 ) Ezz_pwl = Tzz/2.0/eta_pwl;
    if ( dislocation == 1 ) Exz_pwl = Txz/2.0/eta_pwl;
    if ( diffusion   == 1 ) Exx_lin = Txx/2.0/eta_lin;
    if ( diffusion   == 1 ) Ezz_lin = Tzz/2.0/eta_lin;
    if ( diffusion   == 1 ) Exz_lin = Txz/2.0/eta_lin;
    if ( peierls     == 1 ) Exx_exp = Txx/2.0/eta_exp;
    if ( peierls     == 1 ) Ezz_exp = Tzz/2.0/eta_exp;
    if ( peierls     == 1 ) Exz_exp = Txz/2.0/eta_exp;
    if ( constant    == 1 ) Exx_cst = Txx/2.0/eta_cst;
    if ( constant    == 1 ) Ezz_cst = Tzz/2.0/eta_cst;
    if ( constant    == 1 ) Exz_cst = Txz/2.0/eta_cst;

//    printf("exx_tot = %2.2e exx_el = %2.2e exx_pwl = %2.2e exx_exp = %2.2e exx_cst = %2.2e exx_exp = %2.2e exx_pl = %2.2e, exx_net = %2.2e\n", exx*scaling->E, *Exx_el*scaling->E, Exx_pwl*scaling->E, Exx_lin*scaling->E, Exx_cst*scaling->E, Exx_exp*scaling->E, Exx_pl*scaling->E, (exx - (*Exx_el+Exx_pwl+Exx_lin+Exx_cst+Exx_exp+Exx_pl))*scaling->E);
//    printf("ezz_tot = %2.2e ezz_el = %2.2e exx_pwl = %2.2e exx_exp = %2.2e exx_cst = %2.2e exx_exp = %2.2e ezz_pl = %2.2e, ezz_net = %2.2e\n", ezz*scaling->E, *Ezz_el*scaling->E, Ezz_pwl*scaling->E, Ezz_lin*scaling->E, Ezz_cst*scaling->E, Ezz_exp*scaling->E, Ezz_pl*scaling->E, (ezz - (*Ezz_el+Ezz_pwl+Ezz_lin+Ezz_cst+Ezz_exp+Ezz_pl))*scaling->E);
//      //       printf("exx_el = %2.2e exx_vis = %2.2e exx_pl = %2.2e, exx_net = %2.2e\n", *Exx_el*scaling->E, (Exx_pwl+Exx_lin+Exx_cst+Exx_exp)*scaling->E, Exx_pl*scaling->E, (exx - (*Exx_el+Exx_pwl+Exx_lin+Exx_cst+Exx_exp+Exx_pl))*scaling->E);
//    printf("exz_tot = %2.2e exz_el = %2.2e exz_pwl = %2.2e exz_exp = %2.2e exz_cst = %2.2e exz_exp = %2.2e exz_pl = %2.2e, exz_net = %2.2e\n", exz*scaling->E, *Exz_el*scaling->E, Exz_pwl*scaling->E, Exz_lin*scaling->E, Exz_cst*scaling->E, Exz_exp*scaling->E, Exz_pl*scaling->E, (exz - (*Exz_el+Exz_pwl+Exz_lin+Exz_cst+Exz_exp+Exz_pl))*scaling->E);
//    exit(1);

    //    printf("%2.10e %2.10e\n", Tii, sqrt(0.5*(Txx*Txx + Tyy*Tyy + Tzz*Tzz) +Txz*Txz) );

    // Partial derivatives VE --> GLOBAL NEWTON ITERATION
    deta     = 0.0;
    if (dislocation == 1) deta += -C_pwl * pow(Tii, n_pwl    - 3.0)*(n_pwl    - 1.0);
    if (peierls     == 1) deta += -C_exp * pow(Tii, n_exp+ST - 3.0)*(n_exp+ST - 1.0);
    if (diffusion   == 1) deta += -C_lin * pow(Tii, n_lin    - 3.0)*(n_lin    - 1.0) * pow(*d1,-m_lin);
//    if (gbs         == 1)
    detadTxx     =     deta*(2.0*Txx + Tzz)*pow(eta_ve,2.0); // Out-of-plane deviatoric stress was substituted here
    detadTzz     =     deta*(2.0*Tzz + Txx)*pow(eta_ve,2.0); // ... and there
    detadTxz     = 2.0*deta*Txz*pow(eta_ve,2.0);
    deta_ve_dExx = detadTxx * 2.0*eta_ve / (1.0 - 2.0*(detadTxx*Exx + detadTzz*Ezz + detadTxz*Exz));
    deta_ve_dEzz = detadTzz * 2.0*eta_ve / (1.0 - 2.0*(detadTxx*Exx + detadTzz*Ezz + detadTxz*Exz));
    deta_ve_dExz = detadTxz * 2.0*eta_ve / (1.0 - 2.0*(detadTxx*Exx + detadTzz*Ezz + detadTxz*Exz));

    dummy = 0.0;
    Jii = Tii*Tii;
    if (dislocation == 1) dummy += - 2.0*pow(Jii, (1.0/2.0)*n_pwl - 1.0/2.0)*dCpwldP;
    if (diffusion   == 1) dummy += - 2.0*pow(Jii, (1.0/2.0)*n_lin - 1.0/2.0)*pow(d, -m_lin)*dClindP;
//    if (gbs         == 1) dummy += - 2*pow(Jii, (1.0/2.0)*n_gbs - 1.0/2.0)*dCgbsdP;
    deta_ve_dP   = pow(eta_ve, 2)*dummy;

    if (  ProgressiveReaction == 1 ) { // if activation volume is activated: This will not work !!!!!!!!!
        deta_ve_dP   = pow(B_pwl, n_pwl)*Eii*eta_ve*pow(eta_el, 2)*pow(Eii*eta_ve, n_pwl)*(-A_pwl*F_pwl*Ea_pwl*dndP + A_pwl*F_pwl*R*T*dndP*n_pwl*log(B_pwl) - A_pwl*F_pwl*R*T*dndP*n_pwl*log(Eii*eta_ve) + A_pwl*F_pwl*R*T*dndP*log(A_pwl) + A_pwl*F_pwl*dQdP*n_pwl + A_pwl*R*T*dfdP*pow(n_pwl, 2) - F_pwl*R*T*dAdP*n_pwl)/(A_pwl*F_pwl*R*T*n_pwl*(pow(B_pwl, 2*n_pwl)*pow(Eii, 2)*pow(eta_ve, 2) + 2*pow(B_pwl, n_pwl)*Eii*eta_ve*eta_el*pow(Eii*eta_ve, n_pwl) + pow(B_pwl, n_pwl)*Eii*pow(eta_el, 2)*n_pwl*pow(Eii*eta_ve, n_pwl) - pow(B_pwl, n_pwl)*Eii*pow(eta_el, 2)*pow(Eii*eta_ve, n_pwl) + pow(eta_el, 2)*pow(Eii*eta_ve, 2*n_pwl)));
    }

    //------------------------------------------------------------------------//

    // Check yield stress
    F_trial = Tii - Tyield;

    double Tyield_trial = Tyield;
    double F_corr1=F_trial, F_corr2=F_trial, Tc, Pc_chk, Tc_chk;
    
    // Select appropriate dilation angle for tensile domain, see SimpleYields.m
    if (tens == 1) {
        eta_vp   = eta_vp0 * pow(eII, 1.0/n_vp - 1);
        Pc_chk    = -(C*cos_fric) / (- Tii/P + sin_fric);
        Tc_chk    = Tii/P*Pc_chk;
        
//        sin_dil = (C*cos_fric*eta_ve + C*cos_fric*eta_vp + P*eta_ve*sin_fric + P*eta_vp*sin_fric - Tc_chk*eta_ve - Tc_chk*eta_vp)/(K*dt*sin_fric*(Tc_chk - Tii));
        sin_dil = (P*eta_ve + P*eta_vp - Pc_chk*eta_ve - Pc_chk*eta_vp)/(K*dt*(C*cos_fric + Pc_chk*sin_fric - Tii));
//        sin_psi2 = (P.*eta_ve + P.*eta_vp - Pc.*eta_ve - Pc.*eta_vp)./(K.*dt.*(C.*cos_phi + Pc.*sin_phi - Tii));
    }

    if (F_trial > 1e-17) {

        // Initial guess - eta_vp = 0
        is_pl    = 1;
        eta_vp   = eta_vp0 * pow(eII, 1.0/n_vp - 1);
        gdot     = F_trial / ( eta_ve + eta_vp + K*dt*sin_fric*sin_dil);
        dQdtxx   = Txx/2.0/Tii;
        dQdtyy   = Tyy/2.0/Tii;
        dQdtzz   = Tzz/2.0/Tii;
        dQdtxz   = Txz/1.0/Tii;
        dQdP     = -sin_dil; //printf("%2.2e %2.2e\n", dQdP, K*scaling->S);
        res0     = F_trial;
        F_trial0 = F_trial;

        // Return mapping --> find plastic multiplier rate (gdot)
        for (it=0; it<nitmax; it++) {

            dQdP    = -sin_dil;
            divp    = -gdot*dQdP;
            Pc      = P + K*dt*divp; // P0 - k*dt*(div-divp) = P + k*dt*divp
            if (noisy>0 && tens==1) { printf("Pc = %2.4e Pc_chk=%2.4e sin_dil = %2.2e\n",Pc, Pc_chk, sin_dil);  };
            eta_vp  = eta_vp0 * pow(fabs(gdot), 1.0/n_vp - 1.0);
//            if (tens==1) sin_dil = (F_trial0 - eta_ve*gdot - eta_vp*gdot)/(K*dt*gdot*sin_fric);
            Tyield  = C*cos_fric + Pc*sin_fric +  gdot*eta_vp;
            F_trial = Tii - eta_ve*gdot - Tyield;

            // Residual check
            res = fabs(F_trial);
            if (noisy>0 ) printf("%02d Visco-Plastic iterations It., tens = %d F = %2.2e Frel = %2.2e --- n_vp = %2.2e, eta_vp = %2.2e\n", it, tens, res, res/F_trial0, n_vp, eta_vp*scaling->eta);
            if ( res < tol || res/F_trial0 < tol ) break;
            dFdgdot  = - eta_ve - eta_vp/n_vp - K*dt*sin_fric*sin_dil;
            gdot    -= F_trial / dFdgdot;
            
        }
        if ( it==nitmax-1 && (res > tol || res/F_trial0 > tol)  ) { printf("Visco-Plastic iterations failed!\n"); exit(0);}

        double Tii_corr= Tii- eta_ve*gdot;
        double Tii_trial=Tii;

        
        Txx     = 2.0*eta_ve*(Exx - gdot*dQdtxx    );
//        Tyy     = 2.0*eta_ve*(Eyy - gdot*dQdtyy    ); // Tyy     = -(Txx+Tzz);
        Tzz     = 2.0*eta_ve*(Ezz - gdot*dQdtzz    );
        Tyy     = -(Txx+Tzz);

        //        printf("%2.6e %2.6e\n", Tyy,  2.0*eta_ve*(Eyy - gdot*dQdtyy    ) );
        Txz     = 2.0*eta_ve*(Exz - gdot*dQdtxz/2.0);
        Tii     = sqrt( 0.5*( pow(Txx,2.0) + pow(Tzz,2.0) + pow(Tyy,2.0) ) + pow(Txz,2.0)  );
        eta_vep = Tii / (2.0*Eii);
        F_corr1  = Tii - Tyield;
        *Eii_pl = gdot/2.0;
        Tii     = 2.0*eta_vep*Eii;
        F_corr2  = Tii - Tyield;
        // Compute over stress
        *OverS = eta_vp*gdot;
        *Pcorr = Pc;
  
    
    if (is_pl ==1 && (fabs(F_corr1) > 1e-7) && (fabs(F_corr1)/res0 > 1e-7) ) {
        printf("Problem with F_corr in phase = %d\n", phase);
        printf("is_pl %d --- it = %03d --- tens = %d --- ty = %2.2e F_trial_0 = %2.2e --- F_trial = %2.2e --- F_corr = %2.2e --- F_corr1 = %2.2e --- F_corr2 = %2.2e\n", is_pl, it, tens, Tyield*scaling->S, res0, res, F_corr, F_corr1, F_corr2);
        printf("Pc = %2.2e P0 = %2.2e\n", Pc*scaling->S, P*scaling->S);
        printf("Tii_trial = %2.2e Tii_corr1 = %2.2e Tii_corr = %2.2e\n", Tii_trial*scaling->S, Tii*scaling->S, Tii_corr*scaling->S);
//        exit(90);
        printf("WARNING: Non linear convergence is expected to be poor\n");
    }
    }
    

    //------------------------------------------------------------------------//
    double dsin_dildExx = 0.0, dsin_dildEzz = 0.0, dsin_dildExz = 0.0, dsin_dildP = 0.0;

    
    // Partial derivatives VEP
    if ( is_pl == 1 ) {
        
        if (tens == 1 ) {
            double dEii_dExx = (2.0*Exx + Ezz)/2.0/Eii;
            double dEii_dEzz = (2.0*Ezz + Exx)/2.0/Eii;
            double dEii_dExz = Exz/1.0/Eii;
            double dEii_dP   = 0.0;
            dsin_dildExx = (2*K*dt*sin_dil*(Eii*deta_ve_dExx + dEii_dExx*eta_ve)*(C*P*cos_fric*sin_fric + pow(2*Eii*eta_ve - P*sin_fric, 2)) + P*(-C*cos_fric*deta_ve_dExx*(2*Eii*eta_ve - P*sin_fric) + 2*C*cos_fric*(eta_ve + eta_vp)*(Eii*deta_ve_dExx + dEii_dExx*eta_ve) + deta_ve_dExx*pow(2*Eii*eta_ve - P*sin_fric, 2)))/(K*dt*(2*Eii*eta_ve - P*sin_fric)*(C*P*cos_fric*sin_fric + (C*cos_fric - 2*Eii*eta_ve)*(2*Eii*eta_ve - P*sin_fric)));
            dsin_dildEzz = (2*K*dt*sin_dil*(Eii*deta_ve_dEzz + dEii_dEzz*eta_ve)*(C*P*cos_fric*sin_fric + pow(2*Eii*eta_ve - P*sin_fric, 2)) + P*(-C*cos_fric*deta_ve_dEzz*(2*Eii*eta_ve - P*sin_fric) + 2*C*cos_fric*(eta_ve + eta_vp)*(Eii*deta_ve_dEzz + dEii_dEzz*eta_ve) + deta_ve_dEzz*pow(2*Eii*eta_ve - P*sin_fric, 2)))/(K*dt*(2*Eii*eta_ve - P*sin_fric)*(C*P*cos_fric*sin_fric + (C*cos_fric - 2*Eii*eta_ve)*(2*Eii*eta_ve - P*sin_fric)));
            dsin_dildExz = (2*K*dt*sin_dil*(Eii*deta_ve_dExz + dEii_dExz*eta_ve)*(C*P*cos_fric*sin_fric + pow(2*Eii*eta_ve - P*sin_fric, 2)) + P*(-C*cos_fric*deta_ve_dExz*(2*Eii*eta_ve - P*sin_fric) + 2*C*cos_fric*(eta_ve + eta_vp)*(Eii*deta_ve_dExz + dEii_dExz*eta_ve) + deta_ve_dExz*pow(2*Eii*eta_ve - P*sin_fric, 2)))/(K*dt*(2*Eii*eta_ve - P*sin_fric)*(C*P*cos_fric*sin_fric + (C*cos_fric - 2*Eii*eta_ve)*(2*Eii*eta_ve - P*sin_fric)));
            dsin_dildP   = (2*C*Eii*cos_fric*(eta_ve + eta_vp)*(P*deta_ve_dP - eta_ve) - C*P*cos_fric*deta_ve_dP*(2*Eii*eta_ve - P*sin_fric) + 2*Eii*K*dt*sin_dil*(C*cos_fric*sin_fric*(P*deta_ve_dP - eta_ve) + deta_ve_dP*pow(2*Eii*eta_ve - P*sin_fric, 2)) + pow(2*Eii*eta_ve - P*sin_fric, 2)*(P*deta_ve_dP + eta_ve + eta_vp))/(K*dt*(2*Eii*eta_ve - P*sin_fric)*(C*P*cos_fric*sin_fric + (C*cos_fric - 2*Eii*eta_ve)*(2*Eii*eta_ve - P*sin_fric)));
        }
        
        dFdExx        =     (2.0*Exx + Ezz)*eta_ve/Eii + 2.0*Eii*deta_ve_dExx;
        dFdEzz        =     (2.0*Ezz + Exx)*eta_ve/Eii + 2.0*Eii*deta_ve_dEzz;
        dFdExz        =            2.0*Exz *eta_ve/Eii + 2.0*Eii*deta_ve_dExz;
        dFdP          =                      -sin_fric + 2.0*Eii*deta_ve_dP;
        g             = 1.0 / ( eta_ve + eta_vp/n_vp + K*dt*sin_fric*sin_dil );
        dlamdExx      = g * (dFdExx - gdot * (deta_ve_dExx + K*dt*sin_fric*dsin_dildExx ));
        dlamdEzz      = g * (dFdEzz - gdot * (deta_ve_dEzz + K*dt*sin_fric*dsin_dildEzz ));
        dlamdExz      = g * (dFdExz - gdot * (deta_ve_dExz + K*dt*sin_fric*dsin_dildExz ));
        dlamdP        = g * (dFdP   - gdot * (deta_ve_dP   + K*dt*sin_fric*dsin_dildP   ));
        deta_vp_dExx  = eta_vp/gdot*dlamdExx*(1.0/n_vp - 1.0);
        deta_vp_dEzz  = eta_vp/gdot*dlamdEzz*(1.0/n_vp - 1.0);
        deta_vp_dExz  = eta_vp/gdot*dlamdExz*(1.0/n_vp - 1.0);
        deta_vp_dP    = eta_vp/gdot*dlamdP  *(1.0/n_vp - 1.0);
        a             =  eta_vp +  K*dt*sin_fric*sin_dil;
        deta_vep_dExx = -0.5*(2.0*Exx + Ezz)*Tyield/(2.0*pow(Eii,3.0)) + (a*dlamdExx + gdot*deta_vp_dExx + gdot*K*dt*sin_fric*dsin_dildExx)/(2.0*Eii);
        deta_vep_dEzz = -0.5*(2.0*Ezz + Exx)*Tyield/(2.0*pow(Eii,3.0)) + (a*dlamdEzz + gdot*deta_vp_dEzz + gdot*K*dt*sin_fric*dsin_dildEzz)/(2.0*Eii);
        deta_vep_dExz =                -Exz *Tyield/(2.0*pow(Eii,3.0)) + (a*dlamdExz + gdot*deta_vp_dExz + gdot*K*dt*sin_fric*dsin_dildExz)/(2.0*Eii);
        deta_vep_dP   =                                      (sin_fric +  a*dlamdP   + gdot*deta_vp_dP   + gdot*K*dt*sin_fric*dsin_dildP  )/(2.0*Eii);
        
        
    }

    // ----------------- Reaction volume changes
    
    if (ProgressiveReaction == 1) {

        rho_ref      = (1.0-X)*rho1 + X*rho2;
        drhodX       = rho2 - rho1;
        drho_ref_dP  = drhodX * dXdP;
        *rho         = rho_ref * exp(1.0/K * Pc - alpha * T);
        *drhodp      = (*rho)/K + exp(1.0/K * Pc - alpha * T) * drho_ref_dP;
//        printf("ph=%d ; X= %2.2e; rho1 = %2.2e, rho2= %2.2e; *rho= %2.2e\n", phase, X, rho1*scaling->rho, rho2*scaling->rho, *rho*scaling->rho);

    }

    else {
        rho_ref      = rho1;
        *rho         = rho_ref * exp(1.0/K * Pc - alpha * T );
        *drhodp      = (*rho)/K;
//        printf("ph=%d ; rho_ref = %2.2e; *rho= %2.2e\n", phase, rho_ref*scaling->rho, *rho*scaling->rho);

    }

    
//    // Activate volume changes only if reaction is taking place
//    if ( VolChangeReac == 1 && fabs(X-X0)>0.0 ) {
//        rho       = rho1 * (1-X) + rho2 * X;
//        drhodX    = rho2 - rho1;
//        drhodP    = drhodX * dXdP;
//        d2rhodP2  = drhodX * d2XdP2;
//        divr      = -1.0/rho * drhodP;
////        printf("%2.2lf %2.2lf divr=%2.2e, rho=%2.2lf x-x0=%2.2lf drhodP = %2.2e dXdP = %2.2e \n", rho1*scaling->rho, rho2*scaling->rho, divr*scaling->E, rho*scaling->rho, X-X0, drhodP*scaling->rho/scaling->S, dXdP/scaling->S);
//
//        ddivrdpc  = -d2rhodP2/rho + pow(drhodP/rho, 2.0);
//        Pc        = Pc + K*dt*divr;
//        *Pcorr    = Pc;
//    }

    // ----------------- Reaction volume changes
    
//    if (P>4.6e9/scaling->S && P<4.7e9/scaling->S) {
//        printf("%2.2e %2.2e\n", eta_ve*scaling->eta, eta_lin*scaling->eta);
//    }

    double inv_eta_diss = 0.0;
    if (peierls    == 1)  inv_eta_diss += (1.0/eta_exp);
    if (dislocation== 1)  inv_eta_diss += (1.0/eta_pwl);
    if (diffusion  == 1)  inv_eta_diss += (1.0/eta_lin);
    if (constant   == 1)  inv_eta_diss += (1.0/eta_cst);
    if (is_pl == 0) {
        (*detadexx) = deta_ve_dExx;
        (*detadezz) = deta_ve_dEzz;
        (*detadexz) = deta_ve_dExz;
        (*detadp)   = deta_ve_dP;
        (*etaVE)    = eta_ve;
        (*ddivpdexx)  = 0.0;
        (*ddivpdezz)  = 0.0;
        (*ddivpdexz)  = 0.0;
        (*ddivpdp)    = ddivrdpc;
    }
    else {
        (*detadexx)   = deta_vep_dExx;
        (*detadezz)   = deta_vep_dEzz;
        (*detadexz)   = deta_vep_dExz;
        (*detadp)     = deta_vep_dP;
        (*etaVE)      = eta_vep;
        (*ddivpdexx)  = -dQdP*dlamdExx + gdot*dsin_dildExx;
        (*ddivpdezz)  = -dQdP*dlamdEzz + gdot*dsin_dildEzz;
        (*ddivpdexz)  = -dQdP*dlamdExz + gdot*dsin_dildExz;
        (*ddivpdp)    = -dQdP*dlamdP   + gdot*dsin_dildP    + ddivrdpc;
        inv_eta_diss += 1.0/eta_vep;
    }
    
//       if (phase==0) printf("phase = %d, eta_cst = %2.2e, eta_ve = %2.2e, eta_vep = %2.2e, eta_up = %2.2e, eta_lo = %2.2e, eta = %2.2e *etaVE=%2.2e elastic=%d constant=%d pl=%d\n",   phase, eta_cst*scaling->eta, eta_ve*scaling->eta, eta_vep*scaling->eta, eta_up*scaling->eta, eta_lo*scaling->eta, eta*scaling->eta, *etaVE*scaling->eta, elastic, constant, is_pl);
    
    

    //    if ( fabs(*detadexx) >1e10 ) {
    //        printf("Viscosity: deta_ve_dExx %2.2e deta_vep_dExx %2.2e\n", deta_ve_dExx, deta_vep_dExx);
    //        printf("exx = %2.2e Exx = %2.2e Tyield = %2.2e, ccos = %2.2e, psin=%2.2e, overS = %2.2e, fric = %2.2e, C = %2.2e\n", exx, Exx, Tyield, C*cos(fric)*scaling->S, Pc*sin(fric)*scaling->S, gdot*eta_vp*scaling->S, fric, C*scaling->S);
    //        printf("Ft =%2.2e Tyt = %2.2e\n", F_trial_trial*scaling->S, Tyield_trial*scaling->S);
    //        exit(1);
    //    }

    /*----------------------------------------------------*/
    /*----------------------------------------------------*/
    /*----------------------------------------------------*/

    eta_pwl  = pow(2.0*C_pwl,-1.0) * pow(Tii, 1.0-n_pwl);
    *Exx_el =  (double)elastic*(Txx-Txx0)/2.0/eta_el;
    *Ezz_el =  (double)elastic*(Tzz-Tzz0)/2.0/eta_el;
    *Exz_el =  (double)elastic*(Txz-Txz0)/2.0/eta_el;
    Exx_pl = gdot*dQdtxx;
    Ezz_pl = gdot*dQdtzz;
    Exz_pl = gdot*dQdtxz/2.0;
    Exx_pwl = Txx/2.0/eta_pwl;
    Ezz_pwl = Tzz/2.0/eta_pwl;
    Exz_pwl = Txz/2.0/eta_pwl;
    
//    printf("%2.8e %2.8e\n", Txz/2.0/eta_pwl, *Eii_pwl * Txz/Tii);

    // Compute dissipative strain rate components
    *Exx_diss =  Exx_pl + Exx_lin +  Exx_pwl + Exx_exp + Exx_gbs + Exx_cst;
    *Ezz_diss =  Ezz_pl + Ezz_lin +  Ezz_pwl + Ezz_exp + Ezz_gbs + Ezz_cst;
    *Exz_diss =  Exz_pl + Exz_lin +  Exz_pwl + Exz_exp + Exz_gbs + Exz_cst;
    
    *Wdiss = Txx*(*Exx_diss) + Tzz*(*Ezz_diss) + Tyy*(-(*Exx_diss + *Ezz_diss)) + 2.0*Txz*(*Exz_diss);
    
//    printf("*Wdiss = %2.2e\n", *Wdiss*scaling->S*scaling->E);
    
    if (*Wdiss<0.0) {
        printf("Negative dissipation = %2.2e // %2.2e %2.2e %2.2e %2.2e\n", *Wdiss*scaling->S*scaling->E,  Txx*(*Exx_diss)*scaling->S*scaling->E, Tzz*(*Ezz_diss)*scaling->S*scaling->E, Tyy*(-(*Exx_diss + *Ezz_diss))*scaling->S*scaling->E, 2.0*Txz*(*Exz_diss)*scaling->S*scaling->E);
        
        printf("is_pl %d --- F_trial = %2.2e --- F_corr = %2.2e --- F_corr1 = %2.2e --- F_corr2 = %2.2e\n", is_pl, F_trial, F_corr, F_corr1, F_corr2);
        
          printf("exx_tot = %2.2e exx_el = %2.2e exx_vis = %2.2e exx_pl = %2.2e, exx_net = %2.2e\n", exx*scaling->E, *Exx_el*scaling->E, (Exx_pwl+Exx_lin+Exx_cst+Exx_exp)*scaling->E, Exx_pl*scaling->E, (exx - (*Exx_el+Exx_pwl+Exx_lin+Exx_cst+Exx_exp+Exx_pl))*scaling->E);
                    printf("ezz_tot = %2.2e ezz_el = %2.2e ezz_vis = %2.2e ezz_pl = %2.2e, ezz_net = %2.2e\n", ezz*scaling->E, *Ezz_el*scaling->E, (Ezz_pwl+Ezz_lin+Ezz_cst+Ezz_exp)*scaling->E, Ezz_pl*scaling->E, (ezz - (*Ezz_el+Ezz_pwl+Ezz_lin+Ezz_cst+Ezz_exp+Ezz_pl))*scaling->E);
            //       printf("exx_el = %2.2e exx_vis = %2.2e exx_pl = %2.2e, exx_net = %2.2e\n", *Exx_el*scaling->E, (Exx_pwl+Exx_lin+Exx_cst+Exx_exp)*scaling->E, Exx_pl*scaling->E, (exx - (*Exx_el+Exx_pwl+Exx_lin+Exx_cst+Exx_exp+Exx_pl))*scaling->E);
                    printf("exz_tot = %2.2e exz_el = %2.2e exz_vis = %2.2e exz_pl = %2.2e, exz_net = %2.2e\n", exz*scaling->E, *Exz_el*scaling->E, (Exz_pwl+Exz_lin+Exz_cst+Exz_exp)*scaling->E, Exz_pl*scaling->E, (exz - (*Exz_el+Exz_pwl+Exz_lin+Exz_cst+Exz_exp+Exz_pl))*scaling->E);

            if (Txx<0.0 && *Exx_diss>0.0) {
                printf("Error 11\n");
            }
            if (Txx>0.0 && *Exx_diss<0.0) {
                printf("Error 12\n");
            }
        
            if (Txz<0.0 && *Exz_diss>0.0) {
                printf("Error 13\n");
            }
        
//        exit(9);
    }
    *Wel   = Txx*(*Exx_el)   + Tzz*(*Ezz_el)   + Tyy*(-(*Exx_el   + *Ezz_el))   + 2.0*Txz*(*Exz_el);
    *Wtot  = Txx*exx         + Tzz*ezz         + Tyy*eyy                        + 2.0*Txz*exz;

    *Eii_el  = (double)elastic*fabs(Tii-Tii0)/2.0/eta_el;
    *Eii_pwl =  Tii/2.0/eta_pwl;

    // Partitioning of volumetric strain
    *div_pl  = divp;
    *div_el  = - (Pc - P0) / (K*dt);
    *div_r   = divr;

    // Viscosity for dissipative processes (no elasticity)
    eta        = 1.0/(inv_eta_diss);//Tii/2.0/Eii_vis;
    *VEcoeff   = eta_ve/eta_el;
    if (elastic==0) *VEcoeff = 0.0;

    // Override viscosity at step 0 (100% visco-plastic)
    if ( model->step == 0 ) *etaVE = eta;
    *txxn = Txx;
    *tzzn = Tzz;
    *txzn = Txz;
    
//    printf("exx_tot = %2.2e exx_el = %2.2e exx_vis = %2.2e exx_pl = %2.2e, exx_net = %2.2e\n", exx*scaling->E, *Exx_el*scaling->E, (Exx_pwl+Exx_lin+Exx_cst+Exx_exp)*scaling->E, Exx_pl*scaling->E, (exx - (*Exx_el+Exx_pwl+Exx_lin+Exx_cst+Exx_exp+Exx_pl))*scaling->E);
//            printf("ezz_tot = %2.2e ezz_el = %2.2e ezz_vis = %2.2e ezz_pl = %2.2e, ezz_net = %2.2e\n", ezz*scaling->E, *Ezz_el*scaling->E, (Ezz_pwl+Ezz_lin+Ezz_cst+Ezz_exp)*scaling->E, Ezz_pl*scaling->E, (ezz - (*Ezz_el+Ezz_pwl+Ezz_lin+Ezz_cst+Ezz_exp+Ezz_pl))*scaling->E);
//    //       printf("exx_el = %2.2e exx_vis = %2.2e exx_pl = %2.2e, exx_net = %2.2e\n", *Exx_el*scaling->E, (Exx_pwl+Exx_lin+Exx_cst+Exx_exp)*scaling->E, Exx_pl*scaling->E, (exx - (*Exx_el+Exx_pwl+Exx_lin+Exx_cst+Exx_exp+Exx_pl))*scaling->E);
//            printf("exz_tot = %2.2e exz_el = %2.2e exz_vis = %2.2e exz_pl = %2.2e, exz_net = %2.2e\n", exz*scaling->E, *Exz_el*scaling->E, (Exz_pwl+Exz_lin+Exz_cst+Exz_exp)*scaling->E, Exz_pl*scaling->E, (exz - (*Exz_el+Exz_pwl+Exz_lin+Exz_cst+Exz_exp+Exz_pl))*scaling->E);
//
//    if (Txx<0.0 && *Exx_diss>0.0) {
//        printf("Error 11\n");
//    }
//    if (Txx>0.0 && *Exx_diss<0.0) {
//        printf("Error 12\n");
//    }
////
////    if (Txz<0.0 && *Exz_diss>0.0) {
////        printf("Error 1\n");
//        printf("Txz = %2.2e, exz = %2.2e exz_diss = %2.2e\n",   Txz*scaling->S, exz*scaling->E, (Exz_pwl+Exz_lin+Exz_cst+Exz_exp)*scaling->E);
//             printf("exx_tot = %2.2e exx_el = %2.2e exx_pwl = %2.2e exx_exp = %2.2e exx_cst = %2.2e exx_exp = %2.2e exx_pl = %2.2e, exx_net = %2.2e\n", exx*scaling->E, *Exx_el*scaling->E, Exx_pwl*scaling->E, Exx_lin*scaling->E, Exx_cst*scaling->E, Exx_exp*scaling->E, Exx_pl*scaling->E, (exx - (*Exx_el+Exx_pwl+Exx_lin+Exx_cst+Exx_exp+Exx_pl))*scaling->E);
//         printf("ezz_tot = %2.2e ezz_el = %2.2e exx_pwl = %2.2e exx_exp = %2.2e exx_cst = %2.2e exx_exp = %2.2e ezz_pl = %2.2e, ezz_net = %2.2e\n", ezz*scaling->E, *Ezz_el*scaling->E, Ezz_pwl*scaling->E, Ezz_lin*scaling->E, Ezz_cst*scaling->E, Ezz_exp*scaling->E, Ezz_pl*scaling->E, (ezz - (*Ezz_el+Ezz_pwl+Ezz_lin+Ezz_cst+Ezz_exp+Ezz_pl))*scaling->E);
//           //       printf("exx_el = %2.2e exx_vis = %2.2e exx_pl = %2.2e, exx_net = %2.2e\n", *Exx_el*scaling->E, (Exx_pwl+Exx_lin+Exx_cst+Exx_exp)*scaling->E, Exx_pl*scaling->E, (exx - (*Exx_el+Exx_pwl+Exx_lin+Exx_cst+Exx_exp+Exx_pl))*scaling->E);
//         printf("exz_tot = %2.2e exz_el = %2.2e exz_pwl = %2.2e exz_exp = %2.2e exz_cst = %2.2e exz_exp = %2.2e exz_pl = %2.2e, exz_net = %2.2e\n", exz*scaling->E, *Exz_el*scaling->E, Exz_pwl*scaling->E, Exz_lin*scaling->E, Exz_cst*scaling->E, Exz_exp*scaling->E, Exz_pl*scaling->E, (exz - (*Exz_el+Exz_pwl+Exz_lin+Exz_cst+Exz_exp+Exz_pl))*scaling->E);
//
//        exit(1); }
//    if (Txz>0.0 && *Exz_diss<0.0) {
//         printf("Error 2\n");
//             printf("exx_tot = %2.2e exx_el = %2.2e exx_pwl = %2.2e exx_exp = %2.2e exx_cst = %2.2e exx_exp = %2.2e exx_pl = %2.2e, exx_net = %2.2e\n", exx*scaling->E, *Exx_el*scaling->E, Exx_pwl*scaling->E, Exx_lin*scaling->E, Exx_cst*scaling->E, Exx_exp*scaling->E, Exx_pl*scaling->E, (exx - (*Exx_el+Exx_pwl+Exx_lin+Exx_cst+Exx_exp+Exx_pl))*scaling->E);
//         printf("ezz_tot = %2.2e ezz_el = %2.2e exx_pwl = %2.2e exx_exp = %2.2e exx_cst = %2.2e exx_exp = %2.2e ezz_pl = %2.2e, ezz_net = %2.2e\n", ezz*scaling->E, *Ezz_el*scaling->E, Ezz_pwl*scaling->E, Ezz_lin*scaling->E, Ezz_cst*scaling->E, Ezz_exp*scaling->E, Ezz_pl*scaling->E, (ezz - (*Ezz_el+Ezz_pwl+Ezz_lin+Ezz_cst+Ezz_exp+Ezz_pl))*scaling->E);
//           //       printf("exx_el = %2.2e exx_vis = %2.2e exx_pl = %2.2e, exx_net = %2.2e\n", *Exx_el*scaling->E, (Exx_pwl+Exx_lin+Exx_cst+Exx_exp)*scaling->E, Exx_pl*scaling->E, (exx - (*Exx_el+Exx_pwl+Exx_lin+Exx_cst+Exx_exp+Exx_pl))*scaling->E);
//         printf("exz_tot = %2.2e exz_el = %2.2e exz_pwl = %2.2e exz_exp = %2.2e exz_cst = %2.2e exz_exp = %2.2e exz_pl = %2.2e, exz_net = %2.2e\n", exz*scaling->E, *Exz_el*scaling->E, Exz_pwl*scaling->E, Exz_lin*scaling->E, Exz_cst*scaling->E, Exz_exp*scaling->E, Exz_pl*scaling->E, (exz - (*Exz_el+Exz_pwl+Exz_lin+Exz_cst+Exz_exp+Exz_pl))*scaling->E);
//    }//exit(1); }

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
        
//    if (P>4.6e9/scaling->S && P<4.7e9/scaling->S) {
//        printf("%2.2e %2.2e\n", eta_ve*scaling->eta, eta_lin*scaling->eta);
//        printf("%2.2e %2.2e\n", *etaVE*scaling->eta, eta_lin*scaling->eta);
//    }

    return eta;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void NonNewtonianViscosityGrid( grid *mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling ) {

    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1, k1, cond;
    double eta, txx1, tzz1, txz1, Pn, Tn, etaVE, VEcoeff=0.0, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1, div_el, div_pl, div_r;
    double exx_pwl, exz_pwl, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, Wtot, Wel, Wdiss;
    int average = model->eta_avg, UnsplitDiffReac = model->UnsplitDiffReac;
    double detadexx, detadezz, detadexz, detadp;
    double Xreac;
    double OverS;
    double ddivpdexx, ddivpdezz, ddivpdexz, ddivpdp, Pcorr, drhodp, rho;

    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;

    InterpCentroidsToVerticesDouble( mesh->div_u,   mesh->div_u_s, mesh, model );
    InterpCentroidsToVerticesDouble( mesh->T,       mesh->T_s,     mesh, model );
    InterpCentroidsToVerticesDouble( mesh->p_in,    mesh->P_s,     mesh, model );
    InterpCentroidsToVerticesDouble( mesh->d0_n,    mesh->d0_s,    mesh, model );
    InterpCentroidsToVerticesDouble( mesh->phi0_n,  mesh->phi0_s,  mesh, model ); // ACHTUNG NOT FRICTION ANGLE

    // Evaluate cell center viscosities
#pragma omp parallel for shared( mesh  ) private( cond, k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, detadexx, detadezz, detadexz, detadp, Xreac, OverS, ddivpdexx, ddivpdezz, ddivpdexz, ddivpdp, Pcorr, drhodp, rho, div_el, div_pl, div_r, Wel, Wdiss, Wtot ) firstprivate( UnsplitDiffReac, materials, scaling, average, model, Ncx, Ncz )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {

        //    for ( l=0; l<Ncz; l++ ) {
        //        for ( k=0; k<Ncx; k++ ) {

        k      = mesh->kp[k1];
        l      = mesh->lp[k1];
        c0     = k  + l*(Ncx);

        mesh->eta_n[c0]       = 0.0;
        mesh->eta_phys_n[c0]  = 0.0;
        mesh->VE_n[c0]        = 0.0;
        mesh->sxxd[c0]        = 0.0;
        mesh->szzd[c0]        = 0.0;
        mesh->sxz_n[c0]       = 0.0;
        mesh->eII_el[c0]      = 0.0;
        mesh->eII_pl[c0]      = 0.0;
        mesh->eII_pwl[c0]     = 0.0;
        mesh->eII_exp[c0]     = 0.0;
        mesh->eII_lin[c0]     = 0.0;
        mesh->eII_gbs[c0]     = 0.0;
        mesh->eII_cst[c0]     = 0.0;
        mesh->d_n[c0]         = 0.0;
        mesh->exx_el[c0]      = 0.0;
        mesh->exx_diss[c0]    = 0.0;
        mesh->ezz_el[c0]      = 0.0;
        mesh->ezz_diss[c0]    = 0.0;
        mesh->detadexx_n[c0]  = 0.0;
        mesh->detadezz_n[c0]  = 0.0;
        mesh->detadgxz_n[c0]  = 0.0;
        mesh->detadp_n[c0]    = 0.0;
        mesh->ddivpdexx_n[c0] = 0.0;
        mesh->ddivpdezz_n[c0] = 0.0;
        mesh->ddivpdgxz_n[c0] = 0.0;
        mesh->ddivpdp_n[c0]   = 0.0;
        mesh->p_corr[c0]      = 0.0;                // Achtung baby
        mesh->div_u_el[c0]    = 0.0;
        mesh->div_u_pl[c0]    = 0.0;
        mesh->div_u_r[c0]     = 0.0;
        mesh->Wtot[c0]        = 0.0;
        mesh->Wel[c0]         = 0.0;
        mesh->Wdiss[c0]       = 0.0;
        //        X                     =  mesh->Xreac_n[c0]; // Save X first
        //        if (model->ProgReac==1) mesh->Xreac_n[c0]    = 0.0;
        if ( UnsplitDiffReac == 0 ) mesh->X_n[c0]        = 0.0;
        mesh->OverS_n[c0]    = 0.0;
        
        
        if ( model->VolChangeReac == 1 ) {
            mesh->rho_n[c0]  = 0.0;
            mesh->drhodp_n[c0] = 0.0;
        }
        
        

        // Loop on grid nodes
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {

            // Loop on phases
            for ( p=0; p<model->Nb_phases; p++) {

                cond =  fabs(mesh->phase_perc_n[p][c0])>1.0e-13;

                if ( cond == 1 ) {
                    eta =  Viscosity( p, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi0_n[c0], mesh->X0_n[c0], mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1, mesh->strain_n[c0], mesh->dil_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_n[c0], mesh->div_u[c0], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel );
                }

                // ARITHMETIC AVERAGE
                if (average == 0) {
                    if ( cond == 1 ) mesh->eta_n[c0]       += mesh->phase_perc_n[p][c0] * etaVE;
                    if ( cond == 1 ) mesh->eta_phys_n[c0]  += mesh->phase_perc_n[p][c0] * eta;
                    if ( cond == 1 ) mesh->detadexx_n[c0]  += mesh->phase_perc_n[p][c0] * detadexx;
                    if ( cond == 1 ) mesh->detadezz_n[c0]  += mesh->phase_perc_n[p][c0] * detadezz;
                    if ( cond == 1 ) mesh->detadgxz_n[c0]  += mesh->phase_perc_n[p][c0] * detadexz/2.0;
                    if ( cond == 1 ) mesh->detadp_n[c0]    += mesh->phase_perc_n[p][c0] * detadp;
                }
                if (average == 0 || average == 2 ) {
                    if ( cond == 1 ) mesh->sxxd[c0]   += mesh->phase_perc_n[p][c0] * txx1;
                    if ( cond == 1 ) mesh->szzd[c0]   += mesh->phase_perc_n[p][c0] * tzz1;
                    if ( cond == 1 ) mesh->sxz_n[c0]  += mesh->phase_perc_n[p][c0] * txz1;
                }
                if ( cond == 1 ) mesh->VE_n[c0]       += mesh->phase_perc_n[p][c0] * VEcoeff;
                if ( cond == 1 ) mesh->eII_el[c0]     += mesh->phase_perc_n[p][c0] * eII_el;
                if ( cond == 1 ) mesh->eII_pl[c0]     += mesh->phase_perc_n[p][c0] * eII_pl;
                if ( cond == 1 ) mesh->eII_pwl[c0]    += mesh->phase_perc_n[p][c0] * eII_pwl;
                if ( cond == 1 ) mesh->eII_exp[c0]    += mesh->phase_perc_n[p][c0] * eII_exp;
                if ( cond == 1 ) mesh->eII_lin[c0]    += mesh->phase_perc_n[p][c0] * eII_lin;
                if ( cond == 1 ) mesh->eII_gbs[c0]    += mesh->phase_perc_n[p][c0] * eII_gbs;
                if ( cond == 1 ) mesh->eII_cst[c0]    += mesh->phase_perc_n[p][c0] * eII_cst;
                if ( cond == 1 ) mesh->d_n[c0]        += mesh->phase_perc_n[p][c0] * 1.0/d1;

                if ( cond == 1 ) mesh->exx_el[c0]     += mesh->phase_perc_n[p][c0] * exx_el;
                if ( cond == 1 ) mesh->exx_diss[c0]   += mesh->phase_perc_n[p][c0] * exx_diss;
                if ( cond == 1 ) mesh->ezz_el[c0]     += mesh->phase_perc_n[p][c0] * ezz_el;
                if ( cond == 1 ) mesh->ezz_diss[c0]   += mesh->phase_perc_n[p][c0] * ezz_diss;
                if ( cond == 1 ) mesh->Wtot[c0]       += mesh->phase_perc_n[p][c0] * Wtot;
                if ( cond == 1 ) mesh->Wdiss[c0]      += mesh->phase_perc_n[p][c0] * Wdiss;
                if ( cond == 1 ) mesh->Wel[c0]        += mesh->phase_perc_n[p][c0] * Wel;
                
                if (cond == 1 && Wdiss<0.0) {printf("negative dissipation: you crazy! --> Wdiss = %2.2e\n", Wdiss*scaling->S*scaling->E); }

                if ( cond == 1 ) mesh->ddivpdexx_n[c0] += mesh->phase_perc_n[p][c0] * ddivpdexx;
                if ( cond == 1 ) mesh->ddivpdezz_n[c0] += mesh->phase_perc_n[p][c0] * ddivpdezz;
                if ( cond == 1 ) mesh->ddivpdgxz_n[c0] += mesh->phase_perc_n[p][c0] * ddivpdexz/2.0;
                if ( cond == 1 ) mesh->ddivpdp_n[c0]   += mesh->phase_perc_n[p][c0] * ddivpdp;
                if ( cond == 1 ) mesh->p_corr[c0]      += mesh->phase_perc_n[p][c0] * Pcorr;
                if ( cond == 1 ) mesh->div_u_el[c0]    += mesh->phase_perc_n[p][c0] * div_el;
                if ( cond == 1 ) mesh->div_u_pl[c0]    += mesh->phase_perc_n[p][c0] * div_pl;
                if ( cond == 1 ) mesh->div_u_r[c0]     += mesh->phase_perc_n[p][c0] * div_r;

                if ( cond == 1 && UnsplitDiffReac == 0) mesh->X_n[c0]         += mesh->phase_perc_n[p][c0] * Xreac;
                if ( cond == 1 ) mesh->OverS_n[c0]     += mesh->phase_perc_n[p][c0] * OverS;
//                if ()

                // HARMONIC AVERAGE
                if (average == 1) {
                    if ( cond == 1 ) mesh->sxxd[c0]       += mesh->phase_perc_n[p][c0] * 1.0/txx1;
                    if ( cond == 1 ) mesh->szzd[c0]       += mesh->phase_perc_n[p][c0] * 1.0/tzz1;
                    if ( cond == 1 ) mesh->sxz_n[c0]      += mesh->phase_perc_n[p][c0] * 1.0/txz1;
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
                }
                
                // Volume changes
                if ( model->VolChangeReac == 1 ) {
                    if ( cond == 1 ) mesh->rho_n[c0]       += mesh->phase_perc_n[p][c0] * rho;
                    if ( cond == 1 ) mesh->drhodp_n[c0]    += mesh->phase_perc_n[p][c0] * drhodp;
                }
                
            }
//             eta =  Viscosity( p, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi0_n[c0], mesh->X0_n[c0], mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1, mesh->strain_n[c0], mesh->dil_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0], &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_n[c0], mesh->div_u[c0], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel );
//            CheckSym( mesh->mu_n, scaling->S, mesh->Nx-1, mesh->Nz-1, "mu_n bfore" );
//            CheckSym( mesh->p_in, scaling->S, mesh->Nx-1, mesh->Nz-1, "p_in bfore" );
//            CheckSym( mesh->d0_n, scaling->L, mesh->Nx-1, mesh->Nz-1, "d bfore" );
//            CheckSym( mesh->phi0_n, 1.0, mesh->Nx-1, mesh->Nz-1, "phi0_n bfore" );
//            CheckSym( mesh->X0_n, 1.0, mesh->Nx-1, mesh->Nz-1, "X0_n bfore" );
//            CheckSym( mesh->exxd, scaling->E, mesh->Nx-1, mesh->Nz-1, "exxd bfore" );
//            CheckSym( mesh->ezzd, scaling->E, mesh->Nx-1, mesh->Nz-1, "ezzd bfore" );
//            CheckSym( mesh->exz_n, scaling->E, mesh->Nx-1, mesh->Nz-1, "exz_n bfore", 1, 1 );
            
//            CheckSym( mesh->sxxd0, 1.0, mesh->Nx-1, mesh->Nz-1, "sxxd0 bfore", 0, 1  );
//            CheckSym( mesh->szzd0, 1.0, mesh->Nx-1, mesh->Nz-1, "szzd0 bfore", 0, 1  );
//            CheckSym( mesh->sxz0_n, 1.0, mesh->Nx-1, mesh->Nz-1, "sxz0_n bfore", 1, 1 );
            
//            CheckSym( mesh->strain_n, 1.0, mesh->Nx-1, mesh->Nz-1, "strain_n bfore", 0, 0  );
//            CheckSym( mesh->dil_n, 1.0, mesh->Nx-1, mesh->Nz-1, "dil_n bfore", 0, 0  );
//            CheckSym( mesh->fric_n, 1.0, mesh->Nx-1, mesh->Nz-1, "fric_n bfore", 0,0  );
//            CheckSym( mesh->C_n, 1.0, mesh->Nx-1, mesh->Nz-1, "C_n bfore", 0, 0  );
            
//            CheckSym( mesh->p0_n, 1.0, mesh->Nx-1, mesh->Nz-1, "p0_n bfore", 0, 0  );
//            CheckSym( mesh->bet_n, 1.0, mesh->Nx-1, mesh->Nz-1, "bet_n bfore", 0, 0  );
//            CheckSym( mesh->div_u, 1.0, mesh->Nx-1, mesh->Nz-1, "div_u bfore", 0, 0  );
//
//            CheckSym( mesh->T, scaling->T, mesh->Nx-1, mesh->Nz-1, "T bfore", 0, 0 );


            mesh->d_n[c0]          = 1.0/mesh->d_n[c0];

            // HARMONIC AVERAGE
            if (average == 1) {
                mesh->sxxd[c0]       = 1.0/mesh->sxxd[c0];
                mesh->szzd[c0]       = 1.0/mesh->szzd[c0];
                mesh->sxz_n[c0]      = 1.0/mesh->sxz_n[c0];
                mesh->eta_n[c0]      = 1.0/mesh->eta_n[c0];
                mesh->eta_phys_n[c0] = 1.0/mesh->eta_phys_n[c0];
                mesh->detadexx_n[c0] *= pow(mesh->eta_n[c0],2.0);
                mesh->detadezz_n[c0] *= pow(mesh->eta_n[c0],2.0);
                mesh->detadgxz_n[c0] *= pow(mesh->eta_n[c0],2.0);
                mesh->detadp_n[c0]   *= pow(mesh->eta_n[c0],2.0);
                if (isinf (mesh->eta_phys_n[c0]) ) {
                    printf("Inf: Problem on cell centers:\n");
                    for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_n[p][c0]);
                    printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n", eta, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0]);
                    printf("flag %d nb part cell = %d cell index = %d\n", mesh->BCp.type[c0],mesh->nb_part_cell[c0], c0);
                    printf("x=%2.2e z=%2.2e\n", mesh->xc_coord[k]*scaling->L/1000.0, mesh->zc_coord[l]*scaling->L/1000.0);
                    exit(1);
                }
                if (isnan (mesh->eta_phys_n[c0]) ) {
                    printf("NaN: Problem on cell centers:\n");
                    printf("ProgReac %d\n", model->ProgReac);
                    for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_n[p][c0]);
                    printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n", eta, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0]);
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

            // ACHTUNG!!!! THIS IS HARD-CODED
            // Anisotropy
            if (model->aniso==1) {
                mesh->sxxd[c0] =  mesh->D11_n[c0]*mesh->exxd[c0] + mesh->D12_n[c0]*mesh->ezzd[c0] +  2.0*mesh->D13_n[c0]*mesh->exz_n[c0];
                mesh->szzd[c0] =  mesh->D21_n[c0]*mesh->exxd[c0] + mesh->D22_n[c0]*mesh->ezzd[c0] +  2.0*mesh->D23_n[c0]*mesh->exz_n[c0];
            }

        }
    }
   // printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
   // printf("Cell centrer rheology updated\n");
   // printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");

    // Calculate vertices viscosity

#pragma omp parallel for shared( mesh  ) private( cond, k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, detadexx, detadezz, detadexz, detadp, Xreac, OverS, ddivpdexx, ddivpdezz, ddivpdexz, ddivpdp, Pcorr, drhodp, rho, div_el, div_pl, div_r, Wtot, Wdiss, Wel ) firstprivate( UnsplitDiffReac, materials, scaling, average, model, Nx, Nz )
    for ( k1=0; k1<Nx*Nz; k1++ ) {

        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;

        mesh->VE_s[c1]       = 0.0;
        mesh->sxz[c1]        = 0.0;
        mesh->eta_phys_s[c1] = 0.0;
        mesh->eta_s[c1]      = 0.0;
        mesh->exz_el[c1]     = 0.0;
        mesh->exz_diss[c1]   = 0.0;
        mesh->detadexx_s[c1] = 0.0;
        mesh->detadezz_s[c1] = 0.0;
        mesh->detadgxz_s[c1] = 0.0;
        mesh->detadp_s[c1]   = 0.0;
        if (UnsplitDiffReac == 0) mesh->X_s[c1]        = 0.0;
        //        X                    = mesh->Xreac_s[c1];
        mesh->OverS_s[c1]    = 0.0;

        if ( mesh->BCg.type[c1] != 30 ) {


            for ( p=0; p<model->Nb_phases; p++) {

                cond = fabs(mesh->phase_perc_s[p][c1])>1.0e-13;

                if ( cond == 1 ) {

                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi0_s[c1], mesh->X0_s[c1], mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1, mesh->strain_s[c1], mesh->dil_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_s[c1], &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_s[c1], mesh->div_u_s[c1], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel );
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
                if ( cond == 1 ) mesh->exz_el[c1]     += mesh->phase_perc_s[p][c1] * exz_el;
                if ( cond == 1 ) mesh->exz_diss[c1]   += mesh->phase_perc_s[p][c1] * exz_diss;
                if ( cond == 1 ) mesh->OverS_s[c1]    += mesh->phase_perc_s[p][c1] * OverS;
                if ( cond == 1 && UnsplitDiffReac == 0) mesh->X_s[c1]        += mesh->phase_perc_s[p][c1] * Xreac;

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
                    printf("Inf: Problem on cell vertices:\n");
                    for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_s[p][c1]);
                    printf("%2.2e %2.2e %2.2e %2.2e %2.2e \n", mesh->mu_s[c1], mesh->exxd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->sxz0[c1]);
                    printf("x=%2.2e z=%2.2e\n", mesh->xg_coord[k]*scaling->L/1000, mesh->zg_coord[l]*scaling->L/1000);
                    exit(1);
                }
                if (isnan (mesh->eta_phys_s[c1]) ) {
                    printf("Nan: Problem on cell vertices:\n");
                    for ( p=0; p<model->Nb_phases; p++) printf("phase %d vol=%2.2e\n", p, mesh->phase_perc_s[p][c1]);
                    printf("%2.2e %2.2e %2.2e %2.2e %2.2e \n", mesh->mu_s[c1],  mesh->exxd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->sxz0[c1]);
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

            // ACHTUNG!!!! THIS IS HARD-CODED
            // Anisotropy
            if (model->aniso==1) {
                //                printf("Stress computed from anisotropy");
                //                printf("%2.2e %2.2e\n",2.0*mesh->D33_s[c1]*mesh->exz[c1], mesh->sxz[c1]);
                mesh->sxz[c1] =  mesh->D31_s[c1]*mesh->exxd_s[c1] + mesh->D32_s[c1]*mesh->ezzd_s[c1] + 2.0*mesh->D33_s[c1]*mesh->exz[c1];
            }

            // if (mesh->exz_diss[c1]>0.0 && mesh->sxz[c1]<0.0) exit(21);
            // if (mesh->exz_diss[c1]<0.0 && mesh->sxz[c1]>0.0) exit(22);


        }
    }
    // printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
    // printf("Cell vert rheology updated\n");
    // printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
    //
    // exit(1);

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

void CohesionFrictionDilationGrid( grid* mesh, markers* particles, mat_prop materials, params model, scale scaling ) {

    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1;
    int average = 0;
    double fric, dil, C, strain_acc, fric0, dil0, C0, mu_strain;
    double dstrain, dfric, dcoh, ddil, *strain_pl;
    int SmoothSoftening = 1;
    int    cent=1, vert=0, prop=1, interp=0;

    Nx = mesh->Nx;
    Nz = mesh->Nz;
    Ncx = Nx-1;
    Ncz = Nz-1;

    // Plastic strain
    strain_pl  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
    P2Mastah( &model, *particles,  particles->strain_pl,     mesh, strain_pl,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);

    // Calculate cell centers cohesion and friction
#pragma omp parallel for shared( mesh, strain_pl ) private( p, c0, strain_acc, fric, dil, C, fric0, dil0, C0, dstrain, ddil, dcoh, dfric, mu_strain ) firstprivate( model, materials, average, Ncx, Ncz )
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {

        // First - initialize to 0
        mesh->fric_n[c0] = 0.0;
        mesh->dil_n[c0]  = 0.0;
        mesh->C_n[c0]    = 0.0;

        // Compute only if below free surface
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {

            // Retrieve accumulated plastic strain
            strain_acc = strain_pl[c0];

            // Loop on phases
            for ( p=0; p<model.Nb_phases; p++) {

                fric = materials.phi[p];
                dil  = materials.psi[p];
                C    = materials.C[p];

                // Apply strain softening
                dstrain   = materials.pls_end[p] - materials.pls_start[p];
                
                if ( SmoothSoftening == 1) {
                
                    mu_strain = 0.5*(materials.pls_end[p] + materials.pls_start[p]);

                    if (materials.phi_soft[p] == 1) {
                        dfric     = materials.phi[p]     - materials.phi_end[p];
                        fric      = materials.phi[p]     - dfric/2.0 *erfc( -(strain_acc - mu_strain) / dstrain ) - materials.phi_end[p];
                        fric0     = materials.phi[p]     - dfric/2.0 *erfc( -(       0.0 - mu_strain) / dstrain ) - materials.phi_end[p];
                        fric      = fric * dfric / fric0 + materials.phi_end[p];
                    }

                    if (materials.psi_soft[p] == 1) {
                        ddil      = materials.psi[p]     - materials.psi_end[p];
                        dil       = materials.psi[p]     - ddil /2.0 *erfc( -(strain_acc - mu_strain) / dstrain ) - materials.psi_end[p];
                        dil0      = materials.psi[p]     - ddil /2.0 *erfc( -(       0.0 - mu_strain) / dstrain ) - materials.psi_end[p];
                        dil       = dil * ddil  / dil0  + materials.psi_end[p];
                    }

                    if (materials.coh_soft[p] == 1) {
                        dcoh      = materials.C[p]       - materials.C_end[p];
                        C         = materials.C[p]       - dcoh /2.0 *erfc( -(strain_acc - mu_strain) / dstrain ) - materials.C_end[p];
                        C0        = materials.C[p]       - dcoh /2.0 *erfc( -(       0.0 - mu_strain) / dstrain ) - materials.C_end[p];
                        C         = C * dcoh / C0 + materials.C_end[p];
                    }
                }
                // Pieciewise linear function
                else {
                    
                    // If we are below the lower strain limit
                    if (strain_acc < materials.pls_start[p]) {
                        if (materials.phi_soft[p] == 1) fric = materials.phi[p];
                        if (materials.psi_soft[p] == 1) dil  = materials.psi[p];
                        if (materials.coh_soft[p] == 1) C    = materials.C[p];
                    }
                    // If we are above the upper strain limit
                    if (strain_acc >= materials.pls_end[p]) {
                        if (materials.phi_soft[p] == 1) fric = materials.phi_end[p];
                        if (materials.psi_soft[p] == 1) dil  = materials.psi_end[p];
                        if (materials.coh_soft[p] == 1) C    = materials.C_end[p];
                    }
                    // If we are in the softening strain range
                    if (strain_acc >= materials.pls_start[p] && strain_acc < materials.pls_end[p] ) {
                        if (materials.phi_soft[p] == 1) fric = materials.phi[p] + (materials.phi_end[p] - materials.phi[p]) / (materials.pls_end[p] - materials.pls_start[p]) *  strain_acc;
                        if (materials.psi_soft[p] == 1) dil  = materials.psi[p] + (materials.psi_end[p] - materials.psi[p]) / (materials.pls_end[p] - materials.pls_start[p]) *  strain_acc;
                        if (materials.coh_soft[p] == 1) C    = materials.C[p]   + (materials.C_end[p]   - materials.C[p]  ) / (materials.pls_end[p] - materials.pls_start[p]) *  strain_acc;
                    }
                    
                }

                // Arithmetic
                if (average ==0) {
                    mesh->fric_n[c0] += mesh->phase_perc_n[p][c0] * fric;
                    mesh->dil_n[c0]  += mesh->phase_perc_n[p][c0] * dil;
                    mesh->C_n[c0]    += mesh->phase_perc_n[p][c0] * C;
                }
                // Harmonic
                if (average == 1) {
                    mesh->fric_n[c0] += mesh->phase_perc_n[p][c0] *  1.0/fric;
                    mesh->dil_n[c0]  += mesh->phase_perc_n[p][c0] *  1.0/dil;
                    mesh->C_n[c0]    += mesh->phase_perc_n[p][c0] *  1.0/C;
                }
                // Geometric
                if (average == 2) {
                    mesh->fric_n[c0] += mesh->phase_perc_n[p][c0] *  log(fric);
                    mesh->dil_n[c0]  += mesh->phase_perc_n[p][c0] *  log(dil);
                    mesh->C_n[c0]    += mesh->phase_perc_n[p][c0] *  log(C);

                }
            }
            // Post-process for geometric/harmonic averages
            if ( average==1 ) mesh->fric_n[c0] = 1.0/mesh->fric_n[c0];
            if ( average==2 ) mesh->fric_n[c0] = exp(mesh->fric_n[c0]);
            if ( average==1 ) mesh->dil_n[c0]  = 1.0/mesh->dil_n[c0];
            if ( average==2 ) mesh->dil_n[c0]  = exp(mesh->dil_n[c0]);
            if ( average==1 ) mesh->C_n[c0]    = 1.0/mesh->C_n[c0];
            if ( average==2 ) mesh->C_n[c0]    = exp(mesh->C_n[c0]);

        }
    }

    // Freedom
    DoodzFree( strain_pl );

    // Plastic strain
    strain_pl  = DoodzCalloc((model.Nx-0)*(model.Nz-0), sizeof(double));
    Interp_P2N ( *particles,  particles->strain_pl, mesh, strain_pl, mesh->xg_coord, mesh->zg_coord, 1, 0, &model );

#pragma omp parallel for shared( mesh, strain_pl ) private( p, c1, strain_acc, fric, dil, C, fric0, dil0, C0, dstrain, ddil, dcoh, dfric, mu_strain ) firstprivate( model, materials, average, Nx, Nz )
    // Calculate vertices cohesion and friction
    for ( c1=0; c1<Nx*Nz; c1++ ) {

        // First - initialize to 0
        mesh->fric_s[c1] = 0.0;
        mesh->dil_s[c1]  = 0.0;
        mesh->C_s[c1]    = 0.0;

        // Compute only if below free surface
        if ( mesh->BCg.type[c1] != 30 ) {

            // Retrieve accumulated plastic strain
            strain_acc = strain_pl[c1];

            // Loop on phases
            for ( p=0; p<model.Nb_phases; p++) {

                fric = materials.phi[p];
                dil  = materials.psi[p];
                C    = materials.C[p];

                // Apply strain softening
                dstrain   = materials.pls_end[p] - materials.pls_start[p];
                
                // Smooth function
                if ( SmoothSoftening == 1) {
                    
                    mu_strain = 0.5*(materials.pls_end[p] + materials.pls_start[p]);

                    if (materials.phi_soft[p] == 1) {
                        dfric     = materials.phi[p]     - materials.phi_end[p];
                        fric      = materials.phi[p]     - dfric/2.0 *erfc( -(strain_acc - mu_strain) / dstrain ) - materials.phi_end[p];
                        fric0     = materials.phi[p]     - dfric/2.0 *erfc( -(       0.0 - mu_strain) / dstrain ) - materials.phi_end[p];
                        fric      = fric * dfric / fric0 + materials.phi_end[p];
                    }

                    if (materials.psi_soft[p] == 1) {
                        ddil      = materials.psi[p]     - materials.psi_end[p];
                        dil       = materials.psi[p]     - ddil /2.0 *erfc( -(strain_acc - mu_strain) / dstrain ) - materials.psi_end[p];
                        dil0      = materials.psi[p]     - ddil /2.0 *erfc( -(       0.0 - mu_strain) / dstrain ) - materials.psi_end[p];
                        dil       = dil * ddil  / dil0  + materials.psi_end[p];
                    }

                    if (materials.coh_soft[p] == 1) {
                        dcoh      = materials.C[p]       - materials.C_end[p];
                        C         = materials.C[p]       - dcoh /2.0 *erfc( -(strain_acc - mu_strain) / dstrain ) - materials.C_end[p];
                        C0        = materials.C[p]       - dcoh /2.0 *erfc( -(       0.0 - mu_strain) / dstrain ) - materials.C_end[p];
                        C         = C * dcoh / C0 + materials.C_end[p];
                    }
                }
                // Pieciewise linear function
                else {
                    
                    // If we are below the lower strain limit
                    if (strain_acc < materials.pls_start[p]) {
                        if (materials.phi_soft[p] == 1) fric = materials.phi[p];
                        if (materials.psi_soft[p] == 1) dil  = materials.psi[p];
                        if (materials.coh_soft[p] == 1) C    = materials.C[p];
                    }
                    // If we are above the upper strain limit
                    if (strain_acc >= materials.pls_end[p]) {
                        if (materials.phi_soft[p] == 1) fric = materials.phi_end[p];
                        if (materials.psi_soft[p] == 1) dil  = materials.psi_end[p];
                        if (materials.coh_soft[p] == 1) C    = materials.C_end[p];
                    }
                    // If we are in the softening strain range
                    if (strain_acc >= materials.pls_start[p] && strain_acc < materials.pls_end[p] ) {
                        if (materials.phi_soft[p] == 1) fric = materials.phi[p] + (materials.phi_end[p] - materials.phi[p]) / (materials.pls_end[p] - materials.pls_start[p]) *  strain_acc;
                        if (materials.psi_soft[p] == 1) dil  = materials.psi[p] + (materials.psi_end[p] - materials.psi[p]) / (materials.pls_end[p] - materials.pls_start[p]) *  strain_acc;
                        if (materials.coh_soft[p] == 1) C    = materials.C[p]   + (materials.C_end[p]   - materials.C[p]  ) / (materials.pls_end[p] - materials.pls_start[p]) *  strain_acc;
                    }
                    
                }
            

                // Arithmetic
                if (average == 0) {
                    mesh->fric_s[c1] += mesh->phase_perc_s[p][c1] * fric;
                    mesh->dil_s[c1]  += mesh->phase_perc_s[p][c1] * dil;
                    mesh->C_s[c1]    += mesh->phase_perc_s[p][c1] * C;
                }
                // Harmonic
                if (average == 1) {
                    mesh->fric_s[c1] += mesh->phase_perc_s[p][c1] *  1.0/fric;
                    mesh->dil_s[c1]  += mesh->phase_perc_s[p][c1] *  1.0/dil;
                    mesh->C_s[c1]    += mesh->phase_perc_s[p][c1] *  1.0/C;
                }
                // Geometric
                if (average == 2) {
                    mesh->fric_s[c1] += mesh->phase_perc_s[p][c1] *  log(fric);
                    mesh->dil_s[c1]  += mesh->phase_perc_s[p][c1] *  log(dil);
                    mesh->C_s[c1]    += mesh->phase_perc_s[p][c1] *  log(C);
                }
            }
            // Post-process for geometric/harmonic averages
            if ( average==1 ) mesh->fric_s[c1] = 1.0/mesh->fric_s[c1];
            if ( average==2 ) mesh->fric_s[c1] = exp(mesh->fric_s[c1]);
            if ( average==1 ) mesh->dil_s[c1]  = 1.0/mesh->dil_s[c1];
            if ( average==2 ) mesh->dil_s[c1]  = exp(mesh->dil_s[c1]);
            if ( average==1 ) mesh->C_s[c1]    = 1.0/mesh->C_s[c1];
            if ( average==2 ) mesh->C_s[c1]    = exp(mesh->C_s[c1]);
        }
    }

    // Freedom
    DoodzFree( strain_pl );

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ShearModCompExpGrid( grid* mesh, mat_prop materials, params model, scale scaling ) {

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
            mesh->mu_n[c0]  = 0.0;
            mesh->bet_n[c0] = 0.0;
            mesh->alp[c0]   = 0.0;

            // Compute only if below free surface
            if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {

                // Loop on phases
                for ( p=0; p<model.Nb_phases; p++) {

                    // Arithmetic
                    if (average == 0) {
                        mesh->mu_n[c0]  += mesh->phase_perc_n[p][c0] * materials.mu[p];
                        mesh->bet_n[c0] += mesh->phase_perc_n[p][c0] * materials.bet[p];
                    }
                    // Harmonic
                    if (average == 1) {
                        mesh->mu_n[c0]  += mesh->phase_perc_n[p][c0] *  1.0/materials.mu[p];
                        mesh->bet_n[c0] += mesh->phase_perc_n[p][c0] *  1.0/materials.bet[p];
                    }
                    // Geometric
                    if (average == 2) {
                        mesh->mu_n[c0]  += mesh->phase_perc_n[p][c0] *  log(materials.mu[p]);
                        mesh->bet_n[c0] += mesh->phase_perc_n[p][c0] *  log(materials.bet[p]);
                    }

                    // Standard arithmetic interpolation
                    mesh->alp  [c0] += mesh->phase_perc_n[p][c0] * materials.alp[p];

                }
                // Post-process for geometric/harmonic averages
                if ( average==1 ) mesh->mu_n[c0] = 1.0/mesh->mu_n[c0];
                if ( average==2 ) mesh->mu_n[c0] = exp(mesh->mu_n[c0]);
                if ( average==1 ) mesh->bet_n[c0] = 1.0/mesh->bet_n[c0];
                if ( average==2 ) mesh->bet_n[c0] = exp(mesh->bet_n[c0]);

            }
        }
    }


    // Calculate vertices shear modulus
    for ( l=0; l<Nz; l++ ) {
        for ( k=0; k<Nx; k++ ) {

            // Vertex index
            c1 = k + l*Nx;

            // First - initialize to 0
            mesh->mu_s[c1]  = 0.0;
            mesh->bet_s[c1] = 0.0;

            // Compute only if below free surface
            if ( mesh->BCg.type[c1] != 30 ) {

                // Loop on phases
                for ( p=0; p<model.Nb_phases; p++) {

                    // Arithmetic
                    if (average == 0) {
                        mesh->mu_s[c1]  += mesh->phase_perc_s[p][c1] * materials.mu[p];
                        mesh->bet_s[c1] += mesh->phase_perc_s[p][c1] * materials.bet[p];
                    }
                    // Harmonic
                    if (average == 1) {
                        mesh->mu_s[c1]  += mesh->phase_perc_s[p][c1] *  1.0/materials.mu[p];
                        mesh->bet_s[c1] += mesh->phase_perc_s[p][c1] *  1.0/materials.bet[p];
                    }
                    // Geometric
                    if (average == 2) {
                        mesh->mu_s[c1]  += mesh->phase_perc_s[p][c1] *  log(materials.mu[p]);
                        mesh->bet_s[c1] += mesh->phase_perc_s[p][c1] *  log(materials.bet[p]);
                    }

                }

                if ( isinf(1.0/mesh->mu_s[c1]) ) {
                    printf("Aaaaargh...!! %2.2e %2.2e ----> ShearModulusCompressibilityExpansivityGrid\n", mesh->phase_perc_s[0][c1], mesh->phase_perc_s[1][c1]);
                }

                // Post-process for geometric/harmonic averages
                if ( average==1 ) mesh->mu_s[c1]  = 1.0/mesh->mu_s[c1];
                if ( average==2 ) mesh->mu_s[c1]  = exp(mesh->mu_s[c1]);
                if ( average==1 ) mesh->bet_s[c1] = 1.0/mesh->bet_s[c1];
                if ( average==2 ) mesh->bet_s[c1] = exp(mesh->bet_s[c1]);
            }
        }
    }

    // Periodic
    double av;
    if (model.isperiodic_x==1) {
        for( l=0; l<Nz; l++) {
            c1 = l*Nx + Nx-1;
            av = 0.5*(mesh->mu_s[c1] + mesh->mu_s[l*Nx]);
            mesh->mu_s[c1] = av; mesh->mu_s[l*Nx] = av;
            av = 0.5*(mesh->bet_s[c1] + mesh->bet_s[l*Nx]);
            mesh->bet_s[c1] = av; mesh->bet_s[l*Nx] = av;
        }
    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// MD6
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

                // T and P dependent density based on EOS
                if ( materials->density_model[p] == 0 ) {
                    rho0   = materials->rho [p];
                    rhop   = rho0;
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
                    rhop   = ((1.0-mesh->X_n[c0])*rho0 + mesh->X_n[c0]*(rho0+drho))*rhop; // Average density based on X
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

                // P-T dependent density
                if ( materials->density_model[p] == 3 ) {

                    rho0   = materials->rho[p];
                    beta   = materials->bet[p];
                    alpha  = materials->alp[p];
                    rhop   = rho0*exp(beta * mesh->p_in[c0] - alpha * mesh->T[c0]);
                }

                // Average density base on phase density and phase volume fraction
                if ( mesh->BCp.type[c0] != 30 ) rhonew += mesh->phase_perc_n[p][c0] * rhop;
            }
        }
        mesh->rho_n[c0]   = rhonew;
    }

    InterpCentroidsToVerticesDouble( mesh->rho_n, mesh->rho_s, mesh, model );

    printf("Updated density fields:\n");
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
    InterpCentroidsToVerticesDouble( mesh->exxd, mesh->exxd_s, mesh, model );
    InterpCentroidsToVerticesDouble( mesh->ezzd, mesh->ezzd_s, mesh, model );
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
    double exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, div_el, div_pl, div_r;
    double Pn = model->Pn , eta;
    int    *dom_mech, ind, mech, loud=0, k;
    double *stress, *visco, t_omp;
    double detadexx, detadezz, detadexz, detadp;
    double Xreac;
    double OverS;
    double ddivpdexx, ddivpdezz, ddivpdexz, ddivpdp, Pcorr, drhodp, rho, Wel, Wdiss, Wtot;

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
                    eta =  Viscosity( k, 0.0, T[ix], Pn, d[iz], 0.0, 0.0, E[iy], E[iy], 0.0, 0.0, 0.0, 0.0, materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss,  &d1, 0.0, materials->psi[k], materials->phi[k], materials->C[k], &detadexx, &detadezz, &detadexz, &detadp, 0.0, &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, 0.0, 0.0, &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel  );

                    //                    if(k==0 && eII_exp*scaling->E>1e-17)printf("eII_exp = %2.2e eII_pwl = %2.2e \n", eII_exp*scaling->E, eII_pwl*scaling->E);

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
    double txx1, tzz1, txz1, etaVE, VEcoeff = 0.0, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1, eta;
    double exx_pwl, exz_pwl, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, exx_pl, exz_pl, div_el, div_pl , div_r;
    int average = model->eta_avg;
    double detadexx, detadezz, detadexz, detadp;
    double Xreac;
    double OverS;
    double ddivpdexx, ddivpdezz, ddivpdexz, ddivpdp, Pcorr, drhodp, rho, Wtot, Wdiss, Wel;


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
#pragma omp parallel for shared( mesh  ) private( eii, cond, k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, etaVE, etaVE_exx, etaVE_ezz, etaVE_exz, etaVE_p, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1, exx_el, exz_el, exx_diss, exz_diss, eta_exx, eta_ezz, eta_exz, eta_p, pert_xx, pert_zz, pert_xz , pert_p, detadexx, detadezz, detadexz, detadp, Xreac, OverS, ddivpdexx, ddivpdezz, ddivpdexz, ddivpdp, Pcorr, drhodp, rho, div_el, div_pl, div_r, Wtot, Wdiss, Wel  ) firstprivate( materials, scaling, average, model, Ncx, Ncz, eps, eps1 )
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

            eii     = sqrt(0.5*pow(mesh->exxd[c0],2) + 0.5*pow(mesh->ezzd[c0],2) + pow(mesh->exz_n[c0],2));
            pert_xx = eps*eii;
            pert_zz = eps*eii;
            pert_xz = eps*eii;
            pert_p  = eps*mesh->p_in[c0];
            pert_p  = eps1;

            // Loop on phases
            for ( p=0; p<model->Nb_phases; p++) {

                cond =  fabs(mesh->phase_perc_n[p][c0])>1.0e-13;

                if ( cond == 1 ) {

                    eta =  Viscosity( p, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi0_n[c0], mesh->X0_n[c0], mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1, mesh->strain_n[c0], mesh->dil_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0],  &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_n[c0], mesh->div_u[c0], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel    );

                    eta =  Viscosity( p, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi0_n[c0], mesh->X0_n[c0], mesh->exxd[c0]+pert_xx, mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling, &txx1, &tzz1, &txz1, &etaVE_exx, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst,  &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1, mesh->strain_n[c0], mesh->dil_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0],  &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_n[c0], mesh->div_u[c0], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel  );

                    eta =  Viscosity( p, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi0_n[c0], mesh->X0_n[c0], mesh->exxd[c0], mesh->ezzd[c0]+pert_zz, mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling,  &txx1, &tzz1, &txz1, &etaVE_ezz, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1, mesh->strain_n[c0], mesh->dil_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0],  &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_n[c0], mesh->div_u[c0], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel  );
                    //
                    eta =  Viscosity( p, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi0_n[c0], mesh->X0_n[c0], mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0]+pert_xz, mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling,  &txx1, &tzz1, &txz1, &etaVE_exz, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1, mesh->strain_n[c0], mesh->dil_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0],  &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_n[c0], mesh->div_u[c0], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel  );
                    //
                    eta =  Viscosity( p, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0]+pert_p, mesh->d0_n[c0], mesh->phi0_n[c0], mesh->X0_n[c0], mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling,  &txx1, &tzz1, &txz1, &etaVE_p  , &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1, mesh->strain_n[c0], mesh->dil_n[c0], mesh->fric_n[c0], mesh->C_n[c0], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_n[c0],  &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_n[c0], mesh->div_u[c0], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel  );

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
#pragma omp parallel for shared( mesh, model ) private( eii, cond, k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, etaVE, etaVE_exx, etaVE_ezz, etaVE_exz, etaVE_p, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, d1s, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, eta_exx, eta_ezz, eta_exz, eta_p, pert_xx, pert_zz, pert_xz , pert_p, detadexx, detadezz, detadexz, detadp, Xreac, OverS, ddivpdexx, ddivpdezz, ddivpdexz, ddivpdp, Pcorr, drhodp, rho, div_el, div_pl, div_r, Wtot, Wdiss, Wel  ) firstprivate( materials, scaling, average, Nx, Nz, eps, eps1  )
    for ( k1=0; k1<Nx*Nz; k1++ ) {

        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;

        mesh->detadexx_s[c1] = 0.0;
        mesh->detadezz_s[c1] = 0.0;
        mesh->detadgxz_s[c1] = 0.0;
        mesh->detadp_s[c1]   = 0.0;

        eta_exx = 0.0;
        eta_ezz = 0.0;
        eta_exz = 0.0;
        eta_p   = 0.0;

        if ( mesh->BCg.type[c1] != 30 ) {

            eii     = sqrt(pow(0.5*mesh->exxd_s[c1],2)+0.5*pow(mesh->ezzd_s[c1],2)+pow(mesh->exz[c1],2));
            pert_xx = eps*eii;
            pert_zz = eps*eii;
            pert_xz = eps*eii;
            pert_p  = eps*mesh->P_s[c1];
            pert_p  = eps1;

            for ( p=0; p<model->Nb_phases; p++) {

                cond = fabs(mesh->phase_perc_s[p][c1])>1.0e-13;

                if ( cond == 1 ) {

                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi0_s[c1], mesh->X0_s[c1], mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1s, mesh->strain_s[c1], mesh->dil_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_s[c1],  &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_s[c1], mesh->div_u_s[c1], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel  );

                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi0_s[c1], mesh->X0_s[c1], mesh->exxd_s[c1]+pert_xx, mesh->ezzd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE_exx, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1s, mesh->strain_s[c1], mesh->dil_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_s[c1],  &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_s[c1], mesh->div_u_s[c1], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel  );
                    //
                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi0_s[c1], mesh->X0_s[c1], mesh->exxd_s[c1], mesh->ezzd_s[c1]+pert_zz, mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE_ezz, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1s, mesh->strain_s[c1], mesh->dil_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_s[c1],  &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_s[c1], mesh->div_u_s[c1], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel  );
                    //
                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi0_s[c1], mesh->X0_s[c1], mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1]+pert_xz, mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE_exz, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1s, mesh->strain_s[c1], mesh->dil_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_s[c1],  &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_s[c1], mesh->div_u_s[c1], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel  );
                    //
                    eta =  Viscosity( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1]+pert_p, mesh->d0_s[c1], mesh->phi0_s[c1], mesh->X0_s[c1], mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE_p, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &d1s, mesh->strain_s[c1], mesh->dil_s[c1], mesh->fric_s[c1], mesh->C_s[c1], &detadexx, &detadezz, &detadexz, &detadp, mesh->p0_s[c1],  &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, mesh->bet_s[c1], mesh->div_u_s[c1], &div_el, &div_pl, &div_r, &Wtot, &Wdiss, &Wel  );
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

void InitialiseDirectorVector (grid* mesh, markers* particles, params* model, mat_prop* materials ) {
    
    int    cent=1, vert=0, prop=1, interp=0;

    //    int Nx, Nz, k, p;
    //    double nx, nz, norm, angle;
    //    Nx = model->Nx;
    //    Nz = model->Nz;
    //
    //#pragma omp parallel for shared ( mesh ) \
    //private ( k, nx, nz, norm, angle )              \
    //firstprivate( model )
    //        for ( k=0; k<Nx*Nz; k++ ) {
    //
    //            mesh->nx_s[k] = 0.0;
    //            mesh->nz_s[k] = 0.0;
    //
    //            if (mesh->BCg.type[k] != 30) {
    //
    //                angle = 0.0;
    //                for ( p=0; p<model->Nb_phases; p++) {
    //                    angle += mesh->phase_perc_s[p][k] * materials->aniso_angle[p];
    //                }
    //                nx            = cos(angle);
    //                nz            = sin(angle);
    //                norm          = sqrt(nx*nx + nz*nz);
    //                nx            = nx/norm;
    //                nz            = nz/norm;
    //                mesh->nx_s[k] = nx;
    //                mesh->nz_s[k] = nz;
    //            }
    //        }
    //
    //    #pragma omp parallel for shared ( mesh ) \
    //    private ( k, nx, nz, norm, angle, p )              \
    //    firstprivate( model )
    //        for ( k=0; k<(Nx-1)*(Nz-1); k++ ) {
    //
    //            mesh->nx_n[k] = 0.0;
    //            mesh->nz_n[k] = 0.0;
    //
    //            if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) {
    //                angle = 0.0;
    //                for ( p=0; p<model->Nb_phases; p++) {
    //                    angle += mesh->phase_perc_n[p][k] * materials->aniso_angle[p];
    //                }
    //                nx            = cos(angle);
    //                nz            = sin(angle);
    //                norm          = sqrt(nx*nx + nz*nz);
    //                nx            = nx/norm;
    //                nz            = nz/norm;
    //                mesh->nx_n[k] = nx;
    //                mesh->nz_n[k] = nz;
    //            }
    //        }

    int k;
    double angle, norm;

#pragma omp parallel for shared( particles ) private( angle, norm )
    for (k=0; k<particles->Nb_part; k++) {

        if ( particles->phase[k] != -1 ) {

            // Set up director vector
            angle             = materials->aniso_angle[particles->phase[k]];
            particles->nx[k]  = cos(angle);
            particles->nz[k]  = sin(angle);
            norm              = sqrt(particles->nx[k]*particles->nx[k] + particles->nz[k]*particles->nz[k]);
            particles->nx[k] /= norm;
            particles->nz[k] /= norm;
        }
    }

    P2Mastah( model, *particles, particles->nx,     mesh, mesh->nx0_n,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
    P2Mastah( model, *particles, particles->nz,     mesh, mesh->nz0_n,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
    P2Mastah( model, *particles, particles->nx,     mesh, mesh->nx0_s,   mesh->BCg.type,  1, 0, interp, vert, model->itp_stencil);
    P2Mastah( model, *particles, particles->nz,     mesh, mesh->nz0_s,   mesh->BCg.type,  1, 0, interp, vert, model->itp_stencil);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void NormalizeDirector( grid* mesh, DoodzFP* nx_n, DoodzFP* nz_n, DoodzFP* nx_s, DoodzFP* nz_s, params *model ) {

    int Nx, Nz, k;
    double nx, nz, norm;
    Nx = model->Nx;
    Nz = model->Nz;

#pragma omp parallel for shared ( mesh, nx_n, nz_n ) \
private ( k, nx, nz, norm )              \
firstprivate( model )
    for ( k=0; k<(Nx-1)*(Nz-1); k++ ) {
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) {
            nx            = nx_n[k];
            nz            = nz_n[k];
            norm          = sqrt(nx*nx + nz*nz);
            nx            = nx/norm;
            nz            = nz/norm;
            nx_n[k]       = nx;
            nz_n[k]       = nz;
        }
    }

#pragma omp parallel for shared ( mesh, nx_s, nz_s ) \
private ( k, nx, nz, norm )              \
firstprivate( model )
    for ( k=0; k<Nx*Nz; k++ ) {
        if (mesh->BCg.type[k] != 30) {
            nx            = nx_s[k];
            nz            = nz_s[k];
            norm          = sqrt(nx*nx + nz*nz);
            nx            = nx/norm;
            nz            = nz/norm;
            nx_s[k]       = nx;
            nz_s[k]       = nz;
        }
    }

}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

//void ComputeIncrementsOnParticles( grid* mesh, markers* particles, params* model, mat_prop* materials, scale *scaling ) {
//
//    int k;
//    double *txxm0, *tzzm0, *txzm0, *nxm0, *nzm0, *dm0, *Xm0, *phim0, *Pm0, *Tm0, *divthm0, *rhom0;
//    int Nx = mesh->Nx;
//    int Nz = mesh->Nz;
//
////    CheckSym( mesh->p0_n, 1.0, mesh->Nx-1, mesh->Nz-1, "mesh->p0_n ComputeIncrementsOnParticles", 0, 1  );
////    CheckSym( mesh->p0_n, 1.0, mesh->Nx-1, mesh->Nz-1, "mesh->p0_n ComputeIncrementsOnParticles", 0, 1  );
//
//
//
//    txxm0   = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//    tzzm0   = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//    txzm0   = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//    dm0     = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//    Xm0     = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//    phim0   = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//    Pm0     = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//    Tm0     = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//    divthm0 = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//    rhom0   = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//
//    // Interpolate old fields to new marker locations
//    Interp_Grid2P_centroids2( *particles, txxm0, mesh, mesh->sxxd0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model );
//    Interp_Grid2P_centroids2( *particles, tzzm0, mesh, mesh->szzd0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model );
//    Interp_Grid2P(           *particles, txzm0, mesh, mesh->sxz0 , mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type        );
//
//    if ( model->aniso == 1 ) {
//        nxm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//        nzm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//        Interp_Grid2P_centroids2( *particles, nxm0, mesh, mesh->nx0_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
//        Interp_Grid2P_centroids2( *particles, nzm0, mesh, mesh->nz0_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
//    }
//
//    Interp_Grid2P_centroids2( *particles, dm0, mesh, mesh->d0_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model );
//    Interp_Grid2P_centroids2( *particles, Xm0, mesh, mesh->X0_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model );
//    Interp_Grid2P_centroids2( *particles, phim0, mesh, mesh->phi0_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model );
//    Interp_Grid2P_centroids2( *particles, rhom0, mesh, mesh->rho0_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model );
//    Interp_Grid2P_centroids2( *particles, divthm0, mesh, mesh->divth0_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model );
//    Interp_Grid2P_centroids2( *particles, Pm0, mesh, mesh->p0_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model );
//    Interp_Grid2P_centroids2( *particles, Tm0, mesh, mesh->T0_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model );
//
//    // Compute total changes on markers
//#pragma omp parallel for shared ( particles, txxm0, tzzm0, txzm0, nxm0, nzm0, dm0, Xm0, phim0, Pm0, Tm0, divthm0 ) \
//private ( k )              \
//firstprivate( model, materials )
//    for ( k=0; k<particles->Nb_part; k++ ) {
//
//        if (particles->phase[k] != -1) {
//
//
//            if (model->StressUpdate==0) {
//            particles->dsxxd[k]  =  particles->sxxd[k] - txxm0[k];
//            particles->dszzd[k]  =  particles->szzd[k] - tzzm0[k];
//            }
//            if (model->StressUpdate==1) {
//                particles->dsxxd[k]  =  particles->sxxd[k] - (txxm0[k]-Pm0[k]);
//                particles->dszzd[k]  =  particles->szzd[k] - (tzzm0[k]-Pm0[k]);
//                particles->dsyy[k]   =  particles->syy[k] - ( -(txxm0[k]+tzzm0[k])-Pm0[k]);
//            }
//            particles->dsxz[k]   =  particles->sxz[k]  - txzm0[k];
//            if ( model->aniso == 1 ) particles->dnx[k] = particles->nx[k] - nxm0[k];
//            if ( model->aniso == 1 ) particles->dnz[k] = particles->nz[k] - nzm0[k];
//            particles->dd[k]     =  particles->d[k]     - dm0[k];
//            particles->dX[k]     =  particles->X[k]     - Xm0[k];
//            particles->dphi[k]   =  particles->phi[k]   - phim0[k];
//            particles->dP[k]     =  particles->P[k]     - Pm0[k];
//            particles->dT[k]     =  particles->T[k]     - Tm0[k];
//            // On the fly - compute thermal divergence rate using total dT
//            particles->divth[k]  = materials->alp[particles->phase[k]] * particles->dT[k] / model->dt;
//            particles->ddivth[k] =  particles->divth[k] - divthm0[k];
//            particles->drho[k]   =  particles->rho[k]   - rhom0[k];
//
//        }
//    }
//
//    DoodzFree(txxm0);
//    DoodzFree(tzzm0);
//    DoodzFree(txzm0);
//    DoodzFree(dm0);
//    DoodzFree(Xm0);
//    DoodzFree(phim0);
//    DoodzFree(Tm0);
//    DoodzFree(divthm0);
//    DoodzFree(Pm0);
//    DoodzFree(rhom0);
//
//    if ( model->aniso == 1 ) {
//        DoodzFree(nxm0);
//        DoodzFree(nzm0);
//    }
//
//}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

//void UpdateGridFields( grid* mesh, markers* particles, params* model, mat_prop* materials, scale *scaling ) {
//
//    int Nx = mesh->Nx;
//    int Nz = mesh->Nz;
//    int Ncx = mesh->Nx-1;
//    int Ncz = mesh->Nz-1;
//    int c0, k1, l, k;
//    int    cent=1, vert=0, prop=1, interp=0;
//
//    double *dP     = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//    P2Mastah( model, *particles, particles->dP,     mesh, dP,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//    ArrayPlusArray(  mesh->p0_n,     dP,     Ncx*Ncz );
//
//    // Elasticity - interpolate advected/rotated stresses
//    if  ( model->iselastic == 1 ) {
//
//        if (model->StressUpdate==0) {
//            double *dsxxd = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//            double *dszzd = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//            double *dsxz  = DoodzCalloc ( Nx*Nz,   sizeof(double));
//            P2Mastah( model, *particles, particles->dsxxd,     mesh, dsxxd,  mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//            P2Mastah( model, *particles, particles->dszzd,     mesh, dszzd,  mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//            P2Mastah( model, *particles, particles->dsxz,      mesh, dsxz,   mesh->BCg.type,  1, 0, interp, vert, model->itp_stencil);
//
//            ArrayPlusArray(  mesh->sxxd0,  dsxxd, Ncx*Ncz );
//            ArrayPlusArray(  mesh->szzd0,  dszzd, Ncx*Ncz );
//            ArrayPlusArray(  mesh->sxz0,   dsxz,   Nx*Nz );
//            DoodzFree(dsxxd);
//            DoodzFree(dszzd);
//            DoodzFree(dsxz );
//        }
//        if (model->StressUpdate==1) {
//            double *dsxxd = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//            double *dszzd = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//            double *dsxz  = DoodzCalloc ( Nx*Nz,   sizeof(double));
//            double *dsyy = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//
//            P2Mastah( model, *particles, particles->dsxxd,     mesh, dsxxd,  mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//            P2Mastah( model, *particles, particles->dszzd,     mesh, dszzd,  mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//            P2Mastah( model, *particles, particles->dsxz,      mesh, dsxz,   mesh->BCg.type,  1, 0, interp, vert, model->itp_stencil);
//            P2Mastah( model, *particles, particles->dsyy,      mesh, dsyy,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//
//                for ( k1=0; k1<Ncx*Ncz; k1++ ) {
//                    k  = mesh->kp[k1];
//                    l  = mesh->lp[k1];
//                    c0 = k  + l*(Nx-1);
//
//                    if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
//                        dP[c0]         = -1.0/3.0*(dsxxd[c0] + dsyy[c0] + dszzd[c0] );
//                        dsxxd[c0]      =  dP[c0] + dsxxd[c0];
//                        dszzd[c0]      =  dP[c0] + dszzd[c0];
//                    }
//                }
//
//            ArrayPlusArray(  mesh->sxxd0,  dsxxd, Ncx*Ncz );
//            ArrayPlusArray(  mesh->szzd0,  dszzd, Ncx*Ncz );
//            ArrayPlusArray(  mesh->sxz0,   dsxz,   Nx*Nz );
//            DoodzFree(dsxxd);
//            DoodzFree(dszzd);
//            DoodzFree(dsxz );
//            DoodzFree(dsyy );
//        }
//    }
//
//    // Interpolate to necessary locations for rheological computations...
//    InterpCentroidsToVerticesDouble( mesh->sxxd0, mesh->sxxd0_s, mesh, model );
//    InterpCentroidsToVerticesDouble( mesh->szzd0, mesh->szzd0_s, mesh, model );
//    InterpVerticesToCentroidsDouble( mesh->sxz0_n,  mesh->sxz0,  mesh, model );
//
//    // Director vector
//    if (model->aniso == 1 ) {
//        double *dnx_n = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//        double *dnz_n = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//        double *dnx_s = DoodzCalloc ( Nx*Nz,   sizeof(double));
//        double *dnz_s = DoodzCalloc ( Nx*Nz,   sizeof(double));
//        P2Mastah( model, *particles, particles->dnx,     mesh, dnx_n,  mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//        P2Mastah( model, *particles, particles->dnz,     mesh, dnz_n,  mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//        P2Mastah( model, *particles, particles->dnx,     mesh, dnx_s,  mesh->BCg.type,  1, 0, interp, vert, model->itp_stencil);
//        P2Mastah( model, *particles, particles->dnz,     mesh, dnz_s,  mesh->BCg.type,  1, 0, interp, vert, model->itp_stencil);
//        ArrayPlusArray(  mesh->nx0_n,  dnx_n, Ncx*Ncz );
//        ArrayPlusArray(  mesh->nx0_s,  dnx_s,   Nx*Nz );
//        ArrayPlusArray(  mesh->nz0_n,  dnz_n, Ncx*Ncz );
//        ArrayPlusArray(  mesh->nz0_s,  dnz_s,   Nx*Nz );
//        NormalizeDirector( mesh, mesh->nx0_n, mesh->nz0_n, mesh->nx0_s, mesh->nz0_s, model );
//        DoodzFree(dnx_n);
//        DoodzFree(dnz_n);
//        DoodzFree(dnx_s);
//        DoodzFree(dnz_s);
//    }
//
//    double *dT     = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//    double *ddivth = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//    double *dX_n   = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//    double *dX_s   = DoodzCalloc (   Nx*Nz, sizeof(double));
//    double *dd     = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//    double *dphi   = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//    double *drho   = DoodzCalloc ( Ncx*Ncz, sizeof(double));
//
//    P2Mastah( model, *particles, particles->dT,     mesh, dT,     mesh->BCp.type,  1, 0, interp, cent, 1);
//    P2Mastah( model, *particles, particles->ddivth, mesh, ddivth, mesh->BCp.type,  1, 0, interp, cent, 1);
//    P2Mastah( model, *particles, particles->dd,     mesh, dd,     mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//    P2Mastah( model, *particles, particles->dphi,   mesh, dphi,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//    P2Mastah( model, *particles, particles->drho,   mesh, drho,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//    P2Mastah( model, *particles, particles->dX,     mesh, dX_n,   mesh->BCp.type,  1, 0, interp, cent, model->itp_stencil);
//    P2Mastah( model, *particles, particles->dX,     mesh, dX_s,   mesh->BCg.type,  1, 0, interp, vert, model->itp_stencil);
//
//    ArrayPlusArray(  mesh->T0_n,     dT,     Ncx*Ncz );
//    ArrayPlusArray(  mesh->divth0_n, ddivth, Ncx*Ncz );
//    ArrayPlusArray(  mesh->X0_n,     dX_n,   Ncx*Ncz );
//    ArrayPlusArray(  mesh->X0_s,     dX_s,   Nx*Nz );
//    ArrayPlusArray(  mesh->d0_n,     dd,     Ncx*Ncz );
//    ArrayPlusArray(  mesh->phi0_n,   dphi,   Ncx*Ncz );
//    ArrayPlusArray(  mesh->rho0_n,   drho,   Ncx*Ncz );
//
//    if ( model->Plith_trick == 1 ) ArrayPlusArray(  mesh->p0_n,   mesh->p_lith0, Ncx*Ncz ); // Add back lithostatic component
//
//    DoodzFree( dT );
//    DoodzFree( ddivth );
//    DoodzFree( dX_n );
//    DoodzFree( dX_s );
//    DoodzFree( dd );
//    DoodzFree( dphi );
//    DoodzFree( dP );
//    DoodzFree( drho );
//
////    // Get T and dTdt from previous step from particles
////    Interp_P2C ( particles, particles.T,     &mesh, mesh.T0_n,    mesh.xg_coord, mesh.zg_coord,  1, 0 );
////    Interp_P2C ( particles, particles.divth, &mesh, mesh.divth0_n, mesh.xg_coord, mesh.zg_coord,  1, 0 );
////
////    Interp_P2C ( particles, particles.X, &mesh, mesh.X0_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
////    Interp_P2N ( particles, particles.X, &mesh, mesh.X0_s, mesh.xg_coord, mesh.zg_coord, 1, 0, &model );
////
////    // Interpolate Grain size
////    Interp_P2C ( particles,   particles.d, &mesh, mesh.d0_n,   mesh.xg_coord, mesh.zg_coord, 1, 0 );
////    ArrayEqualArray(  mesh.d_n,  mesh.d0_n, Ncx*Ncz );
////
////    // Interpolate Melt fraction
////    Interp_P2C ( particles, particles.phi, &mesh, mesh.phi0_n,  mesh.xg_coord, mesh.zg_coord, 1, 0 );
////
////    //-----------------------------------------------------------------------------------------------------------
////    // Interp P --> p0_n , p0_s
////    Interp_P2C ( particles, particles.P, &mesh, mesh.p0_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
////    ArrayPlusArray(  mesh.p0_n,   mesh.p_lith0, Ncx*Ncz ); // Add back lithostatic component
//
//
//}
