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

/* 'isout' deactivates particles that are outside of the domain, phase changes to '-1' */
void isout( markers *particles, params model ) {

    int k, count=0;

    // x periodic case
    if ( model.isperiodic_x) {
#pragma omp parallel for shared ( particles )    \
private ( k )                        \
firstprivate ( model ) reduction (+:count) //schedule( static )

        for(k=0; k<particles->Nb_part; k++) {
            if ( particles->x[k] < model.xmin ) {
                // Correct position in x
                particles->x[k] = model.xmax - ABSV(model.xmin - particles->x[k]);
            }
            if ( particles->x[k] > model.xmax ) {
                // Correct position in x
                particles->x[k] = model.xmin + ABSV(model.xmax - particles->x[k]);
            }
            if ( particles->z[k] < model.zmin || particles->z[k] > model.zmax ) {
                particles->phase[k] = -1;
                count++;
            }
        }
    }

    else {
        // General statement
#pragma omp parallel for shared ( particles )    \
private ( k )                        \
firstprivate ( model ) reduction (+:count) //schedule( static )

        for(k=0; k<particles->Nb_part; k++) {
            if (particles->x[k] < model.xmin || particles->x[k] > model.xmax || particles->z[k] < model.zmin || particles->z[k] > model.zmax) {
                particles->phase[k] = -1;
                count++;
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void isoutPart( markers *particles, params *model, int k ) {

     if ( model->isperiodic_x) {
         if ( particles->x[k] < model->xmin ) {
             // Correct position in x
             particles->x[k] = model->xmax - ABSV(model->xmin - particles->x[k]);
         }
         if ( particles->x[k] > model->xmax ) {
             // Correct position in x
             particles->x[k] = model->xmin + ABSV(model->xmax - particles->x[k]);
         }
         if ( particles->z[k] < model->zmin || particles->z[k] > model->zmax ) {
             particles->phase[k] = -1;
         }
     }
     else {
         if (particles->x[k] < model->xmin || particles->x[k] > model->xmax || particles->z[k] < model->zmin || particles->z[k] > model->zmax) {
             particles->phase[k] = -1;
         }
     }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void RogerGunther( markers *particles, params model, grid mesh, int precise, scale scaling ) {

    DoodzFP *VxA, *VxB, *VxC, *VxD;
    DoodzFP *VzA, *VzB, *VzC, *VzD;
    DoodzFP *OmA, *OmB, *OmC, *OmD, *om_n;
    DoodzFP *xA, *zA;
    int k, k1, l, c1, c3, Nb_part = particles->Nb_part;
    clock_t t_omp = (double)omp_get_wtime();


    VxA = DoodzCalloc( Nb_part,sizeof(DoodzFP));
    VzA = DoodzCalloc( Nb_part,sizeof(DoodzFP));

    VxB = DoodzCalloc( Nb_part,sizeof(DoodzFP));
    VzB = DoodzCalloc( Nb_part,sizeof(DoodzFP));

    VxC = DoodzCalloc( Nb_part,sizeof(DoodzFP));
    VzC = DoodzCalloc( Nb_part,sizeof(DoodzFP));

    VxD = DoodzCalloc( Nb_part,sizeof(DoodzFP));
    VzD = DoodzCalloc( Nb_part,sizeof(DoodzFP));

    xA  = DoodzCalloc( Nb_part,sizeof(DoodzFP));
    zA  = DoodzCalloc( Nb_part,sizeof(DoodzFP));

    // Caculate rotation rate of the stress tensor
    if ( model.iselastic == 1 ) {

        om_n = DoodzMalloc ((model.Nx)*(model.Nz)*sizeof(double));
        OmA  = DoodzCalloc( Nb_part,sizeof(DoodzFP));
        OmB  = DoodzCalloc( Nb_part,sizeof(DoodzFP));
        OmC  = DoodzCalloc( Nb_part,sizeof(DoodzFP));
        OmD  = DoodzCalloc( Nb_part,sizeof(DoodzFP));

#pragma omp parallel for shared ( mesh, om_n ) \
private ( k, l, k1, c1, c3 )                            \
firstprivate( model ) // schedule( static )
        for ( k1=0; k1<model.Nx*model.Nz; k1++ ) {
            k  = mesh.kn[k1];
            l  = mesh.ln[k1];

            //        for (k=0; k<model.Nx; k++) {
            //            for (l=0; l<model.Nz; l++) {
            c1 = k + l*model.Nx;
            c3 = k + l*(model.Nx+1);
            om_n[c1] = -(mesh.v_in[c3+1] - mesh.v_in[c3])/model.dx + (mesh.u_in[c1+model.Nx] - mesh.u_in[c1])/model.dz;
            om_n[c1] = 0.5*om_n[c1];
            //            om_n[c1] = model.EpsBG/2.0;
            //}
        }
        Interp_Grid2P( *particles, OmA, &mesh, om_n, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
    }

    // Print2DArrayDouble( om_n, model.Nx, model.Nz, scaling.t );

    // Initial position save
    ArrayEqualArray( xA, particles->x, particles->Nb_part );
    ArrayEqualArray( zA, particles->z, particles->Nb_part );

    // Initial velocity save
    ArrayEqualArray( VxA, particles->Vx, particles->Nb_part );
    ArrayEqualArray( VzA, particles->Vz, particles->Nb_part );

    // Calculate Runge-Kutta velocity (4th order)
    if ( precise == 1 ) {

#pragma omp parallel for shared ( particles, xA, zA, VxA, VzA ) \
private ( k )                            \
firstprivate( Nb_part, model ) // schedule( static )
        for(k=0; k<Nb_part; k++) {

            // Marker shoot #1
            if (particles->phase[k] != -1) {
                particles->x[k] = xA[k] + 0.5 * model.dt * VxA[k];
                particles->z[k] = zA[k] + 0.5 * model.dt * VzA[k];
            }
        }

        // Check if particles are outside of the box
        isout( particles, model );

        //-----------------------------------------------------------------------------------------------------------------------

        // Get the velocity after dt/2
        //        Interp_Grid2P( *particles, VxB, &mesh, mesh.u_in, mesh.xg_coord,  mesh.zvx_coord, mesh.Nx,   mesh.Nz+1, mesh.BCu.type );
        //        Interp_Grid2P( *particles, VzB, &mesh, mesh.v_in, mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz, mesh.BCv.type   );

//        VelocitiesToParticles( &mesh, particles, VxB, VzB, model, scaling );
        if (model.iselastic == 1) Interp_Grid2P( *particles, OmB, &mesh, om_n, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );

        // Get the velocity after dt/2
        Interp_Grid2P( *particles, VxB, &mesh, mesh.u_in, mesh.xg_coord,  mesh.zvx_coord, mesh.Nx,   mesh.Nz+1, mesh.BCu.type );
        Interp_Grid2P( *particles, VzB, &mesh, mesh.v_in, mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz, mesh.BCv.type   );

        if ( model.RK == 4 ) {

            // Get xB and zB
#pragma omp parallel for shared ( particles, xA, zA, VxB, VzB ) \
private ( k )                            \
firstprivate( Nb_part, model )  //schedule( static )

            for(k=0; k<Nb_part; k++) {
                // Marker shoot #2
                if (particles->phase[k] != -1) {
                    particles->x[k] = xA[k] + 0.5 * model.dt * VxB[k];
                    particles->z[k] = zA[k] + 0.5 * model.dt * VzB[k];
                }
            }

            // Check if particles are outside of the box
            isout( particles, model );

            //-----------------------------------------------------------------------------------------------------------------------

            // Get the velocity after dt/2
            //            Interp_Grid2P( *particles, VxC, &mesh, mesh.u_in, mesh.xg_coord,  mesh.zvx_coord, mesh.Nx,   mesh.Nz+1,  mesh.BCu.type );
            //            Interp_Grid2P( *particles, VzC, &mesh, mesh.v_in, mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz,  mesh.BCv.type   );
//            VelocitiesToParticles( &mesh, particles, VxC, VzC, model, scaling );

            // Get the velocity after dt/2
            Interp_Grid2P( *particles, VxC, &mesh, mesh.u_in, mesh.xg_coord,  mesh.zvx_coord, mesh.Nx,   mesh.Nz+1,  mesh.BCu.type );
            Interp_Grid2P( *particles, VzC, &mesh, mesh.v_in, mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz,  mesh.BCv.type   );

            if (model.iselastic == 1) Interp_Grid2P( *particles, OmC, &mesh, om_n, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );

            // Get xC and zC
#pragma omp parallel for shared ( particles, xA, zA, VxC, VzC ) \
private ( k )                            \
firstprivate( Nb_part, model )  //schedule( static )

            for(k=0; k<Nb_part; k++) {
                // Marker shoot #3
                if (particles->phase[k] != -1) {
                    particles->x[k] = xA[k] + 1.0 * model.dt * VxC[k];
                    particles->z[k] = zA[k] + 1.0 * model.dt * VzC[k];
                }
            }

            // Check if particles are outside of the box
            isout( particles, model );

            //-----------------------------------------------------------------------------------------------------------------------

            // Get the velocity after dt
            //            Interp_Grid2P( *particles, VxD, &mesh, mesh.u_in, mesh.xg_coord,  mesh.zvx_coord, mesh.Nx,   mesh.Nz+1, mesh.BCu.type );
            //            Interp_Grid2P( *particles, VzD, &mesh, mesh.v_in, mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz,  mesh.BCv.type   );
//            VelocitiesToParticles( &mesh, particles, VxD, VzD, model, scaling );

            // Get the velocity after dt
            Interp_Grid2P( *particles, VxD, &mesh, mesh.u_in, mesh.xg_coord,  mesh.zvx_coord, mesh.Nx,   mesh.Nz+1, mesh.BCu.type );
            Interp_Grid2P( *particles, VzD, &mesh, mesh.v_in, mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz,  mesh.BCv.type   );

            if (model.iselastic == 1) Interp_Grid2P( *particles, OmD, &mesh, om_n, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
            //-----------------------------------------------------------------------------------------------------------------------

        }

        // Calculate Roger-Gunther velocity
#pragma omp parallel for shared ( particles, VxA, VzA ,VxB, VzB, VxC, VzC, VxD, VzD, OmA, OmB, OmC, OmD) \
private ( k )                              \
firstprivate( Nb_part, model ) //schedule( static )

        for(k=0; k<Nb_part; k++) {

            // RK2
            if ( model.RK == 2 ) {
                if (particles->phase[k] != -1) {
                    VxA[k] = 0.5 * (VxA[k] +  VxB[k]);
                    VzA[k] = 0.5 * (VzA[k] +  VzB[k]);
                    if ( model.iselastic == 1 ) OmA[k] = 0.5 * (OmA[k] +  OmB[k]);
                }
            }

            // RK4
            if ( model.RK == 4 ) {
                if (particles->phase[k] != -1) {
                    VxA[k] = (1.0/6.0) * ( VxA[k] + 2.0 * VxB[k] + 2.0 * VxC[k] + VxD[k]);
                    VzA[k] = (1.0/6.0) * ( VzA[k] + 2.0 * VzB[k] + 2.0 * VzC[k] + VzD[k]);
                    if ( model.iselastic == 1 ) OmA[k] = (1.0/6.0) * ( OmA[k] + 2.0 * OmB[k] + 2.0 * OmC[k] + OmD[k]);
                }
            }
        }

    }

    // Regular first order in time
#pragma omp parallel for shared ( particles, VxA, VzA, xA, zA )    \
private ( k )                              \
firstprivate( Nb_part, model ) //schedule ( static )
    for(k=0; k<Nb_part; k++) {
        if (particles->phase[k] != -1) {
            particles->x[k]    = xA[k] + model.dt * VxA[k];
            particles->z[k]    = zA[k] + model.dt * VzA[k];
        }
    }
    //    }

    // Check if particles are outside of the box
    isout( particles, model );

    DoodzFree(VxA);
    DoodzFree(VzA);
    DoodzFree(VzB);
    DoodzFree(VxB);
    DoodzFree(VzC);
    DoodzFree(VxC);
    DoodzFree(VzD);
    DoodzFree(VxD);
    DoodzFree(xA);
    DoodzFree(zA);

    if ( model.iselastic == 1 ) {
        DoodzFree(om_n);
        DoodzFree(OmA);
        DoodzFree(OmB);
        DoodzFree(OmC);
        DoodzFree(OmD);
    }

    printf("** Time for Roger Gunther = %lf sec\n",  (double)((double)omp_get_wtime() - t_omp) );

}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void RogerGuntherII( markers *particles, params model, grid mesh, int precise, scale scaling ) {

    DoodzFP VxA, VxB, VxC, VxD;
    DoodzFP VzA, VzB, VzC, VzD;
    DoodzFP OmA, OmB, OmC, OmD, *om_s;
    DoodzFP xA, zA;
    int k, k1, l, c0, c1, c2, c3, Nb_part = particles->Nb_part;
    clock_t t_omp = (double)omp_get_wtime();
    double dx, dz;
    double txx, tzz, txz, angle;
    double *dudx_n, *dvdz_n, *dudz_s, *dvdx_s;
    double dudxA, dvdzA, dudzA, dvdxA, dudxB, dvdzB, dudzB, dvdxB, dudxC, dvdzC, dudzC, dvdxC, dudxD, dvdzD, dudzD, dvdxD, VEA,VEB,VEC,VED;
    double nx, nz, ndotx, ndotz, w12, exxA, exzA, exxB, exzB, exxC, exzC, exxD, exzD, exx, exz, norm;

    int new = 1;
    dx = mesh.dx;
    dz = mesh.dz;

    // Caculate rotation rate of the stress tensor
    if ( model.iselastic == 1 || model.aniso == 1 ) {

        om_s   = DoodzCalloc ((model.Nx-0)*(model.Nz-0),sizeof(double));
        dudx_n = DoodzCalloc ((model.Nx-1)*(model.Nz-1),sizeof(double));
        dvdz_n = DoodzCalloc ((model.Nx-1)*(model.Nz-1),sizeof(double));
        dudz_s = DoodzCalloc ((model.Nx-0)*(model.Nz-0),sizeof(double));
        dvdx_s = DoodzCalloc ((model.Nx-0)*(model.Nz-0),sizeof(double));

#pragma omp parallel for shared ( mesh, om_s, dudz_s, dvdx_s ) \
private ( k, l, k1, c1, c3 )                                   \
firstprivate( model )
        for ( k1=0; k1<model.Nx*model.Nz; k1++ ) {
            k  = mesh.kn[k1];
            l  = mesh.ln[k1];
            c1 = k + l*model.Nx;
            c3 = k + l*(model.Nx+1);
            if ( mesh.BCg.type[c1] != 30 ) {
                om_s[c1]   = 0.5*( (mesh.v_in[c3+1] - mesh.v_in[c3])/model.dx - (mesh.u_in[c1+model.Nx] - mesh.u_in[c1])/model.dz);
                dudz_s[c1] = (mesh.u_in[c1+model.Nx] - mesh.u_in[c1])/model.dz;
                dvdx_s[c1] = (mesh.v_in[c3+1       ] - mesh.v_in[c3])/model.dx;
            }
        }

#pragma omp parallel for shared ( mesh, dudx_n, dvdz_n ) \
private ( k, l, k1, c0, c1, c2 )                         \
firstprivate( model )
        for ( k1=0; k1<(model.Nx-1)*(model.Nz-1); k1++ ) {
            k  = mesh.kp[k1];
            l  = mesh.lp[k1];
            c0 = k  + l*(model.Nx-1);
            c1 = k  + l*(model.Nx);
            c2 = k  + l*(model.Nx+1);
            if ( mesh.BCp.type[c0] != 30 && mesh.BCp.type[c0] != 31) {
                dudx_n[c0]  = (mesh.u_in[c1+1+model.Nx]     - mesh.u_in[c1+model.Nx] )/model.dx;
                dvdz_n[c0]  = (mesh.v_in[c2+1+(model.Nx+1)] - mesh.v_in[c2+1]        )/model.dz;
            }
        }
    }

#pragma omp parallel for shared ( particles, mesh, om_s ) \
private ( k, xA, zA, VxA, VzA, VxB, VzB, VxC, VzC, VxD, VzD, OmA, OmB, OmC, OmD, txx, tzz, txz, angle, dudxA, dvdzA, dudzA, dvdxA, dudxB, dvdzB, dudzB, dvdxB, dudxC, dvdzC, dudzC, dvdxC, dudxD, dvdzD, dudzD, dvdxD, VEA,VEB,VEC,VED, nx, nz, ndotx, ndotz, w12, exxA, exzA, exxB, exzB, exxC, exzC, exxD, exzD, exx, exz, norm ) \
firstprivate( model, dx, dz, new )
    for (k=0;k<Nb_part;k++) {

        if ( model.iselastic == 1 || model.aniso == 1 ) OmA = Grid2P( particles, om_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dudzA = Grid2P( particles, dudz_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dvdxA = Grid2P( particles, dvdx_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dudxA = Grid2P( particles, dudx_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dvdzA = Grid2P( particles, dvdz_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) VEA = Grid2P( particles, mesh.VE_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) exxA = Grid2P( particles, mesh.exxd, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) exzA = Grid2P( particles, mesh.exz , mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        xA           = particles->x[k];
        zA           = particles->z[k];
//        VxA          = particles->Vx[k];
//        VzA          = particles->Vz[k];
        if (new==0) VxA = Grid2P( particles, mesh.u_in, mesh.xg_coord,  mesh.zvx_coord, mesh.Nx,   mesh.Nz+1, mesh.BCu.type, dx, dz, k );
        if (new==0) VzA = Grid2P( particles, mesh.v_in, mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz, mesh.BCv.type, dx, dz, k );
        if (new==1) V2P( &VxA, &VzA, particles, mesh.u_in,  mesh.v_in, mesh.xg_coord, mesh.zg_coord, mesh.zvx_coord, mesh.xvz_coord, mesh.Nx, mesh.Nz, mesh.Nz+1, mesh.Nx+1, mesh.BCu.type, mesh.BCv.type, dx, dz, k );
        if (particles->phase[k] != -1) {
            particles->x[k] = xA + 0.5 * model.dt * VxA;
            particles->z[k] = zA + 0.5 * model.dt * VzA;
        }
        isoutPart( particles, &model, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) OmB = Grid2P( particles, om_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dudzB = Grid2P( particles, dudz_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dvdxB = Grid2P( particles, dvdx_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dudxB = Grid2P( particles, dudx_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dvdzB = Grid2P( particles, dvdz_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) VEB = Grid2P( particles, mesh.VE_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) exxB = Grid2P( particles, mesh.exxd, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) exzB = Grid2P( particles, mesh.exz , mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if (new==0) VxB = Grid2P( particles, mesh.u_in, mesh.xg_coord,  mesh.zvx_coord, mesh.Nx,   mesh.Nz+1, mesh.BCu.type, dx, dz, k );
        if (new==0) VzB = Grid2P( particles, mesh.v_in, mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz, mesh.BCv.type, dx, dz, k );
        if (new==1) V2P( &VxB, &VzB, particles, mesh.u_in,  mesh.v_in, mesh.xg_coord, mesh.zg_coord, mesh.zvx_coord, mesh.xvz_coord, mesh.Nx, mesh.Nz, mesh.Nz+1, mesh.Nx+1, mesh.BCu.type, mesh.BCv.type, dx, dz, k );
        if (particles->phase[k] != -1) {
            particles->x[k] = xA + 0.5 * model.dt * VxB;
            particles->z[k] = zA + 0.5 * model.dt * VzB;
        }
        isoutPart( particles, &model, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) OmC = Grid2P( particles, om_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dudzC = Grid2P( particles, dudz_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dvdxC = Grid2P( particles, dvdx_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dudxC = Grid2P( particles, dudx_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dvdzC = Grid2P( particles, dvdz_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) VEC = Grid2P( particles, mesh.VE_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) exxC = Grid2P( particles, mesh.exxd, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) exzC = Grid2P( particles, mesh.exz , mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if (new==0) VxC = Grid2P( particles, mesh.u_in, mesh.xg_coord,  mesh.zvx_coord, mesh.Nx,   mesh.Nz+1, mesh.BCu.type, dx, dz, k );
        if (new==0) VzC = Grid2P( particles, mesh.v_in, mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz, mesh.BCv.type, dx, dz, k );
        if (new==1) V2P( &VxC, &VzC, particles, mesh.u_in,  mesh.v_in, mesh.xg_coord, mesh.zg_coord, mesh.zvx_coord, mesh.xvz_coord, mesh.Nx, mesh.Nz, mesh.Nz+1, mesh.Nx+1, mesh.BCu.type, mesh.BCv.type, dx, dz, k );
        if (particles->phase[k] != -1) {
            particles->x[k] = xA + 1.0 * model.dt * VxC;
            particles->z[k] = zA + 1.0 * model.dt * VzC;
        }
        isoutPart( particles, &model, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) OmD = Grid2P( particles, om_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dudzD = Grid2P( particles, dudz_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dvdxD = Grid2P( particles, dvdx_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dudxD = Grid2P( particles, dudx_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) dvdzD = Grid2P( particles, dvdz_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) VED = Grid2P( particles, mesh.VE_s, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) exxD = Grid2P( particles, mesh.exxd, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1,   mesh.Nz-1, mesh.BCp.type, dx, dz, k );
        if ( model.iselastic == 1 || model.aniso == 1 ) exzD = Grid2P( particles, mesh.exz , mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type, dx, dz, k );
        if (new==0) VxD = Grid2P( particles, mesh.u_in, mesh.xg_coord,  mesh.zvx_coord, mesh.Nx,   mesh.Nz+1, mesh.BCu.type, dx, dz, k );
        if (new==0) VzD = Grid2P( particles, mesh.v_in, mesh.xvz_coord, mesh.zg_coord,  mesh.Nx+1, mesh.Nz,   mesh.BCv.type, dx, dz, k );
        if (new==1) V2P( &VxD, &VzD, particles, mesh.u_in,  mesh.v_in, mesh.xg_coord, mesh.zg_coord, mesh.zvx_coord, mesh.xvz_coord, mesh.Nx, mesh.Nz, mesh.Nz+1, mesh.Nx+1, mesh.BCu.type, mesh.BCv.type, dx, dz, k );
        VxA = (1.0/6.0) * ( VxA + 2.0 * VxB + 2.0 * VxC + VxD);
        VzA = (1.0/6.0) * ( VzA + 2.0 * VzB + 2.0 * VzC + VzD);
        if ( model.iselastic == 1 || model.aniso == 1 ) OmA   = (1.0/6.0) * ( OmA   + 2.0 * OmB   + 2.0 * OmC   + OmD);
        if ( model.iselastic == 1 || model.aniso == 1 ) dudzA = (1.0/6.0) * ( dudzA + 2.0 * dudzB + 2.0 * dudzC + dudzD);
        if ( model.iselastic == 1 || model.aniso == 1 ) dvdxA = (1.0/6.0) * ( dvdxA + 2.0 * dvdxB + 2.0 * dvdxC + dvdxD);
        if ( model.iselastic == 1 || model.aniso == 1 ) dudxA = (1.0/6.0) * ( dudxA + 2.0 * dudxB + 2.0 * dudxC + dudxD);
        if ( model.iselastic == 1 || model.aniso == 1 ) dvdzA = (1.0/6.0) * ( dvdzA + 2.0 * dvdzB + 2.0 * dvdzC + dvdzD);
        if ( model.iselastic == 1 || model.aniso == 1 ) VEA   = (1.0/6.0) * ( VEA   + 2.0 * VEB   + 2.0 * VEC   + VED);
        if ( model.iselastic == 1 || model.aniso == 1 ) exxA  = (1.0/6.0) * ( exxA   + 2.0 * exxB   + 2.0 * exxC   + exxD);
        if ( model.iselastic == 1 || model.aniso == 1 ) exzA  = (1.0/6.0) * ( exzA   + 2.0 * exzB   + 2.0 * exzC   + exzD);

        if (particles->phase[k] != -1) {
            particles->x[k]    = xA + model.dt * VxA;
            particles->z[k]    = zA + model.dt * VzA;

            // Stress rotation
            if ( model.iselastic == 1 || model.aniso == 1 ) {
                txx   = particles->sxxd[k];
                tzz   = particles->szzd[k];
                txz   = particles->sxz[k];
                if (model.StressRotation==1) {
                    angle = model.dt*OmA;
                    particles->sxxd[k] = (txx*cos(angle) - txz*sin(angle))*cos(angle) - (txz*cos(angle) - tzz*sin(angle))*sin(angle);
                    particles->szzd[k] = (txx*sin(angle) + txz*cos(angle))*sin(angle) + (txz*sin(angle) + tzz*cos(angle))*cos(angle);
                    particles->sxz[k]  = (txx*cos(angle) - txz*sin(angle))*sin(angle) + (txz*cos(angle) - tzz*sin(angle))*cos(angle);
                }
                if (model.StressRotation==2) {
                    particles->sxxd[k] -= model.dt * VEA * ( -2.0*txx*dudxA - 2.0*txz*dudzA);
                    particles->szzd[k]  = -particles->sxxd[k];
                    particles->sxz[k]  -= model.dt * VEA * (      txx*dudzA -     txx*dvdxA - txz*(dudxA + dvdzA) );
                }
            }
            
            if ( model.aniso == 1 ) {
            // Director vector rotation
            nx = particles->nx[k];// = 0.0;
            nz = particles->nz[k];// = 1.0;
            w12 = OmA; //0.5*(dudzA - dvdxA);//
            exx = exxA;
            exz = exzA;
            
            ndotx = (-(dudxA-dvdzA)*nx*nz - dvdxA*nz*nz + dudzA*nx*nx)*nz;
            ndotz = ( (dudxA-dvdzA)*nx*nz + dvdxA*nz*nz - dudzA*nx*nx)*nx;
            
//            printf("ndotx = %2.2e\n", ndotx);
//            printf("ndotz = %2.2e\n", ndotz);
            
//            ndotx = nx*(nx*(exx*nx + exz*nz) + nz*(-exx*nz + exz*nx)) + nz*w12 - (pow(nx, 2) + pow(nz, 2))*(exx*nx + exz*nz);
//            ndotz = -nx*w12 + nz*(nx*(exx*nx + exz*nz) + nz*(-exx*nz + exz*nx)) - (pow(nx, 2) + pow(nz, 2))*(-exx*nz + exz*nx);
            nx   += ndotx*model.dt;
            nz   += ndotz*model.dt;
            norm  = sqrt(nx*nx + nz*nz);
            nx   /= norm;
            nz   /= norm;
            particles->nx[k] = nx;
            particles->nz[k] = nz;
        }
        }
        isoutPart( particles, &model, k );
    }

    if ( model.iselastic == 1 || model.aniso == 1 ) {
        DoodzFree(om_s);
        DoodzFree(dudx_n);
        DoodzFree(dvdz_n);
        DoodzFree(dvdx_s);
        DoodzFree(dudz_s);
    }
    printf("** Time for Roger Gunther = %lf sec\n",  (double)((double)omp_get_wtime() - t_omp) );

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void PureShearALE( params *model, grid *Mmesh, markers *topo_chain, scale scaling ) {

    double VxW=0.0, VxE=0.0, VzS=0.0, VzN=0.0;
    int i, j, icount1, icount2, c1, c2;
    double MinFreeSurfaceUpliftRate = 0.0, MaxFreeSurfaceUpliftRate = 0.0;
    double MinFreeSurfaceAltitude   = 0.0, MaxFreeSurfaceAltitude   = 0.0;

    // Get max. surface uplift velocity and max. altitude
    if (model->ispureshear_ale == 2) {
        MinMaxArrayVal( topo_chain->Vz, topo_chain->Nb_part, &MinFreeSurfaceUpliftRate, &MaxFreeSurfaceUpliftRate );
        MinMaxArrayVal( topo_chain->z, topo_chain->Nb_part, &MinFreeSurfaceAltitude, &MaxFreeSurfaceAltitude );
    }

    // Loop on E-W sides: calculate average boundary velocity
    icount1 = 0; icount2 = 0;
    for (j=0; j<Mmesh->Nz; j++) {
        c1 = 0 +j*(Mmesh->Nx);
        c2 = (Mmesh->Nx-1) +j*(Mmesh->Nx);

        if (Mmesh->BCu.type[c1] == 0) {
            VxW += Mmesh->u_in[c1];
            icount1++;
        }
        if (Mmesh->BCu.type[c2] == 0) {
            VxE += Mmesh->u_in[c2];
            icount2++;
        }
    }
    if (icount1>0) VxW /= icount1;
    if (icount2>0) VxE /= icount2;

    // Loop on N-S sides: calculate average boundary velocity
    icount1 = 0; icount2 = 0;
    for (i=0; i<Mmesh->Nx; i++) {
        c1 = i + 0*(Mmesh->Nx+1);
        c2 = i + (Mmesh->Nz-1)*(Mmesh->Nx+1);

        if (Mmesh->BCv.type[c1] == 0) {
            if (model->free_surf == 0) {
                VzS += Mmesh->v_in[c1];
            }
            else {
                VzS += Mmesh->BCv.val[c1];
            }

            icount1++;
        }
        if (Mmesh->BCv.type[c2] == 0 || Mmesh->BCv.type[c2] == 30) {
            if (model->free_surf == 0) {
                VzN += Mmesh->v_in[c2];
            }
            else {
                // Force the N boundary velocity in case ispureshear_ale == 2
                if (model->ispureshear_ale == 2) VzN += MaxFreeSurfaceUpliftRate;
                if (Mmesh->BCv.type[c2] == 0) VzN += Mmesh->BCv.val[c2];
                //                if (model->ispureshear_ale == 2) VzN += Mmesh->BCv.val[c2];
                //                if (Mmesh->BCv.type[c2] == 0) VzN += Mmesh->BCv.val[c2];
            }
            icount2++;
        }
    }
    if (icount1>0) VzS /= icount1;
    if (icount2>0) VzN /= icount2;

    //    printf("xmin = %lf, xmax = %lf\n", model->xmin*scaling.L, model->xmax*scaling.L);
    //    printf("zmin = %lf, zmax = %lf\n", model->zmin*scaling.L, model->zmax*scaling.L);
    //    printf("VzN=%2.2e VzS=%2.2e\n",VzN*scaling.V,VzS*scaling.V);

    // Individual displacement of each boundary
    double dxmin = VxW * model->dt;
    double dxmax = VxE * model->dt;

    double dzmin = VzS * model->dt;
    double dzmax = VzN * model->dt;

    // Update min/max x-position of the mesh
    model->xmin +=  dxmin;
    model->xmax +=  dxmax;

    // Update min/max z-position of the mesh
    model->zmin +=  dzmin;
    if (model->ispureshear_ale == 2)  model->zmax  = MaxFreeSurfaceAltitude + 10e3/scaling.L;// dzmax + 0.25*dzmax;//model->zmax = MaxFreeSurfaceAltitude + 0.25*MaxFreeSurfaceAltitude;
    else                              model->zmax +=  dzmax;

    //    printf ("Max alt. %2.2e Max new %2.2e\n", (MaxFreeSurfaceAltitude + 0.25*MaxFreeSurfaceAltitude)*scaling.L, model->zmax*scaling.L);

    //    if (model->ispureshear_ale == 2)  model->zmax +=  dzmax;
    //    else                              model->zmax +=  dzmax;


    printf("Adjusting the mesh: Epsilon_xx = %2.2e, Volume = %2.2e\n", model->EpsBG * model->dt, (model->xmax - model->xmin) *(model->zmax - model->zmin) * scaling.L*scaling.L);
    printf("xmin = %lf, xmax = %lf\n", model->xmin*scaling.L, model->xmax*scaling.L);
    printf("zmin = %lf, zmax = %lf\n", model->zmin*scaling.L, model->zmax*scaling.L);

    // Remesh
    SetGridCoordinates( Mmesh, model, model->Nx, model->Nz);

    // Re-generate homogeneous pure shear fields
    printf("Re-generate homogeneous pure shear fields\n");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
void VelocitiesToParticles( grid *mesh, markers *particles, DoodzFP *Vx, DoodzFP *Vz, params model, scale scaling ) {

    DoodzFP *VxFromCenters, *VzFromCenters, *VxCenters, *VzCenters;
    int k;

    VxFromCenters = DoodzCalloc( particles->Nb_part, sizeof(DoodzFP) );
    VzFromCenters = DoodzCalloc( particles->Nb_part, sizeof(DoodzFP) );
    VxCenters     = DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1), sizeof(DoodzFP) );
    VzCenters     = DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1), sizeof(DoodzFP) );

    // Interp Vx, Vz on centroids
    VelocitiesOnCenters(  mesh->u_in, mesh->v_in, VxCenters, VzCenters, mesh->Nx, mesh->Nz, scaling );
    Interp_Grid2P( *particles, VxFromCenters, mesh, VxCenters, mesh->xc_coord,  mesh->zc_coord, mesh->Nx-1,   mesh->Nz-1, mesh->BCp.type ); //
    Interp_Grid2P( *particles, VzFromCenters, mesh, VzCenters, mesh->xc_coord,  mesh->zc_coord, mesh->Nx-1,   mesh->Nz-1, mesh->BCp.type );

    // Vx  <-- u_in
    Interp_Grid2P( *particles, Vx, mesh, mesh->u_in, mesh->xg_coord,  mesh->zvx_coord, mesh->Nx,   mesh->Nz+1, mesh->BCu.type );
    Interp_Grid2P( *particles, Vz, mesh, mesh->v_in, mesh->xvz_coord, mesh->zg_coord,  mesh->Nx+1, mesh->Nz,   mesh->BCv.type );


#pragma omp parallel for shared ( particles, VxFromCenters, VzFromCenters ) \
private ( k )
    for ( k=0; k<particles->Nb_part; k++ ) {
        Vx[k] *= 0.666666666666;
        Vz[k] *= 0.666666666666;
        Vx[k] += 0.333333333333*VxFromCenters[k];
        Vz[k] += 0.333333333333*VzFromCenters[k];
    }

    DoodzFree( VxFromCenters );
    DoodzFree( VzFromCenters );
    DoodzFree( VxCenters );
    DoodzFree( VzCenters );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DefineInitialTimestep( params *model, grid *Mmesh, markers particles, mat_prop materials, scale scaling ) {
    int k;
    double dt_maxwell, minMaxwell = 1e100, maxMaxwell = 0.0;

// Define initial timestep is elasticity is turned on
if ( model->iselastic == 1 && model->dt_constant != 1 ) {

    Interp_P2C ( particles, materials.mu,  Mmesh, Mmesh->mu_n,      Mmesh->xg_coord, Mmesh->zg_coord, 0, 1 );
    Interp_P2N ( particles, materials.mu,  Mmesh, Mmesh->mu_s,      Mmesh->xg_coord, Mmesh->zg_coord, 0, 1, model );

    for ( k=0; k<Mmesh->Nx*Mmesh->Nz; k++) {

        if ( Mmesh->BCg.type[k] != 30 ) {
            if ( Mmesh->eta_s[k]/Mmesh->mu_s[k] < minMaxwell ) {
                minMaxwell = Mmesh->eta_s[k]/Mmesh->mu_s[k];
            }
            if ( Mmesh->eta_s[k]/Mmesh->mu_s[k] > maxMaxwell ) {
                maxMaxwell = Mmesh->eta_s[k]/Mmesh->mu_s[k];
            }
        }
    }
    dt_maxwell = minMaxwell;// exp((log(minMaxwell)+log(maxMaxwell))/2); // Good for stress loading

    if (dt_maxwell < model->dt) {
        model->dt = dt_maxwell;
        model->dt0 = model->dt;
        printf("Setting initial dt to minimum Maxwell time: %2.2e\n", dt_maxwell*scaling.t);
    }
}

//// Initial Courant limit
//    if ( model->dt_constant != 1 && fabs(model->EpsBG>0.0) && model->dt > 0.1*model->dx/(fabs(model->EpsBG) * (model->xmax - model->xmin)) ) {
//    printf("Initial Courant check:\n");
//    printf("Suggested dt = %2.2e s (Vchar= %2.2e m/s)\n", model->dt*scaling.t, model->EpsBG * (model->xmax - model->xmin)*scaling.V);
//    model->dt  = 0.1*model->dx/(fabs(model->EpsBG) * (model->xmax - model->xmin));
//    model->dt0 = model->dt;
//    printf("new dt = %2.2e s\n", model->dt*scaling.t);
//}

if (model->iselastic == 1) printf("min. Maxwell = %2.2e s, max. Maxwell = %2.2e s\n", minMaxwell*scaling.t, maxMaxwell*scaling.t);
if (model->iselastic == 1)printf("Suggested dt = %2.2e s, VE dt = %2.2e s\n", model->dt*scaling.t, exp((log(minMaxwell)+log(maxMaxwell))/2.0)*scaling.t );
printf("Initial timestep = %2.2e s\n", model->dt*scaling.t);

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void EvaluateCourantCriterion( double* Vx, double* Vz, params *model, scale scaling, grid *mesh, int quiet ) {

    int k, l, c;
    double minVx=0.0, minVz=0.0, maxVx=0.0, maxVz=0.0, dmin, dtc=0.0, vmax, vmin;
    double C = model->Courant;
//    double time_reaction = 3.1558e11;
//    int reaction_in_progress;
    
    model->dt0 = model->dt;

    for (k=0; k<model->Nx; k++) {
        for (l=0; l<model->Nz+1; l++) {
            c = k + l*model->Nx;
            maxVx = MAXV(maxVx, (Vx[c]));
            minVx = MINV(minVx, (Vx[c]));
        }
    }

    for (k=0; k<model->Nx+1; k++) {
        for (l=0; l<model->Nz; l++) {
            c = k + l*(model->Nx+1);
            maxVz = MAXV(maxVz, (Vz[c]));
            minVz = MINV(minVz, (Vz[c]));
        }
    }
    if (quiet==0) printf("Min Vxm = %2.2e m/s / Max Vxm = %2.2e m/s\n", minVx * scaling.V, maxVx * scaling.V);
    if (quiet==0) printf("Min Vzm = %2.2e m/s / Max Vzm = %2.2e m/s\n", minVz * scaling.V, maxVz * scaling.V);

    dmin = MINV(model->dx, model->dz);
    vmax = MAXV(fabs(maxVx), fabs(maxVz));
    vmin = MAXV(fabs(minVx), fabs(minVz));
    vmax = MAXV(fabs(vmax),  fabs(vmin));

    if (model->dt_constant == 0) {

        double fact = 1.0;

        // Timestep increase factor
        if ( model->iselastic == 1 ) {
            fact = 1.25;
        }
        else {
            fact = 2.00;
        }

        // Courant dt
        dtc = C * dmin / fabs(vmax);

        printf("Courant number = %2.2e --- dtc = %2.2e\n", C, dtc*scaling.t);


        // Timestep cutoff : Do not allow for very large timestep increase
        if (dtc > fact*model->dt0 ) {
            dtc = fact*model->dt0;
        }

        // If timestep is adaptive
        if ( model->dt_constant != 1 ) {
            if (dtc<model->dt) printf("Timestep limited by Courant\n");
            model->dt = dtc;
        }


//        // REACTION DT:
//        // is there a reaction somewhere?
//        reaction_in_progress = 0;
//        for (k=0; k<model->Nx; k++) {
//            for (l=0; l<model->Nz; l++) {
//                c = k + l*(model->Nx);
//                if (mesh->ttrans0_s[c] > 0.0 && mesh->ttrans0_s[c]< time_reaction) reaction_in_progress += 1;
//            }
//        }
//
//        if (reaction_in_progress > 0){
//        printf("Reaction in progress in %d mesh(es)\n", reaction_in_progress);
//        if (model->dt>=0.1*time_reaction/scaling.t){
//            model->dt = 0.1*time_reaction/scaling.t;
//            printf("!!! => Timestep potentially limited by Chemical Reaction Time\n");
//            printf("!! model dt = 0.1 x Reaction Time = %2.2e\n", model->dt*scaling.t);
//            }
//        }
//        else{
//            printf("No Reaction in progress \n");
//            printf("No timestep limitation due to Chemical Reaction \n");
//        }

        // If there is no motion, then the timestep becomes huge: cut off the motion.
        if( model->dt>1.0e30 || vmax<1.0e-30) {
            dtc = 0.0;
            model->dt = model->dt_start;
        }

        if (quiet==0) printf("Current dt = %2.2e s / Courant dt = %2.2e s\n", model->dt * scaling.t, dtc * scaling.t );
    }
    else {
        model->dt = model->dt_start;
        if (quiet==0) printf("Fixed timestep dt = %2.2e s\n", model->dt * scaling.t );
    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Check_dt_for_advection( double* Vx, double* Vz, params *model, scale scaling, grid *mesh, int quiet ) {
    
    int k, l, c;
    double minVx=0.0, minVz=0.0, maxVx=0.0, maxVz=0.0, dmin, dtc=0.0, vmax, vmin;
    double C = model->Courant;
    double dt_solve =0.0;
    
    if (model->dt_constant == 0) {
    // Get current dt value;
    dt_solve = model->dt;
    
    // Compute dt_Courant value;
    for (k=0; k<model->Nx; k++) {
        for (l=0; l<model->Nz+1; l++) {
            c = k + l*model->Nx;
            maxVx = MAXV(maxVx, (Vx[c]));
            minVx = MINV(minVx, (Vx[c]));
        }
    }
    
    for (k=0; k<model->Nx+1; k++) {
        for (l=0; l<model->Nz; l++) {
            c = k + l*(model->Nx+1);
            maxVz = MAXV(maxVz, (Vz[c]));
            minVz = MINV(minVz, (Vz[c]));
        }
    }
    if (quiet==0) printf("Min Vxm = %2.2e m/s / Max Vxm = %2.2e m/s\n", minVx * scaling.V, maxVx * scaling.V);
    if (quiet==0) printf("Min Vzm = %2.2e m/s / Max Vzm = %2.2e m/s\n", minVz * scaling.V, maxVz * scaling.V);
    
    dmin = MINV(model->dx, model->dz);
    vmax = MAXV(fabs(maxVx), fabs(maxVz));
    vmin = MAXV(fabs(minVx), fabs(minVz));
    vmax = MAXV(fabs(vmax),  fabs(vmin));
        
    // Courant dt
    dtc = C * dmin / fabs(vmax);
    
    // If timestep is adaptive
        printf("dt_Courant = %2.2e\n", dtc*scaling.t);
        printf("dt_Solve   = %2.2e\n", dt_solve*scaling.t);
        
    model->dt = MINV(dtc,dt_solve);
        
        printf("dt selected for advection = %2.2e\n",  model->dt*scaling.t);
        
        // If there is no motion, then the timestep becomes huge: cut off the motion.
        if( model->dt>1.0e30 || vmax<1.0e-30) {
            dtc = 0.0;
            model->dt = model->dt_start;
        }
        
        //if (quiet==0) printf("Current dt = %2.2e s / Courant dt = %2.2e s\n", model->dt * scaling.t, dtc * scaling.t );
    }
    else {
        model->dt = model->dt_start;
        if (quiet==0) printf("Fixed timestep dt = %2.2e s\n", model->dt * scaling.t );
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

//void EvaluateCourantCriterionParticles( markers particles, params *model, scale scaling ) {
//
//    int k;
//    double minVx=0.0, minVz=0.0, maxVx=0.0, maxVz=0.0, dmin, dtc=0.0, vmax;
//    double C = 0.5;
//
//    model->dt0 = model->dt;
//
//    if (model->dt_constant == 0) {
//
//        #pragma omp parallel
//        {
//            double pminVx=0.0, pminVz=0.0, pmaxVx=0.0, pmaxVz=0.0;
//            #pragma omp for schedule( static )
//            for (k=0; k<particles.Nb_part; k++) {
//                if (particles.Vx[k]>pmaxVx) pmaxVx = particles.Vx[k];
//                if (particles.Vx[k]<pminVx) pminVx = particles.Vx[k];
//                if (particles.Vz[k]>pmaxVz) pmaxVz = particles.Vz[k];
//                if (particles.Vz[k]<pminVz) pminVz = particles.Vz[k];
//            }
//            #pragma omp flush (maxVx)
//            if (pmaxVx>maxVx) {
//                #pragma omp critical
//                {
//                    if (pmaxVx>maxVx) maxVx = pmaxVx;
//                }
//            }
//
//            #pragma omp flush (maxVz)
//            if (pmaxVz>maxVz) {
//                #pragma omp critical
//                {
//                    if (pmaxVz>maxVz) maxVz = pmaxVz;
//                }
//            }
//
//            #pragma omp flush (minVx)
//            if (pminVx<minVx) {
//                #pragma omp critical
//                {
//                    if (pminVx<minVx) minVx = pminVx;
//                }
//            }
//
//            #pragma omp flush (minVz)
//            if (pminVz<minVz) {
//                #pragma omp critical
//                {
//                    if (pminVz<minVz) minVz = pminVz;
//                }
//            }
//        }
//
//        printf("Min( Vx ) = %2.2e m/s / Min( Vz ) = %2.2e m/s\n", minVx * scaling.V, minVz * scaling.V);
//        printf("Max( Vx ) = %2.2e m/s / Max( Vz ) = %2.2e m/s\n", maxVx * scaling.V, maxVz * scaling.V);
//
//        dmin = MINV(model->dx, model->dz);
//        vmax = MAXV(maxVx, maxVz);
//
//        double fact;
//
//        // Timestep increase factor
//        if ( model->iselastic == 1 ) {
//            fact = 1.05;
//        }
//        else {
//            fact = 2.0;
//        }
//
//        // Courant dt
//        dtc = C * dmin / fabs(vmax);
//
//        // Timestep cutoff : Do not allow for very large timestep increase
//        if (dtc > fact*model->dt0 ) {
//            dtc = fact*model->dt0;
//        }
//
//        // If timestep is adaptive
//        if ( model->dt_constant != 1 ) {
//            model->dt = dtc;
//        }
//
//        // If there is no motion, then the timestep becomes huge: cut off the motion.
//        if( model->dt>1.0e30 ) {
//            dtc = 0;
//            model->dt = 0.0;
//        }
//
//        // Cutoff infinitely slow motion
//        if (vmax<1.0e-30) {
//            dtc = 0.0;
//            model->dt = 0.0;
//        }
//    }
//
//    printf("Running with Courant dt = %2.2e s / previous dt = %2.2e s\n", model->dt * scaling.t, model->dt0 * scaling.t );
//    //printf("Current dt = %2.2e s / Courant dt = %2.2e s\n", model->dt * scaling.t, dtc * scaling.t );
//
//}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FirstOrderUpwindAdvection( double* Vx, double* Vz, double* Field0, double* Field, grid *mesh, int Nx, int Nz, params model, scale scaling, int quiet ) {

    double dx = mesh->dx;
    double dz = mesh->dz;

    double qx, qz;
    int    it, k, l, cc;
    int    nit = 2;
    double dt=model.dt/nit;

    for ( it=0; it<nit; it++ ) {

    for ( l=1; l<Nz-1; l++ ) {
        for ( k=1; k<Nx-1; k++ ) {

            cc = l*Nx + k;

            if ( Vx[cc] > 0.0 ) {
                qx = Vx[cc]/dx * ( Field0[cc] - Field0[cc-1] );
            }
            else {
                qx = Vx[cc]/dx * ( Field0[cc+1] - Field0[cc] );
            }

            if ( Vz[cc] > 0.0 ) {
                qz = Vz[cc]/dz * ( Field0[cc] - Field0[cc-Nx] );
            }
            else {
                qz = Vz[cc]/dz * ( Field0[cc+Nx] - Field0[cc] );
            }
            Field[cc] = Field0[cc] - dt*qx - dt*qz;
        }
    }
    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void VelocitiesOnCenters( double *VxG, double *VzG, double *VxC, double *VzC, int Nx, int Nz, scale scaling ) {

    int    k, l, cC, cGx, cGz;

#pragma omp parallel for shared ( VxC, VxG, VzC, VzG ) \
private ( k, l, cC, cGx, cGz  )	\
firstprivate( Nx, Nz )
    for ( l=0; l<Nz-1; l++ ) {
        for ( k=0; k<Nx-1; k++ ) {
            cC  = l*(Nx-1) + k;
            cGx = l*(Nx)   + k + Nx;
            cGz = l*(Nx+1) + k + 1;
            VxC[cC] = 0.5* ( VxG[cGx] + VxG[cGx+1]    );
            VzC[cC] = 0.5* ( VzG[cGz] + VzG[cGz+Nx+1] );
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void VelocitiesOnVertices( double *VxG, double *VzG, double *VxN, double *VzN, int Nx, int Nz, scale scaling ) {

    int    k, l, cN, cGx, cGz;

#pragma omp parallel for shared ( VxN, VzN, VxG, VzG ) \
private ( k, l, cN, cGx, cGz  )	\
firstprivate( Nx, Nz )
    for ( l=0; l<Nz; l++ ) {
        for ( k=0; k<Nx; k++ ) {
            cN  = l*(Nx)   + k;
            cGx = l*(Nx)   + k;
            cGz = l*(Nx+1) + k;
            VxN[cN] = 0.5* ( VxG[cGx] + VxG[cGx+Nx] );
            VzN[cN] = 0.5* ( VzG[cGz] + VzG[cGz+1]  );
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void PureShearALE_X( params *model,  grid *Mmesh, markers *topo_chain, scale scaling ) {

    double VxW=0.0, VxE=0.0, VzS=0.0, VzN=0.0;
    int i, j, icount1, icount2, c1, c2;
    double MinFreeSurfaceUpliftRate = 0.0, MaxFreeSurfaceUpliftRate = 0.0;
    double MinFreeSurfaceAltitude   = 0.0, MaxFreeSurfaceAltitude   = 0.0;
    double dxmin, dxmax, dzmin, dzmax;

    // Get max. surface uplift velocity and max. altitude
    if (model->ispureshear_ale == 2) {
        MinMaxArrayVal( topo_chain->Vz, topo_chain->Nb_part, &MinFreeSurfaceUpliftRate, &MaxFreeSurfaceUpliftRate );
        MinMaxArrayVal( topo_chain->z, topo_chain->Nb_part, &MinFreeSurfaceAltitude, &MaxFreeSurfaceAltitude );
    }

    // Loop on E-W sides: calculate average boundary velocity
    icount1 = 0; icount2 = 0;
    for (j=0; j<Mmesh->Nz; j++) {
        c1 = 0 +j*(Mmesh->Nx);
        c2 = (Mmesh->Nx-1) +j*(Mmesh->Nx);

        if (Mmesh->BCu.type[c1] == 0) {
            VxW += Mmesh->u_in[c1];
            icount1++;
        }
        if (Mmesh->BCu.type[c2] == 0) {
            VxE += Mmesh->u_in[c2];
            icount2++;
        }
    }
    if (icount1>0) VxW /= icount1;
    if (icount2>0) VxE /= icount2;

    // Loop on N-S sides: calculate average boundary velocity
    icount1 = 0; icount2 = 0;
    for (i=0; i<Mmesh->Nx; i++) {
        c1 = i + 0*(Mmesh->Nx+1);
        c2 = i + (Mmesh->Nz-1)*(Mmesh->Nx+1);

        if (Mmesh->BCv.type[c1] == 0) {
            if (model->free_surf == 0) {
                VzS += Mmesh->v_in[c1];
            }
            else {
                VzS += Mmesh->BCv.val[c1];
            }

            icount1++;
        }
        if (Mmesh->BCv.type[c2] == 0 || Mmesh->BCv.type[c2] == 30) {
            if (model->free_surf == 0) {
                VzN += Mmesh->v_in[c2];
            }
            else {
                // Force the N boundary velocity in case ispureshear_ale == 2
                if (model->ispureshear_ale == 2) VzN += MaxFreeSurfaceUpliftRate;
                if (Mmesh->BCv.type[c2] == 0) VzN += Mmesh->BCv.val[c2];
            }
            icount2++;
        }
    }
    if (icount1>0) VzS /= icount1;
    if (icount2>0) VzN /= icount2;

    // Individual displacement of each boundary
    dxmin = VxW * model->dt;
    dxmax = VxE * model->dt;

    dzmin = 0.0*VzS * model->dt;
    dzmax = 0.0*VzN * model->dt;

    // Update min/max x-position of the mesh
    model->xmin +=  dxmin;
    model->xmax +=  dxmax;

    // Update min/max z-position of the mesh
    model->zmin +=  dzmin;
    if (model->ispureshear_ale == 2) model->zmax = MaxFreeSurfaceAltitude + 10.0e3/scaling.L;
    else         {
       model->zmax +=  dzmax;
    }

    printf("Adjusting the mesh: Epsilon_xx = %2.2e, Volume = %2.2e\n", model->EpsBG * model->dt, (model->xmax - model->xmin) *(model->zmax - model->zmin) * scaling.L*scaling.L);
    printf("xmin = %lf, xmax = %lf\n", model->xmin*scaling.L, model->xmax*scaling.L);
    printf("zmin = %lf, zmax = %lf\n", model->zmin*scaling.L, model->zmax*scaling.L);

    // Remesh
    SetGridCoordinates( Mmesh, model, model->Nx, model->Nz);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void PureShearALE_Z( params *model,  grid *Mmesh, markers *topo_chain, scale scaling ) {

    double VxW=0.0, VxE=0.0, VzS=0.0, VzN=0.0;
    int i, j, icount1, icount2, c1, c2;
    double MinFreeSurfaceUpliftRate = 0.0, MaxFreeSurfaceUpliftRate = 0.0;
    double MinFreeSurfaceAltitude   = 0.0, MaxFreeSurfaceAltitude   = 0.0;

    // Get max. surface uplift velocity and max. altitude
    if (model->ispureshear_ale == 2) {
        MinMaxArrayVal( topo_chain->Vz, topo_chain->Nb_part, &MinFreeSurfaceUpliftRate, &MaxFreeSurfaceUpliftRate );
        MinMaxArrayVal( topo_chain->z, topo_chain->Nb_part, &MinFreeSurfaceAltitude, &MaxFreeSurfaceAltitude );
    }

    // Loop on E-W sides: calculate average boundary velocity
    icount1 = 0; icount2 = 0;
    for (j=0; j<Mmesh->Nz; j++) {
        c1 = 0 +j*(Mmesh->Nx);
        c2 = (Mmesh->Nx-1) +j*(Mmesh->Nx);

        if (Mmesh->BCu.type[c1] == 0) {
            VxW += Mmesh->u_in[c1];
            icount1++;
        }
        if (Mmesh->BCu.type[c2] == 0) {
            VxE += Mmesh->u_in[c2];
            icount2++;
        }
    }
    if (icount1>0) VxW /= icount1;
    if (icount2>0) VxE /= icount2;

    // Loop on N-S sides: calculate average boundary velocity
    icount1 = 0; icount2 = 0;
    for (i=0; i<Mmesh->Nx; i++) {
        c1 = i + 0*(Mmesh->Nx+1);
        c2 = i + (Mmesh->Nz-1)*(Mmesh->Nx+1);

        if (Mmesh->BCv.type[c1] == 0) {
            if (model->free_surf == 0) {
                VzS += Mmesh->v_in[c1];
            }
            else {
                VzS += Mmesh->BCv.val[c1];
            }

            icount1++;
        }
        if (Mmesh->BCv.type[c2] == 0 || Mmesh->BCv.type[c2] == 30) {
            if (model->free_surf == 0) {
                VzN += Mmesh->v_in[c2];
            }
            else {
                // Force the N boundary velocity in case ispureshear_ale == 2
                if (model->ispureshear_ale == 2) VzN += MaxFreeSurfaceUpliftRate;
                if (Mmesh->BCv.type[c2] == 0) VzN += Mmesh->BCv.val[c2];
            }
            icount2++;
        }
    }
    if (icount1>0) VzS /= icount1;
    if (icount2>0) VzN /= icount2;


    // Individual displacement of each boundary
    double dxmin = 0.0*VxW * model->dt;
    double dxmax = 0.0*VxE * model->dt;

    double dzmin = VzS * model->dt;
    double dzmax = VzN * model->dt;

    // Update min/max x-position of the mesh
    model->xmin +=  dxmin;
    model->xmax +=  dxmax;

    // Update min/max z-position of the mesh
    model->zmin +=  dzmin;
    if (model->ispureshear_ale == 2)  model->zmax = MaxFreeSurfaceAltitude + 10.0e3/scaling.L;
    else                              model->zmax +=  dzmax;

    printf("Adjusting the mesh: Epsilon_xx = %2.2e, Volume = %2.2e\n", model->EpsBG * model->dt, (model->xmax - model->xmin) *(model->zmax - model->zmin) * scaling.L*scaling.L);
    printf("xmin = %lf, xmax = %lf\n", model->xmin*scaling.L, model->xmax*scaling.L);
    printf("zmin = %lf, zmax = %lf\n", model->zmin*scaling.L, model->zmax*scaling.L);

    // Remesh
    SetGridCoordinates( Mmesh, model, model->Nx, model->Nz);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
