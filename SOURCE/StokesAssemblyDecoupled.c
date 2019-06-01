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

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddCoeff3( int* J, double*A, int eqn, int jeq, int *nnzc, double coeff, int NODE_TYPE, double NODE_VAL, double* RHS ) {
    
    if (NODE_TYPE == 0 || NODE_TYPE == 31 || NODE_TYPE == 11) { // || NODE_TYPE == 13
        RHS[eqn]  -= coeff * NODE_VAL;
    }
    else {
        if (jeq <= eqn ) {
            J[(*nnzc)] = jeq;
            A[(*nnzc)] = coeff;
            (*nnzc)++;
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Continuity_InnerNodesDecoupled( SparseMat *Stokes, SparseMat *StokesC, SparseMat *StokesD, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempC, double **AtempC, int *nnzc2C, int **JtempD, double **AtempD, int *nnzc2D, int i, int j ) {
    
    double pc=0.0, uW=0.0, uE=0.0, vN=0.0, vS=0.0;
    
    // Compressibility
    if ( comp == 0 ) {
        pc = 0.0;
    }
    else {
        pc = mesh->bet[c2]/model.dt;
    }
    
    // div u
    if ( mesh->BCu.type[c1     ] != 13 ) uW = -one_dx;
    if ( mesh->BCu.type[c1+1   ] != 13 ) uE =  one_dx;
    if ( mesh->BCv.type[c3     ] != 13 ) vS = -one_dz;
    if ( mesh->BCv.type[c3+nxvz] != 13 ) vN =  one_dz;
    
    // Stencil assembly / residual
    if ( Assemble == 1 ) {
        StokesC->b[eqn] *= celvol;
        StokesD->b[eqn] *= celvol;
        if ( mesh->BCu.type[c1     ] != 13 )    AddCoeff2( JtempC[ith], AtempC[ith], eqn, Stokes->eqn_u[c1],      &(nnzc2C[ith]), uW*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      StokesC->bbc );
        if ( mesh->BCu.type[c1+1   ] != 13 )    AddCoeff2( JtempC[ith], AtempC[ith], eqn, Stokes->eqn_u[c1+1],    &(nnzc2C[ith]), uE*celvol, mesh->BCu.type[c1+1],    mesh->BCu.val[c1+1],    StokesC->bbc );
        if ( mesh->BCv.type[c3     ] != 13 )    AddCoeff2( JtempC[ith], AtempC[ith], eqn, Stokes->eqn_v[c3],      &(nnzc2C[ith]), vS*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      StokesC->bbc );
        if ( mesh->BCv.type[c3+nxvz] != 13 )    AddCoeff2( JtempC[ith], AtempC[ith], eqn, Stokes->eqn_v[c3+nxvz], &(nnzc2C[ith]), vN*celvol, mesh->BCv.type[c3+nxvz], mesh->BCv.val[c3+nxvz], StokesC->bbc );

        if (comp == 0)  AddCoeff2( JtempD[ith], AtempD[ith], eqn, eqn, &(nnzc2D[ith]), 1.0, mesh->BCp.type[c2],      mesh->BCp.val[c2],      StokesC->bbc );
        if (comp == 1)  AddCoeff2( JtempD[ith], AtempD[ith], eqn, eqn, &(nnzc2D[ith]), pc*celvol, mesh->BCp.type[c2],      mesh->BCp.val[c2],      StokesC->bbc );
        
//        AddCoeff2( JtempD[ith], AtempD[ith], eqn, eqn                   , &(nnzc2D[ith]), 1.0/celvol, mesh->BCp.type[c2],      mesh->BCp.val[c2],      StokesC->bbc );
    }
    else {
        StokesC->F[eqn] = pc*p[c2] + uW*u[c1] + uE*u[c1+1] + vS*v[c3] + vN*v[c3+nxvz];
        StokesC->F[eqn] -= StokesC->b[eqn];
        StokesC->F[eqn] *= celvol;
        StokesD->F[eqn] = StokesC->F[eqn] ;

    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_InnerNodesDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B  ) {
    
    double etaE = mesh->eta_n[c2+1];
    double etaW = mesh->eta_n[c2];
    double etaN = mesh->eta_s[c1];
    double etaS = mesh->eta_s[c1-nx];
//
    double uS=0.0, uN=0.0, uW=0.0, uE=0.0, uC=0.0, vSW=0.0, vSE=0.0, vNW=0.0, vNE=0.0, pE=0.0, pW=0.0;

    double inE=0.0, inW=0.0, inSu=0.0, inSv=0.0, inNu=0.0, inNv=0.0;

    // dsxx/dx
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1+1] != 30 ) { inW = 1.0; inE = 1.0; }
    if ( mesh->BCu.type[c1-1] == 30 && mesh->BCu.type[c1+1] != 30 ) inE = 1.0;
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1+1] == 30 ) inW = 1.0;

    // eta*dvx/dz
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) inSu = 1.0;
        if ( mesh->BCu.type[c1+nx] != 13 ) inNu = 1.0;
    }
    if ( mesh->BCu.type[c1-nx] == 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1+nx] != 13 ) inNu = 1.0;
    }
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] == 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) inSu = 1.0;
    }

    // eta*dvz/dx
    if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
        inSv = 1.0;
    }

    if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3] != 30  ) {
        inNv = 1.0;
    }

    double dx = mesh->dx;
    double dz = mesh->dz;

    uW = -2*etaW*inW*((1.0L/3.0L)*comp/dx - 1/dx)/dx;
    uC = (-etaN*inNu/dz - etaS*inSu/dz)/dz + (2*etaE*inE*((1.0L/3.0L)*comp/dx - 1/dx) - 2*etaW*inW*(-1.0L/3.0L*comp/dx + 1.0/dx))/dx;
    uE = 2*etaE*inE*(-1.0L/3.0L*comp/dx + 1.0/dx)/dx;
    uS = etaS*inSu/pow(dz, 2);
    uN = etaN*inNu/pow(dz, 2);
    vSW = -2.0L/3.0L*comp*etaW*inW/(dx*dz) + etaS*inSv/(dx*dz);
    vSE = (2.0L/3.0L)*comp*etaE*inE/(dx*dz) - etaS*inSv/(dx*dz);
    vNW = (2.0L/3.0L)*comp*etaW*inW/(dx*dz) - etaN*inNv/(dx*dz);
    vNE = -2.0L/3.0L*comp*etaE*inE/(dx*dz) + etaN*inNv/(dx*dz);

    // Add contribulition from non-conforming Dirichlets
    if ( mesh->BCu.type[c1-nx] == 11 ) uC  +=  -one_dz_dz * etaS;
    if ( mesh->BCu.type[c1+nx] == 11 ) uC  +=  -one_dz_dz * etaN;

    // Pressure gradient
    if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[c2] != 30 ) {
        pW  =   one_dx;
        pE  =  -pW;
    }
    
    // Pressure gradient
    //    if ( mesh->BCp.type[c2+1] != 30 ) pW  =   one_dx;
    //    if ( mesh->BCp.type[c2]   != 30 ) pE  =  -one_dx;
    
    //    // Inertia
    //    if ( model.isinertial == 1 || model.isinertial == 2 ) {
    //        uC -= sign * rhoVx/model.dt;
    //    }
    //    if ( model.isinertial == 2 ) {
    //        uN -= rhoVx/2*one_dz*mesh->VzVx[c1];
    //        uS += rhoVx/2*one_dz*mesh->VzVx[c1];
    //        uW += rhoVx/2*one_dx*mesh->u_in[c1];
    //        uE -= rhoVx/2*one_dx*mesh->u_in[c1];
    //    }
    
    //    // Stabilisation with density gradients
    //    if ( stab == 1 ) {
    ////        double correction = - om*0.5*model.dt * model.gx * (mesh->rho_app_n[c2+1] - mesh->rho_app_n[c2]) * one_dx;
    ////        uC += correction;
    ////        correction = - om*0.5*model.dt * model.gx * (mesh->rho_app_s[c1] - mesh->rho_app_s[c1-nx]) * one_dz;
    ////        vSW += 0.25*correction;
    ////        vNE += 0.25*correction;
    ////        vNW += 0.25*correction;
    ////        vSE += 0.25*correction;
    //    }
    
    uS=-uS, uN=-uN, uW=-uW, uE=-uE, uC=-uC, vSW=-vSW, vSE=-vSE, vNW=-vNW, vNE=-vNE, pE=-pE, pW=-pW;
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        // dsxx/dx - normal stencil
        if (mesh->BCu.type[c1-nx] != 30) {
            if (mesh->BCu.type[c1-nx] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], mesh->BCu.val[c1-nx], StokesA->bbc);
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], 2.0*mesh->BCu.val[c1-nx], StokesA->bbc);
        }
        if (mesh->BCu.type[c1-1]  != 30)     AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],  &(nnzc2A[ith]), uW*celvol, mesh->BCu.type[c1-1],  mesh->BCu.val[c1-1],  StokesA->bbc);
        AddCoeff2( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    StokesA->bbc);
        if (mesh->BCu.type[c1+1]  != 30)     AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+1],  &(nnzc2A[ith]), uE*celvol, mesh->BCu.type[c1+1],  mesh->BCu.val[c1+1],  StokesA->bbc);
        if (mesh->BCu.type[c1+nx] != 30) {
            if (mesh->BCu.type[c1+nx] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], mesh->BCu.val[c1+nx], StokesA->bbc);
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], 2*mesh->BCu.val[c1+nx], StokesA->bbc);
        }
        //--------------------
        
        if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz],   &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz],   mesh->BCv.val[c3-nxvz],   StokesA->bbc);
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz+1], &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1], mesh->BCv.val[c3-nxvz+1], StokesA->bbc);
        }
        
        if ( mesh->BCv.type[c3] != 30 && mesh->BCv.type[c3+1] != 30 )  {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3],        &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[c3],        mesh->BCv.val[c3],        StokesA->bbc);
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],      &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[c3+1],      mesh->BCv.val[c3+1],      StokesA->bbc);
        }
        
        //--------------------
        if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[c2] != 30 ) {
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2] - Stokes->neq_mom,   &(nnzc2B[ith]), pW*celvol, mesh->BCp.type[c2],   mesh->BCp.val[c2],   StokesB->bbc);
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2+1] - Stokes->neq_mom, &(nnzc2B[ith]), pE*celvol, mesh->BCp.type[c2+1], mesh->BCp.val[c2+1], StokesB->bbc);
        }
        //        if ( mesh->BCp.type[c2+1] != 30 ) AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2] - Stokes->neq_mom,   &(nnzc2B[ith]), pW*celvol, mesh->BCp.type[c2],   mesh->BCp.val[c2],   StokesB->bbc);
        //        if ( mesh->BCp.type[c2]   != 30 ) AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2+1] - Stokes->neq_mom, &(nnzc2B[ith]), pE*celvol, mesh->BCp.type[c2+1], mesh->BCp.val[c2+1], StokesB->bbc);
        
        //        if (eqn==35931) {
        //
        //            printf("x momentum --  etaW = %2.2e etaE = %2.2e etaS = %2.2e etaN = %2.2e\n", etaW, etaE, etaS, etaN);
        //            printf("uC = %02d %2.2e %2.12e\n",     mesh->BCu.type[c1+0],  uC, u[c1] );
        //            printf("uW = %02d %2.2e %2.12e\n",     mesh->BCu.type[c1-1],  uW, u[c1-1] );
        //            printf("uE = %02d %2.2e %2.2e\n",      mesh->BCu.type[c1+1],  uE, u[c1+1] );
        //            printf("uS = %02d %2.2e %2.2e \n",     mesh->BCu.type[c1-nx], uS, u[c1-nx] );
        //            printf("uN = %02d %2.2e %2.2e \n",     mesh->BCu.type[c1+nx], uN, u[c1+nx] );
        //            printf("vSE= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3-nxvz+1], vSE, v[c3-nxvz+1] );
        //            printf("vNE= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3+1], vNE, v[c3+1] );
        //            printf("vSW= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3-nxvz], vSW, v[c3-nxvz] );
        //            printf("vNW= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3], vNW, v[c3] );
        //            printf("pW = %02d %2.2e %2.2e \n",     mesh->BCp.type[c2], pW, p[c2] );
        //            printf("pE = %02d %2.2e %2.2e \n",     mesh->BCp.type[c2+1], pE, p[c2+1] );
        //
        //        }
        
        
        
        /*
         if (eqn==7656) {
         printf("\nFx = %2.2e %d\n",StokesA->F[eqn]*4e-5, mesh->BCu.type[c1] );
         
         
         ////        double F1 = pW*p[c2] + pE*p[c2+1] + vSW*v[c3-nxvz] + vSE*v[c3-nxvz+1] + vNW*v[c3] + vNE*v[c3+1] + uS*u[c1-nx] + uN*u[c1+nx] + uW*u[c1-1] + uE*u[c1+1] + uC*u[c1] + (StokesA->b[eqn] - StokesA->bbc[eqn]);
         ////        if (fabs(F1)>1e-6) {
         ////        printf("F0 = %2.2e %2.2e\n", StokesA->F[eqn], F1*celvol);
         printf("x momentum --  etaW = %2.2e etaE = %2.2e etaS = %2.2e etaN = %2.2e\n", etaW, etaE, etaS, etaN);
         printf("uC = %02d %2.2e %2.12e\n",      mesh->BCu.type[c1+0],  uC, u[c1] );
         printf("uW = %02d %2.2e %2.12e\n",      mesh->BCu.type[c1-1],  uW, u[c1-1] );
         printf("uE = %02d %2.2e %2.2e\n",      mesh->BCu.type[c1+1],  uE, u[c1+1] );
         printf("uS = %02d %2.2e %2.2e \n",     mesh->BCu.type[c1-nx], uS, u[c1-nx] );
         printf("uN = %02d %2.2e %2.2e \n",     mesh->BCu.type[c1+nx], uN, u[c1+nx] );
         printf("vSE= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3-nxvz+1], vSE, v[c3-nxvz+1] );
         printf("vNE= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3+1], vNE, v[c3+1] );
         printf("vSW= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3-nxvz], vSW, v[c3-nxvz] );
         printf("vNW= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3], vNW, v[c3] );
         printf("pW = %02d %2.2e %2.2e \n",     mesh->BCp.type[c2], pW, p[c2] );
         printf("pE = %02d %2.2e %2.2e %2.2e %2.2e bbc = %2.2e b = %2.2e\n",mesh->BCp.type[c2+1], pE, p[c2+1], StokesA->b[eqn] + StokesA->bbc[eqn], mesh->BCu.val[c1], StokesA->bbc[eqn] , StokesA->b[eqn]);
         printf("%2.2e EQN=%d\n", (uC*u[c1] + uW*u[c1-1])*4e-5* celvol, eqn);
         }
         */
        
    }
    else {
        // Residual function
        StokesA->F[eqn] = uC*u[c1];
        if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[c2] != 30 ) {
            StokesA->F[eqn]  += pW*p[c2] + pE*p[c2+1];
        }
        if ( mesh->BCv.type[c3-nxvz+1] != 30 && mesh->BCv.type[c3-nxvz] != 30  ) {
            StokesA->F[eqn] += vSW*v[c3-nxvz] + vSE*v[c3-nxvz+1];
        }
        if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3] != 30 ) {
            StokesA->F[eqn] += vNW*v[c3] + vNE*v[c3+1];
        }
        if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1-nx] != 11  ) StokesA->F[eqn] += uS*u[c1-nx];
        if ( mesh->BCu.type[c1+nx] != 30 && mesh->BCu.type[c1+nx] != 11  ) StokesA->F[eqn] += uN*u[c1+nx];
        if ( mesh->BCu.type[c1-1]  != 30 ) StokesA->F[eqn] += uW*u[c1-1];
        if ( mesh->BCu.type[c1+1]  != 30 ) StokesA->F[eqn] += uE*u[c1+1];
        if ( mesh->BCu.type[c1-nx] == 11 ) StokesA->F[eqn] += -2.0*etaS*one_dz_dz*mesh->BCu.val[c1-nx];
        if ( mesh->BCu.type[c1+nx] == 11 ) StokesA->F[eqn] += -2.0*etaN*one_dz_dz*mesh->BCu.val[c1+nx];
        StokesA->F[eqn] -= (StokesA->b[eqn]);// + Stokes->bbc[eqn];
        StokesA->F[eqn] *= celvol;
        
        /*
         if (eqn==7656) {
         printf("\nFx = %2.2e %d\n",StokesA->F[eqn]*4e-5, mesh->BCu.type[c1] );
         
         ////        double F1 = pW*p[c2] + pE*p[c2+1] + vSW*v[c3-nxvz] + vSE*v[c3-nxvz+1] + vNW*v[c3] + vNE*v[c3+1] + uS*u[c1-nx] + uN*u[c1+nx] + uW*u[c1-1] + uE*u[c1+1] + uC*u[c1] + (StokesA->b[eqn] - StokesA->bbc[eqn]);
         ////        if (fabs(F1)>1e-6) {
         ////        printf("F0 = %2.2e %2.2e\n", StokesA->F[eqn], F1*celvol);
         printf("x momentum --  etaW = %2.2e etaE = %2.2e etaS = %2.2e etaN = %2.2e\n", etaW, etaE, etaS, etaN);
         printf("uC = %02d %2.2e %2.12e\n",      mesh->BCu.type[c1+0],  uC, u[c1] );
         printf("uW = %02d %2.2e %2.12e\n",      mesh->BCu.type[c1-1],  uW, u[c1-1] );
         printf("uE = %02d %2.2e %2.2e\n",      mesh->BCu.type[c1+1],  uE, u[c1+1] );
         printf("uS = %02d %2.2e %2.2e \n",     mesh->BCu.type[c1-nx], uS, u[c1-nx] );
         printf("uN = %02d %2.2e %2.2e \n",     mesh->BCu.type[c1+nx], uN, u[c1+nx] );
         printf("vSE= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3-nxvz+1], vSE, v[c3-nxvz+1] );
         printf("vNE= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3+1], vNE, v[c3+1] );
         printf("vSW= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3-nxvz], vSW, v[c3-nxvz] );
         printf("vNW= %02d %2.2e %2.2e \n",     mesh->BCv.type[c3], vNW, v[c3] );
         printf("pW = %02d %2.2e %2.2e \n",     mesh->BCp.type[c2], pW, p[c2] );
         printf("pE = %02d %2.2e %2.2e %2.2e %2.2e bbc = %2.2e b = %2.2e\n",mesh->BCp.type[c2+1], pE, p[c2+1], StokesA->b[eqn] + StokesA->bbc[eqn], mesh->BCu.val[c1], StokesA->bbc[eqn] , StokesA->b[eqn]);
         printf("%2.2e EQN=%d\n", (uC*u[c1] + uW*u[c1-1])*4e-5* celvol, eqn);
         }
         */
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_WestPeriodicDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B  ) {
    
    int iVxW = c1+nx-2, iPrW = c2+ncx, iVzSW = c3-2, iVzNW = c3+nxvz-2;
    double etaE = mesh->eta_n[c2+1];
    double etaW = mesh->eta_n[iPrW];
    double etaN = mesh->eta_s[c1];
    double etaS = mesh->eta_s[c1-nx];
    //    double rhoVx = 0.5*(mesh->rho_app_s[c1] + mesh->rho_app_s[c1-nx]);
    double uS=0.0, uN=0.0, uW=0.0, uE=0.0, uC=0.0, vSW=0.0, vSE=0.0, vNW=0.0, vNE=0.0, pE=0.0, pW=0.0;
    
    // dsxx/dx
    if ( mesh->BCu.type[iVxW] != 30 && mesh->BCu.type[c1+1] != 30 ) {
        uW  =  2.0*one_dx_dx * etaW;//         - comp*2.0/3.0*one_dx_dx * etaW;
        uC  = -2.0*one_dx_dx * (etaE + etaW);//  + comp*2.0/3.0*one_dx_dx * etaE + comp*2.0/3.0*one_dx_dx * etaW;
        uE  =  2.0*one_dx_dx * etaE;//         - comp*2.0/3.0*one_dx_dx * etaE;
    }
    if ( mesh->BCu.type[iVxW] == 30 && mesh->BCu.type[c1+1] != 30 ) {
        uC  = -2.0*one_dx_dx * (etaE) + comp*2.0/3.0*one_dx_dx * etaE;
        uE  =  2.0*one_dx_dx * etaE - comp*2.0/3.0*one_dx_dx * etaE;
    }
    if ( mesh->BCu.type[iVxW] != 30 && mesh->BCu.type[c1+1] == 30 ) {
        uW  =  2.0*one_dx_dx * etaW - comp*2.0/3.0*one_dx_dx * etaW;
        uC  = -2.0*one_dx_dx * (etaW) + comp*2.0/3.0*one_dx_dx * etaW;
    }
    
    // eta*dvx/dz
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) {uS   =  one_dz_dz * etaS; uC  +=  -one_dz_dz * etaS;}
        if ( mesh->BCu.type[c1+nx] != 13 ) {uN   =  one_dz_dz * etaN; uC  +=  -one_dz_dz * etaN;}
    }
    if ( mesh->BCu.type[c1-nx] == 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1+nx] != 13 ) uC  +=  -one_dz_dz * etaN;
        if ( mesh->BCu.type[c1+nx] != 13 ) uN   =   one_dz_dz * etaN;
    }
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] == 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) uC  +=  -one_dz_dz * etaS;
        if ( mesh->BCu.type[c1-nx] != 13 ) uS   =   one_dz_dz * etaS;
    }
    
    if ( mesh->BCu.type[c1-nx] == 11 ) uC  +=  -one_dz_dz * etaS;
    if ( mesh->BCu.type[c1+nx] == 11 ) uC  +=  -one_dz_dz * etaN;
    
    // eta*dvz/dx
    if ( mesh->BCv.type[iVzSW] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
        vSE = -one_dx_dz * etaS + comp*2.0/3.0*one_dx_dz * etaE;
        vSW =  one_dx_dz * etaS - comp*2.0/3.0*one_dx_dz * etaW;
    }
    
    if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[iVzNW] != 30  ) {
        vNE =  one_dx_dz * etaN - comp*2.0/3.0*one_dx_dz * etaE;
        vNW = -one_dx_dz * etaN + comp*2.0/3.0*one_dx_dz * etaW;
    }
    
    // Pressure gradient
    if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[iPrW] != 30 ) {
        pW  =   one_dx;
        pE  =  -pW;
    }
    
    // Pressure gradient
    //    if ( mesh->BCp.type[c2+1] != 30 ) pW  =   one_dx;
    //    if ( mesh->BCp.type[c2]   != 30 ) pE  =  -one_dx;
    
    //    // Inertia
    //    if ( model.isinertial == 1 || model.isinertial == 2 ) {
    //        uC -= sign * rhoVx/model.dt;
    //    }
    //    if ( model.isinertial == 2 ) {
    //        uN -= rhoVx/2*one_dz*mesh->VzVx[c1];
    //        uS += rhoVx/2*one_dz*mesh->VzVx[c1];
    //        uW += rhoVx/2*one_dx*mesh->u_in[c1];
    //        uE -= rhoVx/2*one_dx*mesh->u_in[c1];
    //    }
    
    //    // Stabilisation with density gradients
    //    if ( stab == 1 ) {
    ////        double correction = - om*0.5*model.dt * model.gx * (mesh->rho_app_n[c2+1] - mesh->rho_app_n[c2]) * one_dx;
    ////        uC += correction;
    ////        correction = - om*0.5*model.dt * model.gx * (mesh->rho_app_s[c1] - mesh->rho_app_s[c1-nx]) * one_dz;
    ////        vSW += 0.25*correction;
    ////        vNE += 0.25*correction;
    ////        vNW += 0.25*correction;
    ////        vSE += 0.25*correction;
    //    }
    
    uS=-uS, uN=-uN, uW=-uW, uE=-uE, uC=-uC, vSW=-vSW, vSE=-vSE, vNW=-vNW, vNE=-vNE, pE=-pE, pW=-pW;
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        // dsxx/dx - normal stencil
        if (mesh->BCu.type[c1-nx] != 30) {
            if (mesh->BCu.type[c1-nx] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], mesh->BCu.val[c1-nx], StokesA->bbc);
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], 2.0*mesh->BCu.val[c1-nx], StokesA->bbc);
        }
        if (mesh->BCu.type[iVxW]  != 30)     AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[iVxW],  &(nnzc2A[ith]), uW*celvol, mesh->BCu.type[iVxW],  mesh->BCu.val[iVxW],  StokesA->bbc);
        AddCoeff2( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    StokesA->bbc);
        if (mesh->BCu.type[c1+1]  != 30)     AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+1],  &(nnzc2A[ith]), uE*celvol, mesh->BCu.type[c1+1],  mesh->BCu.val[c1+1],  StokesA->bbc);
        if (mesh->BCu.type[c1+nx] != 30) {
            if (mesh->BCu.type[c1+nx] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], mesh->BCu.val[c1+nx], StokesA->bbc);
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], 2*mesh->BCu.val[c1+nx], StokesA->bbc);
        }
        //--------------------
        
        if ( mesh->BCv.type[iVzSW] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzSW],   &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[iVzSW],   mesh->BCv.val[iVzSW],   StokesA->bbc);
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz+1], &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1], mesh->BCv.val[c3-nxvz+1], StokesA->bbc);
        }
        
        if ( mesh->BCv.type[iVzNW] != 30 && mesh->BCv.type[c3+1] != 30 )  {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[iVzNW],        &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[iVzNW],        mesh->BCv.val[iVzNW],        StokesA->bbc);
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],      &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[c3+1],      mesh->BCv.val[c3+1],      StokesA->bbc);
        }
        
        //--------------------
        if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[iPrW] != 30 ) {
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[iPrW] - Stokes->neq_mom,   &(nnzc2B[ith]), pW*celvol, mesh->BCp.type[iPrW],   mesh->BCp.val[iPrW],   StokesB->bbc);
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2+1] - Stokes->neq_mom, &(nnzc2B[ith]), pE*celvol, mesh->BCp.type[c2+1], mesh->BCp.val[c2+1], StokesB->bbc);
        }
        
    }
    else {
        // Residual function
        StokesA->F[eqn] = uC*u[c1];
        if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[iPrW] != 30 ) {
            StokesA->F[eqn]  += pW*p[iPrW] + pE*p[c2+1];
        }
        if ( mesh->BCv.type[c3-nxvz+1] != 30 && mesh->BCv.type[iVzSW] != 30  ) {
            StokesA->F[eqn] += vSW*v[iVzSW] + vSE*v[c3-nxvz+1];
        }
        if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[iVzNW] != 30 ) {
            StokesA->F[eqn] += vNW*v[iVzNW] + vNE*v[c3+1];
        }
        if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1-nx] != 11  ) StokesA->F[eqn] += uS*u[c1-nx];
        if ( mesh->BCu.type[c1+nx] != 30 && mesh->BCu.type[c1+nx] != 11  ) StokesA->F[eqn] += uN*u[c1+nx];
        if ( mesh->BCu.type[iVxW]  != 30 ) StokesA->F[eqn] += uW*u[iVxW];
        if ( mesh->BCu.type[c1+1]  != 30 ) StokesA->F[eqn] += uE*u[c1+1];
        if ( mesh->BCu.type[c1-nx] == 11 ) StokesA->F[eqn] += -2.0*etaS*one_dz_dz*mesh->BCu.val[c1-nx];
        if ( mesh->BCu.type[c1+nx] == 11 ) StokesA->F[eqn] += -2.0*etaN*one_dz_dz*mesh->BCu.val[c1+nx];
        StokesA->F[eqn] -= (StokesA->b[eqn]);// + Stokes->bbc[eqn];
        StokesA->F[eqn] *= celvol;
    }
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_WestNeumannDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B  ) {
    
    double etaE, etaW, etaN, etaS;
    etaE = mesh->eta_n[c2+1];
    etaW = etaE;
    etaN = mesh->eta_s[c1];
    etaS = mesh->eta_s[c1-nx];
    
    
    //    double rhoVx = 0.5*(mesh->rho_app_s[c1] + mesh->rho_app_s[c1-nx]);
    
    double uS=0.0, uN=0.0, uW=0.0, uE=0.0, uC=0.0, vSW=0.0, vSE=0.0, vNW=0.0, vNE=0.0, pE=0.0, pW=0.0;
    
    StokesA->b[eqn] -= mesh->BCu.val[c1]*one_dx;
    
    // dsxx/dx
    if ( mesh->BCu.type[c1+1] != 30 ) {
        uC  = -2.0*one_dx_dx * (etaE)  + comp*2.0/3.0*one_dx_dx * etaE;
        uE  =  2.0*one_dx_dx * etaE    - comp*2.0/3.0*one_dx_dx * etaE;
    }
    
    // eta*dvx/dz
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) {uS   =  one_dz_dz * etaS; uC  +=  -one_dz_dz * etaS;}
        if ( mesh->BCu.type[c1+nx] != 13 ) {uN   =  one_dz_dz * etaN; uC  +=  -one_dz_dz * etaN;}
    }
    if ( mesh->BCu.type[c1-nx] == 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1+nx] != 13 ) uC  +=  -one_dz_dz * etaN;
        if ( mesh->BCu.type[c1+nx] != 13 ) uN   =   one_dz_dz * etaN;
    }
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] == 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) uC  +=  -one_dz_dz * etaS;
        if ( mesh->BCu.type[c1-nx] != 13 ) uS   =   one_dz_dz * etaS;
    }
    
    if ( mesh->BCu.type[c1-nx] == 11 ) uC  +=  -one_dz_dz * etaS;
    if ( mesh->BCu.type[c1+nx] == 11 ) uC  +=  -one_dz_dz * etaN;
    
    // eta*dvz/dx
    if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
        vSE = -one_dx_dz * etaS + comp*2.0/3.0*one_dx_dz * etaE;
        vSW =  one_dx_dz * etaS - comp*2.0/3.0*one_dx_dz * etaW;
    }
    
    if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3] != 30  ) {
        vNE =  one_dx_dz * etaN - comp*2.0/3.0*one_dx_dz * etaE;
        vNW = -one_dx_dz * etaN + comp*2.0/3.0*one_dx_dz * etaW;
    }
    
    // Pressure gradient
    if ( mesh->BCp.type[c2+1] != 30 ) {
        pW  =   one_dx;
        pE  =  -pW;
    }
    //
    //    // Inertia
    //    if ( model.isinertial == 1 || model.isinertial == 2 ) {
    //        uC -= sign * rhoVx/model.dt;
    //    }
    //    if ( model.isinertial == 2 ) {
    //        uN -= rhoVx/2*one_dz*mesh->VzVx[c1];
    //        uS += rhoVx/2*one_dz*mesh->VzVx[c1];
    //        uW += rhoVx/2*one_dx*mesh->u_in[c1];
    //        uE -= rhoVx/2*one_dx*mesh->u_in[c1];
    //    }
    
    uS=-uS, uN=-uN, uW=-uW, uE=-uE, uC=-uC, vSW=-vSW, vSE=-vSE, vNW=-vNW, vNE=-vNE, pE=-pE, pW=-pW;
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        // dsxx/dx - normal stencil
        if (mesh->BCu.type[c1-nx] != 30) {
            if (mesh->BCu.type[c1-nx] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], mesh->BCu.val[c1-nx], StokesA->bbc);
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], 2.0*mesh->BCu.val[c1-nx], StokesA->bbc);
        }
        AddCoeff2( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    StokesA->bbc);
        if (mesh->BCu.type[c1+1]  != 30)     AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+1],  &(nnzc2A[ith]), uE*celvol, mesh->BCu.type[c1+1],  mesh->BCu.val[c1+1],  StokesA->bbc);
        if (mesh->BCu.type[c1+nx] != 30) {
            if (mesh->BCu.type[c1+nx] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], mesh->BCu.val[c1+nx], StokesA->bbc);
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], 2*mesh->BCu.val[c1+nx], StokesA->bbc);
        }
        //--------------------
        
        if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz],   &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz],   mesh->BCv.val[c3-nxvz],   StokesA->bbc);
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz+1], &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1], mesh->BCv.val[c3-nxvz+1], StokesA->bbc);
        }
        
        if ( mesh->BCv.type[c3] != 30 && mesh->BCv.type[c3+1] != 30 )  {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3],        &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[c3],        mesh->BCv.val[c3],        StokesA->bbc);
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],      &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[c3+1],      mesh->BCv.val[c3+1],      StokesA->bbc);
        }
        //--------------------
        if ( mesh->BCp.type[c2+1] != 30  ) {
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2+1] - Stokes->neq_mom, &(nnzc2B[ith]), pE*celvol, mesh->BCp.type[c2+1], mesh->BCp.val[c2+1], StokesB->bbc);
        }
    }
    else {
        
        // RHS
        StokesA->F[eqn] = -(StokesA->b[eqn] + StokesA->bbc[eqn]/celvol);
        StokesA->F[eqn] += uC*u[c1];
        if (mesh->BCu.type[c1-nx] != 30) {
            if (mesh->BCu.type[c1-nx] != 11) StokesA->F[eqn] += uS*u[c1-nx];
        }
        if (mesh->BCu.type[c1+1]  != 30)     StokesA->F[eqn] += uE*u[c1+1];
        if (mesh->BCu.type[c1+nx] != 30) {
            if (mesh->BCu.type[c1+nx] != 11) StokesA->F[eqn] += uN*u[c1+nx];
        }
        //--------------------
        if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
            if (mesh->BCv.type[c3-nxvz+1] != 11 && mesh->BCv.type[c3-nxvz+1] != 0) StokesA->F[eqn] += vSE*v[c3-nxvz+1];
            if (mesh->BCv.type[c3-nxvz  ] != 11 && mesh->BCv.type[c3-nxvz  ] != 0) StokesA->F[eqn] += vSW*v[c3-nxvz];
        }
        if ( mesh->BCv.type[c3] != 30 && mesh->BCv.type[c3+1] != 30 )  {
            if (mesh->BCv.type[c3+1] != 11 && mesh->BCv.type[c3+1] != 0) StokesA->F[eqn] += vNE*v[c3+1];
            if (mesh->BCv.type[c3  ] != 11 && mesh->BCv.type[c3  ] != 0) StokesA->F[eqn] += vNW*v[c3];
        }
        //--------------------
        if ( mesh->BCp.type[c2+1] != 30  ) StokesA->F[eqn]  +=  pE*p[c2+1];
        StokesA->F[eqn] *= celvol;
        
        //        printf("%2.4f %2.4f\n", u[c1-nx], uS);
        //        printf("%2.4f %2.4f\n", u[c1], uC);
        //        printf("%2.4f %2.4f\n", u[c1+1],uE);
        ////        printf("%2.2e\n", u[c1+nx]);
        ////        printf("%2.2e\n", v[c3-nxvz]);
        //        printf("%2.4f %2.4f\n", v[c3-nxvz+1], vSE);
        ////        printf("%2.2e\n", v[c3]);
        ////        printf("%2.2e\n", v[c3+1]);
        //        printf("%2.4f %2.4f\n", p[c2+1], pE);
        //        printf("F1 = %2.4f\n", u[c1-nx]*uS+u[c1]*uC+u[c1+1]*uE+p[c2+1]*pE+v[c3-nxvz+1]*vSE);
        
        
        //        // Residual function
        //        StokesA->F[eqn] = uC*u[c1];
        //        if ( mesh->BCp.type[c2+1] != 30  ) {
        //            StokesA->F[eqn]  +=  pE*p[c2+1];
        //        }
        //        if ( mesh->BCv.type[c3-nxvz+1] == -1 || mesh->BCv.type[c3-nxvz+1] == 2 ) StokesA->F[eqn] += vSE*v[c3-nxvz+1];
        //        if ( mesh->BCv.type[c3-nxvz]   == -1 || mesh->BCv.type[c3-nxvz]   == 2 ) StokesA->F[eqn] += vSW*v[c3-nxvz];
        //        if ( mesh->BCv.type[c3+1]      == -1 || mesh->BCv.type[c3+1]      == 2 ) StokesA->F[eqn] += vNE*v[c3+1];
        //        if ( mesh->BCv.type[c3]        == -1 || mesh->BCv.type[c3]        == 2 ) StokesA->F[eqn] += vNW*v[c3];
        //
        //        if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1-nx] != 11  ) StokesA->F[eqn] += uS*u[c1-nx];
        //        if ( mesh->BCu.type[c1+nx] != 30 && mesh->BCu.type[c1+nx] != 11  ) StokesA->F[eqn] += uN*u[c1+nx];
        //        if ( mesh->BCu.type[c1+1]  != 30 ) StokesA->F[eqn] += uE*u[c1+1];
        //        if ( mesh->BCu.type[c1-nx] == 11 ) StokesA->F[eqn] += -2*etaS*one_dz_dz*mesh->BCu.val[c1-nx];
        //        if ( mesh->BCu.type[c1+nx] == 11 ) StokesA->F[eqn] += -2*etaN*one_dz_dz*mesh->BCu.val[c1+nx];
        //
        //        StokesA->F[eqn] -= (StokesA->b[eqn]) + StokesA->bbc[eqn];
        //        StokesA->F[eqn] *= celvol;
        
        
        //        printf("\n F = %2.4e %2.2e %2.2e %2.2e\n\n", StokesA->F[eqn], StokesA->b[eqn], StokesA->bbc[eqn],  StokesA->b[eqn] + StokesA->bbc[eqn]);
        //        printf("%2.4e %2.4e  %2.4e %2.4e %2.4e\n", u[c1-nx], u[c1], u[c1+1], v[c3-nxvz+1], p[c2+1] );
        //        printf("F = %2.2e\n", (StokesA->b[eqn]+ StokesA->bbc[eqn] ) - (u[c1-nx]*uS + u[c1]*uC + u[c1+1]*uE + v[c3-nxvz+1]*vSE + p[c2+1]*pE));
    }
    
    ////    if (mesh->BCu.type[c1+nx]==30) {
    //    printf("x momentum W --  etaE = %2.2e etaS = %2.2e etaN = %2.2e\n", etaE, etaS, etaN);
    //    printf("uC = %02d %2.2e\n",      mesh->BCu.type[c1+0],  uC );
    //    //                    printf("uW = %02d %2.2e\n",      mesh->BCu.type[c1-1],  uW );
    //    printf("uE = %02d %2.2e\n",      mesh->BCu.type[c1+1],  uE );
    //    printf("uS = %02d %2.2e \n",     mesh->BCu.type[c1-nx], uS );
    //    printf("uN = %02d %2.2e \n",     mesh->BCu.type[c1+nx], uN );
    //    printf("vSE= %02d %2.2e \n",     mesh->BCv.type[c3-nxvz+1], vSE );
    //    printf("vNE= %02d %2.2e \n",     mesh->BCv.type[c3+1], vNE );
    //    printf("vSW= %02d %2.2e \n",     mesh->BCv.type[c3-nxvz], vSW );
    //    printf("vNW= %02d %2.2e \n",     mesh->BCv.type[c3], vNW );
    ////    //                    printf("pW = %02d %2.2e \n",     mesh->BCp.type[c2], pW );
    //    printf("pE = %02d %2.2e \n",mesh->BCp.type[c2+1], pE);
    ////    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xmomentum_EastNeumannDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B  ) {
    
    double etaW, etaN, etaS, etaE;
    etaW = mesh->eta_n[c2];
    etaE = etaW;
    etaN = mesh->eta_s[c1];
    etaS = mesh->eta_s[c1-nx];
    
    //    double rhoVx = 0.5*(mesh->rho_app_s[c1] + mesh->rho_app_s[c1-nx]);
    double uS=0.0, uN=0.0, uW=0.0, uE=0.0, uC=0.0, vSW=0.0, vSE=0.0, vNW=0.0, vNE=0.0, pE=0.0, pW=0.0;
    
    StokesA->b[eqn] += mesh->BCu.val[c1]*one_dx;
    
    // dsxx/dx
    if ( mesh->BCu.type[c1-1] != 30  ) {
        uW  =  2.0*one_dx_dx * etaW         - comp*2.0/3.0*one_dx_dx * etaW;
        uC  = -2.0*one_dx_dx * (etaW )      + comp*2.0/3.0*one_dx_dx * etaW; //!!!!!!!!
        //        uE  =  2*one_dx_dx * etaE         - comp*2.0/3.0*one_dx_dx * etaE;
    }
    
    // eta*dvx/dz
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) {uS   =  one_dz_dz * etaS; uC  +=  -one_dz_dz * etaS;}
        if ( mesh->BCu.type[c1+nx] != 13 ) {uN   =  one_dz_dz * etaN; uC  +=  -one_dz_dz * etaN;}
    }
    if ( mesh->BCu.type[c1-nx] == 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1+nx] != 13 ) uC  +=  -one_dz_dz * etaN;
        if ( mesh->BCu.type[c1+nx] != 13 ) uN   =   one_dz_dz * etaN;
    }
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] == 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) uC  +=  -one_dz_dz * etaS;
        if ( mesh->BCu.type[c1-nx] != 13 ) uS   =   one_dz_dz * etaS;
    }
    
    if ( mesh->BCu.type[c1-nx] == 11 ) uC  +=  -one_dz_dz * etaS;
    if ( mesh->BCu.type[c1+nx] == 11 ) uC  +=  -one_dz_dz * etaN;
    
    // eta*dvz/dx
    if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
        vSE = -one_dx_dz * etaS + comp*2.0/3.0*one_dx_dz * etaE;
        vSW =  one_dx_dz * etaS - comp*2.0/3.0*one_dx_dz * etaW;
    }
    
    if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3] != 30  ) {
        vNE =  one_dx_dz * etaN - comp*2.0/3.0*one_dx_dz * etaE;
        vNW = -one_dx_dz * etaN + comp*2.0/3.0*one_dx_dz * etaW;
    }
    // Pressure gradient
    if ( mesh->BCp.type[c2] != 30) {
        pW  =   one_dx;
    }
    
    //    // Inertia
    //    if ( model.isinertial == 1 || model.isinertial == 2 ) {
    //        uC -= sign * rhoVx/model.dt;
    //    }
    //    if ( model.isinertial == 2 ) {
    //        uN -= rhoVx/2*one_dz*mesh->VzVx[c1];
    //        uS += rhoVx/2*one_dz*mesh->VzVx[c1];
    //        //        uW += rhoVx/2*one_dx*mesh->u_in[c1];
    //        uE -= rhoVx/2*one_dx*mesh->u_in[c1];
    //    }
    
    // Stabilisation with density gradients
    if ( stab == 1 ) {
        //        double correction = - om*0.5*model.dt * model.gx * (mesh->rho_app_n[c2+1] - mesh->rho_app_n[c2]) * one_dx;
        //        uC += correction;
        //        correction = - om*0.5*model.dt * model.gx * (mesh->rho_app_s[c1] - mesh->rho_app_s[c1-nx]) * one_dz;
        //        vSW += 0.25*correction;
        //        vNE += 0.25*correction;
        //        vNW += 0.25*correction;
        //        vSE += 0.25*correction;
    }
    
    uS=-uS, uN=-uN, uW=-uW, uE=-uE, uC=-uC, vSW=-vSW, vSE=-vSE, vNW=-vNW, vNE=-vNE, pE=-pE, pW=-pW;
    
    if ( Assemble == 1 ) {
        
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        // dsxx/dx - normal stencil
        if (mesh->BCu.type[c1-nx] != 30) {
            if (mesh->BCu.type[c1-nx] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], mesh->BCu.val[c1-nx], StokesA->bbc);
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], 2.0*mesh->BCu.val[c1-nx], StokesA->bbc);
        }
        if (mesh->BCu.type[c1-1]  != 30)     AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],  &(nnzc2A[ith]), uW*celvol, mesh->BCu.type[c1-1],  mesh->BCu.val[c1-1],  StokesA->bbc);
        AddCoeff2( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    StokesA->bbc);
        if (mesh->BCu.type[c1+nx] != 30) {
            if (mesh->BCu.type[c1+nx] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], mesh->BCu.val[c1+nx], StokesA->bbc);
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], 2.0*mesh->BCu.val[c1+nx], StokesA->bbc);
        }
        //--------------------
        
        if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz],   &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz],   mesh->BCv.val[c3-nxvz],   StokesA->bbc);
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz+1], &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1], mesh->BCv.val[c3-nxvz+1], StokesA->bbc);
        }
        
        if ( mesh->BCv.type[c3] != 30 && mesh->BCv.type[c3+1] != 30 )  {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3],        &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[c3],        mesh->BCv.val[c3],        StokesA->bbc);
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],      &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[c3+1],      mesh->BCv.val[c3+1],      StokesA->bbc);
        }
        
        //--------------------
        if ( mesh->BCp.type[c2] != 30  ) {
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2] - Stokes->neq_mom,   &(nnzc2B[ith]), pW*celvol, mesh->BCp.type[c2],   mesh->BCp.val[c2],   StokesB->bbc);
        }
    }
    else {
        // RHS
        StokesA->F[eqn] = -(StokesA->b[eqn] + StokesA->bbc[eqn]/celvol);
        StokesA->F[eqn] += uC*u[c1];
        if (mesh->BCu.type[c1-nx] != 30) {
            if (mesh->BCu.type[c1-nx] != 11) StokesA->F[eqn] += uS*u[c1-nx];
        }
        if (mesh->BCu.type[c1-1]  != 30)     StokesA->F[eqn] += uW*u[c1-1];
        if (mesh->BCu.type[c1+nx] != 30) {
            if (mesh->BCu.type[c1+nx] != 11) StokesA->F[eqn] += uN*u[c1+nx];
        }
        //--------------------
        if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
            if (mesh->BCv.type[c3-nxvz+1] != 11 && mesh->BCv.type[c3-nxvz+1] != 0) StokesA->F[eqn] += vSE*v[c3-nxvz+1];
            if (mesh->BCv.type[c3-nxvz  ] != 11 && mesh->BCv.type[c3-nxvz  ] != 0) StokesA->F[eqn] += vSW*v[c3-nxvz];
        }
        if ( mesh->BCv.type[c3] != 30 && mesh->BCv.type[c3+1] != 30 )  {
            if (mesh->BCv.type[c3+1] != 11 && mesh->BCv.type[c3+1] != 0) StokesA->F[eqn] += vNE*v[c3+1];
            if (mesh->BCv.type[c3  ] != 11 && mesh->BCv.type[c3  ] != 0) StokesA->F[eqn] += vNW*v[c3];
        }
        //--------------------
        if ( mesh->BCp.type[c2] != 30  ) StokesA->F[eqn]  +=  pW*p[c2];
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_InnerNodesDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B ) {
    
    double uSW=0.0, uSE=0.0, uNW=0.0, uNE=0.0, vS=0.0, vW=0.0, vC=0.0, vE=0.0, vN=0.0, pN=0.0, pS=0.0;
    
    // Coefficients
    double etaS  = mesh->eta_n[c2];
    double etaN  = mesh->eta_n[c2+ncx];
    double etaE  = mesh->eta_s[c1];
    double etaW  = mesh->eta_s[c1-1];

    if (mesh->BCp.type[c2+ncx] == 31)     {
        StokesA->b[eqn] += mesh->BCv.val[c3]*one_dz;
        //        printf("Syy = %2.2e\n", mesh->BCv.val[c3]);
    }
    
    
//    // (eta du/dx) S
//    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
//        uSW =  one_dx_dz * etaW - comp*2.0/3.0*etaS*one_dx_dz;
//        uSE = -one_dx_dz * etaE + comp*2.0/3.0*etaS*one_dx_dz;
//    }
//
//    // (eta du/dx) N
//    if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
//        uNW = -one_dx_dz * etaW + comp*2.0/3.0*etaN*one_dx_dz;
//        uNE =  one_dx_dz * etaE - comp*2.0/3.0*etaN*one_dx_dz;
//    }
//
//    //    // (eta du/dx) S
//    //    if ( mesh->BCu.type[c1+nx] != 30 && mesh->BCu.type[c1] != 30 ) {
//    //        uNE =  one_dx_dz * etaE - comp*2.0/3.0*etaN*one_dx_dz;
//    //        uSE = -one_dx_dz * etaE + comp*2.0/3.0*etaS*one_dx_dz;
//    //    }
//    //
//    //    // (eta du/dx) N
//    //    if ( mesh->BCu.type[c1+nx-1] != 30 &&  mesh->BCu.type[c1-1] != 30 ) {
//    //        uNW = -one_dx_dz * etaW + comp*2.0/3.0*etaN*one_dx_dz;
//    //        uSW =  one_dx_dz * etaW - comp*2.0/3.0*etaS*one_dx_dz;
//    //    }
//
//    // dsyy/dz
//    if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3+nxvz] !=30 ) {
//        vS  =  2.0*one_dz_dz * etaS;//        - comp*2.0/3.0*etaS*one_dz_dz;
//        vC  = -2.0*one_dz_dz * (etaN + etaS);// + comp*2.0/3.0*etaN*one_dz_dz + comp*2.0/3.0*etaS*one_dz_dz;
//        vN  =  2.0*one_dz_dz * etaN;//        - comp*2.0/3.0*etaN*one_dz_dz;
//    }
//    if ( mesh->BCv.type[c3-nxvz] == 30 && mesh->BCv.type[c3+nxvz] != 30 ) {
//        vC  = -2.0*one_dz_dz * (etaN) + comp*2.0/3.0*etaN*one_dz_dz;
//        vN  =  2.0*one_dz_dz * etaN   - comp*2.0/3.0*etaN*one_dz_dz;
//    }
//    if ( mesh->BCv.type[c3-nxvz] != 30  && mesh->BCv.type[c3+nxvz] == 30 ) {
//        vS  =  2.0*one_dz_dz * etaS   - comp*2.0/3.0*etaS*one_dz_dz;
//        vC  = -2.0*one_dz_dz * (etaS) + comp*2.0/3.0*etaS*one_dz_dz;
//    }
//
//    // d/dx eta*dvzdx
//    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] != 30 ) {
//        if ( mesh->BCv.type[c3-1] != 13 ) {vW  =  one_dx_dx * etaW; vC += -one_dx_dx * etaW;}
//        if ( mesh->BCv.type[c3+1] != 13 ) {vE  =  one_dx_dx * etaE; vC += -one_dx_dx * etaE;}
//    }
//    // d/dx eta*dvzdx
//    if ( mesh->BCv.type[c3-1] == 30 && mesh->BCv.type[c3+1] != 30 ) {
//        if ( mesh->BCv.type[c3+1] != 13 ) vE  =  one_dx_dx * etaE;
//        if ( mesh->BCv.type[c3+1] != 13 ) vC += -one_dx_dx * etaE;
//    }
//    // d/dx eta*dvzdx
//    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] == 30 ) {
//        if ( mesh->BCv.type[c3-1] != 13 ) vW  =  one_dx_dx * etaW;
//        if ( mesh->BCv.type[c3-1] != 13 ) vC += -one_dx_dx * etaW;
//    }
    
    
    double inS=0.0,inN=0.0,inWu=0.0,inWv=0.0,inEu=0.0,inEv=0.0;
    
//    // (eta du/dx) S
//    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
//        uSW =  one_dx_dz * etaW - comp*2.0/3.0*etaS*one_dx_dz;
//        uSE = -one_dx_dz * etaE + comp*2.0/3.0*etaS*one_dx_dz;
//    }
//
//    // (eta du/dx) N
//    if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
//        uNW = -one_dx_dz * etaW + comp*2.0/3.0*etaN*one_dx_dz;
//        uNE =  one_dx_dz * etaE - comp*2.0/3.0*etaN*one_dx_dz;
//    }
//
        // (eta du/dx) S
        if ( mesh->BCu.type[c1+nx] != 30 && mesh->BCu.type[c1] != 30 ) {
            inEu = 1.0;
        }
    
        // (eta du/dx) N
        if ( mesh->BCu.type[c1+nx-1] != 30 &&  mesh->BCu.type[c1-1] != 30 ) {
            inWu = 1.0;
        }
    
    // dsyy/dz
    if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3+nxvz] !=30 ) {
        inS = 1.0;
        inN = 1.0;
    }
    if ( mesh->BCv.type[c3-nxvz] == 30 && mesh->BCv.type[c3+nxvz] != 30 ) {
        inN = 1.0;
    }
    if ( mesh->BCv.type[c3-nxvz] != 30  && mesh->BCv.type[c3+nxvz] == 30 ) {
        inS = 1.0;
    }
    
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) inWv = 1.0;
        if ( mesh->BCv.type[c3+1] != 13 ) inEv = 1.0;
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] == 30 && mesh->BCv.type[c3+1] != 30 ) {
        inEv = 1.0;
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] == 30 ) {
        inWv = 1.0;
    }

    double dx = mesh->dx;
    double dz = mesh->dz;
    
    vW = etaW*inWv/pow(dx, 2);
    vC = (2*etaN*inN*((1.0L/3.0L)*comp/dz - 1/dz) - 2*etaS*inS*(-1.0L/3.0L*comp/dz + 1.0/dz))/dz + (-etaE*inEv/dx - etaW*inWv/dx)/dx;
    vE = etaE*inEv/pow(dx, 2);
    vS = -2*etaS*inS*((1.0L/3.0L)*comp/dz - 1/dz)/dz;
    vN = 2*etaN*inN*(-1.0L/3.0L*comp/dz + 1.0/dz)/dz;
    uSW = -2.0L/3.0L*comp*etaS*inS/(dx*dz) + etaW*inWu/(dx*dz);
    uSE = (2.0L/3.0L)*comp*etaS*inS/(dx*dz) - etaE*inEu/(dx*dz);
    uNW = (2.0L/3.0L)*comp*etaN*inN/(dx*dz) - etaW*inWu/(dx*dz);
    uNE = -2.0L/3.0L*comp*etaN*inN/(dx*dz) + etaE*inEu/(dx*dz);
    
    if ( mesh->BCv.type[c3-1] == 11 ) vC  +=  -one_dx_dx * etaW;
    if ( mesh->BCv.type[c3+1] == 11 ) vC  +=  -one_dx_dx * etaE;
    
    // Stabilisation with density gradients
    if (stab==1) {
        double drhodz  = (mesh->rho_app_n[c2+ncx] - mesh->rho_app_n[c2]);
//        double vC_corr = - 1.00 * om * model.dt * model.gz * drhodz;
        double correction = - 1.0/om * 0.5 * model.dt * model.gz * drhodz * one_dz;
        vC += correction;

//        om = 0.15;
//        double vC_corr = -  1.0/om * 0.5 * model.dt * model.gz * drhodz;
        // Importante trique, voire meme gigantesque!

//        if (vC+vC_corr<0.0)  vC += vC_corr;
        //        if (vW+.25*vC_corr>0)  vW += .25*vC_corr;
        //        if (vE+.25*vC_corr>0)  vE += .25*vC_corr;
        //        if (vS+.25*vC_corr>0)  vS += .25*vC_corr;
        //        if (vN+.25*vC_corr>0)  vN += .25*vC_corr;

        // Non-symmetric contibution to the system of equation
        //        double drhodx = (mesh->rho_app_s[c1]   - mesh->rho_app_s[c1-1])*one_dx;
        //        uNW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uNE += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSE += - 0.25 * om * model.dt * model.gz * drhodx;
    }

    // Pressure gradient
    if ( mesh->BCp.type[c2] != 30 && mesh->BCp.type[c2+ncx] != 30 ) {
        pS  =   one_dz;
        pN  =  -one_dz;
    }
    
    uSW=-uSW, uSE=-uSE, uNW=-uNW, uNE=-uNE, vS=-vS, vW=-vW, vC=-vC, vE=-vE, vN=-vN, pN=-pN, pS=-pS;
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
//        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
//            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    StokesA->bbc );
//            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      StokesA->bbc );
//        }
//        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30) {
//            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], StokesA->bbc );
//            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   StokesA->bbc );
//        }
        
                if ( mesh->BCu.type[c1+nx] != 30 && mesh->BCu.type[c1] != 30 ) {
                    AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   StokesA->bbc );
                    AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      StokesA->bbc );
                }
                if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1-1] != 30) {
                    AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], StokesA->bbc );
                    AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    StokesA->bbc );
                }
    
        //--------------------
        if ( mesh->BCv.type[c3-nxvz] != 30  ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz],  &(nnzc2A[ith]), vS*celvol, mesh->BCv.type[c3-nxvz], mesh->BCv.val[c3-nxvz], StokesA->bbc );
        if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3-1] != -12  ) {
            if( mesh->BCv.type[c3-1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    mesh->BCv.val[c3-1],    StokesA->bbc );
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    2.0*mesh->BCv.val[c3-1],    StokesA->bbc );
        }
        // Periodic
        if ( mesh->BCv.type[c3-1] == -12 ) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+nxvz-3],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3+nxvz-3],    mesh->BCv.val[c3+nxvz-3],    StokesA->bbc );
        }
        AddCoeff2( JtempA[ith], AtempA[ith], eqn, eqn,                     &(nnzc2A[ith]), vC*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      StokesA->bbc );
        if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3+1] != -12  ) {
            if ( mesh->BCv.type[c3+1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    mesh->BCv.val[c3+1],    StokesA->bbc );
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    2.0*mesh->BCv.val[c3+1],    StokesA->bbc );
        }
        // Periodic
        if ( mesh->BCv.type[c3+1] == -12 ) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz+3],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3-nxvz+3],    mesh->BCv.val[c3-nxvz+3],    StokesA->bbc );
        }
        
        if ( mesh->BCv.type[c3+nxvz] != 30  ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+nxvz],  &(nnzc2A[ith]), vN*celvol, mesh->BCv.type[c3+nxvz], mesh->BCv.val[c3+nxvz], StokesA->bbc );
        //--------------------
        if ( mesh->BCp.type[c2] != 30 && mesh->BCp.type[c2+ncx] != 30) {
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2] - Stokes->neq_mom,      &(nnzc2B[ith]), pS*celvol, mesh->BCp.type[c2],     mesh->BCp.val[c2],     StokesB->bbc );
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2+ncx] - Stokes->neq_mom,  &(nnzc2B[ith]), pN*celvol, mesh->BCp.type[c2+ncx], mesh->BCp.val[c2+ncx], StokesB->bbc );
        }
        
    }
    else {
        
        // Residual
        StokesA->F[eqn] = vC*v[c3];
        if ( mesh->BCp.type[c2] != 30 && mesh->BCp.type[c2+ncx] != 30 ) {
            StokesA->F[eqn]  += pS*p[c2] + pN*p[c2+ncx];
        }
        if ( mesh->BCv.type[c3-nxvz] != 30 ) StokesA->F[eqn] += vS*v[c3-nxvz];
        if ( mesh->BCv.type[c3-1]    != 30 && mesh->BCv.type[c3-1] != 11 && mesh->BCv.type[c3-1] != -12 ) StokesA->F[eqn] += vW*v[c3-1];
        if ( mesh->BCv.type[c3+1]    != 30 && mesh->BCv.type[c3+1] != 11 && mesh->BCv.type[c3+1] != -12 ) StokesA->F[eqn] += vE*v[c3+1];
        if ( mesh->BCv.type[c3+nxvz] != 30 ) StokesA->F[eqn] += vN*v[c3+nxvz];
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            StokesA->F[eqn] += uSW*u[c1-1] + uSE*u[c1];
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
            StokesA->F[eqn] += uNW*u[c1+nx-1] + uNE*u[c1+nx];
        }
        if ( mesh->BCv.type[c3-1] == 11 )   StokesA->F[eqn] += -2*etaW*one_dx_dx*mesh->BCv.val[c3-1];
        if ( mesh->BCv.type[c3+1] == 11 )   StokesA->F[eqn] += -2*etaE*one_dx_dx*mesh->BCv.val[c3+1];
        if ( mesh->BCv.type[c3-1] == -12 ) StokesA->F[eqn] += vW*v[c3+nxvz-3];
        if ( mesh->BCv.type[c3+1] == -12 ) StokesA->F[eqn] += vE*v[c3-nxvz+3];
        StokesA->F[eqn] -= (StokesA->b[eqn]);
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_SouthNeumannDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B ) {
    
    double uSW=0.0, uSE=0.0, uNW=0.0, uNE=0.0, vS=0.0, vW=0.0, vC=0.0, vE=0.0, vN=0.0, pN=0.0, pS=0.0;
    
    // Coefficients
    double etaS, etaN, etaW, etaE;
    etaN  = mesh->eta_n[c2+ncx];
    etaS  = etaN;
    etaE  = mesh->eta_s[c1];
    etaW  = mesh->eta_s[c1-1];
    
    StokesA->b[eqn] -= mesh->BCv.val[c3]*one_dz;
    
    // (eta du/dx) S
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
        uSW =  one_dx_dz * etaW - comp*2.0/3.0*etaS*one_dx_dz;
        uSE = -one_dx_dz * etaE + comp*2.0/3.0*etaS*one_dx_dz;
    }
    
    // (eta du/dx) N
    if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        uNW = -one_dx_dz * etaW + comp*2.0/3.0*etaN*one_dx_dz;
        uNE =  one_dx_dz * etaE - comp*2.0/3.0*etaN*one_dx_dz;
    }
    
    // dsyy/dz
    if (  mesh->BCv.type[c3+nxvz] !=30 ) { //mesh->BCv.type[c3-nxvz] != 30 &&
        vC  = -2.0*one_dz_dz * ( etaN ) + comp*2.0/3.0*etaS*one_dz_dz;
        vN  =  2.0*one_dz_dz *   etaN   - comp*2.0/3.0*etaN*one_dz_dz;
    }
    
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) {vW  =  one_dx_dx * etaW; vC += -one_dx_dx * etaW;}
        if ( mesh->BCv.type[c3+1] != 13 ) {vE  =  one_dx_dx * etaE; vC += -one_dx_dx * etaE;}
        
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] == 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3+1] != 13 ) vE  =  one_dx_dx * etaE;
        if ( mesh->BCv.type[c3+1] != 13 ) vC += -one_dx_dx * etaE;
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] == 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) vW  =  one_dx_dx * etaW;
        if ( mesh->BCv.type[c3-1] != 13 ) vC += -one_dx_dx * etaW;
    }
    
    if ( mesh->BCv.type[c3-1] == 11 ) vC  +=  -one_dx_dx * etaW;
    if ( mesh->BCv.type[c3+1] == 11 ) vC  +=  -one_dx_dx * etaE;
    
    // Stabilisation with density gradients
    if (stab==1) {
        double drhodz = (mesh->rho_app_n[c2+ncx] - mesh->rho_app_n[c2])*one_dz;
        //        double drhodx = (mesh->rho_app_s[c1]   - mesh->rho_app_s[c1-1])*one_dx;
        vC  +=  1.00 * om * model.dt * model.gz * drhodz;
        //        // Non-symmetric contibution to the system of equation
        //        uNW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uNE += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSE += - 0.25 * om * model.dt * model.gz * drhodx;
    }
    
    // Pressure gradient
    if ( mesh->BCp.type[c2+ncx] != 30 ) {
        pN  =  -one_dz;
    }
    
    uSW=-uSW, uSE=-uSE, uNW=-uNW, uNE=-uNE, vS=-vS, vW=-vW, vC=-vC, vE=-vE, vN=-vN, pN=-pN, pS=-pS;
    
    //    printf("%1.1e %1.1e %1.1e %1.1e %1.1e %1.1e %1.1e %1.1e %1.1e\n", vC, vW,  vE, vN, uSW, uSE, uNW, uNE, pN);
    
    
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    StokesA->bbc );
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      StokesA->bbc );
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], StokesA->bbc );
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   StokesA->bbc );
        }
        //--------------------
        if ( mesh->BCv.type[c3-1] != 30 ) {
            if( mesh->BCv.type[c3-1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    mesh->BCv.val[c3-1],    StokesA->bbc );
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    2*mesh->BCv.val[c3-1],    StokesA->bbc );
        }
        AddCoeff2( JtempA[ith], AtempA[ith], eqn, eqn,                     &(nnzc2A[ith]), vC*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      StokesA->bbc );
        if ( mesh->BCv.type[c3+1] != 30 ) {
            if ( mesh->BCv.type[c3+1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    mesh->BCv.val[c3+1],    StokesA->bbc );
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    2*mesh->BCv.val[c3+1],    StokesA->bbc );
        }
        if ( mesh->BCv.type[c3+nxvz] != 30  ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+nxvz],  &(nnzc2A[ith]), vN*celvol, mesh->BCv.type[c3+nxvz], mesh->BCv.val[c3+nxvz], StokesA->bbc );
        //--------------------
        if ( mesh->BCp.type[c2+ncx] != 30) {
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2+ncx] - Stokes->neq_mom,  &(nnzc2B[ith]), pN*celvol, mesh->BCp.type[c2+ncx], mesh->BCp.val[c2+ncx], StokesB->bbc );
        }
        //        printf("Vz South residual = %2.2e b= %2.2e sum.c=%2.2e\n", StokesA->F[eqn], StokesA->b[eqn], uSW+uSE+uNW+uNE+vS+vW+vC+vE+vN+pN+pS);
        
    }
    else {
        // RHS
        StokesA->F[eqn] = -(StokesA->b[eqn] + StokesA->bbc[eqn]/celvol);
        StokesA->F[eqn] += vC*v[c3];
        if (mesh->BCv.type[c3-1] != 30) {
            if (mesh->BCv.type[c3-1] != 11) StokesA->F[eqn] += vW*v[c3-1];
        }
        if (mesh->BCv.type[c3+nxvz]  != 30) StokesA->F[eqn] += vN*v[c3+nxvz];
        if (mesh->BCv.type[c3+1] != 30) {
            if (mesh->BCv.type[c3+1] != 11) StokesA->F[eqn] += vE*v[c3+1];
        }
        //--------------------
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            if (mesh->BCu.type[c1-1] != 11 && mesh->BCu.type[c1-1] != 0) StokesA->F[eqn] += uSW*u[c1-1];
            if (mesh->BCu.type[c1  ] != 11 && mesh->BCu.type[c1  ] != 0) StokesA->F[eqn] += uSE*u[c1  ];
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
            if (mesh->BCu.type[c1+nx-1] != 11 && mesh->BCu.type[c1+nx-1] != 0) StokesA->F[eqn] += uNW*u[c1+nx-1];
            if (mesh->BCu.type[c1+nx  ] != 11 && mesh->BCu.type[c1+nx  ] != 0) StokesA->F[eqn] += uNE*u[c1+nx];
        }
        //--------------------
        if ( mesh->BCp.type[c2+ncx] != 30  ) StokesA->F[eqn]  +=  pN*p[c2+ncx];
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zmomentum_NorthNeumannDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B ) {
    
    double uSW=0.0, uSE=0.0, uNW=0.0, uNE=0.0, vS=0.0, vW=0.0, vC=0.0, vE=0.0, vN=0.0, pN=0.0, pS=0.0;
    
    // Coefficients
    double etaS, etaW, etaE, etaN;
    etaS  = mesh->eta_n[c2];
    etaE  = mesh->eta_s[c1];
    etaW  = mesh->eta_s[c1-1];
    etaN  = etaS;
    
    StokesA->b[eqn] += mesh->BCv.val[c3]*one_dz;
    
    // (eta du/dx) S
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
        uSW =  one_dx_dz * etaW - comp*2.0/3.0*etaS*one_dx_dz;
        uSE = -one_dx_dz * etaE + comp*2.0/3.0*etaS*one_dx_dz;
    }
    
    // (eta du/dx) N
    if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        uNW = -one_dx_dz * etaW + comp*2.0/3.0*etaN*one_dx_dz;
        uNE =  one_dx_dz * etaE - comp*2.0/3.0*etaN*one_dx_dz;
    }
    
    // dsyy/dz
    if ( mesh->BCv.type[c3-nxvz] != 30 ) {
        vS  =  2.0*one_dz_dz * etaS   - comp*2.0/3.0*etaS*one_dz_dz;
        vC  = -2.0*one_dz_dz * (etaS) + comp*2.0/3.0*etaS*one_dz_dz;
    }
    
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) {vW  =  one_dx_dx * etaW; vC += -one_dx_dx * etaW;}
        if ( mesh->BCv.type[c3+1] != 13 ) {vE  =  one_dx_dx * etaE; vC += -one_dx_dx * etaE;}
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] == 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3+1] != 13 ) vE  =  one_dx_dx * etaE;
        if ( mesh->BCv.type[c3+1] != 13 ) vC += -one_dx_dx * etaE;
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] == 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) vW  =  one_dx_dx * etaW;
        if ( mesh->BCv.type[c3-1] != 13 ) vC += -one_dx_dx * etaW;
    }
    
    if ( mesh->BCv.type[c3-1] == 11 ) vC  +=  -one_dx_dx * etaW;
    if ( mesh->BCv.type[c3+1] == 11 ) vC  +=  -one_dx_dx * etaE;
    
    // Stabilisation with density gradients
    if (stab==1) {
        double drhodz = (mesh->rho_app_n[c2+ncx] - mesh->rho_app_n[c2])*one_dz;
        //        double drhodx = (mesh->rho_app_s[c1]   - mesh->rho_app_s[c1-1])*one_dx;
        vC  += - 1.00 * om * model.dt * model.gz * drhodz;
        //        // Non-symmetric contibution to the system of equation
        //        uNW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uNE += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSE += - 0.25 * om * model.dt * model.gz * drhodx;
    }
    
    // Pressure gradient
    if (  mesh->BCp.type[c2] != 30 ) {
        pS  =   one_dz;
    }
    
    uSW=-uSW, uSE=-uSE, uNW=-uNW, uNE=-uNE, vS=-vS, vW=-vW, vC=-vC, vE=-vE, vN=-vN, pN=-pN, pS=-pS;
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    StokesA->bbc );
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      StokesA->bbc );
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], StokesA->bbc );
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   StokesA->bbc );
        }
        //--------------------
        if ( mesh->BCv.type[c3-nxvz] != 30  ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz],  &(nnzc2A[ith]), vS*celvol, mesh->BCv.type[c3-nxvz], mesh->BCv.val[c3-nxvz], StokesA->bbc );
        if ( mesh->BCv.type[c3-1] != 30 ) {
            if( mesh->BCv.type[c3-1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    mesh->BCv.val[c3-1],    StokesA->bbc );
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    2.0*mesh->BCv.val[c3-1],    StokesA->bbc );
        }
        AddCoeff2( JtempA[ith], AtempA[ith], eqn, eqn,                     &(nnzc2A[ith]), vC*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      StokesA->bbc );
        if ( mesh->BCv.type[c3+1] != 30 ) {
            if ( mesh->BCv.type[c3+1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    mesh->BCv.val[c3+1],    StokesA->bbc );
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    2.0*mesh->BCv.val[c3+1],    StokesA->bbc );
        }
        //--------------------
        if (  mesh->BCp.type[c2] != 30) {
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2] - Stokes->neq_mom,      &(nnzc2B[ith]), pS*celvol, mesh->BCp.type[c2],     mesh->BCp.val[c2],     StokesB->bbc );
        }
    }
    else {
        // RHS
        StokesA->F[eqn] = -(StokesA->b[eqn] + StokesA->bbc[eqn]/celvol);
        StokesA->F[eqn] += vC*v[c3];
        if (mesh->BCv.type[c3-1] != 30) {
            if (mesh->BCv.type[c3-1] != 11) StokesA->F[eqn] += vW*v[c3-1];
        }
        if (mesh->BCv.type[c3-nxvz]  != 30) StokesA->F[eqn] += vS*v[c3-nxvz];
        if (mesh->BCv.type[c3+1] != 30) {
            if (mesh->BCv.type[c3+1] != 11) StokesA->F[eqn] += vE*v[c3+1];
        }
        //--------------------
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            if (mesh->BCu.type[c1-1] != 11 && mesh->BCu.type[c1-1] != 0) StokesA->F[eqn] += uSW*u[c1-1];
            if (mesh->BCu.type[c1  ] != 11 && mesh->BCu.type[c1  ] != 0) StokesA->F[eqn] += uSE*u[c1  ];
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
            if (mesh->BCu.type[c1+nx-1] != 11 && mesh->BCu.type[c1+nx-1] != 0) StokesA->F[eqn] += uNW*u[c1+nx-1];
            if (mesh->BCu.type[c1+nx  ] != 11 && mesh->BCu.type[c1+nx  ] != 0) StokesA->F[eqn] += uNE*u[c1+nx];
        }
        //--------------------
        if ( mesh->BCp.type[c2-ncx] != 30  ) StokesA->F[eqn]  +=  pS*p[c2-ncx];
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MergeParallelMatrixDecoupled( SparseMat *Stokes, double **Atemp, int **Jtemp, int **Itemp, grid* mesh, int* estart, int* eend, int *nnzc, int *nnzc2, int *last_eqn, int n_th, char* BC_type, int *eqn_number ) {
    
    int eqn, ith, begin[n_th], end[n_th], c, last_eqn_serial=0, k, l;
    
    for (ith=0; ith<n_th; ith++) {
        
        // Get last equation number
        if (last_eqn[ith]>last_eqn_serial) {
            last_eqn_serial = last_eqn[ith];
        }
        
        if (ith==0) {
            begin[ith] = *(nnzc);
            end[ith]   = *(nnzc) + nnzc2[ith]-1;
        }
        else {
            begin[ith] = end[ith-1]+1;
            end[ith]   = end[ith-1]+1+nnzc2[ith]-1;
        }
    }
    
    // Recombine A and J
#pragma omp parallel private(ith, eqn, c, k, l ) shared( Stokes, Atemp, Jtemp, Itemp, begin, BC_type, eqn_number )
    {
        ith = omp_get_thread_num();
        
        for( c=estart[ith]; c<eend[ith]+1; c++) {
            
            // Get equation numbering
            //            if ( BC_type[c] != 0 && BC_type[c] != 30 && BC_type[c] != 31 && BC_type[c] != 13 && BC_type[c] != 11 ) {
            if ( BC_type[c] != 0 && BC_type[c] != 30 && BC_type[c] != 31 && BC_type[c] != 13 && BC_type[c] != 11 && BC_type[c] != -12  ) {
                eqn = eqn_number[c] -  Stokes->neq_mom;
                Stokes->Ic[eqn] = Itemp[ith][eqn] + begin[ith];
            }
        }
        
        for (k=0; k<nnzc2[ith]; k++) {
            l = begin[ith] + k;
            Stokes->A[l] = Atemp[ith][k];
            Stokes->J[l] = Jtemp[ith][k];
        }
    }
    
    // Update total non-zero number
    for (ith=0; ith<n_th; ith++) {
        *(nnzc) += nnzc2[ith];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FreeDecompositionDecoupled( int *estart, int *eend, int *DD, int *last_eqn  ) {
    DoodzFree( DD );
    DoodzFree( estart );
    DoodzFree( eend );
    DoodzFree( last_eqn );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AllocateDecompositionDecoupled( int n_th, int **estart, int **eend, int **DD, int **last_eqn ) {
    
    /* Array containing the size of the block for each thread */
    *DD       = DoodzMalloc( sizeof(int) * n_th );
    *estart   = DoodzMalloc( sizeof(int) * n_th );
    *eend     = DoodzMalloc( sizeof(int) * n_th );
    *last_eqn = DoodzCalloc( sizeof(int), n_th );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DomainDecompositionDecoupled( int N, int n_th, int *estart, int *eend, int *DD ) {
    
    int ith, eps, block_size, x0 = 0, x1 = 0;
    eps = N%n_th;
    block_size  = (N - eps) / n_th;
    
    // Divide domain into blocks
    for (ith=0; ith<n_th; ith++) {
        if (ith<n_th-1)
            DD[ith] = block_size;
        else {
            DD[ith] = block_size + eps;
        }
        
        x1 = x0 + DD[ith];
        estart[ith] = x0;
        eend[ith]   = x1-1;
        x0          = x1;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AllocateTempMatArraysDecoupled(double ***Atemp, int ***Itemp, int ***Jtemp, int n_th, int nnz, int neq, int* DD, int **nnzc2 ) {
    
    int ith;
    
    // Temporary parallel Matrix arrays
    *Atemp = DoodzMalloc( sizeof(double*) * n_th );
    *Itemp = DoodzMalloc( sizeof(int*) * n_th );
    *Jtemp = DoodzMalloc( sizeof(int*) * n_th );
    *nnzc2 = DoodzCalloc( sizeof(int), n_th );
    
    for (ith=0; ith<n_th; ith++) {
        (*Atemp)[ith] = DoodzMalloc( sizeof(double) * (int)(nnz/n_th) );
        (*Itemp)[ith] = DoodzCalloc( sizeof(int), (neq+1) );
        (*Jtemp)[ith] = DoodzMalloc( sizeof(int) * (int)(nnz/n_th) );
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FreeTempMatArraysDecoupled(double **Atemp, int **Itemp, int **Jtemp, int n_th, int *nnzc2 ) {
    
    int ith;
    
    // Temporary parallel Matrix arrays
    for (ith=0; ith<n_th; ith++) {
        DoodzFree( Atemp[ith] );
        DoodzFree( Itemp[ith] );
        DoodzFree( Jtemp[ith] );
    }
    DoodzFree( Atemp );
    DoodzFree( Itemp );
    DoodzFree( Jtemp );
    DoodzFree( nnzc2 );
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void BuildStokesOperatorDecoupled( grid *mesh, params model, int lev, double *p, double *u, double *v, SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, int Assemble ) {
    
    // FLAG
    // Assemble=1: Build the linear system of discrete equations
    // Assemble=0: Evaluate Stokes residuals
    
    int    cc, k, l, c1, c2, c3, nx=mesh->Nx, nz=mesh->Nz, nxvz=nx+1, nzvx=nz+1, ncx=nx-1, ncz=nz-1;
    
    // Pre-calculate FD coefs
    double celvol    = mesh->dx*mesh->dz;
    double one_dx    = 1.0/mesh->dx;
    double one_dz    = 1.0/mesh->dz;
    double one_dx_dx = 1.0/mesh->dx/mesh->dx;
    double one_dz_dz = 1.0/mesh->dz/mesh->dz;
    double one_dx_dz = 1.0/mesh->dx/mesh->dz;
    
    // Switches
    int eqn,  sign = 1, comp = model.compressible, stab = 0;
    if (model.free_surf_stab>0) stab = 1;
    double theta = model.free_surf_stab;
    
    //    if ( stab==1 ) printf("It tastes like Bo'\n");
    
    // Decompose domain
    int n_th, N, ith;
    int *DD, *estart, *eend, *last_eqn;
    
    // Temporary parallel matrix
    double **AtempA;
    int **JtempA, **ItempA;
    double **AtempB;
    int **JtempB, **ItempB;
    double **AtempC;
    int **JtempC, **ItempC;
    double **AtempD;
    int **JtempD, **ItempD;
    
    int nnzA, nnzB, nnzC, nnzD, nnzcA=0, nnzcB=0, nnzcC=0, nnzcD=0;
    int *nnzc2A, *nnzc2B, *nnzc2C, *nnzc2D;
    
#pragma omp parallel shared(n_th)
    {
        n_th = omp_get_num_threads();
    }
#pragma omp barrier
    
    // Matrix initialisation
    if ( Assemble == 1 ) {
        
        nnzA  = 9*((mesh->Nx-1) * mesh->Nz + (mesh->Nz-1) * mesh->Nx);
        nnzB  = 5*((mesh->Nx-1) * (mesh->Nz-1));
        nnzC  = nnzB;
        nnzD  = 1*((mesh->Nx-1) * (mesh->Nz-1));
        StokesA->neq = Stokes->neq_mom;
        StokesB->neq = Stokes->neq_mom;
        StokesC->neq = Stokes->neq_cont; StokesC->neq_mom = Stokes->neq_mom;
        StokesD->neq = Stokes->neq_cont; StokesD->neq_mom = Stokes->neq_mom;
        printf("Assembling  decoupled Stokes matrix...\n");
        AllocMat( StokesA, nnzA );
        AllocMat( StokesB, nnzB );
        AllocMat( StokesC, nnzC );
        AllocMat( StokesD, nnzD );
        
    }
    
    // Build velocity block RHS
    int inc=0;
    for( l=0; l<nzvx; l++) {
        for( k=0; k<nx; k++) {
            cc = k + l*nx;
            if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 && mesh->BCu.type[cc] != -12 ) {
                StokesA->b[inc] = mesh->roger_x[cc];
                inc++;
            }
        }
    }
    for( l=0; l<nz; l++) {
        for( k=0; k<nxvz; k++) {
            cc = k + l*nxvz;
            if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13  && mesh->BCv.type[cc] != -12  ) {
                StokesA->b[inc] = mesh->roger_z[cc];
                inc++;
            }
        }
    }
    
    // Build pressure block RHS
    inc=0;
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {
            cc = k + l*ncx;
            if ( mesh->BCp.type[cc] != 31 && mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30 ) {
                StokesC->b[inc] = mesh->rhs_p[cc];
                inc++;
            }
        }
    }
    
    //------------------------------------------------------------------//
    //------------------------- U-momentum -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = nx*nzvx;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempA, &ItempA, &JtempA, n_th, nnzA, Stokes->neq_mom, DD, &nnzc2A  );
        AllocateTempMatArraysDecoupled( &AtempB, &ItempB, &JtempB, n_th, nnzB, Stokes->neq_mom, DD, &nnzc2B  );
    }
    
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesA, StokesB, u, v, p, nx, ncx, nzvx, nnzc2A, AtempA, JtempA, ItempA, nnzc2B, AtempB, JtempB, ItempB, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, celvol, comp )
    {
        ith = omp_get_thread_num();
        
        for( c1=estart[ith]; c1<eend[ith]+1; c1++) {
            
            k      = mesh->kvx[c1];
            l      = mesh->lvx[c1];
            c2     = k-1 + (l-1)*ncx;
            c3     = k   + l*nxvz;
            
            // Get equation numbering
            if ( mesh->BCu.type[c1]==-1 || mesh->BCu.type[c1]==2 || mesh->BCu.type[c1]==-2 ) {
                eqn = Stokes->eqn_u[c1];
                last_eqn[ith]   = eqn;
                
                if ( Assemble == 1 ) {
                    ItempA[ith][eqn] = nnzc2A[ith];
                    ItempB[ith][eqn] = nnzc2B[ith];
                }
                
                //--------------------- INNER NODES ---------------------//
                if ( l>0 && l<nzvx-1 && k>0 && k<nx-1 ) {
                    Xmomentum_InnerNodesDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                if ( k==0    && mesh->BCu.type[c1]==2 ) {
                    //                    printf("%2.2e\n", mesh->eta_n[c2+1]);
                    Xmomentum_WestNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                // Periodic
                if ( k==0    && mesh->BCu.type[c1]==-2 ) {
                    Xmomentum_WestPeriodicDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                if ( k==nx-1 && mesh->BCu.type[c1]==2 ) {
                    Xmomentum_EastNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
            }
        }
    }
    
    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrix( StokesA, AtempA, JtempA, ItempA, mesh, estart, eend, &nnzcA, nnzc2A, last_eqn, n_th, mesh->BCu.type, Stokes->eqn_u );
        MergeParallelMatrix( StokesB, AtempB, JtempB, ItempB, mesh, estart, eend, &nnzcB, nnzc2B, last_eqn, n_th, mesh->BCu.type, Stokes->eqn_u );
    }
    
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempA, ItempA, JtempA, n_th, nnzc2A );
        FreeTempMatArraysDecoupled( AtempB, ItempB, JtempB, n_th, nnzc2B );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //------------------------- V-momentum -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = nxvz*nz;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempA, &ItempA, &JtempA, n_th, nnzA, Stokes->neq_mom, DD, &nnzc2A  );
        AllocateTempMatArraysDecoupled( &AtempB, &ItempB, &JtempB, n_th, nnzB, Stokes->neq_mom, DD, &nnzc2B  );
    }
    
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesA, StokesB, u, v, p, nx, ncx, nzvx, nnzc2A, AtempA, JtempA, ItempA, nnzc2B, AtempB, JtempB, ItempB, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, comp, celvol )
    {
        ith = omp_get_thread_num();
        
        for( c3=estart[ith]; c3<eend[ith]+1; c3++) {
            
            k  = mesh->kvz[c3];
            l  = mesh->lvz[c3];
            c1 = k   + l*nx;
            c2 = k-1 + (l-1)*ncx;
            
            // Get equation numbering
            if ( mesh->BCv.type[c3] == -1 || mesh->BCv.type[c3]==2 ) {
                eqn = Stokes->eqn_v[c3];
                last_eqn[ith]   = eqn;
                
                if ( Assemble == 1 ) {
                    ItempA[ith][eqn] = nnzc2A[ith];
                    ItempB[ith][eqn] = nnzc2B[ith];
                }
                
                //--------------------- INNER NODES ---------------------//
                if ( k>0 && k<nxvz-1 && l>0 && l<nz-1 ) {
                    
                    Zmomentum_InnerNodesDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                //--------------------- INNER NODES ---------------------//
                
                if ( l==0 && mesh->BCv.type[c3] == 2 ) {
                    Zmomentum_SouthNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                //--------------------- INNER NODES ---------------------//
                
                if ( l==nz-1 && mesh->BCv.type[c3] == 2 ) {
                    Zmomentum_NorthNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
            }
        }
    }
    
    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrix( StokesA, AtempA, JtempA, ItempA, mesh, estart, eend, &nnzcA, nnzc2A, last_eqn, n_th, mesh->BCv.type, Stokes->eqn_v );
        MergeParallelMatrix( StokesB, AtempB, JtempB, ItempB, mesh, estart, eend, &nnzcB, nnzc2B, last_eqn, n_th, mesh->BCv.type, Stokes->eqn_v );
    }
    
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempA, ItempA, JtempA, n_th, nnzc2A );
        FreeTempMatArraysDecoupled( AtempB, ItempB, JtempB, n_th, nnzc2B );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //------------------------- Continuity -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = ncx*ncz;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempC, &ItempC, &JtempC, n_th, nnzC, Stokes->neq_cont, DD, &nnzc2C  );
        AllocateTempMatArraysDecoupled( &AtempD, &ItempD, &JtempD, n_th, nnzD, Stokes->neq_cont, DD, &nnzc2D  );
    }
    
    //    printf("NC = %d %d %d %d %d\n", Stokes->neq_cont, estart[0], eend[0], DD[0], last_eqn[0]);
    
    //#pragma omp parallel shared( eend, estart, mesh, Stokes, u, v, p, nx, ncx, nzvx, nnzc2C, AtempC, JtempC, ItempC, nnzc2D, AtempD, JtempD, ItempD, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, comp, celvol )
    //    {
    //
    //        ith = omp_get_thread_num();
    //
    //        for( c2=estart[ith]; c2<eend[ith]+1; c2++) {
    //
    //            k   = mesh->kp[c2];
    //            l   = mesh->lp[c2];
    //            c1  = k   + (l+1)*nx;
    //            c3  = k   + l*nxvz + 1;
    //
    //            //--------------------- INNER NODES ---------------------//
    //            if ( mesh->BCp.type[c2] == -1) {
    //
    //
    //                eqn = Stokes->eqn_p[c2]  - Stokes->neq_mom;
    //                last_eqn[ith]   = eqn ;
    //
    //                if ( Assemble == 1 ) {
    //                    ItempC[ith][eqn] = nnzc2C[ith];
    //                    ItempD[ith][eqn] = nnzc2D[ith]; //printf("%d ",  nnzc2D[ith]);
    //                }
    //                //
    //                Continuity_InnerNodesDecoupled( Stokes, StokesC, StokesD, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempC, AtempC, nnzc2C, JtempD, AtempD, nnzc2D, k, l );
    //            }
    //        }
    //    }
    
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesC, StokesD, u, v, p, nnzc2C, AtempC, JtempC, ItempC, nnzc2D, AtempD, JtempD, ItempD, last_eqn )  private( ith, l, k, c1, c2, c3, eqn, comp ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, celvol, nx, ncx, nxvz, nzvx )
    {
        
        ith = omp_get_thread_num();
        
        for( c2=estart[ith]; c2<eend[ith]+1; c2++) {
            
            k   = mesh->kp[c2];
            l   = mesh->lp[c2];
            c1  = k   + (l+1)*nx;
            c3  = k   + l*nxvz + 1;
            
            //--------------------- INNER NODES ---------------------//
            if ( mesh->BCp.type[c2] == -1) {
                
                comp = mesh->comp_cells[c2];

                eqn = Stokes->eqn_p[c2]  - Stokes->neq_mom;
                last_eqn[ith]   = eqn ;
                
                if ( Assemble == 1 ) {
                    ItempC[ith][eqn] = nnzc2C[ith];
                    ItempD[ith][eqn] = nnzc2D[ith]; //printf("%d ",  nnzc2D[ith]);
                }
                //
                Continuity_InnerNodesDecoupled( Stokes, StokesC, StokesD, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempC, AtempC, nnzc2C, JtempD, AtempD, nnzc2D, k, l );
            }
        }
    }
    
    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrixDecoupled( StokesC, AtempC, JtempC, ItempC, mesh, estart, eend, &nnzcC, nnzc2C, last_eqn, n_th, mesh->BCp.type, Stokes->eqn_p );
        MergeParallelMatrixDecoupled( StokesD, AtempD, JtempD, ItempD, mesh, estart, eend, &nnzcD, nnzc2D, last_eqn, n_th, mesh->BCp.type, Stokes->eqn_p );
    }
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempC, ItempC, JtempC, n_th, nnzc2C );
        FreeTempMatArraysDecoupled( AtempD, ItempD, JtempD, n_th, nnzc2D );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //----------------------------- End --------------------------------//
    //------------------------------------------------------------------//
    
    
    
    if ( Assemble == 1 ) {
        
        // Add contribution from the BC's
        for (k=0; k<StokesA->neq; k++) {
            StokesA->b[k] += StokesA->bbc[k];
            StokesB->b[k] = StokesA->b[k];
        }
        
        //            MinMaxArray(StokesA->b, 1, StokesA->neq, "rhs_mom_here" );
        //            MinMaxArray(StokesC->b, 1, StokesC->neq, "rhs_cont_here" );
        
        
        // Add contribution from the BC's
        for (k=0; k<StokesC->neq; k++) {
            StokesC->b[k] += StokesC->bbc[k];
            StokesD->b[k] = StokesC->b[k];
        }
        
        
        
        // Final index
        StokesA->Ic[StokesA->neq] = nnzcA;
        StokesB->Ic[StokesB->neq] = nnzcB;
        StokesC->Ic[StokesC->neq] = nnzcC;
        StokesD->Ic[StokesD->neq] = nnzcD;
        
        //        for (ieq=0; ieq<Stokes->neq_cont+1; ieq++) printf( "%d ", StokesC->I[ieq] );
        
        StokesA->nnz = nnzcA;
        StokesB->nnz = nnzcB;
        StokesC->nnz = nnzcC;
        StokesD->nnz = nnzcD;
        
        // Resize arrays to proper number of non-zeros
        double *bufd;
        int *bufi;
        bufi      = DoodzRealloc(StokesA->J, nnzcA*sizeof(int));
        bufd      = DoodzRealloc(StokesA->A, nnzcA*sizeof(double));
        StokesA->J = bufi;
        StokesA->A = bufd;
        
        bufi      = DoodzRealloc(StokesB->J, nnzcB*sizeof(int));
        bufd      = DoodzRealloc(StokesB->A, nnzcB*sizeof(double));
        StokesB->J = bufi;
        StokesB->A = bufd;
        
        bufi      = DoodzRealloc(StokesC->J, nnzcC*sizeof(int));
        bufd      = DoodzRealloc(StokesC->A, nnzcC*sizeof(double));
        StokesC->J = bufi;
        StokesC->A = bufd;
        
        bufi      = DoodzRealloc(StokesD->J, nnzcD*sizeof(int));
        bufd      = DoodzRealloc(StokesD->A, nnzcD*sizeof(double));
        StokesD->J = bufi;
        StokesD->A = bufd;
        
        
        printf("System size: ndof = %d, nzA = %d nzB = %d nzC = %d nzD = %d\n", Stokes->neq, nnzcA, nnzcB, nnzcC, nnzcD);
        
        printf("Number of momentum equations: %d\n", StokesA->neq);
        
//        // Extract Diagonal of A - Viscous block
//        int i, j, locNNZ;
//        int I1, J1;
//
//        if (model.diag_scaling == 1) {
//
//#pragma omp parallel for shared(StokesA) private(I1,J1, i, j, locNNZ )
//        for (i=0;i<StokesA->neq; i++) {
//            I1     = StokesA->Ic[i];
//            locNNZ = StokesA->Ic[i+1] - StokesA->Ic[i];
//            for (J1=0;J1<locNNZ; J1++) {
//                 j = StokesA->J[I1 + J1];
//                if (i==j) StokesA->d[i] = StokesA->A[I1 + J1];
//             }
//        }
//        // Extract Diagonal of D - Pressure block
//#pragma omp parallel for shared(StokesD) private(i)
//        for (i=0;i<StokesD->neq; i++) {
//            StokesC->d[i] = 1.0;
//        }
//
//        // Scale A
//#pragma omp parallel for shared(StokesA) private(I1,J1, i, j, locNNZ )
//        for (i=0;i<StokesA->neq; i++) {
//
//            StokesA->b[i] *= StokesA->d[i];  // scale RHS
//
//            I1     = StokesA->Ic[i];
//            locNNZ = StokesA->Ic[i+1] - StokesA->Ic[i];
//            for (J1=0;J1<locNNZ; J1++) {
//                j = StokesA->J[I1 + J1];
//                StokesA->A[I1 + J1] *= StokesA->d[i]*StokesA->d[j];
//            }
//        }
//
//        // Scale B
//#pragma omp parallel for shared(StokesB,StokesA) private(I1,J1, i, j, locNNZ )
//        for (i=0;i<StokesB->neq; i++) {
//            I1     = StokesB->Ic[i];
//            locNNZ = StokesB->Ic[i+1] - StokesB->Ic[i];
//            for (J1=0;J1<locNNZ; J1++) {
//                j = StokesB->J[I1 + J1];
//                StokesB->A[I1 + J1] *= StokesA->d[i]*StokesC->d[j];
//            }
//        }
//
//        // Scale C
//#pragma omp parallel for shared(StokesC,StokesA) private(I1,J1, i, j, locNNZ )
//        for (i=0;i<StokesC->neq; i++) {
//
//            StokesC->b[i] *= StokesC->d[i]; // scale RHS
//
//            I1     = StokesC->Ic[i];
//            locNNZ = StokesC->Ic[i+1] - StokesC->Ic[i];
//            for (J1=0;J1<locNNZ; J1++) {
//                j = StokesC->J[I1 + J1];
//                StokesC->A[I1 + J1] *= StokesC->d[i]*StokesA->d[j];
//            }
//        }
//
//    }
//
//        MinMaxArray(StokesA->d, 1, StokesA->neq, "diag. A" );
        
        //--------------------------------------//

#ifndef _VG_
        if ( model.write_debug == 1 ) {
            
            char *filename;
            asprintf( &filename, "MatrixDecoupled.gzip_%dcpu.h5", n_th );
            printf("Writing Matrix file to disk...\n");
            
            // Fill in DD data structure
            OutputSparseMatrix OutputDDA, OutputDDB, OutputDDC, OutputDDD;
            
            OutputDDA.V = StokesA->A;
            OutputDDA.Ic = StokesA->Ic;
            OutputDDA.J = StokesA->J;
            OutputDDA.b = StokesA->b;
            OutputDDB.V = StokesB->A;
            OutputDDB.Ic = StokesB->Ic;
            OutputDDB.J = StokesB->J;
            OutputDDB.b = StokesB->b;
            OutputDDC.V = StokesC->A;
            OutputDDC.Ic = StokesC->Ic;
            OutputDDC.J = StokesC->J;
            OutputDDC.b = StokesC->b;
            OutputDDD.V = StokesD->A;
            OutputDDD.Ic = StokesD->Ic;
            OutputDDD.J = StokesD->J;
            OutputDDD.b = StokesD->b;
            OutputDDA.eta_cell = mesh->eta_n;
            OutputDDA.params[0] = nx;
            OutputDDA.params[1] = nz;
            OutputDDA.params[2] = mesh->dx;
            OutputDDA.params[3] = mesh->dz;
            OutputDDA.eqn_u = DoodzMalloc(nx*nzvx*sizeof(int));
            OutputDDA.eqn_v = DoodzMalloc(nxvz*nz*sizeof(int));
            OutputDDA.eqn_p = DoodzMalloc(ncx*ncz*sizeof(int));
            
            for( l=0; l<nzvx; l++) {
                for( k=0; k<nx; k++) {
                    cc = k + l*nx;
                    c1 = k + l*nx;
                    //                OutputDD.eqn_u[cc]=c1;
                    OutputDDA.eqn_u[cc]=Stokes->eqn_u[cc];
                }
            }
            for( l=0; l<nz; l++) {
                for( k=0; k<nxvz; k++) {
                    cc = k + l*nxvz;
                    c3 = nx*nzvx + k + l*nxvz;
                    //                OutputDD.eqn_v[cc]=c3;
                    OutputDDA.eqn_v[cc]=Stokes->eqn_v[cc];
                }
            }
            for( l=0; l<ncz; l++) {
                for( k=0; k<ncx; k++) {
                    cc = k + l*ncx;
                    c2 = nx*nzvx + nxvz*nz + k + l*ncx;
                    //                OutputDD.eqn_p[cc]=c2;
                    OutputDDA.eqn_p[cc]=Stokes->eqn_p[cc];
                }
            }
            
            // Send data to file
            create_output_hdf5( filename );
            AddGroup_to_hdf5( filename, "model" );
            AddGroup_to_hdf5( filename, "matrix" );
            AddGroup_to_hdf5( filename, "numbering" );
            AddGroup_to_hdf5( filename, "fields" );
            AddFieldToGroup_generic( _TRUE_, filename, "model", "params" , 'd', 4, OutputDDA.params,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "numbering", "eqn_u" , 'i', nx*nzvx, OutputDDA.eqn_u,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "numbering", "eqn_v" , 'i', nxvz*nz, OutputDDA.eqn_v,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "numbering", "eqn_p" , 'i', ncx*ncz, OutputDDA.eqn_p,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "eta_n" , 'd', ncx*ncz, mesh->eta_n,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "eta_s" , 'd', nx*nz, mesh->eta_s,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_s" , 'c', nx*nz, mesh->BCg.type,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_n" , 'c', ncx*ncz, mesh->BCp.type,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_u" , 'c', nx*nzvx, mesh->BCu.type,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_v" , 'c', nxvz*nz, mesh->BCv.type,  1 );
            
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "IA" , 'i', StokesA->neq+1, OutputDDA.Ic,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "JA" , 'i', nnzcA,  OutputDDA.J,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "VA" , 'd', nnzcA,  OutputDDA.V,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "IB" , 'i', StokesB->neq+1, OutputDDB.Ic,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "JB" , 'i', nnzcB,  OutputDDB.J,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "VB" , 'd', nnzcB,  OutputDDB.V,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "IC" , 'i', StokesC->neq+1, OutputDDC.Ic,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "JC" , 'i', nnzcC,  OutputDDC.J,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "VC" , 'd', nnzcC,  OutputDDC.V,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "ID" , 'i', StokesD->neq+1, OutputDDD.Ic,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "JD" , 'i', nnzcD,  OutputDDD.J,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "VD" , 'd', nnzcD,  OutputDDD.V,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "eta_cell", 'd', ncx*ncz, OutputDDA.eta_cell,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "rhs_mom" , 'd', StokesA->neq, OutputDDA.b,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "rhs_cont", 'd', StokesC->neq, OutputDDC.b,  1 );
            //            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "rhs", 'd', Stokes->neq, OutputDD.b,  1 );
            DoodzFree(OutputDDA.eqn_u);
            DoodzFree(OutputDDA.eqn_v);
            DoodzFree(OutputDDA.eqn_p);
            free(filename);
        }
#endif
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Xjacobian_InnerNodesDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B  ) {
    
    double  uC=0.0;
    double  uS=0.0,  uN=0.0,  uW=0.0,  uE=0.0,  vSW=0.0,  vSE=0.0,  vNW=0.0,  vNE=0.0, pE=0.0, pW=0.0;
    double uSW=0.0, uSE=0.0, uNW=0.0, uNE=0.0, vSWW=0.0, vSEE=0.0, vNWW=0.0, vNEE=0.0;
    
    double etaE  = mesh->eta_n[c2+1], etaW  = mesh->eta_n[c2], etaN  = mesh->eta_s[c1], etaS  = mesh->eta_s[c1-nx];
    double ExxE_pwl = mesh->exx_pwl_n[c2+1], ExxW_pwl = mesh->exx_pwl_n[c2], ExxN_pwl = mesh->exx_pwl_s[c1], ExxS_pwl = mesh->exx_pwl_s[c1-nx];
    double ExzE_pwl = mesh->exz_pwl_n[c2+1], ExzW_pwl = mesh->exz_pwl_n[c2], ExzN_pwl = mesh->exz_pwl_s[c1], ExzS_pwl = mesh->exz_pwl_s[c1-nx];
    double A2E_pwl  = mesh->A2_pwl_n [c2+1], A2W_pwl  = mesh->A2_pwl_n [c2], A2N_pwl  = mesh->A2_pwl_s [c1], A2S_pwl  = mesh->A2_pwl_s [c1-nx];
    
    // Total strain rate components
    double ExxW = one_dx*( mesh->u_in[c1  ] -  mesh->u_in[c1-1] );
    double ExxE = one_dx*( mesh->u_in[c1+1] -  mesh->u_in[c1  ] );
    double ExzS = 0.5*one_dz*( mesh->u_in[c1   ] -  mesh->u_in[c1-nx] ) + 0.5*one_dx*( mesh->v_in[c3-nxvz+1] - mesh->v_in[c3-nxvz]);
    double ExzN = 0.5*one_dz*( mesh->u_in[c1+nx] -  mesh->u_in[c1   ] ) + 0.5*one_dx*( mesh->v_in[c3+1     ] - mesh->v_in[c3     ]);
    
    //    printf("ExxE_pwl=%2.2e ExxW_pwl=%2.2e ExxN_pwl=%2.2e ExxS_pwl=%2.2e ExxW=%2.2e ExxE=%2.2e\n", ExxE_pwl, ExxW_pwl, ExxN_pwl, ExxS_pwl, ExxW, ExxE);
    //    printf("ExzE_pwl=%2.2e ExzW_pwl=%2.2e ExzN_pwl=%2.2e ExzN=%2.2e ExzS_pwl=%2.2e ExzS=%2.2e\n", ExzE_pwl, ExzW_pwl, ExzN_pwl, ExzN, ExzS_pwl, ExzS);
    
    
    // West derivatives
    double detaWdVxC   =  A2W_pwl * one_dx * ExxW_pwl;
    double detaWdVxW   = -detaWdVxC;
    double detaWdVxS   = -0.125*A2W_pwl*ExzW_pwl*one_dz;
    double detaWdVxN   = -detaWdVxS;
    double detaWdVxSW  =  detaWdVxS;
    double detaWdVxNW  = -detaWdVxS;
    double detaWdVzSW  =  0.0;
    double detaWdVzSE  =  0.125*A2W_pwl*ExzW_pwl*one_dx;
    double detaWdVzNW  =  0.0;
    double detaWdVzNE  =  detaWdVzSE;
    double detaWdVzSWW = -detaWdVzSE;
    double detaWdVzNWW = -detaWdVzSE;
    
    //    printf("VxC=%2.2e VxW=%2.2e VxS=%2.2e  VxN=%2.2e VxSW=%2.2e VxNW=%2.2e VzSW=%2.2e VzSE=%2.2e VzNW=%2.2e VzNE=%2.2e VzSWW=%2.2e VzNWW=%2.2e\n", detaWdVxC, detaWdVxW, detaWdVxS, detaWdVxN, detaWdVxSW, detaWdVxNW, detaWdVzSW, detaWdVzSE, detaWdVzNW, detaWdVzNE, detaWdVzSWW, detaWdVzNWW);
    
    // East derivatives
    double detaEdVxC   = -A2E_pwl * one_dx * ExxE_pwl;
    double detaEdVxE   = -detaEdVxC;
    double detaEdVxS   = -0.125*A2E_pwl*ExzE_pwl*one_dz;
    double detaEdVxN   = -detaEdVxS;
    double detaEdVxSE  =  detaEdVxS;
    double detaEdVxNE  = -detaEdVxS;
    double detaEdVzSW  = -0.125*A2E_pwl*ExzE_pwl*one_dx;
    double detaEdVzSE  =  0.0;
    double detaEdVzNW  =  detaEdVzSW;
    double detaEdVzNE  =  0.0;
    double detaEdVzSEE = -detaEdVzSW;
    double detaEdVzNEE = -detaEdVzSW;
    
    // South derivatives
    double detaSdVxC   =  0.5*A2S_pwl*ExzS_pwl*one_dz;
    double detaSdVxS   = -detaSdVxC;
    double detaSdVxW   = -0.25*A2S_pwl*ExxS_pwl*one_dx;
    double detaSdVxE   = -detaSdVxW;
    double detaSdVxSW  =  detaSdVxW;
    double detaSdVxSE  = -detaSdVxW;
    double detaSdVzSW  = -0.5*A2S_pwl*ExzS_pwl*one_dx;
    double detaSdVzSE  = -detaSdVzSW;
    
    // North derivatives
    double detaNdVxC   = -0.5*A2N_pwl*ExzN_pwl*one_dz;
    double detaNdVxN   = -detaNdVxC;
    double detaNdVxW   = -0.25*A2N_pwl*ExxN_pwl*one_dx;
    double detaNdVxE   = -detaNdVxW;
    double detaNdVxNW  =  detaNdVxW;
    double detaNdVxNE  = -detaNdVxW;
    double detaNdVzNW  = -0.5*A2N_pwl*ExzN_pwl*one_dx;
    double detaNdVzNE  = -detaNdVzNW;
    
    
    // dsxx/dx
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1+1] != 30 ) {
        uW  =  2.0*one_dx_dx * etaW         - comp*2.0/3.0*one_dx_dx * etaW;
        uW   +=  one_dx*(                  -2.0*detaWdVxW*ExxW)  + one_dz*( 2.0*detaNdVxW *ExzN -2.0*detaSdVxW*ExzS);
        uC  = -2.0*one_dx_dx * (etaE + etaW)  + comp*2.0/3.0*one_dx_dx * etaE + comp*2.0/3.0*one_dx_dx * etaW;
        uC +=  one_dx*( 2.0*detaEdVxC*ExxE   -2.0*detaWdVxC*ExxW);
        uE  =  2.0*one_dx_dx * etaE         - comp*2.0/3.0*one_dx_dx * etaE;
        uE   +=  one_dx*( 2.0*detaEdVxE*ExxE                  )  + one_dz*( 2.0*detaNdVxE *ExzN -2.0*detaSdVxE*ExzS);
    }
    
    if ( mesh->BCu.type[c1-1] == 30 && mesh->BCu.type[c1+1] != 30 ) {
        uC  = -2.0*one_dx_dx * (etaE) + comp*2.0/3.0*one_dx_dx * etaE;
        uC +=  one_dx*( 2.0*detaEdVxC*ExxE );
        uE  =  2.0*one_dx_dx *  etaE  - comp*2.0/3.0*one_dx_dx * etaE;
        uE +=  one_dx*( 2.0*detaEdVxE*ExxE                  )  + one_dz*( 2.0*detaNdVxE *ExzN -2.0*detaSdVxE*ExzS);
    }
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1+1] == 30 ) {
        uW  =  2.0*one_dx_dx *  etaW  - comp*2.0/3.0*one_dx_dx * etaW;
        uW +=  one_dx*(                  -2.0*detaWdVxW*ExxW)  + one_dz*( 2.0*detaNdVxW *ExzN -2.0*detaSdVxW*ExzS);
        uC  = -2.0*one_dx_dx * (etaW) + comp*2.0/3.0*one_dx_dx * etaW;
        uC +=  one_dx*( -2.0*detaWdVxC*ExxW);
    }
    
    if ( mesh->BCu.type[c1-nx-1] != 30) uSW  +=  one_dx*(-2.0*detaWdVxSW*ExxW)                   - one_dz*( 2.0*detaSdVxSW*ExzS);
    if ( mesh->BCu.type[c1-nx+1] != 30) uSE  +=  one_dx*( 2.0*detaEdVxSE*ExxE)                   - one_dz*( 2.0*detaSdVxSE*ExzS);
    if ( mesh->BCu.type[c1+nx-1] != 30) uNW  +=  one_dx*(-2.0*detaWdVxNW*ExxW)                   + one_dz*( 2.0*detaNdVxNW*ExzN);
    if ( mesh->BCu.type[c1+nx+1] != 30) uNE  +=  one_dx*( 2.0*detaEdVxNE*ExxE)                   + one_dz*( 2.0*detaNdVxNE*ExzN);
    
    // eta*dvx/dz
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) {
            uS   =  one_dz_dz * etaS;
            uS  +=  one_dx*( 2.0*detaEdVxS*ExxE   -2.0*detaWdVxS*ExxW)  + one_dz*(-2.0*detaSdVxS *ExzS);
            uC  +=  -one_dz_dz * etaS;
            uC  +=   one_dz*( -2.0*detaSdVxC*ExzS);
            
        }
        if ( mesh->BCu.type[c1+nx] != 13 ) {
            uN   =   one_dz_dz * etaN;
            uN  +=   one_dx*( 2.0*detaEdVxN*ExxE   -2.0*detaWdVxN*ExxW)  + one_dz*( 2.0*detaNdVxN *ExzN);
            uC  +=  -one_dz_dz * etaN;
            uC  +=   one_dz*( 2.0*detaNdVxC *ExzN );
        }
    }
    if ( mesh->BCu.type[c1-nx] == 30 && mesh->BCu.type[c1+nx] != 30 ) {
        if ( mesh->BCu.type[c1+nx] != 13 ) { uC  +=  -one_dz_dz * etaN; uC   +=  one_dz*( 2.0*detaNdVxC *ExzN); }
        if ( mesh->BCu.type[c1+nx] != 13 ) { uN   =   one_dz_dz * etaN; uN   +=  one_dz*( 2.0*detaNdVxN *ExzN); }
        
    }
    if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1+nx] == 30 ) {
        if ( mesh->BCu.type[c1-nx] != 13 ) { uC  +=  -one_dz_dz * etaS; uC   +=  one_dz*( -2.0*detaSdVxC*ExzS); }
        if ( mesh->BCu.type[c1-nx] != 13 ) { uS   =   one_dz_dz * etaS; uS   +=  one_dz*(-2.0*detaSdVxS *ExzS); }
    }
    
    if ( mesh->BCu.type[c1-nx] == 11 ) { uC  +=  -one_dz_dz * etaS; uC   +=  one_dz*( -2.0*detaSdVxC*ExzS); }
    if ( mesh->BCu.type[c1+nx] == 11 ) { uC  +=  -one_dz_dz * etaN; uC   +=  one_dz*( 2.0*detaNdVxC *ExzN); }
    
    // eta*dvz/dx WW
    if ( mesh->BCv.type[c3-nxvz-1] != 30 && mesh->BCv.type[c3-1] != 30 )  {
        vSWW +=  one_dx*(-2.0*detaWdVzSWW*ExxW);
        vNWW +=  one_dx*(-2.0*detaWdVzNWW*ExxW);
    }
    
    // eta*dvz/dx EE
    if ( mesh->BCv.type[c3-nxvz+2] != 30 && mesh->BCv.type[c3+2] != 30 )  {
        vSEE +=  one_dx*( 2.0*detaEdVzSEE*ExxE);
        vNEE +=  one_dx*( 2.0*detaEdVzNEE*ExxE);
    }
    
    // eta*dvz/dx S
    if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
        vSE  = -one_dx_dz * etaS + comp*2.0/3.0*one_dx_dz * etaE;
        vSW  =  one_dx_dz * etaS - comp*2.0/3.0*one_dx_dz * etaW;
        vSW +=  one_dx*( 2.0*detaEdVzSW*ExxE  -2.0*detaWdVzSW*ExxW) + one_dz*(-2.0*detaSdVzSW*ExzS);
        vSE +=  one_dx*( 2.0*detaEdVzSE*ExxE  -2.0*detaWdVzSE*ExxW) + one_dz*(-2.0*detaSdVzSE*ExzS);
    }
    
    // eta*dvz/dx N
    if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3] != 30  ) {
        vNE  =  one_dx_dz * etaN - comp*2.0/3.0*one_dx_dz * etaE;
        vNW  = -one_dx_dz * etaN + comp*2.0/3.0*one_dx_dz * etaW;
        vNW +=  one_dx*( 2.0*detaEdVzNW*ExxE  -2.0*detaWdVzNW*ExxW) + one_dz*( 2.0*detaNdVzNW*ExzN);
        vNE +=  one_dx*( 2.0*detaEdVzNE*ExxE  -2.0*detaWdVzNE*ExxW) + one_dz*( 2.0*detaNdVzNE*ExzN);
    }
    //    // Newton
    //    double A;
    //    double dx = 1.0/one_dx;
    //    double dy = 1.0/one_dz;
    //    double ExyN=ExzN_pwl, ExyS=ExzS_pwl;
    //    double ExyW=ExzW_pwl, ExyE=ExzE_pwl;
    //    double ExxS=ExxS_pwl, ExxN=ExxN_pwl;
    
    //    A         =  A2W_pwl;
    //    detaWdVxC   =  A2W_pwl * one_dx * ExxW_pwl;
    //    detaWdVxW   = -detaWdVxC;
    //    detaWdVxS   = -0.1250*ExyW*A/dy; //-(A2W_pwl/8.0)*ExzW_pwl*one_dz
    //    detaWdVxN   = -detaWdVxS;
    //    detaWdVxSW  =  detaWdVxS;
    //    detaWdVxNW  = -detaWdVxS;
    //    detaWdVzSW  =  0.0;
    //    detaWdVzSE  =  0.1250*ExyW*A/dx;
    //    detaWdVzNW  =  0.0;
    //    detaWdVzNE  =  detaWdVzSE;
    //    detaWdVzSWW = -detaWdVzSE;
    //    detaWdVzNWW = -detaWdVzSE;
    //
    //    A         = A2E_pwl;
    //    detaEdVxC   = -ExxE*A/dx;
    //    detaEdVxE   = -detaEdVxC;
    //    detaEdVxS   = -0.1250*ExyE*A/dy;
    //    detaEdVxN   = -detaEdVxS;
    //    detaEdVxSE  =  detaEdVxS;
    //    detaEdVxNE  = -detaEdVxS;
    //    detaEdVzSW  = -0.1250*ExyE*A/dx;
    //    detaEdVzSE  =  0.0;
    //    detaEdVzNW  =  detaEdVzSW;
    //    detaEdVzNE  =  0.0;
    //    detaEdVzSEE = -detaEdVzSW;
    //    detaEdVzNEE = -detaEdVzSW;
    //
    //    A         = A2S_pwl;
    //    detaSdVxC   =  0.5*ExyS*A/dy;
    //    detaSdVxS   = -detaSdVxC;
    //    detaSdVxW   = -0.25*ExxS*A/dx;
    //    detaSdVxE   = -detaSdVxW;
    //    detaSdVxSW  =  detaSdVxW;
    //    detaSdVxSE  = -detaSdVxW;
    //    detaSdVzSW  = -0.5*ExyS*A/dx;
    //    detaSdVzSE  = -detaSdVzSW;
    //
    //    A         = A2N_pwl;
    //    detaNdVxC   = -0.5*ExyN*A/dy;
    //    detaNdVxN   = -detaNdVxC;
    //    detaNdVxW   = -0.25*ExxN*A/dx;
    //    detaNdVxE   = -detaNdVxW;
    //    detaNdVxNW  =  detaNdVxW;
    //    detaNdVxNE  = -detaNdVxW;
    //    detaNdVzNW  = -0.5*ExyN*A/dx;
    //    detaNdVzNE  = -detaNdVzNW;
    
    //    uC   -= (-(2.0*detaEdVxC*ExxE - 2.0*detaWdVxC*ExxW )/dx -(2.0*detaNdVxC*ExyN -2.0*detaSdVxC*ExyS)/dy);
    //    uW   -= ( (2.0*detaWdVxW*ExxW   )/dx - (2.0*detaNdVxW*ExyN-2.0*detaSdVxW*ExyS)/dy);
    //    uE   -= (-(2.0*detaEdVxE*ExxE   )/dx - (2.0*detaNdVxE*ExyN-2.0*detaSdVxE*ExyS)/dy);
    //    uS   -= (-(2.0*detaEdVxS*ExxE - 2.0*detaWdVxS*ExxW)/dx -(-2.0*detaSdVxS*ExyS)/dy);
    //    uN   -= (-(2.0*detaEdVxN*ExxE - 2.0*detaWdVxN*ExxW)/dx -( 2.0*detaNdVxN*ExyN)/dy);
    //    uSW  -= (  2.0*detaWdVxSW*ExxW/dx + 2.0*detaSdVxSW*ExyS/dy);
    //    uSE  -= ( -2.0*detaEdVxSE*ExxE/dx + 2.0*detaSdVxSE*ExyS/dy);
    //    uNW  -= (  2.0*detaWdVxNW*ExxW/dx - 2.0*detaNdVxNW*ExyN/dy);
    //    uNE  -= ( -2.0*detaEdVxNE*ExxE/dx - 2.0*detaNdVxNE*ExyN/dy);
    //    vSW  -= (-(2.0*detaEdVzSW*ExxE    - 2.0*detaWdVzSW*ExxW)/dx - (-2.0*detaSdVzSW*ExyS)/dy);
    //    vSE  -= (-(2.0*detaEdVzSE*ExxE    - 2.0*detaWdVzSE*ExxW)/dx - (-2.0*detaSdVzSE*ExyS)/dy);
    //    vNW  -= (-(2.0*detaEdVzNW*ExxE    - 2.0*detaWdVzNW*ExxW)/dx - ( 2.0*detaNdVzNW*ExyN)/dy);
    //    vNE  -= (-(2.0*detaEdVzNE*ExxE    - 2.0*detaWdVzNE*ExxW)/dx - ( 2.0*detaNdVzNE*ExyN)/dy);
    //    vSWW -= ( 2.0*detaWdVzSWW*ExxW/dx);
    //    vSEE -= (-2.0*detaEdVzSEE*ExxE/dx);
    //    vNWW -= ( 2.0*detaWdVzNWW*ExxW/dx);
    //    vNEE -= (-2.0*detaEdVzNEE*ExxE/dx);
    
    //    uC   +=  one_dx*( 2*detaEdVxC*ExxE   -2*detaWdVxC*ExxW)  + one_dz*( 2*detaNdVxC *ExzN -2*detaSdVxC*ExzS);
    //    uW   +=  one_dx*(                  -2*detaWdVxW*ExxW)  + one_dz*( 2*detaNdVxW *ExzN -2*detaSdVxW*ExzS);
    //    uE   +=  one_dx*( 2*detaEdVxE*ExxE                  )  + one_dz*( 2*detaNdVxE *ExzN -2*detaSdVxE*ExzS);
    //    uN   +=  one_dx*( 2*detaEdVxN*ExxE   -2*detaWdVxN*ExxW)  + one_dz*( 2*detaNdVxN *ExzN);
    //    uS   +=  one_dx*( 2*detaEdVxS*ExxE   -2*detaWdVxS*ExxW)  + one_dz*(-2*detaSdVxS *ExzS);
    //    uSW  +=  one_dx*(-2*detaWdVxSW*ExxW)                   - one_dz*( 2*detaSdVxSW*ExzS);
    //    uSE  +=  one_dx*( 2*detaEdVxSE*ExxE)                   - one_dz*( 2*detaSdVxSE*ExzS);
    //    uNW  +=  one_dx*(-2*detaWdVxNW*ExxW)                   + one_dz*( 2*detaNdVxNW*ExzN);
    //    uNE  +=  one_dx*( 2*detaEdVxNE*ExxE)                   + one_dz*( 2*detaNdVxNE*ExzN);
    //    vSW +=  one_dx*( 2*detaEdVzSW*ExxE  -2*detaWdVzSW*ExxW) + one_dz*(-2*detaSdVzSW*ExzS);
    //    vSE +=  one_dx*( 2*detaEdVzSE*ExxE  -2*detaWdVzSE*ExxW) + one_dz*(-2*detaSdVzSE*ExzS);
    //    vNW +=  one_dx*( 2*detaEdVzNW*ExxE  -2*detaWdVzNW*ExxW) + one_dz*( 2*detaNdVzNW*ExzN);
    //    vNE +=  one_dx*( 2*detaEdVzNE*ExxE  -2*detaWdVzNE*ExxW) + one_dz*( 2*detaNdVzNE*ExzN);
    //    vSWW +=  one_dx*(-2*detaWdVzSWW*ExxW);
    //    vSEE +=  one_dx*( 2*detaEdVzSEE*ExxE);
    //    vNWW +=  one_dx*(-2*detaWdVzNWW*ExxW);
    //    vNEE +=  one_dx*( 2*detaEdVzNEE*ExxE);
    //
    // Pressure gradient
    if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[c2] != 30 ) {
        pW  =   one_dx;
        pE  =  -pW;
    }
    
    // Minus sign of all coefficients
    uS=-uS, uN=-uN, uW=-uW, uE=-uE, uC=-uC, vSW=-vSW, vSE=-vSE, vNW=-vNW, vNE=-vNE, pE=-pE, pW=-pW;
    uSW=-uSW, uSE=-uSE, uNW=-uNW, uNE=-uNE, vSWW=-vSWW, vSEE=-vSEE, vNWW=-vNWW, vNEE=-vNEE;
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        // dsxx/dx - normal stencil
        
        // uSW (Newton)
        if (mesh->BCu.type[c1-nx-1] != 30) {
            if (mesh->BCu.type[c1-nx-1] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx-1], &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[c1-nx-1], mesh->BCu.val[c1-nx-1], StokesA->bbc);
            else                                    AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx-1], &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[c1-nx-1], 2.0*mesh->BCu.val[c1-nx-1], StokesA->bbc);
        }
        
        // uS
        if (mesh->BCu.type[c1-nx] != 30) {
            if (mesh->BCu.type[c1-nx] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], mesh->BCu.val[c1-nx], StokesA->bbc);
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx], &(nnzc2A[ith]), uS*celvol, mesh->BCu.type[c1-nx], 2.0*mesh->BCu.val[c1-nx], StokesA->bbc);
        }
        
        // uSE (Newton)
        if (mesh->BCu.type[c1-nx+1] != 30) {
            if (mesh->BCu.type[c1-nx+1] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx+1], &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[c1-nx+1], mesh->BCu.val[c1-nx+1], StokesA->bbc);
            else                                    AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx+1], &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[c1-nx+1], 2.0*mesh->BCu.val[c1-nx+1], StokesA->bbc);
        }
        
        // uW
        if (mesh->BCu.type[c1-1]  != 30)     AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],  &(nnzc2A[ith]), uW*celvol, mesh->BCu.type[c1-1],  mesh->BCu.val[c1-1],  StokesA->bbc);
        
        // uC
        AddCoeff2( JtempA[ith], AtempA[ith], eqn, eqn,                  &(nnzc2A[ith]), uC*celvol, mesh->BCu.type[c1],    mesh->BCu.val[c1],    StokesA->bbc);
        
        // uE
        if (mesh->BCu.type[c1+1]  != 30)     AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+1],  &(nnzc2A[ith]), uE*celvol, mesh->BCu.type[c1+1],  mesh->BCu.val[c1+1],  StokesA->bbc);
        
        // uNW (Newton)
        if (mesh->BCu.type[c1+nx-1] != 30) {
            if (mesh->BCu.type[c1+nx-1] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx-1], &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1],   mesh->BCu.val[c1+nx-1], StokesA->bbc);
            else                                    AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx-1], &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], 2.0*mesh->BCu.val[c1+nx-1], StokesA->bbc);
        }
        
        // uN
        if (mesh->BCu.type[c1+nx] != 30) {
            if (mesh->BCu.type[c1+nx] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], mesh->BCu.val[c1+nx], StokesA->bbc);
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx], &(nnzc2A[ith]), uN*celvol, mesh->BCu.type[c1+nx], 2*mesh->BCu.val[c1+nx], StokesA->bbc);
        }
        
        // uNE (Newton)
        if (mesh->BCu.type[c1+nx+1] != 30) {
            if (mesh->BCu.type[c1+nx+1] != 11) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx+1], &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[c1+nx+1],   mesh->BCu.val[c1+nx+1], StokesA->bbc);
            else                                    AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx+1], &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[c1+nx+1], 2.0*mesh->BCu.val[c1+nx+1], StokesA->bbc);
        }
        //--------------------
        
        // vSWW (Newton)
        if ( mesh->BCv.type[c3-nxvz-1] != 30 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz-1],   &(nnzc2A[ith]), vSWW*celvol, mesh->BCv.type[c3-nxvz-1],   mesh->BCv.val[c3-nxvz-1],   StokesA->bbc);
        
        // vSW && vSE
        if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3-nxvz+1] != 30 )  {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz],   &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz],   mesh->BCv.val[c3-nxvz],   StokesA->bbc);
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz+1], &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1], mesh->BCv.val[c3-nxvz+1], StokesA->bbc);
        }
        
        // vSEE (Newton)
        if ( mesh->BCv.type[c3-nxvz+2] != 30 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz+2], &(nnzc2A[ith]), vSEE*celvol, mesh->BCv.type[c3-nxvz+2], mesh->BCv.val[c3-nxvz+2], StokesA->bbc);
        
        // vNWW (Newton)
        if ( mesh->BCv.type[c3-1] != 30) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],        &(nnzc2A[ith]), vNWW*celvol, mesh->BCv.type[c3-1],        mesh->BCv.val[c3-1],        StokesA->bbc);
        
        // vNW && vNE
        if ( mesh->BCv.type[c3] != 30 && mesh->BCv.type[c3+1] != 30 )  {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3],        &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[c3],        mesh->BCv.val[c3],        StokesA->bbc);
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],      &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[c3+1],      mesh->BCv.val[c3+1],      StokesA->bbc);
        }
        
        // vNEE (Newton)
        if ( mesh->BCv.type[c3+2] != 30 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+2],      &(nnzc2A[ith]), vNEE*celvol, mesh->BCv.type[c3+2],      mesh->BCv.val[c3+2],      StokesA->bbc);
        
        //--------------------
        
        // pE && pW
        if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[c2] != 30 ) {
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2] - Stokes->neq_mom,   &(nnzc2B[ith]), pW*celvol, mesh->BCp.type[c2],   mesh->BCp.val[c2],   StokesB->bbc);
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2+1] - Stokes->neq_mom, &(nnzc2B[ith]), pE*celvol, mesh->BCp.type[c2+1], mesh->BCp.val[c2+1], StokesB->bbc);
        }
    }
    else {
        // Residual function
        StokesA->F[eqn] = uC*u[c1];
        if ( mesh->BCp.type[c2+1] != 30 && mesh->BCp.type[c2] != 30 ) {
            StokesA->F[eqn]  += pW*p[c2] + pE*p[c2+1];
        }
        if ( mesh->BCv.type[c3-nxvz+1] != 30 && mesh->BCv.type[c3-nxvz] != 30  ) {
            StokesA->F[eqn] += vSW*v[c3-nxvz] + vSE*v[c3-nxvz+1];
        }
        if ( mesh->BCv.type[c3+1] != 30 && mesh->BCv.type[c3] != 30 ) {
            StokesA->F[eqn] += vNW*v[c3] + vNE*v[c3+1];
        }
        if ( mesh->BCu.type[c1-nx] != 30 && mesh->BCu.type[c1-nx] != 11  ) StokesA->F[eqn] += uS*u[c1-nx];
        if ( mesh->BCu.type[c1+nx] != 30 && mesh->BCu.type[c1+nx] != 11  ) StokesA->F[eqn] += uN*u[c1+nx];
        if ( mesh->BCu.type[c1-1]  != 30 ) StokesA->F[eqn] += uW*u[c1-1];
        if ( mesh->BCu.type[c1+1]  != 30 ) StokesA->F[eqn] += uE*u[c1+1];
        if ( mesh->BCu.type[c1-nx] == 11 ) StokesA->F[eqn] += -2*etaS*one_dz_dz*mesh->BCu.val[c1-nx];
        if ( mesh->BCu.type[c1+nx] == 11 ) StokesA->F[eqn] += -2*etaN*one_dz_dz*mesh->BCu.val[c1+nx];
        StokesA->F[eqn] -= (StokesA->b[eqn]);// + Stokes->bbc[eqn];
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Zjacobian_InnerNodesDecoupled( SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, int Assemble, int lev, int stab, int comp, double om, int sign, params model, double one_dx, double one_dz, double one_dx_dx, double one_dz_dz, double one_dx_dz, double celvol, grid* mesh, int ith, int c1, int c2, int c3, int nx, int ncx, int nxvz, int eqn, double* u, double* v, double* p, int **JtempA, double **AtempA, int *nnzc2A, int **JtempB, double **AtempB, int *nnzc2B ) {
    
    double  vC=0.0;
    double uSW=0.0, uSE=0.0, uNW=0.0, uNE=0.0,   vS=0.0,   vW=0.0,   vE=0.0,   vN=0.0, pN=0.0, pS=0.0;
    double vSW=0.0, vSE=0.0, vNW=0.0, vNE=0.0, uSSW=0.0, uSSE=0.0, uNNW=0.0, uNNE=0.0;
    
    // Coefficients
    double etaE  = mesh->eta_s[c1], etaW  = mesh->eta_s[c1-1], etaN  = mesh->eta_n[c2+ncx], etaS  = mesh->eta_n[c2];
    double ExxE_pwl = mesh->exx_pwl_s[c1], ExxW_pwl = mesh->exx_pwl_s[c1-1], ExxN_pwl = mesh->exx_pwl_n[c2+ncx], ExxS_pwl = mesh->exx_pwl_n[c2];
    double ExzE_pwl = mesh->exz_pwl_s[c1], ExzW_pwl = mesh->exz_pwl_s[c1-1], ExzN_pwl = mesh->exz_pwl_n[c2+ncx], ExzS_pwl = mesh->exz_pwl_n[c2];
    double A2E_pwl  = mesh->A2_pwl_s [c1], A2W_pwl  = mesh->A2_pwl_s [c1-1], A2N_pwl  = mesh->A2_pwl_n [c2+ncx], A2S_pwl  = mesh->A2_pwl_n [c2];
    
    // Total strain rate components
    double EzzN = -one_dx*( mesh->u_in[c1+nx] - mesh->u_in[c1+nx-1] );
    double EzzS = -one_dx*( mesh->u_in[c1   ] - mesh->u_in[c1   -1] );
    double ExzW = 0.5*one_dx*( mesh->v_in[c3  ] - mesh->v_in[c3-1]) + 0.5*one_dz*( mesh->u_in[c1+nx-1] - mesh->u_in[c1-1] );
    double ExzE = 0.5*one_dx*( mesh->v_in[c3+1] - mesh->v_in[c3  ]) + 0.5*one_dz*( mesh->u_in[c1+nx]   - mesh->u_in[c1  ] );
    
    //    printf("EzzE_pwl=%2.2e EzzW_pwl=%2.2e EzzN_pwl=%2.2e EzzN=%2.2e EzzS_pwl=%2.2e EzzS=%2.2e\n", -ExxE_pwl, -ExxW_pwl, -ExxN_pwl, EzzN, -ExxS_pwl, EzzS);
    //    printf("ExzE_pwl=%2.2e ExzE=%2.2e ExzW_pwl=%2.2e ExzW=%2.2e ExzN_pwl=%2.2e ExzS_pwl=%2.2e\n", ExzE_pwl, ExzE, ExzW_pwl, ExzW, ExzN_pwl, ExzS_pwl);
    
    // South derivatives
    double detaSdVzC   =  A2S_pwl*(-ExxS_pwl)*one_dz;
    double detaSdVzS   = -detaSdVzC;
    double detaSdVzW   = -0.125*A2S_pwl*ExzS_pwl*one_dx;
    double detaSdVzE   = -detaSdVzW;
    double detaSdVzSW  =  detaSdVzW;
    double detaSdVzSE  = -detaSdVzW;
    double detaSdVxSSW = -0.125*A2S_pwl*ExzS_pwl*one_dz;
    double detaSdVxSSE =  detaSdVxSSW;
    double detaSdVxSW  =  0.0;
    double detaSdVxSE  =  0.0;
    double detaSdVxNW  = -detaSdVxSSW;
    double detaSdVxNE  = -detaSdVxSSW;
    
    // North derivatives
    double detaNdVzC   = -A2N_pwl*(-ExxN_pwl)*one_dz;
    double detaNdVzN   = -detaNdVzC;
    double detaNdVzW   = -0.125*A2N_pwl*ExzN_pwl*one_dx;
    double detaNdVzE   = -detaNdVzW;
    double detaNdVzNW  =  detaNdVzW;
    double detaNdVzNE  = -detaNdVzW;
    double detaNdVxNNW =  0.125*A2N_pwl*ExzN_pwl*one_dz;
    double detaNdVxNNE =  detaNdVxNNW;
    double detaNdVxSW  = -detaNdVxNNW;
    double detaNdVxSE  = -detaNdVxNNW;
    double detaNdVxNW  = 0.0;
    double detaNdVxNE  = 0.0;
    
    // West derivatives
    double detaWdVzC   =  0.5*A2W_pwl*ExzW_pwl*one_dx;
    double detaWdVzW   = -detaWdVzC;
    double detaWdVzS   = -0.25*A2W_pwl*(-ExxW_pwl)*one_dz;
    double detaWdVzN   = -detaWdVzS;
    double detaWdVzSW  =  detaWdVzS;
    double detaWdVzNW  = -detaWdVzS;
    double detaWdVxSW  = -0.5*A2W_pwl*ExzW_pwl*one_dz;
    double detaWdVxNW  = -detaWdVxSW;
    
    // East derivatives
    double detaEdVzC   = -0.5*A2E_pwl*ExzE_pwl*one_dx;
    double detaEdVzE   = -detaEdVzC;
    double detaEdVzS   = -0.25*A2E_pwl*(-ExxE_pwl)*one_dz;
    double detaEdVzN   = -detaEdVzS;
    double detaEdVzSE  =  detaEdVzS;
    double detaEdVzNE  = -detaEdVzS;
    double detaEdVxSE  = -0.5*A2E_pwl*ExzE_pwl*one_dz;
    double detaEdVxNE  = -detaEdVxSE;
    
    // (eta du/dx) SS
    if ( mesh->BCu.type[c1-nx-1] != 30 && mesh->BCu.type[c1-nx] != 30 ) {
        uSSW += one_dz*(-2.0*detaSdVxSSW*EzzS);
        uSSE += one_dz*(-2.0*detaSdVxSSE*EzzS);
    }
    
    // (eta du/dx) S
    if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
        uSW =  one_dx_dz * etaW - comp*2.0/3.0*etaS*one_dx_dz;
        uSE = -one_dx_dz * etaE + comp*2.0/3.0*etaS*one_dx_dz;
        uSW  += one_dz*( 2.0*detaNdVxSW*EzzN-2.0*detaSdVxSW*EzzS) + one_dx*(-2.0*detaWdVxSW*ExzW);
        uSE  += one_dz*( 2.0*detaNdVxSE*EzzN-2.0*detaSdVxSE*EzzS) + one_dx*( 2.0*detaEdVxSE*ExzE);
    }
    
    // (eta du/dx) N
    if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
        uNW = -one_dx_dz * etaW + comp*2.0/3.0*etaN*one_dx_dz;
        uNE =  one_dx_dz * etaE - comp*2.0/3.0*etaN*one_dx_dz;
        uNW  += one_dz*( 2.0*detaNdVxNW*EzzN-2*detaSdVxNW*EzzS) + one_dx*(-2.0*detaWdVxNW*ExzW);
        uNE  += one_dz*( 2.0*detaNdVxNE*EzzN-2*detaSdVxNE*EzzS) + one_dx*( 2.0*detaEdVxNE*ExzE);
    }
    
    // (eta du/dx) NN
    if ( mesh->BCu.type[c1+2*nx-1] != 30 && mesh->BCu.type[c1+2*nx] != 30 ) {
        uNNW += one_dz*( 2.0*detaNdVxNNW*EzzN);
        uNNE += one_dz*( 2.0*detaNdVxNNE*EzzN);
    }
    
    // dsyy/dz
    if ( mesh->BCv.type[c3-nxvz] != 30 && mesh->BCv.type[c3+nxvz] !=30 ) {
        vS  =  2.0*one_dz_dz * etaS;//        - comp*2.0/3.0*etaS*one_dz_dz;
        vS += one_dz*(-2.0*detaSdVzS*EzzS)                  + one_dx*( 2.0*detaEdVzS*ExzE - 2.0*detaWdVzS*ExzW);
        vC  = -2.0*one_dz_dz * (etaN + etaS);// + comp*2.0/3.0*etaN*one_dz_dz + comp*2.0/3.0*etaS*one_dz_dz;
        vC += one_dz*( 2.0*detaNdVzC*EzzN  -2.0*detaSdVzC*EzzS);
        vN  =  2.0*one_dz_dz * etaN;//        - comp*2.0/3.0*etaN*one_dz_dz;
        vN += one_dz*( 2.0*detaNdVzN*EzzN)                  + one_dx*( 2.0*detaEdVzN*ExzE - 2.0*detaWdVzN*ExzW);
        
    }
    if ( mesh->BCv.type[c3-nxvz] == 30 && mesh->BCv.type[c3+nxvz] != 30 ) {
        vC  = -2.0*one_dz_dz * (etaN) + comp*2.0/3.0*etaN*one_dz_dz;
        vC += one_dz*( 2.0*detaNdVzC*EzzN);
        vN  =  2.0*one_dz_dz * etaN   - comp*2.0/3.0*etaN*one_dz_dz;
        vN += one_dz*( 2.0*detaNdVzN*EzzN)                  + one_dx*( 2.0*detaEdVzN*ExzE - 2.0*detaWdVzN*ExzW);
    }
    if ( mesh->BCv.type[c3-nxvz] != 30  && mesh->BCv.type[c3+nxvz] == 30 ) {
        vS  =  2.0*one_dz_dz * etaS   - comp*2.0/3.0*etaS*one_dz_dz;
        vS += one_dz*(-2.0*detaSdVzS*EzzS)                  + one_dx*( 2.0*detaEdVzS*ExzE - 2.0*detaWdVzS*ExzW);
        vC  = -2.0*one_dz_dz * (etaS) + comp*2.0/3.0*etaS*one_dz_dz;
        vC += one_dz*(                 -2.0*detaSdVzC*EzzS);
    }
    
    if ( mesh->BCv.type[c3-nxvz-1] != 30 ) vSW  += one_dz*(-2.0*detaSdVzSW*EzzS)                 - one_dx*( 2.0*detaWdVzSW*ExzW);
    if ( mesh->BCv.type[c3-nxvz+1] != 30 ) vSE  += one_dz*(-2.0*detaSdVzSE*EzzS)                 + one_dx*( 2.0*detaEdVzSE*ExzE);
    if ( mesh->BCv.type[c3+nxvz-1] != 30 ) vNW  += one_dz*( 2.0*detaNdVzNW*EzzN)                 - one_dx*( 2.0*detaWdVzNW*ExzW);
    if ( mesh->BCv.type[c3+nxvz+1] != 30 ) vNE  += one_dz*( 2.0*detaNdVzNE*EzzN)                 + one_dx*( 2.0*detaEdVzNE*ExzE);
    
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) {
            vW  =  one_dx_dx * etaW;
            vW += one_dz*( 2.0*detaNdVzW*EzzN  -2.0*detaSdVzW*EzzS) + one_dx*(                - 2.0*detaWdVzW*ExzW );
            vC += -one_dx_dx * etaW;
            vC +=  one_dx*( - 2.0*detaWdVzC*ExzW );
        }
        if ( mesh->BCv.type[c3+1] != 13 ) {
            vE  =  one_dx_dx * etaE;
            vE += one_dz*( 2.0*detaNdVzE*EzzN  -2.0*detaSdVzE*EzzS) + one_dx*( 2.0*detaEdVzE*ExzE                 );
            vC += -one_dx_dx * etaE;
            vC +=  one_dx*( 2.0*detaEdVzC*ExzE );
        }
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] == 30 && mesh->BCv.type[c3+1] != 30 ) {
        if ( mesh->BCv.type[c3+1] != 13 ) { vE  =  one_dx_dx * etaE; vE += one_dx*( 2*detaEdVzE*ExzE ); }
        if ( mesh->BCv.type[c3+1] != 13 ) { vC += -one_dx_dx * etaE; vC += one_dx*( 2*detaEdVzC*ExzE ); }
    }
    // d/dx eta*dvzdx
    if ( mesh->BCv.type[c3-1] != 30 && mesh->BCv.type[c3+1] == 30 ) {
        if ( mesh->BCv.type[c3-1] != 13 ) { vW  =  one_dx_dx * etaW; vW += one_dx*(- 2*detaWdVzW*ExzW ); }
        if ( mesh->BCv.type[c3-1] != 13 ) { vC += -one_dx_dx * etaW; vC += one_dx*(- 2*detaWdVzC*ExzW ); }
    }
    
    if ( mesh->BCv.type[c3-1] == 11 ) { vC  +=  -one_dx_dx * etaW; vC += one_dx*(- 2*detaWdVzC*ExzW );  }
    if ( mesh->BCv.type[c3+1] == 11 ) { vC  +=  -one_dx_dx * etaE; vC += one_dx*(  2*detaEdVzC*ExzE );  }
    
    // Stabilisation with density gradients
    if (stab==1) {
        double drhodz = (mesh->rho_app_n[c2+ncx] - mesh->rho_app_n[c2])*one_dz;
        //        double drhodx = (mesh->rho_app_s[c1]   - mesh->rho_app_s[c1-1])*one_dx;
        vC  += - 1.00 * om * model.dt * model.gz * drhodz;
        // Non-symmetric contibution to the system of equation
        //        uNW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uNE += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSW += - 0.25 * om * model.dt * model.gz * drhodx;
        //        uSE += - 0.25 * om * model.dt * model.gz * drhodx;
    }
    
    
    
    
    //    // Newton
    //    double A;
    //    double dx = 1.0/one_dx;
    //    double dy = 1.0/one_dz;
    //    double EyyN = EzzN, EyyS = EzzS, ExyW=ExzW, ExyE=ExzE ;
    //    double ExyS=ExzS_pwl, ExyN=ExzN_pwl;
    //    double EyyW=-ExxS_pwl, EyyE=-ExxE_pwl;
    //
    //    A = A2S_pwl;
    //    detaSdVzC   =  EyyS*A/dy;
    //    detaSdVzS   = -detaSdVzC;
    //    detaSdVzW   = -0.1250*ExyS*A/dx;
    //    detaSdVzE   = -detaSdVzW;
    //    detaSdVzSW  =  detaSdVzW;
    //    detaSdVzSE  = -detaSdVzW;
    //    detaSdVxSSW = -0.1250*ExyS*A/dy;
    //    detaSdVxSSE =  detaSdVxSSW;
    //    detaSdVxSW  = 0.0;
    //    detaSdVxSE  = 0.0;
    //    detaSdVxNW  = -detaSdVxSSW;
    //    detaSdVxNE  = -detaSdVxSSW;
    //
    //    A = A2N_pwl;
    //    detaNdVzC   = -EyyN*A/dy;
    //    detaNdVzN   = -detaNdVzC;
    //    detaNdVzW   = -0.1250*ExyN*A/dx;
    //    detaNdVzE   = -detaNdVzW;
    //    detaNdVzNW  =  detaNdVzW;
    //    detaNdVzNE  = -detaNdVzW;
    //    detaNdVxNNW =  0.1250*ExyN*A/dy;
    //    detaNdVxNNE =  detaNdVxNNW;
    //    detaNdVxSW  = -detaNdVxNNW;
    //    detaNdVxSE  = -detaNdVxNNW;
    //    detaNdVxNW  = 0.0;
    //    detaNdVxNE  = 0.0;
    //
    //    A = A2W_pwl;
    //    detaWdVzC   =  0.5*ExyW*A/dx;
    //    detaWdVzW   = -detaWdVzC;
    //    detaWdVzS   = -0.25*EyyW*A/dy;
    //    detaWdVzN   = -detaWdVzS;
    //    detaWdVzSW  =  detaWdVzS;
    //    detaWdVzNW  = -detaWdVzS;
    //    detaWdVxSW  = -0.5*ExyW*A/dy;
    //    detaWdVxNW  = -detaWdVxSW;
    //
    //    A = A2E_pwl;
    //    detaEdVzC   = -0.5*ExyE*A/dx;
    //    detaEdVzE   = -detaEdVzC;
    //    detaEdVzS   = -0.25*EyyE*A/dy;
    //    detaEdVzN   = -detaEdVzS;
    //    detaEdVzSE  =  detaEdVzS;
    //    detaEdVzNE  = -detaEdVzS;
    //    detaEdVxSE  = -0.5*ExyE*A/dy;
    //    detaEdVxNE  = -detaEdVxSE;
    //
    //    vC   -= (-(2.0*detaNdVzC*EyyN-2.0-2.0*detaSdVzC*EyyS-2.0)/dy-(2.0*detaEdVzC*ExyE- 2.0*detaWdVzC*ExyW)/dx);
    //    vW   -= (-(2.0*detaNdVzW *EyyN - 2.0*detaSdVzW*EyyS)/dy - (-2.0*detaWdVzW*ExyW)/dx);
    //    vE   -= (-(2.0*detaNdVzE *EyyN - 2.0*detaSdVzE*EyyS)/dy - ( 2.0*detaEdVzE*ExyE)/dx);
    //    vS   -= ( (2.0*detaSdVzS *EyyS-2.0)/dy - (2.0*detaEdVzS*ExyE-2.0*detaWdVzS*ExyW)/dx);
    //    vN   -= (-(2.0*detaNdVzN *EyyN+2.0)/dy - (2.0*detaEdVzN*ExyE-2.0*detaWdVzN*ExyW)/dx);
    //    vSW  -= (  2.0*detaSdVzSW*EyyS/dy + 2.0*detaWdVzSW*ExyW/dx);
    //    vSE  -= (  2.0*detaSdVzSE*EyyS/dy - 2.0*detaEdVzSE*ExyE/dx);
    //    vNW  -= ( -2.0*detaNdVzNW*EyyN/dy + 2.0*detaWdVzNW*ExyW/dx);
    //    vNE  -= ( -2.0*detaNdVzNE*EyyN/dy - 2.0*detaEdVzNE*ExyE/dx);
    //    uSW  -= (-(2.0*detaNdVxSW*EyyN-2.0*detaSdVxSW*EyyS)/dy - (-2.0*detaWdVxSW*ExyW)/dx);
    //    uSE  -= (-(2.0*detaNdVxSE*EyyN-2.0*detaSdVxSE*EyyS)/dy - ( 2.0*detaEdVxSE*ExyE)/dx);
    //    uNW  -=( -(2.0*detaNdVxNW*EyyN-2.0*detaSdVxNW*EyyS)/dy - (-2.0*detaWdVxNW*ExyW)/dx);
    //    uNE  -=( -(2.0*detaNdVxNE*EyyN-2.0*detaSdVxNE*EyyS)/dy - ( 2.0*detaEdVxNE*ExyE)/dx);
    //    uSSW -=(  2.0*detaSdVxSSW*EyyS/dy);
    //    uSSE -=(  2.0*detaSdVxSSE*EyyS/dy);
    //    uNNW -=( -2.0*detaNdVxNNW*EyyN/dy);
    //    uNNE -=( -2.0*detaNdVxNNE*EyyN/dy);
    
    //    vC   += one_dz*( 2*detaNdVzC*EzzN  -2*detaSdVzC*EzzS) + one_dx*( 2*detaEdVzC*ExzE - 2*detaWdVzC*ExzW );
    //    vW   += one_dz*( 2*detaNdVzW*EzzN  -2*detaSdVzW*EzzS) + one_dx*(                - 2*detaWdVzW*ExzW );
    //    vE   += one_dz*( 2*detaNdVzE*EzzN  -2*detaSdVzE*EzzS) + one_dx*( 2*detaEdVzE*ExzE                 );
    //    vS   += one_dz*(-2*detaSdVzS*EzzS)                  + one_dx*( 2*detaEdVzS*ExzE - 2*detaWdVzS*ExzW);
    //    vN   += one_dz*( 2*detaNdVzN*EzzN)                  + one_dx*( 2*detaEdVzN*ExzE - 2*detaWdVzN*ExzW);
    //    vSW  += one_dz*(-2*detaSdVzSW*EzzS)                 - one_dx*( 2*detaWdVzSW*ExzW);
    //    vSE  += one_dz*(-2*detaSdVzSE*EzzS)                 + one_dx*( 2*detaEdVzSE*ExzE);
    //    vNW  += one_dz*( 2*detaNdVzNW*EzzN)                 - one_dx*( 2*detaWdVzNW*ExzW);
    //    vNE  += one_dz*( 2*detaNdVzNE*EzzN)                 + one_dx*( 2*detaEdVzNE*ExzE);
    //    uSW  += one_dz*( 2*detaNdVxSW*EzzN-2*detaSdVxSW*EzzS) + one_dx*(-2*detaWdVxSW*ExzW);
    //    uSE  += one_dz*( 2*detaNdVxSE*EzzN-2*detaSdVxSE*EzzS) + one_dx*( 2*detaEdVxSE*ExzE);
    //    uNW  += one_dz*( 2*detaNdVxNW*EzzN-2*detaSdVxNW*EzzS) + one_dx*(-2*detaWdVxNW*ExzW);
    //    uNE  += one_dz*( 2*detaNdVxNE*EzzN-2*detaSdVxNE*EzzS) + one_dx*( 2*detaEdVxNE*ExzE);
    //    uSSW += one_dz*(-2*detaSdVxSSW*EzzS);
    //    uSSE += one_dz*(-2*detaSdVxSSE*EzzS);
    //    uNNW += one_dz*( 2*detaNdVxNNW*EzzN);
    //    uNNE += one_dz*( 2*detaNdVxNNE*EzzN);
    
    // Pressure gradient
    if ( mesh->BCp.type[c2] != 30 && mesh->BCp.type[c2+ncx] != 30 ) {
        pS  =   one_dz;
        pN  =  -one_dz;
    }
    
    // Minus sign of all coefficients
    uSW=-uSW, uSE=-uSE, uNW=-uNW, uNE=-uNE, vS=-vS, vW=-vW, vC=-vC, vE=-vE, vN=-vN, pN=-pN, pS=-pS;
    vSW=-vSW, vSE=-vSE, vNW=-vNW, vNE=-vNE, uSSW=-uSSW, uSSE=-uSSE, uNNW=-uNNW, uNNE=-uNNE;
    
    if ( Assemble == 1 ) {
        StokesA->b[eqn] *= celvol;
        StokesB->b[eqn] *= celvol;
        //--------------------
        
        // uSSW (Newton)
        if ( mesh->BCu.type[c1-1-nx] != 30  ) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1-nx],   &(nnzc2A[ith]), uSSW*celvol, mesh->BCu.type[c1-1-nx],    mesh->BCu.val[c1-1-nx],    StokesA->bbc );
        }
        
        // uSSE (Newton)
        if ( mesh->BCu.type[c1-nx] != 30 ) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-nx],     &(nnzc2A[ith]), uSSE*celvol, mesh->BCu.type[c1-nx],    mesh->BCu.val[c1-nx],    StokesA->bbc );
        }
        
        // uSW && uSE
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1-1],     &(nnzc2A[ith]), uSW*celvol, mesh->BCu.type[c1-1],    mesh->BCu.val[c1-1],    StokesA->bbc );
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1],       &(nnzc2A[ith]), uSE*celvol, mesh->BCu.type[c1],      mesh->BCu.val[c1],      StokesA->bbc );
        }
        
        // uNW && uNE
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx-1],  &(nnzc2A[ith]), uNW*celvol, mesh->BCu.type[c1+nx-1], mesh->BCu.val[c1+nx-1], StokesA->bbc );
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+nx],    &(nnzc2A[ith]), uNE*celvol, mesh->BCu.type[c1+nx],   mesh->BCu.val[c1+nx],   StokesA->bbc );
        }
        
        // uNNW (Newton)
        if ( mesh->BCu.type[c1+2*nx-1] != 30 ) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+2*nx-1],  &(nnzc2A[ith]), uNNW*celvol, mesh->BCu.type[c1+2*nx-1], mesh->BCu.val[c1+2*nx-1], StokesA->bbc );
        }
        
        // uNNE (Newton)
        if ( mesh->BCu.type[c1+2*nx] != 30 ) {
            AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_u[c1+2*nx],  &(nnzc2A[ith]), uNNE*celvol, mesh->BCu.type[c1+2*nx], mesh->BCu.val[c1+2*nx], StokesA->bbc );
        }
        
        //--------------------
        
        // vSW (Newton)
        if ( mesh->BCv.type[c3-nxvz-1] != 30  ) {
            if( mesh->BCv.type[c3-nxvz-1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz-1],     &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz-1],      mesh->BCv.val[c3-nxvz-1],    StokesA->bbc );
            else                                       AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz-1],     &(nnzc2A[ith]), vSW*celvol, mesh->BCv.type[c3-nxvz-1],    2*mesh->BCv.val[c3-nxvz-1],    StokesA->bbc );
        }
        
        // vS
        if ( mesh->BCv.type[c3-nxvz] != 30  ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz],  &(nnzc2A[ith]), vS*celvol, mesh->BCv.type[c3-nxvz], mesh->BCv.val[c3-nxvz], StokesA->bbc );
        
        // vSE (Newton)
        if ( mesh->BCv.type[c3-nxvz+1] != 30  ) {
            if ( mesh->BCv.type[c3-nxvz+1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz+1],     &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1],      mesh->BCv.val[c3-nxvz+1],    StokesA->bbc );
            else                                        AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-nxvz+1],     &(nnzc2A[ith]), vSE*celvol, mesh->BCv.type[c3-nxvz+1],    2*mesh->BCv.val[c3-nxvz+1],    StokesA->bbc );
        }
        
        // vW
        if ( mesh->BCv.type[c3-1] != 30 ) {
            if( mesh->BCv.type[c3-1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    mesh->BCv.val[c3-1],    StokesA->bbc );
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3-1],     &(nnzc2A[ith]), vW*celvol, mesh->BCv.type[c3-1],    2*mesh->BCv.val[c3-1],    StokesA->bbc );
        }
        
        // vC
        AddCoeff2( JtempA[ith], AtempA[ith], eqn, eqn,                     &(nnzc2A[ith]), vC*celvol, mesh->BCv.type[c3],      mesh->BCv.val[c3],      StokesA->bbc );
        
        // vE
        if ( mesh->BCv.type[c3+1] != 30 ) {
            if ( mesh->BCv.type[c3+1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    mesh->BCv.val[c3+1],    StokesA->bbc );
            else AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+1],     &(nnzc2A[ith]), vE*celvol, mesh->BCv.type[c3+1],    2*mesh->BCv.val[c3+1],    StokesA->bbc );
        }
        
        // vNW (Newton)
        if ( mesh->BCv.type[c3+nxvz-1] != 30  ) {
            if( mesh->BCv.type[c3+nxvz-1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+nxvz-1],     &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[c3+nxvz-1],      mesh->BCv.val[c3+nxvz-1],    StokesA->bbc );
            else                                       AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+nxvz-1],     &(nnzc2A[ith]), vNW*celvol, mesh->BCv.type[c3+nxvz-1],    2*mesh->BCv.val[c3+nxvz-1],    StokesA->bbc );
        }
        
        // vN
        if ( mesh->BCv.type[c3+nxvz] != 30  ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+nxvz],  &(nnzc2A[ith]), vN*celvol, mesh->BCv.type[c3+nxvz], mesh->BCv.val[c3+nxvz], StokesA->bbc );
        
        // vNE (Newton)
        if ( mesh->BCv.type[c3+nxvz+1] != 30  ) {
            if ( mesh->BCv.type[c3+nxvz+1] != 11 ) AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+nxvz+1],     &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[c3+nxvz+1],      mesh->BCv.val[c3+nxvz+1],    StokesA->bbc );
            else                                        AddCoeff2( JtempA[ith], AtempA[ith], eqn, Stokes->eqn_v[c3+nxvz+1],     &(nnzc2A[ith]), vNE*celvol, mesh->BCv.type[c3+nxvz+1],    2*mesh->BCv.val[c3+nxvz+1],    StokesA->bbc );
        }
        //--------------------
        
        // pS && pN
        if ( mesh->BCp.type[c2] != 30 && mesh->BCp.type[c2+ncx] != 30) {
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2] - Stokes->neq_mom,      &(nnzc2B[ith]), pS*celvol, mesh->BCp.type[c2],     mesh->BCp.val[c2],     StokesB->bbc );
            AddCoeff2( JtempB[ith], AtempB[ith], eqn, Stokes->eqn_p[c2+ncx] - Stokes->neq_mom,  &(nnzc2B[ith]), pN*celvol, mesh->BCp.type[c2+ncx], mesh->BCp.val[c2+ncx], StokesB->bbc );
        }
    }
    else {
        
        // Residual
        StokesA->F[eqn] = vC*v[c3];
        if ( mesh->BCp.type[c2] != 30 && mesh->BCp.type[c2+ncx] != 30 ) {
            StokesA->F[eqn]  += pS*p[c2] + pN*p[c2+ncx];
        }
        if ( mesh->BCv.type[c3-nxvz] != 30 ) StokesA->F[eqn] += vS*v[c3-nxvz];
        if ( mesh->BCv.type[c3-1]    != 30 && mesh->BCv.type[c3-1] != 11 ) StokesA->F[eqn] += vW*v[c3-1];
        if ( mesh->BCv.type[c3+1]    != 30 && mesh->BCv.type[c3+1] != 11 ) StokesA->F[eqn] += vE*v[c3+1];
        if ( mesh->BCv.type[c3+nxvz] != 30 ) StokesA->F[eqn] += vN*v[c3+nxvz];
        if ( mesh->BCu.type[c1-1] != 30 && mesh->BCu.type[c1] != 30 ) {
            StokesA->F[eqn] += uSW*u[c1-1] + uSE*u[c1];
        }
        if ( mesh->BCu.type[c1+nx-1] != 30 && mesh->BCu.type[c1+nx] != 30 ) {
            StokesA->F[eqn] += uNW*u[c1+nx-1] + uNE*u[c1+nx];
        }
        if ( mesh->BCv.type[c3-1] == 11 )   StokesA->F[eqn] += -2*etaW*one_dx_dx*mesh->BCv.val[c3-1];
        if ( mesh->BCv.type[c3+1] == 11 )   StokesA->F[eqn] += -2*etaE*one_dx_dx*mesh->BCv.val[c3+1];
        StokesA->F[eqn] -= (StokesA->b[eqn]);
        StokesA->F[eqn] *= celvol;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void BuildJacobianOperatorDecoupled( grid *mesh, params model, int lev, double *p, double *u, double *v, SparseMat *Stokes, SparseMat *StokesA, SparseMat *StokesB, SparseMat *StokesC, SparseMat *StokesD, int Assemble ) {
    
    // FLAG
    // Assemble=1: Build the linear system of discrete equations
    // Assemble=0: Evaluate Stokes residuals
    
    int    cc, k, l, c1, c2, c3, nx=mesh->Nx, nz=mesh->Nz, nxvz=nx+1, nzvx=nz+1, ncx=nx-1, ncz=nz-1;
    
    // Pre-calculate FD coefs
    double celvol    = mesh->dx*mesh->dz;
    double one_dx    = 1.0/mesh->dx;
    double one_dz    = 1.0/mesh->dz;
    double one_dx_dx = 1.0/mesh->dx/mesh->dx;
    double one_dz_dz = 1.0/mesh->dz/mesh->dz;
    double one_dx_dz = 1.0/mesh->dx/mesh->dz;
    
    // Switches
    int eqn,  sign = 1, comp = model.compressible, stab = 0;
    if (model.free_surf_stab>0) stab = 1;
    double theta = model.free_surf_stab*0.5;
    
    if ( stab==1 ) printf("It tastes like Bo'\n");
    
    // Decompose domain
    int n_th, N, ith;
    int *DD, *estart, *eend, *last_eqn;
    
    // Temporary parallel matrix
    double **AtempA;
    int **JtempA, **ItempA;
    double **AtempB;
    int **JtempB, **ItempB;
    double **AtempC;
    int **JtempC, **ItempC;
    double **AtempD;
    int **JtempD, **ItempD;
    
    int nnzA, nnzB, nnzC, nnzD, nnzcA=0, nnzcB=0, nnzcC=0, nnzcD=0;
    int *nnzc2A, *nnzc2B, *nnzc2C, *nnzc2D;
    
#pragma omp parallel shared(n_th)
    {
        n_th = omp_get_num_threads();
    }
#pragma omp barrier
    
    // Matrix initialisation
    if ( Assemble == 1 ) {
        
        nnzA  = 18*((mesh->Nx-1) * mesh->Nz + (mesh->Nz-1) * mesh->Nx);
        nnzB  = 5*((mesh->Nx-1) * (mesh->Nz-1));
        nnzC  = nnzB;
        nnzD  = 1*((mesh->Nx-1) * (mesh->Nz-1));
        StokesA->neq = Stokes->neq_mom;
        StokesB->neq = Stokes->neq_mom;
        StokesC->neq = Stokes->neq_cont; StokesC->neq_mom = Stokes->neq_mom;
        StokesD->neq = Stokes->neq_cont; StokesD->neq_mom = Stokes->neq_mom;
        printf("Assembling  decoupled Jacobian matrix...\n");
        AllocMat( StokesA, nnzA );
        AllocMat( StokesB, nnzB );
        AllocMat( StokesC, nnzC );
        AllocMat( StokesD, nnzD );
        
    }
    
    // Build velocity block RHS
    int inc=0;
    for( l=0; l<nzvx; l++) {
        for( k=0; k<nx; k++) {
            cc = k + l*nx;
            if ( mesh->BCu.type[cc] != 0 && mesh->BCu.type[cc] != 30 && mesh->BCu.type[cc] != 11 && mesh->BCu.type[cc] != 13 ) {
                StokesA->b[inc] = mesh->roger_x[cc];
                inc++;
            }
        }
    }
    for( l=0; l<nz; l++) {
        for( k=0; k<nxvz; k++) {
            cc = k + l*nxvz;
            if ( mesh->BCv.type[cc] != 0 && mesh->BCv.type[cc] != 30 && mesh->BCv.type[cc] != 11 && mesh->BCv.type[cc] != 13 ) {
                StokesA->b[inc] = mesh->roger_z[cc];
                inc++;
            }
        }
    }
    
    // Build pressure block RHS
    inc=0;
    for( l=0; l<ncz; l++) {
        for( k=0; k<ncx; k++) {
            cc = k + l*ncx;
            if ( mesh->BCp.type[cc] != 31 && mesh->BCp.type[cc] != 0 && mesh->BCp.type[cc] != 30 ) {
                StokesC->b[inc] = mesh->rhs_p[cc] ;
                inc++;
            }
        }
    }
    
    //------------------------------------------------------------------//
    //------------------------- U-momentum -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = nx*nzvx;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempA, &ItempA, &JtempA, n_th, nnzA, Stokes->neq_mom, DD, &nnzc2A  );
        AllocateTempMatArraysDecoupled( &AtempB, &ItempB, &JtempB, n_th, nnzB, Stokes->neq_mom, DD, &nnzc2B  );
    }
    
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesA, StokesB, u, v, p, nx, ncx, nzvx, nnzc2A, AtempA, JtempA, ItempA, nnzc2B, AtempB, JtempB, ItempB, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, comp, celvol )
    {
        ith = omp_get_thread_num();
        
        for( c1=estart[ith]; c1<eend[ith]+1; c1++) {
            
            k      = mesh->kvx[c1];
            l      = mesh->lvx[c1];
            c2     = k-1 + (l-1)*ncx;
            c3     = k   + l*nxvz;
            
            // Get equation numbering
            if ( mesh->BCu.type[c1]==-1 || mesh->BCu.type[c1]==2 ) {
                eqn = Stokes->eqn_u[c1];
                last_eqn[ith]   = eqn;
                
                if ( Assemble == 1 ) {
                    ItempA[ith][eqn] = nnzc2A[ith];
                    ItempB[ith][eqn] = nnzc2B[ith];
                }
                
                //--------------------- INNER NODES ---------------------//
                if ( l>0 && l<nzvx-1 && k>0 && k<nx-1 ) {
                    Xjacobian_InnerNodesDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                if ( k==0    && mesh->BCu.type[c1]==2 ) {
                    //                    printf("%2.2e\n", mesh->eta_n[c2+1]);
                    Xmomentum_WestNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                if ( k==nx-1 && mesh->BCu.type[c1]==2 ) {
                    Xmomentum_EastNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
            }
        }
    }
    
    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrix( StokesA, AtempA, JtempA, ItempA, mesh, estart, eend, &nnzcA, nnzc2A, last_eqn, n_th, mesh->BCu.type, Stokes->eqn_u );
        MergeParallelMatrix( StokesB, AtempB, JtempB, ItempB, mesh, estart, eend, &nnzcB, nnzc2B, last_eqn, n_th, mesh->BCu.type, Stokes->eqn_u );
    }
    
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempA, ItempA, JtempA, n_th, nnzc2A );
        FreeTempMatArraysDecoupled( AtempB, ItempB, JtempB, n_th, nnzc2B );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //------------------------- V-momentum -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = nxvz*nz;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempA, &ItempA, &JtempA, n_th, nnzA, Stokes->neq_mom, DD, &nnzc2A  );
        AllocateTempMatArraysDecoupled( &AtempB, &ItempB, &JtempB, n_th, nnzB, Stokes->neq_mom, DD, &nnzc2B  );
    }
    
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesA, StokesB, u, v, p, nx, ncx, nzvx, nnzc2A, AtempA, JtempA, ItempA, nnzc2B, AtempB, JtempB, ItempB, last_eqn )  private( ith, l, k, c1, c2, c3, eqn ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, comp, celvol )
    {
        ith = omp_get_thread_num();
        
        for( c3=estart[ith]; c3<eend[ith]+1; c3++) {
            
            k  = mesh->kvz[c3];
            l  = mesh->lvz[c3];
            c1 = k   + l*nx;
            c2 = k-1 + (l-1)*ncx;
            
            // Get equation numbering
            if ( mesh->BCv.type[c3] == -1 || mesh->BCv.type[c3]==2 ) {
                eqn = Stokes->eqn_v[c3];
                last_eqn[ith]   = eqn;
                
                if ( Assemble == 1 ) {
                    ItempA[ith][eqn] = nnzc2A[ith];
                    ItempB[ith][eqn] = nnzc2B[ith];
                }
                
                //--------------------- INNER NODES ---------------------//
                if ( k>0 && k<nxvz-1 && l>0 && l<nz-1 ) {
                    
                    Zjacobian_InnerNodesDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                //--------------------- INNER NODES ---------------------//
                
                if ( l==0 && mesh->BCv.type[c3] == 2 ) {
                    Zmomentum_SouthNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
                //--------------------- INNER NODES ---------------------//
                
                if ( l==nz-1 && mesh->BCv.type[c3] == 2 ) {
                    Zmomentum_NorthNeumannDecoupled( Stokes, StokesA, StokesB, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempA, AtempA, nnzc2A, JtempB, AtempB, nnzc2B );
                }
                
            }
        }
    }
    
    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrix( StokesA, AtempA, JtempA, ItempA, mesh, estart, eend, &nnzcA, nnzc2A, last_eqn, n_th, mesh->BCv.type, Stokes->eqn_v );
        MergeParallelMatrix( StokesB, AtempB, JtempB, ItempB, mesh, estart, eend, &nnzcB, nnzc2B, last_eqn, n_th, mesh->BCv.type, Stokes->eqn_v );
    }
    
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempA, ItempA, JtempA, n_th, nnzc2A );
        FreeTempMatArraysDecoupled( AtempB, ItempB, JtempB, n_th, nnzc2B );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //------------------------- Continuity -----------------------------//
    //------------------------------------------------------------------//
    
    // Parallel decomposition
    N = ncx*ncz;
    AllocateDecompositionDecoupled(  n_th, &estart, &eend, &DD, &last_eqn );
    DomainDecompositionDecoupled( N, n_th, estart, eend, DD );
    
    // Allocate parallel matrix
    if ( Assemble == 1 ) {
        AllocateTempMatArraysDecoupled( &AtempC, &ItempC, &JtempC, n_th, nnzC, Stokes->neq_cont, DD, &nnzc2C  );
        AllocateTempMatArraysDecoupled( &AtempD, &ItempD, &JtempD, n_th, nnzD, Stokes->neq_cont, DD, &nnzc2D  );
    }
    
    
#pragma omp parallel shared( eend, estart, mesh, Stokes, StokesC, StokesD, u, v, p, nnzc2C, AtempC, JtempC, ItempC, nnzc2D, AtempD, JtempD, ItempD, last_eqn )  private( ith, l, k, c1, c2, c3, eqn, comp ) firstprivate( model, Assemble, lev, one_dx_dx, one_dz_dz, one_dx_dz, one_dx, one_dz, sign, theta, stab, celvol, nx, ncx, nxvz, nzvx )
    {
        
        ith = omp_get_thread_num();
        
        for( c2=estart[ith]; c2<eend[ith]+1; c2++) {
            
            k   = mesh->kp[c2];
            l   = mesh->lp[c2];
            c1  = k   + (l+1)*nx;
            c3  = k   + l*nxvz + 1;
            
            //--------------------- INNER NODES ---------------------//
            if ( mesh->BCp.type[c2] == -1) {
                
                comp = mesh->comp_cells[c1];
                
                eqn = Stokes->eqn_p[c2]  - Stokes->neq_mom;
                last_eqn[ith]   = eqn ;
                
                if ( Assemble == 1 ) {
                    ItempC[ith][eqn] = nnzc2C[ith];
                    ItempD[ith][eqn] = nnzc2D[ith]; //printf("%d ",  nnzc2D[ith]);
                }
                //
                Continuity_InnerNodesDecoupled( Stokes, StokesC, StokesD, Assemble, lev, stab, comp, theta, sign, model, one_dx, one_dz, one_dx_dx, one_dz_dz, one_dx_dz, celvol, mesh, ith, c1, c2, c3, nx, ncx, nxvz, eqn, u, v, p, JtempC, AtempC, nnzc2C, JtempD, AtempD, nnzc2D, k, l );
            }
        }
    }
    
    //-----------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        MergeParallelMatrixDecoupled( StokesC, AtempC, JtempC, ItempC, mesh, estart, eend, &nnzcC, nnzc2C, last_eqn, n_th, mesh->BCp.type, Stokes->eqn_p );
        MergeParallelMatrixDecoupled( StokesD, AtempD, JtempD, ItempD, mesh, estart, eend, &nnzcD, nnzc2D, last_eqn, n_th, mesh->BCp.type, Stokes->eqn_p );
    }
    if ( Assemble == 1 ) {
        FreeTempMatArraysDecoupled( AtempC, ItempC, JtempC, n_th, nnzc2C );
        FreeTempMatArraysDecoupled( AtempD, ItempD, JtempD, n_th, nnzc2D );
    }
    FreeDecompositionDecoupled( estart, eend, DD, last_eqn );
    
    //------------------------------------------------------------------//
    //----------------------------- End --------------------------------//
    //------------------------------------------------------------------//
    
    if ( Assemble == 1 ) {
        
        // Add contribution from the BC's
        for (k=0; k<StokesA->neq; k++) {
            StokesA->b[k] += StokesA->bbc[k];
            StokesB->b[k] = StokesA->b[k];
        }

        // Add contribution from the BC's
        for (k=0; k<StokesC->neq; k++) {
            StokesC->b[k] += StokesC->bbc[k];
            StokesD->b[k] = StokesC->b[k];
        }
        
        // Final index
        StokesA->Ic[StokesA->neq] = nnzcA;
        StokesB->Ic[StokesB->neq] = nnzcB;
        StokesC->Ic[StokesC->neq] = nnzcC;
        StokesD->Ic[StokesD->neq] = nnzcD;
        
        StokesA->nnz = nnzcA;
        StokesB->nnz = nnzcB;
        StokesC->nnz = nnzcC;
        StokesD->nnz = nnzcD;
        
        // Resize arrays to proper number of non-zeros
        double *bufd;
        int *bufi;
        bufi      = DoodzRealloc(StokesA->J, nnzcA*sizeof(int));
        bufd      = DoodzRealloc(StokesA->A, nnzcA*sizeof(double));
        StokesA->J = bufi;
        StokesA->A = bufd;
        
        bufi      = DoodzRealloc(StokesB->J, nnzcB*sizeof(int));
        bufd      = DoodzRealloc(StokesB->A, nnzcB*sizeof(double));
        StokesB->J = bufi;
        StokesB->A = bufd;
        
        bufi      = DoodzRealloc(StokesC->J, nnzcC*sizeof(int));
        bufd      = DoodzRealloc(StokesC->A, nnzcC*sizeof(double));
        StokesC->J = bufi;
        StokesC->A = bufd;
        
        bufi      = DoodzRealloc(StokesD->J, nnzcD*sizeof(int));
        bufd      = DoodzRealloc(StokesD->A, nnzcD*sizeof(double));
        StokesD->J = bufi;
        StokesD->A = bufd;
    
        printf("System size: ndof = %d, nzA = %d nzB = %d nzC = %d nzD = %d\n", Stokes->neq, nnzcA, nnzcB, nnzcC, nnzcD);
        
//        // Extract Diagonal of A - Viscous block
//        int i, j, locNNZ;
//        int I1, J1;
//
//        if (model.diag_scaling == 1) {
//
//#pragma omp parallel for shared(StokesA) private(I1,J1, i, j, locNNZ )
//            for (i=0;i<StokesA->neq; i++) {
//                I1     = StokesA->Ic[i];
//                locNNZ = StokesA->Ic[i+1] - StokesA->Ic[i];
//                for (J1=0;J1<locNNZ; J1++) {
//                    j = StokesA->J[I1 + J1];
//                    if (i==j) StokesA->d[i] = StokesA->A[I1 + J1];
//                }
//            }
//            // Extract Diagonal of D - Pressure block
//#pragma omp parallel for shared(StokesD) private(i)
//            for (i=0;i<StokesD->neq; i++) {
//                StokesC->d[i] = 1.0;
//            }
//
//            // Scale A
//#pragma omp parallel for shared(StokesA) private(I1,J1, i, j, locNNZ )
//            for (i=0;i<StokesA->neq; i++) {
//
//                StokesA->b[i] *= StokesA->d[i];  // scale RHS
//
//                I1     = StokesA->Ic[i];
//                locNNZ = StokesA->Ic[i+1] - StokesA->Ic[i];
//                for (J1=0;J1<locNNZ; J1++) {
//                    j = StokesA->J[I1 + J1];
//                    StokesA->A[I1 + J1] *= StokesA->d[i]*StokesA->d[j];
//                }
//            }
//
//            // Scale B
//#pragma omp parallel for shared(StokesB,StokesA) private(I1,J1, i, j, locNNZ )
//            for (i=0;i<StokesB->neq; i++) {
//                I1     = StokesB->Ic[i];
//                locNNZ = StokesB->Ic[i+1] - StokesB->Ic[i];
//                for (J1=0;J1<locNNZ; J1++) {
//                    j = StokesB->J[I1 + J1];
//                    StokesB->A[I1 + J1] *= StokesA->d[i]*StokesC->d[j];
//                }
//            }
//
//            // Scale C
//#pragma omp parallel for shared(StokesC,StokesA) private(I1,J1, i, j, locNNZ )
//            for (i=0;i<StokesC->neq; i++) {
//
//                StokesC->b[i] *= StokesC->d[i]; // scale RHS
//
//                I1     = StokesC->Ic[i];
//                locNNZ = StokesC->Ic[i+1] - StokesC->Ic[i];
//                for (J1=0;J1<locNNZ; J1++) {
//                    j = StokesC->J[I1 + J1];
//                    StokesC->A[I1 + J1] *= StokesC->d[i]*StokesA->d[j];
//                }
//            }
//
//        }
//
//        MinMaxArray(StokesA->d, 1, StokesA->neq, "diag. A" );
        
                //--------------------------------------//
#ifdef _HDF5_
        if ( model.write_debug == 1 ) {
            
            char *filename;
            asprintf( &filename, "Jacobian.gzip_%dcpu.h5", n_th );
            
            // Fill in DD data structure
            data4DD OutputDDA, OutputDDB, OutputDDC, OutputDDD;
            
            OutputDDA.V = StokesA->A;
            OutputDDA.Ic = StokesA->Ic;
            OutputDDA.J = StokesA->J;
            OutputDDA.b = StokesA->b;
            OutputDDB.V = StokesB->A;
            OutputDDB.Ic = StokesB->Ic;
            OutputDDB.J = StokesB->J;
            OutputDDB.b = StokesB->b;
            OutputDDC.V = StokesC->A;
            OutputDDC.Ic = StokesC->Ic;
            OutputDDC.J = StokesC->J;
            OutputDDC.b = StokesC->b;
            OutputDDD.V = StokesD->A;
            OutputDDD.Ic = StokesD->Ic;
            OutputDDD.J = StokesD->J;
            OutputDDD.b = StokesD->b;
            OutputDDA.eta_cell = mesh->eta_n;
            OutputDDA.params[0] = nx;
            OutputDDA.params[1] = nz;
            OutputDDA.params[2] = mesh->dx;
            OutputDDA.params[3] = mesh->dz;
            OutputDDA.eqn_u = DoodzMalloc(nx*nzvx*sizeof(int));
            OutputDDA.eqn_v = DoodzMalloc(nxvz*nz*sizeof(int));
            OutputDDA.eqn_p = DoodzMalloc(ncx*ncz*sizeof(int));
            
            for( l=0; l<nzvx; l++) {
                for( k=0; k<nx; k++) {
                    cc = k + l*nx;
                    c1 = k + l*nx;
                    //                OutputDD.eqn_u[cc]=c1;
                    OutputDDA.eqn_u[cc]=Stokes->eqn_u[cc];
                }
            }
            for( l=0; l<nz; l++) {
                for( k=0; k<nxvz; k++) {
                    cc = k + l*nxvz;
                    c3 = nx*nzvx + k + l*nxvz;
                    //                OutputDD.eqn_v[cc]=c3;
                    OutputDDA.eqn_v[cc]=Stokes->eqn_v[cc];
                }
            }
            for( l=0; l<ncz; l++) {
                for( k=0; k<ncx; k++) {
                    cc = k + l*ncx;
                    c2 = nx*nzvx + nxvz*nz + k + l*ncx;
                    //                OutputDD.eqn_p[cc]=c2;
                    OutputDDA.eqn_p[cc]=Stokes->eqn_p[cc];
                }
            }
            
            // Send data to file
            create_output_hdf5( filename );
            AddGroup_to_hdf5( filename, "model" );
            AddGroup_to_hdf5( filename, "matrix" );
            AddGroup_to_hdf5( filename, "numbering" );
            AddGroup_to_hdf5( filename, "fields" );
            AddFieldToGroup_generic( _TRUE_, filename, "model", "params" , 'd', 4, OutputDDA.params,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "numbering", "eqn_u" , 'i', nx*nzvx, OutputDDA.eqn_u,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "numbering", "eqn_v" , 'i', nxvz*nz, OutputDDA.eqn_v,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "numbering", "eqn_p" , 'i', ncx*ncz, OutputDDA.eqn_p,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "eta_n" , 'd', ncx*ncz, mesh->eta_n,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "eta_s" , 'd', nx*nz, mesh->eta_s,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_s" , 'c', nx*nz, mesh->BCg.type,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_n" , 'c', ncx*ncz, mesh->BCp.type,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_u" , 'c', nx*nzvx, mesh->BCu.type,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "fields", "tag_v" , 'c', nxvz*nz, mesh->BCv.type,  1 );
            
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "IA" , 'i', StokesA->neq+1, OutputDDA.Ic,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "JA" , 'i', nnzcA,  OutputDDA.J,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "VA" , 'd', nnzcA,  OutputDDA.V,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "IB" , 'i', StokesB->neq+1, OutputDDB.Ic,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "JB" , 'i', nnzcB,  OutputDDB.J,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "VB" , 'd', nnzcB,  OutputDDB.V,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "IC" , 'i', StokesC->neq+1, OutputDDC.Ic,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "JC" , 'i', nnzcC,  OutputDDC.J,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "VC" , 'd', nnzcC,  OutputDDC.V,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "ID" , 'i', StokesD->neq+1, OutputDDD.Ic,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "JD" , 'i', nnzcD,  OutputDDD.J,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "VD" , 'd', nnzcD,  OutputDDD.V,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "eta_cell", 'd', ncx*ncz, OutputDDA.eta_cell,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "rhs_mom" , 'd', StokesA->neq, OutputDDA.b,  1 );
            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "rhs_cont", 'd', StokesC->neq, OutputDDC.b,  1 );
            //            AddFieldToGroup_generic( _TRUE_, filename, "matrix", "rhs", 'd', Stokes->neq, OutputDD.b,  1 );
            
            DoodzFree(OutputDDA.eqn_u);
            DoodzFree(OutputDDA.eqn_v);
            DoodzFree(OutputDDA.eqn_p);
            free(filename);
        }
#endif
    }
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/



