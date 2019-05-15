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
#include "ctype.h"
#include "header_MDOODZ.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#endif

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DeletePreviousBreakpoint( int step, int writer_step ) {
    char *name, *command;
    int success;
    asprintf(&name, "Breakpoint%05d.dat", step- 2*writer_step);
//    asprintf(&command, "rm -rf %s", name );
//    success = system( command );
        success = remove( name );
    if ( success!=-1 ) {
        printf("File %s was successfully deleted\n", name);
    }
    free(name);
//    free(command);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void LoadBreakpointParticles( markers *particles, grid* mesh, markers *topo_chain, markers *topo_chain_ini, params *model, surface* topo, surface* topo_ini, scale scaling ) {
    
    char *name;
    int s1=0, s2=0, s3=0, s4=0;
    int k, l, Nx=model->Nx, Nz=model->Nz, c;
    int Ncx = Nx-1, Ncz = Nz-1;
    
    //---------------------------------------------------------------------------------------------------------//
    
    FILE *file;
    
    // Filename
//    asprintf(&name, "Breakpoint%05d.dat", model->step);
    if( model->rm_break == -1 ) asprintf(&name, "BreakpointXXXXX.dat");
    else asprintf(&name, "Breakpoint%05d.dat", model->step);
    if (fopen(name, "rb")==NULL) {
        printf("Error opening %s, this file is probably corrupted or non-existant\nExiting...\n", name);
        exit(1);
    }
    else {
        printf("Load setup from %s...\n", name);
        file = fopen(name, "rb");
    }
    fread( &s1, 1, 1, file);
    fread( &s2, 1, 1, file);
    fread( &s3, 1, 1, file);
    fread( &s4, 1, 1, file);
    fread( &model->dt,   s2, 1, file);
    fread( &model->dt0,  s2, 1, file);
    fread( &model->time, s2, 1, file);
    fread( &model->xmin, s2, 1, file);
    fread( &model->xmax, s2, 1, file);
    fread( &model->zmin, s2, 1, file);
    fread( &model->zmax, s2, 1, file);
    
    
    if (model->free_surf == 1) {
        fread(&topo_chain->Nb_part,   s1,                   1, file );
        fread( topo_chain->x,         s3, topo_chain->Nb_part, file );
        fread( topo_chain->z,         s3, topo_chain->Nb_part, file );
        fread( topo_chain->Vx,        s3, topo_chain->Nb_part, file );
        fread( topo_chain->Vz,        s3, topo_chain->Nb_part, file );
        
        fread(&topo_chain_ini->Nb_part,   s1,                       1, file );
        fread( topo_chain_ini->x,         s3, topo_chain_ini->Nb_part, file );
        fread( topo_chain_ini->z,         s3, topo_chain_ini->Nb_part, file );
        fread( topo_chain_ini->Vx,        s3, topo_chain_ini->Nb_part, file );
        fread( topo_chain_ini->Vz,        s3, topo_chain_ini->Nb_part, file );
    }
    
    fread( &particles->Nb_part,  s1, 1, file);
    fread( particles->x,    s3, particles->Nb_part, file);
    fread( particles->z,    s3, particles->Nb_part, file);
    fread( particles->P,    s3, particles->Nb_part, file);
    fread( particles->Vx,   s3, particles->Nb_part, file);
    fread( particles->Vz,   s3, particles->Nb_part, file);
    fread( particles->phi,  s3, particles->Nb_part, file);
    fread( particles->X  ,  s3, particles->Nb_part, file);
    fread( particles->phase, s1, particles->Nb_part, file);
    
    //        if  (model->moving_front == 1 ) {
    //            fread( particles->generation, s1, particles->Nb_part, file);
    //            fread( particles->progress, s1, particles->Nb_part, file);
    //        }
    
    if (model->iselastic == 1) {
        fread( particles->sxxd,   s3, particles->Nb_part, file );
        fread( particles->sxz,    s3, particles->Nb_part, file );
        
        fread( mesh->eta_n, s3, (Nx-1)*(Nz-1), file );
        fread( mesh->VE_n, s3, (Nx-1)*(Nz-1), file );
        fread( mesh->eta_s, s3, (Nx)*(Nz), file );
        fread( mesh->VE_s, s3, (Nx)*(Nz), file );
    }
    
    fread( particles->strain,     s3, particles->Nb_part, file);
    fread( particles->strain_el,  s3, particles->Nb_part, file);
    fread( particles->strain_pl,  s3, particles->Nb_part, file);
    fread( particles->strain_pwl, s3, particles->Nb_part, file);
    fread( particles->strain_exp, s3, particles->Nb_part, file);
    fread( particles->strain_lin, s3, particles->Nb_part, file);
    fread( particles->strain_gbs, s3, particles->Nb_part, file);
    fread( particles->d         , s3, particles->Nb_part, file);
    
    if (model->fstrain == 1) {
        fread( particles->Fxx         , s3, particles->Nb_part, file);
        fread( particles->Fxz         , s3, particles->Nb_part, file);
        fread( particles->Fzx         , s3, particles->Nb_part, file);
        fread( particles->Fzz         , s3, particles->Nb_part, file);
    }
    
    if (model->rec_T_P_x_z == 1) {
        fread( particles->T0         , s3, particles->Nb_part, file);
        fread( particles->P0         , s3, particles->Nb_part, file);
        fread( particles->x0         , s3, particles->Nb_part, file);
        fread( particles->z0         , s3, particles->Nb_part, file);
        fread( particles->Tmax       , s3, particles->Nb_part, file);
        fread( particles->Pmax       , s3, particles->Nb_part, file);
    }
    
    if (model->isthermal == 1) {
        fread( particles->T, s3, particles->Nb_part, file);
    }
    
    if (model->eqn_state > 0) {
        fread( particles->rho, s3, particles->Nb_part, file);
    }
    
    
    fread( mesh->eta_phys_n, s3, (Nx-1)*(Nz-1), file );
    fread( mesh->eta_phys_s, s3, (Nx)*(Nz), file );
    
    //        if (model->isinertial == 2) {
    //            fread( mesh->u_adv, s3, (Nx*(Nz+1)), file );
    //            fread( mesh->v_adv, s3, (Nz*(Nx+1)), file );
    //        }
    
    fread( mesh->u_in, s3, (Nx*(Nz+1)), file );
    fread( mesh->v_in, s3, (Nz*(Nx+1)), file );
    fread( mesh->p_in, s3, ((Nz-1)*(Nx-1)), file );
    
    fread( &mesh->Ut,  s3, 1, file );
    fread( &mesh->W,   s3, 1, file );
    fread( &model->L0, s3, 1, file );
    
    // This is to avoid any problem with restarting - more data needs to be stored

        // Topo related
        if (model->free_surf == 1) {
            fread( topo->height,  s3,  Nx, file );
            fread( topo->height0, s3,  Nx, file );
            fread( topo->a,       s3, Ncx, file );
            fread( topo->a0,      s3, Ncx, file );
            fread( topo->b,       s3, Ncx, file );
            fread( topo->b0,      s3, Ncx, file );
            fread( topo->vx,      s3,  Nx, file );
            fread( topo->vz,      s3,Nx+1, file );
            fread( topo_ini->height,  s3,  Nx, file );
            fread( topo_ini->height0, s3,  Nx, file );
            fread( topo_ini->a,       s3, Ncx, file );
            fread( topo_ini->a0,      s3, Ncx, file );
            fread( topo_ini->b,       s3, Ncx, file );
            fread( topo_ini->b0,      s3, Ncx, file );
            fread( topo_ini->vx,      s3,  Nx, file );
            fread( topo_ini->vz,      s3,Nx+1, file );
        }
        // Cell flags
        fread( mesh->BCt.type,     s4,  Ncx*Ncz,   file );
        fread( mesh->BCp.type,  s4,  Ncx*Ncz,   file );
        fread( mesh->BCu.type,  s4,  Nx*(Nz+1), file );
        fread( mesh->BCv.type,  s4,  (Nx+1)*Nz, file );
        fread( mesh->BCg.type,     s4,  Nx *Nz ,   file );
        fread( mesh->BCp.val,   s3,  Ncx*Ncz,   file );
        fread( mesh->BCu.val,   s3,  Nx*(Nz+1), file );
        fread( mesh->BCv.val,   s3,  (Nx+1)*Nz, file );
        fread( mesh->BCg.val,      s3,  Nx *Nz ,   file );
        printf("Loading phase proportions:\n");
        // Phase proportions
        fread( &model->Nb_phases,  s1,        1,   file );
        for ( k=0; k< model->Nb_phases; k++ ) {
            fread( mesh->phase_perc_n[k], s3,  Ncx*Ncz,   file );
            fread( mesh->phase_perc_s[k], s3,  Nx *Nz ,   file );
            //            printf("phase %02d:\n", k);
            //            MinMaxArray(mesh->phase_perc_n[k], 1.0, Ncx*Ncz, "phase_perc_n");
            //            MinMaxArray(mesh->phase_perc_s[k], 1.0, Nx *Nz, "phase_perc_s");
        }
        fread( mesh->nb_part_cell, s1,  Ncx*Ncz,   file );
        fread( mesh->nb_part_vert, s1,  Nx *Nz ,   file );
    
    
    fclose(file);
    free(name);
    
    //---------------------------------------------------------------------------------------------------------//
    
    // scale
    model->dt   /= scaling.t;
    model->dt0  /= scaling.t;
    model->time /= scaling.t;
    model->xmin /= scaling.L;
    model->xmax /= scaling.L;
    model->zmin /= scaling.L;
    model->zmax /= scaling.L;
    
#pragma omp parallel for shared( particles, model, scaling )
    for ( k=0; k<particles->Nb_part; k++ ) {
        particles->x[k]     /=scaling.L;
        particles->z[k]     /=scaling.L;
        particles->P[k]     /=scaling.S;
        particles->Vx[k]    /=scaling.V;
        particles->Vz[k]    /=scaling.V;
        particles->phi[k]   /=1.0;
        particles->X[k]     /=1.0;
        if ( model->iselastic == 1 ) {
            particles->sxxd[k]    /= scaling.S;
            particles->sxz[k]     /= scaling.S;
        }
        if ( model->isthermal == 1 ) {
            particles->T[k]  /= scaling.T;
        }
        
        if ( model->eqn_state > 0 ) {
            particles->rho[k]  /=scaling.rho;
        }
        
        if ( model->rec_T_P_x_z == 1) {
            particles->T0[k]  /= scaling.T;
            particles->P0[k]  /= scaling.S;
            particles->x0[k]  /= scaling.L;
            particles->z0[k]  /= scaling.L;
            particles->Tmax[k]/= scaling.T;
            particles->Pmax[k]/= scaling.S;
        }
        particles->d[k]  /=scaling.L;
    }
    
    // Grid data
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz+1; l++) {
            c = k +l*Nx;
            mesh->u_adv[c] /= scaling.V;
            mesh->u_in[c]  /= scaling.V;
        }
    }
    
    for (k=0; k<Nx+1; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx+1);
            mesh->v_adv[c] /= scaling.V;
            mesh->v_in[c]  /= scaling.V;
        }
    }
    
    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            c = k +l*(Nx-1);
            mesh->eta_n[c] /= scaling.eta;
            mesh->VE_n[c] /= 1.0;
            mesh->eta_phys_n[c] /= scaling.eta;
            mesh->p_in[c] /= scaling.S;
        }
    }
    
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx);
            mesh->eta_phys_s[c] /= scaling.eta;
            mesh->eta_s[c] /= scaling.eta;
            mesh->VE_s[c] /= 1.0;
        }
    }
    
    if (model->free_surf == 1) {
        for (k=0; k<topo_chain->Nb_part; k++) {
            topo_chain->x[k]  /= scaling.L;
            topo_chain->z[k]  /= scaling.L;
            topo_chain->Vx[k] /= scaling.V;
            topo_chain->Vz[k] /= scaling.V;
        }
        for (k=0; k<topo_chain_ini->Nb_part; k++) {
            topo_chain_ini->x[k] /= scaling.L;
            topo_chain_ini->z[k] /= scaling.L;
            topo_chain_ini->Vx[k] /= scaling.V;
            topo_chain_ini->Vz[k] /= scaling.V;
        }
        //    }
        
        // This is to avoid any problem with restarting - more data needs to be stored
        for (k=0;k<mesh->Nx;k++) {
            topo->height[k]      /= scaling.L;
            topo->height0[k]     /= scaling.L;
            topo->vx[k]          /= scaling.V;
            topo_ini->height[k]  /= scaling.L;
            topo_ini->height0[k] /= scaling.L;
            topo_ini->vx[k]      /= scaling.V;
        }
        
        for (k=0;k<mesh->Nx+1;k++) {
            topo->vz[k]          /= scaling.V;
            topo_ini->vz[k]      /= scaling.V;
        }
        
        for (k=0;k<mesh->Nx-1;k++) {
            topo->b0[k] /= scaling.L;
            topo->b[k]  /= scaling.L;
            topo_ini->b0[k] /= scaling.L;
            topo_ini->b[k]  /= scaling.L;
        }
    }
    
    mesh->Ut /= (scaling.rhoE*scaling.L*scaling.L);
    mesh->W  /= (scaling.rhoE*scaling.L*scaling.L);
    model->L0  /= (scaling.L);
    
    MinMaxArray(particles->T, scaling.T, particles->Nb_part, "T part");
    
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MakeBreakpointParticles( markers *particles,  grid* mesh, markers *topo_chain, markers *topo_chain_ini, params model, surface *topo, surface* topo_ini, scale scaling ) {
    
    char *name;
    int s1=sizeof(int), s2=sizeof(double), s3=sizeof(DoodzFP), s4=sizeof(char);
    int k, l, Nx=model.Nx, Nz=model.Nz, c;
    int Ncx = Nx-1, Ncz = Nz-1;
    
    // Scale such that dimensional data goes to the breakpoint file
    model.dt   *= scaling.t;
    model.dt0  *= scaling.t;
    model.time *= scaling.t;
    model.xmin *= scaling.L;
    model.xmax *= scaling.L;
    model.zmin *= scaling.L;
    model.zmax *= scaling.L;
    
#pragma omp parallel for shared( particles, model, scaling )
    for ( k=0; k<particles->Nb_part; k++ ) {
        particles->x[k]      *= scaling.L;
        particles->z[k]      *= scaling.L;
        particles->P[k]      *= scaling.S;
        particles->Vx[k]     *= scaling.V;
        particles->Vz[k]     *= scaling.V;
        particles->phi[k]    *= 1.0;
        particles->X[k]      *= 1.0;
        if (model.iselastic == 1) {
            particles->sxxd[k]   *= scaling.S;
            particles->sxz[k]    *= scaling.S;
        }
        if (model.isthermal == 1) {
            particles->T[k]  *= scaling.T;
        }
        if (model.eqn_state > 0) {
            particles->rho[k]  *= scaling.rho;
        }
        particles->d[k]  *= scaling.L;
        
        if ( model.rec_T_P_x_z == 1) {
            particles->T0[k]  *= scaling.T;
            particles->P0[k]  *= scaling.S;
            particles->x0[k]  *= scaling.L;
            particles->z0[k]  *= scaling.L;
            particles->Tmax[k]*= scaling.T;
            particles->Pmax[k]*= scaling.S;
        }
    }
    
    // grid data
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz+1; l++) {
            c = k +l*Nx;
            mesh->u_adv[c] *= scaling.V;
            mesh->u_in[c]  *= scaling.V;
        }
    }
    
    for (k=0; k<Nx+1; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx+1);
            mesh->v_adv[c] *= scaling.V;
            mesh->v_in[c]  *= scaling.V;
        }
    }
    
    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            c = k +l*(Nx-1);
            mesh->eta_n[c] *= scaling.eta;
            mesh->VE_n[c] *= 1.0;
            mesh->eta_phys_n[c] *= scaling.eta;
            mesh->p_in[c] *= scaling.S;
        }
    }
    
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx);
            mesh->eta_phys_s[c] *= scaling.eta;
            mesh->eta_s[c] *= scaling.eta;
            mesh->VE_s[c] *= 1.0;
        }
    }
    
    if (model.free_surf == 1) {
        for (k=0; k<topo_chain->Nb_part; k++) {
            topo_chain->x[k] *= scaling.L;
            topo_chain->z[k] *= scaling.L;
            topo_chain->Vx[k] *= scaling.V;
            topo_chain->Vz[k] *= scaling.V;
        }
        for (k=0; k<topo_chain_ini->Nb_part; k++) {
            topo_chain_ini->x[k] *= scaling.L;
            topo_chain_ini->z[k] *= scaling.L;
            topo_chain_ini->Vx[k] *= scaling.V;
            topo_chain_ini->Vz[k] *= scaling.V;
        }
        
        
        // This is to avoid any problem with restarting - more data needs to be stored
        for (k=0;k<mesh->Nx;k++) {
            topo->height[k]      *= scaling.L;
            topo->height0[k]     *= scaling.L;
            topo->vx[k]          *= scaling.V;
            topo_ini->height[k]  *= scaling.L;
            topo_ini->height0[k] *= scaling.L;
            topo_ini->vx[k]      *= scaling.V;
        }
        
        for (k=0;k<mesh->Nx+1;k++) {
            topo->vz[k]          *= scaling.V;
            topo_ini->vz[k]      *= scaling.V;
        }
        
        for (k=0;k<mesh->Nx-1;k++) {
            topo->b0[k] *= scaling.L;
            topo->b[k]  *= scaling.L;
            topo_ini->b0[k] *= scaling.L;
            topo_ini->b[k]  *= scaling.L;
        }
    }
    
    mesh->Ut *= (scaling.rhoE*scaling.L*scaling.L);
    mesh->W  *= (scaling.rhoE*scaling.L*scaling.L);
    model.L0 *= (scaling.L);
    
    
    //---------------------------------------------------------------------------------------------------------//
    
    FILE *file;
    
    // Filename
    if( model.rm_break == -1 ) asprintf(&name, "BreakpointXXXXX.dat");
    else asprintf(&name, "Breakpoint%05d.dat", model.step);
    
    if (fopen(name, "wb")==NULL) {
        printf("Error opening %s, this file is probably corrupted or non-existant\nExiting...\n", name);
        free(name);
        exit(1);
    }
    else {
        printf("Save setup in %s...\n", name);
        file = fopen(name, "wb");
    }
    fwrite( &s1, 1, 1, file);
    fwrite( &s2, 1, 1, file);
    fwrite( &s3, 1, 1, file);
    fwrite( &s4, 1, 1, file);
    fwrite( &model.dt,  s2, 1, file);
    fwrite( &model.dt0, s2, 1, file);
    fwrite( &model.time, s2, 1, file);
    fwrite( &model.xmin, s2, 1, file);
    fwrite( &model.xmax, s2, 1, file);
    fwrite( &model.zmin, s2, 1, file);
    fwrite( &model.zmax, s2, 1, file);
    
    if (model.free_surf == 1) {
        fwrite( &topo_chain->Nb_part,  s1,                   1, file );
        fwrite( topo_chain->x,         s3, topo_chain->Nb_part, file );
        fwrite( topo_chain->z,         s3, topo_chain->Nb_part, file );
        fwrite( topo_chain->Vx,        s3, topo_chain->Nb_part, file );
        fwrite( topo_chain->Vz,        s3, topo_chain->Nb_part, file );
        
        fwrite( &topo_chain_ini->Nb_part,  s1,                       1, file );
        fwrite( topo_chain_ini->x,         s3, topo_chain_ini->Nb_part, file );
        fwrite( topo_chain_ini->z,         s3, topo_chain_ini->Nb_part, file );
        fwrite( topo_chain_ini->Vx,        s3, topo_chain_ini->Nb_part, file );
        fwrite( topo_chain_ini->Vz,        s3, topo_chain_ini->Nb_part, file );
    }
    
    fwrite( &particles->Nb_part,  s1, 1, file);
    fwrite( particles->x,     s3, particles->Nb_part, file);
    fwrite( particles->z,     s3, particles->Nb_part, file);
    fwrite( particles->P,     s3, particles->Nb_part, file);
    fwrite( particles->Vx,    s3, particles->Nb_part, file);
    fwrite( particles->Vz,    s3, particles->Nb_part, file);
    fwrite( particles->phi,   s3, particles->Nb_part, file);
    fwrite( particles->X  ,   s3, particles->Nb_part, file);
    fwrite( particles->phase, s1, particles->Nb_part, file);
    
    //        if  (model.moving_front == 1 ) {
    //            fwrite( particles->generation, s1, particles->Nb_part, file);
    //            fwrite( particles->progress,   s1, particles->Nb_part, file);
    //        }
    
    if (model.iselastic == 1) {
        fwrite( particles->sxxd,   s3, particles->Nb_part, file );
        fwrite( particles->sxz,    s3, particles->Nb_part, file );
        
        fwrite( mesh->eta_n, s3, (Nx-1)*(Nz-1), file );
        fwrite( mesh->VE_n, s3, (Nx-1)*(Nz-1), file );
        fwrite( mesh->eta_s, s3, (Nx)*(Nz), file );
        fwrite( mesh->VE_s, s3, (Nx)*(Nz), file );
    }
    
    fwrite( particles->strain,     s3, particles->Nb_part, file);
    fwrite( particles->strain_el,  s3, particles->Nb_part, file);
    fwrite( particles->strain_pl,  s3, particles->Nb_part, file);
    fwrite( particles->strain_pwl, s3, particles->Nb_part, file);
    fwrite( particles->strain_exp, s3, particles->Nb_part, file);
    fwrite( particles->strain_lin, s3, particles->Nb_part, file);
    fwrite( particles->strain_gbs, s3, particles->Nb_part, file);
    fwrite( particles->d         , s3, particles->Nb_part, file);
    
    if (model.fstrain == 1) {
        fwrite( particles->Fxx         , s3, particles->Nb_part, file);
        fwrite( particles->Fxz         , s3, particles->Nb_part, file);
        fwrite( particles->Fzx         , s3, particles->Nb_part, file);
        fwrite( particles->Fzz         , s3, particles->Nb_part, file);
    }
    
    if (model.rec_T_P_x_z == 1) {
        fwrite( particles->T0         , s3, particles->Nb_part, file);
        fwrite( particles->P0         , s3, particles->Nb_part, file);
        fwrite( particles->x0         , s3, particles->Nb_part, file);
        fwrite( particles->z0         , s3, particles->Nb_part, file);
        fwrite( particles->Tmax       , s3, particles->Nb_part, file);
        fwrite( particles->Pmax       , s3, particles->Nb_part, file);
    }
    
    if (model.isthermal == 1) {
        fwrite( particles->T, s3, particles->Nb_part, file);
    }
    if (model.eqn_state > 0) {
        fwrite( particles->rho, s3, particles->Nb_part, file);
    }
    
    fwrite( mesh->eta_phys_n, s3, (Nx-1)*(Nz-1), file );
    fwrite( mesh->eta_phys_s, s3, (Nx)*(Nz), file );
    
    //        if (model.isinertial == 2) {
    //            fwrite( mesh->u_adv, s3, (Nx*(Nz+1)), file );
    //            fwrite( mesh->v_adv, s3, (Nz*(Nx+1)), file );
    //        }
    fwrite( mesh->u_in, s3, (Nx*(Nz+1)), file );
    fwrite( mesh->v_in, s3, (Nz*(Nx+1)), file );
    fwrite( mesh->p_in, s3, ((Nx-1)*(Nz-1)), file );
    
    fwrite( &mesh->Ut, s3, 1, file );
    fwrite( &mesh->W,  s3, 1, file );
    fwrite( &model.L0, s3, 1, file );
    
    // This is to avoid any problem with restarting - more data needs to be stored
    if (model.free_surf == 1) {
        // Topo related
        fwrite( topo->height,  s3,  Nx, file );
        fwrite( topo->height0, s3,  Nx, file );
        fwrite( topo->a,       s3, Ncx, file );
        fwrite( topo->a0,      s3, Ncx, file );
        fwrite( topo->b,       s3, Ncx, file );
        fwrite( topo->b0,      s3, Ncx, file );
        fwrite( topo->vx,      s3,  Nx, file );
        fwrite( topo->vz,      s3,Nx+1, file );
        fwrite( topo_ini->height,  s3,  Nx, file );
        fwrite( topo_ini->height0, s3,  Nx, file );
        fwrite( topo_ini->a,       s3, Ncx, file );
        fwrite( topo_ini->a0,      s3, Ncx, file );
        fwrite( topo_ini->b,       s3, Ncx, file );
        fwrite( topo_ini->b0,      s3, Ncx, file );
        fwrite( topo_ini->vx,      s3,  Nx, file );
        fwrite( topo_ini->vz,      s3,Nx+1, file );
    }
    // Cell flags
    fwrite( mesh->BCt.type,     s4,  Ncx*Ncz,   file );
    fwrite( mesh->BCp.type,  s4,  Ncx*Ncz,   file );
    fwrite( mesh->BCu.type,  s4,  Nx*(Nz+1), file );
    fwrite( mesh->BCv.type,  s4,  (Nx+1)*Nz, file );
    fwrite( mesh->BCg.type,     s4,  Nx *Nz ,   file );
    fwrite( mesh->BCp.val,   s3,  Ncx*Ncz,   file );
    fwrite( mesh->BCu.val,   s3,  Nx*(Nz+1), file );
    fwrite( mesh->BCv.val,   s3,  (Nx+1)*Nz, file );
    fwrite( mesh->BCg.val,      s3,  Nx *Nz ,   file );
    // Phase proportions
    printf("Writing phase proportions:\n");
    fwrite( &model.Nb_phases, s1,        1,   file );
    for ( k=0; k< model.Nb_phases; k++ ) {
        fwrite( mesh->phase_perc_n[k], s3,  Ncx*Ncz,   file );
        fwrite( mesh->phase_perc_s[k], s3,  Nx *Nz ,   file );
        //            printf("phase %02d:\n", k);
        //            MinMaxArray(mesh->phase_perc_n[k], 1.0, Ncx*Ncz, "phase_perc_n");
        //            MinMaxArray(mesh->phase_perc_s[k], 1.0, Nx *Nz, "phase_perc_s");
    }
    fwrite( mesh->nb_part_cell, s1,  Ncx*Ncz,   file );
    fwrite( mesh->nb_part_vert, s1,  Nx *Nz ,   file );
    //    }
    
    
    fclose(file);
    free(name);
    
    //---------------------------------------------------------------------------------------------------------//
    
    // scale such that dimensions vanish
    model.dt   /= scaling.t;
    model.dt0  /= scaling.t;
    model.time /= scaling.t;
    model.xmin /= scaling.L;
    model.xmax /= scaling.L;
    model.zmin /= scaling.L;
    model.zmax /= scaling.L;
    
#pragma omp parallel for shared( particles, model, scaling )
    for ( k=0; k<particles->Nb_part; k++ ) {
        particles->x[k]     /= scaling.L;
        particles->z[k]     /= scaling.L;
        particles->P[k]     /= scaling.S;
        particles->Vx[k]    /= scaling.V;
        particles->Vz[k]    /= scaling.V;
        particles->phi[k]   /= 1.0;
        particles->X[k]     /= 1.0;
        
        if (model.iselastic == 1) {
            particles->sxxd[k]   /= scaling.S;
            particles->sxz[k]    /= scaling.S;
        }
        if (model.isthermal == 1) {
            particles->T[k]  /= scaling.T;
        }
        if (model.eqn_state > 0) {
            particles->rho[k]  /= scaling.rho;
        }
        particles->d[k]  /= scaling.L;
        
        if (model.rec_T_P_x_z == 1) {
            particles->T0[k]  /= scaling.T;
            particles->P0[k]  /= scaling.S;
            particles->x0[k]  /= scaling.L;
            particles->z0[k]  /= scaling.L;
            particles->Tmax[k]/= scaling.T;
            particles->Pmax[k]/= scaling.S;
        }
    }
    
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz+1; l++) {
            c = k +l*Nx;
            mesh->u_adv[c] /= scaling.V;
            mesh->u_in[c]  /= scaling.V;
        }
    }
    
    for (k=0; k<Nx+1; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx+1);
            mesh->v_adv[c] /= scaling.V;
            mesh->v_in[c]  /= scaling.V;
        }
    }
    
    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            c = k +l*(Nx-1);
            mesh->eta_n[c] /= scaling.eta;
            mesh->VE_n[c] /= 1.0;
            mesh->eta_phys_n[c] /= scaling.eta;
            mesh->p_in[c] /= scaling.S;
        }
    }
    
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx);
            mesh->eta_phys_s[c] /= scaling.eta;
            mesh->eta_s[c] /= scaling.eta;
            mesh->VE_s[c] /= 1.0;
        }
    }
    
    if (model.free_surf == 1) {
        for (k=0; k<topo_chain->Nb_part; k++) {
            topo_chain->x[k]  /= scaling.L;
            topo_chain->z[k]  /= scaling.L;
            topo_chain->Vx[k] /= scaling.V;
            topo_chain->Vz[k] /= scaling.V;
        }
        for (k=0; k<topo_chain_ini->Nb_part; k++) {
            topo_chain_ini->x[k]  /= scaling.L;
            topo_chain_ini->z[k]  /= scaling.L;
            topo_chain_ini->Vx[k] /= scaling.V;
            topo_chain_ini->Vz[k] /= scaling.V;
        }
        //}
        
        // This is to avoid any problem with restarting - more data needs to be stored
        for (k=0;k<mesh->Nx;k++) {
            topo->height[k]      /= scaling.L;
            topo->height0[k]     /= scaling.L;
            topo->vx[k]          /= scaling.V;
            topo_ini->height[k]  /= scaling.L;
            topo_ini->height0[k] /= scaling.L;
            topo_ini->vx[k]      /= scaling.V;
        }
        
        for (k=0;k<mesh->Nx+1;k++) {
            topo->vz[k]          /= scaling.V;
            topo_ini->vz[k]      /= scaling.V;
        }
        
        for (k=0;k<mesh->Nx-1;k++) {
            topo->b0[k] /= scaling.L;
            topo->b[k]  /= scaling.L;
            topo_ini->b0[k] /= scaling.L;
            topo_ini->b[k]  /= scaling.L;
        }
    }
    
    mesh->Ut  /= (scaling.rhoE*scaling.L*scaling.L);
    mesh->W   /= (scaling.rhoE*scaling.L*scaling.L);
    model.L0  /= (scaling.L);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ReadInputFile( char* fin_name, int *istep, int *irestart, int *writer, int *writer_step, params *model, scale *scaling, mat_prop *materials, markers *particles, Nparams *Nmodel ) {
    
    //------------------------------------------------------------------------------------------------------------------------------//
    // MODEL PARAMETERS
    //------------------------------------------------------------------------------------------------------------------------------//
    
    FILE *fin;
    int k, gsel;
    
    if (fopen(fin_name,"rt") == NULL) {
        printf("Setup file '%s' does not exist\nExiting...\n", fin_name);
        exit(1);
    }
    else {
        fin = fopen(fin_name,"rt");
    }
    
    // Simulation start/restart from Breakpoint
    *istep                 = ReadInt2( fin, "istep", 0 );
    *irestart              = ReadInt2( fin, "irestart", 0 );
    
    if ( *istep == 0 ) {
        *irestart = 0; // Override the restart step number written in the text file (istep)
    }
    
    // Output
    *writer                = ReadInt2( fin, "writer",          0 );
    *writer_step           = ReadInt2( fin, "writer_step",     1 );
    model->write_markers   = ReadInt2( fin, "writer_markers",  0 );
    model->write_debug     = ReadInt2( fin, "writer_debug",    0 );
    model->write_energies  = ReadInt2( fin, "writer_energies", 0 );
    
    // Input
    model->input_file      = ReadChar( fin, "input_file", "blah.bin");
    
    // Read scales for non-dimensionalisation
    scaling->eta           = ReadDou2( fin, "eta", 1.0  );
    scaling->L             = ReadDou2( fin, "L",   1.0  );
    scaling->V             = ReadDou2( fin, "V",   1.0  );
    scaling->T             = ReadDou2( fin, "T",   1.0  );
    ScaleMe( scaling );
    
    // Domain size
    model->Nx              = ReadInt2( fin, "Nx",   10 );
    model->Nz              = ReadInt2( fin, "Nz",   10 );
    model->Nt              = ReadInt2( fin, "Nt",    1 );
    model->xmin            = ReadDou2( fin, "xmin",  1.0  ) / scaling->L; model->xmin0 = model->xmin;
    model->zmin            = ReadDou2( fin, "zmin",  1.0  ) / scaling->L; model->zmin0 = model->zmin;
    model->xmax            = ReadDou2( fin, "xmax",  1.0  ) / scaling->L; model->xmax0 = model->xmax;
    model->zmax            = ReadDou2( fin, "zmax",  1.0  ) / scaling->L; model->zmax0 = model->zmax;
    model->dt              = ReadDou2( fin, "dt",    0.0  ) / scaling->t;
    model->Courant         = ReadDou2( fin, "Courant",       0.5 );
    model->penalty         = ReadDou2( fin, "penalty",      1.0e10 );
    model->abs_tol_div     = ReadDou2( fin, "abs_tol_div", 1.0e-14 );
    model->rel_tol_div     = ReadDou2( fin, "rel_tol_div", 1.0e-5 );
    model->auto_penalty    = ReadDou2( fin, "auto_penalty",    0.0  );
    model->decoupled_solve = ReadInt2( fin, "decoupled_solve",    1 );
    
    // Switches
    model->ismechanical    = ReadInt2( fin, "ismechanical",    1 );
    model->dt_constant     = ReadInt2( fin, "dt_constant",     0 );
    model->RK              = ReadInt2( fin, "RK",              4 );
    model->isperiodic_x    = ReadInt2( fin, "isperiodic_x",    0 );
    model->ispureshear_ale = ReadInt2( fin, "ispureshear_ALE", 0 );
    model->isinertial      = ReadInt2( fin, "isinertial",      0 );
    model->iselastic       = ReadInt2( fin, "iselastic",       0 );
    model->isthermal       = ReadInt2( fin, "isthermal",       0 );
    model->line_search     = ReadInt2( fin, "line_search",     0 );
    model->free_surf       = ReadInt2( fin, "free_surf",       0 );
    model->free_surf_stab  = ReadDou2( fin, "free_surf_stab",  0 );
    model->eqn_state       = ReadInt2( fin, "eqn_state",       0 );
    model->thermal_eq      = ReadInt2( fin, "thermal_eq",      0 );
    model->cooling_time    = ReadDou2( fin, "cooling_time", 1000e6*365.25*3600*24/scaling->t);
    model->subgrid_diff    = ReadInt2( fin, "subgrid_diff",    0 );
    model->shear_heat      = ReadInt2( fin, "shear_heat",      1 );
    model->adiab_heat      = ReadInt2( fin, "adiab_heat",      0 );
    model->isPl_soft       = ReadInt2( fin, "isPl_soft",       0 );
    model->surf_processes  = ReadInt2( fin, "surf_processes",  0 );
    model->surf_remesh     = ReadInt2( fin, "surf_remesh",     1 );
    model->cpc             = ReadInt2( fin, "cpc",             1 );
    model->advection       = ReadInt2( fin, "advection",       1 );
    model->loc_iter        = ReadInt2( fin, "loc_iter",        1 );
    model->therm_pert      = ReadInt2( fin, "therm_pert",      0 );
    model->fstrain         = ReadInt2( fin, "fstrain",         0 );
    model->cut_noise       = ReadInt2( fin, "cut_noise",       0 );
    model->accu            = ReadDou2( fin, "accu",         1.0e13 );
    model->rheo_on_cells   = ReadInt2( fin, "rheo_on_cells",   0 );
    model->DefectCorrectionForm = ReadInt2( fin, "DefectCorrectionForm",      0 );
    model->HsOnly          = ReadInt2( fin, "HsOnly",          0 );
    model->HomoFields      = ReadInt2( fin, "HomoFields",      0 );
    model->rec_T_P_x_z     = ReadInt2( fin, "rec_T_P_x_z",     0 );
    model->rm_break        = ReadInt2( fin, "rm_break",        1 );
    materials->eta_VP      = ReadDou2( fin, "eta_VP",        0.0 ) / scaling->S / scaling->t;
    model->topografix      = ReadInt2( fin, "topografix",      0 );
    model->aniso           = ReadInt2( fin, "aniso",           0 );
    
    // Setup dependant
    model->EpsBG           = ReadDou2( fin, "EpsBG",           0.0 ) / scaling->E;
    model->PrBG            = ReadDou2( fin, "PrBG",            0.0 ) / scaling->S;
    // Surface processes
    model->surf_diff       = ReadDou2( fin, "surf_diff",       0.0 ) / (pow(scaling->L,2.0)/scaling->t);
    model->surf_ised1      = ReadInt2( fin, "surf_ised1",      0.0 );
    model->surf_ised2      = ReadInt2( fin, "surf_ised2",      0.0 );
    model->surf_sedirate   = ReadDou2( fin, "surf_sedirate",   0.0 ) / scaling->V;
    model->surf_baselev    = ReadDou2( fin, "surf_baselev",    0.0 ) / scaling->L;
    // Initial thermal perturbation
    model->therm_pert_x0   = ReadDou2( fin, "therm_pert_x0",   0.0 ) / scaling->L;
    model->therm_pert_z0   = ReadDou2( fin, "therm_pert_z0",   0.0 ) / scaling->L;
    model->therm_pert_rad  = ReadDou2( fin, "therm_pert_rad",  0.0 ) / scaling->L;
    model->therm_pert_dT   = ReadDou2( fin, "therm_pert_dT" ,  0.0 ) / scaling->T;
    // For rheological database reasons...
    model->force_act_vol_ast = ReadInt2( fin, "force_act_vol_ast",   0 );
    model->act_vol_dis_ast   = ReadDou2( fin, "act_vol_dis_ast" ,  0.0 );
    model->act_vol_dif_ast   = ReadDou2( fin, "act_vol_dif_ast" ,  0.0 );
    
    // Model user's delights
    model->user0           = ReadDou2( fin, "user0",           0.0 );
    model->user1           = ReadDou2( fin, "user1",           0.0 );
    model->user2           = ReadDou2( fin, "user2",           0.0 );
    model->user3           = ReadDou2( fin, "user3",           0.0 );
    model->user4           = ReadDou2( fin, "user4",           0.0 );
    model->user5           = ReadDou2( fin, "user5",           0.0 );
    model->user6           = ReadDou2( fin, "user6",           0.0 );
    model->user7           = ReadDou2( fin, "user7",           0.0 );
    model->user8           = ReadDou2( fin, "user8",           0.0 );
    
    // Derived quantities
    model->dx              = (model->xmax - model->xmin) / (model->Nx - 1); printf("dx = %2.6e\n", model->dx );
    model->dz              = (model->zmax - model->zmin) / (model->Nz - 1);
    model->dt0             = model->dt;
    model->dt_start        = model->dt;
    model->eta_avg         = ReadInt2( fin, "eta_avg",       0 ); // 0 : arithmetic mean
    model->p_avg           = 0;                                   // 0 : arithmetic mean
    
    // Gravity
    model->gx              = ReadDou2( fin, "gx",  0.0 ) / scaling->a;
    model->gz              = ReadDou2( fin, "gz",  0.0 ) / scaling->a;
    
    // Material properties
    model->Nb_phases = materials->Nb_phases =  ReadInt2( fin, "Nb_phases", 0 );
    for ( k=0; k<materials->Nb_phases; k++) {
        // Read general parameters
        materials->rho[k]  = ReadMatProps( fin, "rho", k,   2700.0 )  / scaling->rho;
        materials->mu[k]   = ReadMatProps( fin, "mu",  k,   1.0e10 )  / scaling->S;
        materials->Cv[k]   = ReadMatProps( fin, "Cv",  k,   1.0e3  )  / scaling->Cv;
        materials->k[k]    = ReadMatProps( fin, "k",   k,   1.0e-6 )  / scaling->k;
        materials->k_eff[k] = materials->k[k];
        materials->Qr[k]   = ReadMatProps( fin, "Qr",  k,   1.0e-30)  / (scaling->W / pow(scaling->L,3.0));
        materials->C[k]    = ReadMatProps( fin, "C",   k,   1.0e7  )  / scaling->S;
        materials->phi[k]  = ReadMatProps( fin, "phi", k,    30.0  )  * M_PI/ 180.0;
        materials->Slim[k] = ReadMatProps( fin, "Slim",k,   1.0e10 )  / scaling->S;
        materials->alp[k]  = ReadMatProps( fin, "alp", k,      0.0)  / (1.0/scaling->T);
        materials->bet[k]  = ReadMatProps( fin, "bet", k,  1.0e-40 )  / (1.0/scaling->S);
        materials->drho[k] = ReadMatProps( fin, "drho",k,      0.0 )  / (scaling->rho);
        materials->T0[k]   = (zeroC) / (scaling->T); // +20
        materials->P0[k]   = 1e5 / (scaling->S);
        // Read flow law settings
        materials->cstv[k]  = ReadMatProps( fin, "cstv",k,    1.0  );
        materials->pwlv[k]  = ReadMatProps( fin, "pwlv",k,    0.0  );
        materials->linv[k]  = ReadMatProps( fin, "linv",k,    0.0  );
        materials->gbsv[k]  = ReadMatProps( fin, "gbsv",k,    0.0  );
        materials->expv[k]  = ReadMatProps( fin, "expv",k,    0.0  );
        gsel                = ReadMatProps( fin, "gsel",k,    0.0  );
        materials->eta0[k]  = ReadMatProps( fin, "eta0",k, 1.0e20  );//  / scaling->eta;
        materials->npwl[k]  = ReadMatProps( fin, "npwl",k,    1.0  );
        materials->Qpwl[k]  = ReadMatProps( fin, "Qpwl",k,    0.0  );//   / scaling->J;
        materials->pref_pwl[k] = ReadMatProps( fin, "pref_pwl",k,    1.0 );    // weakening prefactor for power law
        materials->gs[k]    = ReadMatProps( fin, "gs",    k,    0.0   );
        materials->gs_ref[k]= ReadMatProps( fin, "gsref" ,k,  2.0e-3  ) /scaling->L;
        // Strain softening
        materials->C_end[k]     = ReadMatProps( fin, "Ce",     k,    1.0e7  )  / scaling->S;
        materials->phi_end[k]   = ReadMatProps( fin, "phie",   k,     30.0  )  / 1.0 * M_PI / 180.0;
        materials->pls_start[k] = ReadMatProps( fin, "plss",   k,    1.0e6  );
        materials->pls_end[k]   = ReadMatProps( fin, "plse",   k,    1.0e6  );
        // Density models
        materials->density_model[k]     = (int)ReadMatProps( fin, "density_model",     k,    1  );
        materials->phase_diagram[k]     = (int)ReadMatProps( fin, "phase_diagram",     k,   -1  );
        
        // Check if any flow law is active
        int sum = abs(materials->cstv[k]) + abs(materials->pwlv[k]) + abs(materials->linv[k]) + abs(materials->gbsv[k]) + abs(materials->expv[k]);
        if ( sum == 0 ) {
            printf ("Phase %0d has no determined flow mechanism\n Simulation will end now!\n", k);
            exit(12);
        }
        
        // Print material parameters
        printf("----------------------------------------- MODEL DOMAIN ------------------------------------------\n");
        printf("Xmin   = %2.1lf  km         Xmax   = %2.1lf  km     Nx   = %3d    dx   = %.2lf m\n", (model->xmin*scaling->L)/1e3, (model->xmax*scaling->L)/1e3, model->Nx, model->dx*scaling->L);
        printf("Zmin   = %2.1lf  km         Zmax   = %2.1lf  km      Nz   = %3d    dz   = %.2lf m\n", (model->zmin*scaling->L)/1e3, (model->zmax*scaling->L)/1e3, model->Nz, model->dz*scaling->L );
        printf("-------------------------------------------- PHASE: %d -------------------------------------------\n", k);
        printf("rho    = %2.2e kg/m^3     mu = %2.2e Pa\n", materials->rho[k]*scaling->rho, materials->mu[k]*scaling->S );
        printf("Cv     = %2.2e J/kg/K      k = %2.2e W/m/K      Qr = %2.2e W/m3\n", materials->Cv[k]*scaling->Cv, materials->k[k]*scaling->k, materials->Qr[k]*(scaling->W / pow(scaling->L,3)) );
        printf("C      = %2.2e Pa        phi = %2.2e deg      Slim = %2.2e Pa\n",  materials->C[k]*scaling->S, materials->phi[k]*180/M_PI, materials->Slim[k]*scaling->S );
        printf("alp    = %2.2e 1/T        T0 = %2.2e K         bet = %2.2e 1/Pa       P0 = %2.2e Pa       drho = %2.2e kg/m^3 \n", materials->alp[k]*(1/scaling->T), materials->T0[k]*(scaling->T), materials->bet[k]*(1/scaling->S), materials->P0[k]*(scaling->S), materials->drho[k]*scaling->rho );
        printf("prefactor for power-law: %2.2e\n", materials->pref_pwl[k]);
                 printf("C_end    = %2.2e Pa        Phi_end = %2.2e deg         pls_start = %2.2e        pls_end = %2.2e \n", materials->C_end[k]*scaling->S, materials->phi_end[k]*180/M_PI, materials->pls_start[k],  materials->pls_end[k] );
        
        printf("Flow law settings:\n");
        if ( abs(materials->cstv[k])>0 ) printf("--->    Constant viscosity activated \n");
        if ( abs(materials->pwlv[k])>0 ) printf("--->   Power law viscosity activated \n");
        if ( abs(materials->linv[k])>0 ) printf("--->      Linear viscosity activated \n");
        if ( abs(materials->gbsv[k])>0 ) printf("--->         GBS viscosity activated \n");
        if ( abs(materials->expv[k])>0 ) printf("---> Exponential viscosity activated \n");
        if ( abs(gsel)              >0 ) printf("--->  Grain size evolution activated \n");
        
        // Call flow law data base
        if ( abs(materials->pwlv[k])>0 ) ReadDataPowerLaw   ( materials, model, k, materials->pwlv[k], scaling );
        if ( abs(materials->linv[k])>0 ) ReadDataLinear     ( materials, model, k, materials->linv[k], scaling );
        if ( abs(materials->gbsv[k])>0 ) ReadDataGBS        ( materials, model, k, materials->gbsv[k], scaling );
        if ( abs(materials->expv[k])>0 ) ReadDataExponential( materials, model, k, materials->expv[k], scaling );
        if ( abs(materials->gs[k])  >0 ) ReadDataGSE        ( materials, model, k, materials->gs[k], scaling );
        
        if ( abs(materials->cstv[k])>0 ) {
            materials->eta0[k]  /= scaling->eta;
            printf("eta0 = %2.2e Pa.s\n", materials->eta0[k]*scaling->eta);
        }
    }
    materials->R = Rg / (scaling->J/scaling->T);
    
    //------------------------------------------------------------------------------------------------------------------------------//
    // PHASE DIAGRAM INFO
    //------------------------------------------------------------------------------------------------------------------------------//
    model->isPD = 0;
    
    for ( k=0; k<materials->Nb_phases; k++) {
        printf("Phase %d ---  density model %d \n",k, materials->density_model[k]);
        if ( materials->density_model[k] == 2 ) {
            printf("The density model of phase %0d relies on phase diagrams\n", k);
            model->isPD = 1;
        }
        if ( materials->density_model[k] == 2 && materials->phase_diagram[k] == -1 ) {
            printf("However the phase diagram index was not set (default -1)\n");
            printf("Therefore the simulation will not run, go fix 'phase_diagram' of phase %02d\n", k);
            exit(5);
        }
    }
    
    // Phase diagrams
    if ( model->isPD == 1 ) {
        
        printf("Loading phase_diagrams...\n");
        int pid;
        model->num_PD = 3;
        
        // Allocate
        AllocatePhaseDiagrams( model );
        
        /**** PHASE DIAGRAMS #00 - Mantle (Jenadi_stx.dat)  ****/
        pid                       = 0;         // Kaus & Connolly, 2005: Effect of mineral phase transitions on sedimentary basin subsidence and uplift
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 1000;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 1000;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 473.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 2273.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 100e6/scaling->S;         // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 15e9 /scaling->S;         // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/Hawaiian_Pyrolite_rho_bin.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);
        
        /**** PHASE DIAGRAMS #00 - Mantle (Jenadi_stx_HR.dat)  ****/
        pid                       = 1;         // Kaus & Connolly, 2005: Effect of mineral phase transitions on sedimentary basin subsidence and uplift
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 1500;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 1500;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 273.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 2273.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 1e5/scaling->S;           // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 25e9 /scaling->S;         // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/Hawaiian_Pyrolite_HR_rho_bin.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);
        
        /**** PHASE DIAGRAMS #01 - Basalt (MORB_L.dat)  ****/
        pid                       = 2;  // Water saturated MORB - Bulk composition taken from Schmidt & Poli 1998 EPSL (Table 1)
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 1000;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 1000;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 573.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 1273.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 100e6/scaling->S;         // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 5.1e9 /scaling->S;        // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/MORB_H2Osat_rho_bin.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);
        
    }
    
    //------------------------------------------------------------------------------------------------------------------------------//
    // DEFORMATION MAP PARAMETERS
    //------------------------------------------------------------------------------------------------------------------------------//
    
    // Create deformation maps or not (default no)
    model->def_maps        = ReadInt2(fin, "def_maps", 0);
    
    // Resolution
    model->nT              = ReadInt2(fin, "nT", 11);
    model->nE              = ReadInt2(fin, "nE", 11);
    model->nd              = ReadInt2(fin, "nd", 11);
    
    // Temperature, strain rate, grain size MIN/MAX & pressure
    model->Tmin            = ReadDou2(fin, "Tmin", 100.0);  model->Tmin += zeroC; model->Tmin /= scaling->T;    // C -> K & non-dimensionalization
    model->Tmax            = ReadDou2(fin, "Tmax", 1000.0); model->Tmax += zeroC; model->Tmax /= scaling->T;    // C -> K & non-dimensionalization
    model->Emin            = ReadDou2(fin, "Emin", -30.0);  //model->Emin /= scaling->E;                           // Non-dimensionalization
    model->Emax            = ReadDou2(fin, "Emax", -4.0 );   //model->Emax /= scaling->E;                           // Non-dimensionalization
    model->dmin            = ReadDou2(fin, "dmin", -7.0 );   //model->dmin /= scaling->L;                           // Non-dimensionalization
    model->dmax            = ReadDou2(fin, "dmax", -2.0 );   //model->dmax /= scaling->L;                           // Non-dimensionalization
    model->Pn              = ReadDou2(fin, "Pn",  5.0e8 );      model->Pn   /= scaling->S;                           // Non-dimensionalization
    
    //------------------------------------------------------------------------------------------------------------------------------//
    // NUMERICAL PARAMETERS
    //------------------------------------------------------------------------------------------------------------------------------//
    
    // Particles
    particles->Nx_part       = ReadInt2( fin, "Nx_part", 0 );
    particles->Nz_part       = ReadInt2( fin, "Nz_part", 0 );
    particles->min_part_cell = ReadInt2( fin, "min_part_cell", 0 );
    particles->Nb_part       = (model->Nx-1)*(model->Nz-1) * particles->Nx_part * particles->Nz_part;
    particles->Nb_part_max   = 4.1*particles->Nb_part;
    
    // Nonlinear iteration parameters
    model->Newton           = ReadInt2( fin, "Newton", 0 );
    Nmodel->nit_max         = ReadInt2( fin, "nit_max", 1 );
    Nmodel->tol_u           = ReadDou2( fin, "tol_u", 1.0e-5 );// / (scaling->F/pow(scaling->L,3.0));
    Nmodel->tol_p           = ReadDou2( fin, "tol_p", 1.0e-5 );// / scaling->E;
    model->mineta           = ReadDou2( fin, "mineta", 1.0e18 ) / scaling->eta;
    model->maxeta           = ReadDou2( fin, "maxeta", 1.0e24 ) / scaling->eta;
    Nmodel->stagnated       = 0;
    
    // Direct solver parameters
    model->lsolver          = ReadInt2( fin, "lsolver", 0 );
    if ( model->Newton == 1 ) model->lsolver         = 1;

    // Close input file
    fclose(fin);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ScaleMe( scale* scale) {
    
    scale->t    = scale->L / scale->V;
    scale->a    = scale->V / scale->t;
    scale->E    = 1.0 / scale->t;
    scale->S    = scale->eta / scale->t;
    scale->m    = scale->S * scale->L * pow(scale->t,2.0);
    scale->rho  = scale->m / pow(scale->L,3.0);
    scale->F    = scale->m * scale->L / pow(scale->t, 2.0);
    scale->J    = scale->m * pow(scale->L,2.0) / pow(scale->t,2.0);
    scale->W    = scale->J / scale->t;
    scale->Cv   = scale->J / scale->m / scale->T;
    scale->k    = scale->W / scale->L / scale->T;
    scale->rhoE = scale->J / scale->m * scale->rho;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double* ReadBin( char A_name[], int nx, int ny, double scale ){
    double *A;
    char* bname; size_t nb_elems = nx*ny; FILE* fid; asprintf(&bname, "%s", A_name);
    A = malloc((nb_elems)*sizeof(double));
    fid=fopen(bname, "rb"); // Open file
    if (!fid){
        fprintf(stderr, "\nUnable to open file %s. I will exit here ... \n", bname); exit(2);
    }
    fread(A, sizeof(double), nb_elems, fid); fclose(fid); free(bname);
    ArrayTimesScalar( A, 1.0/scale, nb_elems );
    return A;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

char* ReadChar( FILE *fin, char FieldName[], char Default[] ) {
    int h, find=0, bufmax=50, length;
    char *string1, *string2;
    char line[bufmax];
    
    string1 = malloc(sizeof(char)*bufmax);
    
    // Get the length of the demanded string
    length = strlen( FieldName );
    
    // Buffer array to contain the string to compare with
    char *param1, *param2;
    param2 = malloc( (length+1)*sizeof(char));
    asprintf(&param1, "%s", FieldName);
    
    // Initialise line
    for (h=0;h<bufmax;h++) {
        line[h]='\0';
    }
    
    // Seach for sign 'equal' in the lines
    while ( find == 0 ) {
        
        // Read new line
        fgets ( line, sizeof(line), fin );
        
        if (feof(fin)) {
            printf("Warning : Parameter '%s' not found in the setup file, running with default value %s\n", FieldName, Default);
            rewind (fin);
            return Default;
        }
        
        // Get the first 'length' characters of the line
        for (h=0;h<length;h++) {
            param2[h] = line[h];
        }
        param2[length] = '\0';
        
        // Check if we found the right parameter name
        if ( strcmp(param1, param2) == 0 ) {
            find = 1;
            int str_size = 0;
            int h1;
            // Search for equal sign
            for (h=0;h<bufmax;h++) {
                if(strlen(line)> 0 && line[h]=='=') {
                    
                    for (h1=0; h1<30; h1++) {
                        if ( isspace(line[h+2+h1]) ) {
                            string1[h1] = '\0';
                            str_size++;
                            break;
                        }
                        else {
                            string1[h1] = (line[h+2+h1]);
                            str_size++;
                        }
                    }
                    
                    string2 = malloc((str_size+1)*sizeof(char));
                    for (h1=0; h1<str_size+1; h1++) {
                        string2[h1] = string1[h1];
                    }
   
                    free(param1);
                    free(param2);
                    free(string1);
                    return string2;
                }
            }
        }
    }
    free(param1);
    free(param2);
    return Default;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

char* ReadPhaseDiagram( FILE *fin, char FieldName[] ) {
    int h, find=0, bufmax=50, length;
    char *string1, *string2;
    char line[bufmax];
    
    string1 = malloc(sizeof(char)*bufmax);
    
    // Get the length of the demanded string
    length = strlen( FieldName );
    
    // Buffer array to contain the string to compare with
    char *param1, *param2;
    param2 = malloc( (length+1)*sizeof(char));
    asprintf(&param1, "%s", FieldName);
    
    // Initialise line
    for (h=0;h<bufmax;h++) {
        line[h]='\0';
    }
    
    // Seach for sign 'equal' in the lines
    while ( find == 0 ) {
        
        // Read new line
        fgets ( line, sizeof(line), fin );
        
        if (feof(fin)) {
            printf("Error: The phase diagram '%s' could not be found in the setup file. I will exit here.\n", FieldName);
            rewind (fin    );
            free   (param1 );
            free   (param2 );
            free   (string1);
            exit(2);
        }
        
        // Get the first 'length' characters of the line
        for (h=0;h<length;h++) {
            param2[h] = line[h];
        }
        param2[length] = '\0';
        
        // Check if we found the right parameter name
        if ( strcmp(param1, param2) == 0 ) {
            find = 1;
            int str_size = 0;
            int h1;
            // Search for equal sign
            for (h=0;h<bufmax;h++) {
                if(strlen(line)> 0 && line[h]=='=') {
                    
                    for (h1=0; h1<30; h1++) {
                        if ( isspace(line[h+2+h1]) ) {
                            //printf("found space last char is %c\n", line[h+2+h1-1]);
                            string1[h1] = '\0';
                            str_size++;
                            break;
                        }
                        else {
                            string1[h1] = (line[h+2+h1]);
                            str_size++;
                        }
                    }
                    
                    string2 = malloc((str_size+1)*sizeof(char));
                    for (h1=0; h1<str_size+1; h1++) {
                        string2[h1] = string1[h1];
                    }
                    //printf("returned string is %s %s %d\n", string1, string2, str_size );
                    free(param1);
                    free(param2);
                    free(string1);
                    return string2;
                }
            }
        }
    }
    free(param1);
    free(param2);
    return NULL;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int ReadInt2( FILE *fin, char FieldName[], int Default )
{
    // Some declaration.
    int     bufmax=1000;
    int     h = 0;
    char    line[bufmax];
    int     value = 0;
    
    // Start from beginning of the file.
    rewind(fin);
    
    // Buffer array to contain the string to compare with
    char *param1;
    asprintf(&param1, "%s", FieldName);
    
    // Loop over all lines in input file
    while(value == 0)
    {
        // Read new line
        fgets ( line, sizeof(line), fin );
        if (feof(fin)) {
            printf("Warning : Parameter '%s' not found in the setup file, running with default value %d\n", FieldName, Default);
            rewind (fin);
            free(param1);
            return Default;
        }
        
        // Determine the length of the parameter string in the current line.
        int InPar_length = 0;
        
        while(line[InPar_length] != ' ')
        {
            InPar_length = InPar_length + 1;
        }
        
        // Allocate memory to save the current parameter string.
        char *param2;
        param2 = malloc( (InPar_length + 1)*sizeof(char));
        
        // Save the current parameter string.
        for (h = 0; h < InPar_length; h++)
        {
            param2[h] = line[h];
        }
        param2[InPar_length] = '\0';
        
        // Find match.
        if( (strcmp(param1, param2)) == 0 )
        {
            // Search for equal sign
            for (h = 0; h < bufmax ;h++)
            {
                if(strlen(line)> 0 && line[h]=='=')
                {
                    value = atoi(&line[h+1]);
                    free(param1);
                    free(param2);
                    rewind(fin);
                    return value;
                }
            }
        }
        free(param2);
    }
    rewind(fin);
    free(param1);
    return Default;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double ReadDou2( FILE *fin, char FieldName[], double Default )
{
    // Some declaration.
    int     bufmax=1000;
    int     h = 0;
    char    line[bufmax];
    double  value = 0.0;
    
    // Start from beginning of the file.
    rewind(fin);
    
    // Buffer array to contain the string to compare with
    char *param1;
    asprintf(&param1, "%s", FieldName);
    
    // Loop over all lines in input file
    while(value == 0)
    {
        // Read new line
        fgets ( line, sizeof(line), fin );
        if (feof(fin)) {
            printf("Warning : Parameter '%s' not found in the setup file, running with default value %2.2e\n", FieldName, Default);
            rewind (fin);
            free(param1);
            return Default;
        }
        
        // Determine the length of the parameter string in the current line.
        int InPar_length = 0;
        while(line[InPar_length] != ' ')
        {
            InPar_length = InPar_length + 1;
        }
        
        // Allocate memory to save the current parameter string.
        char *param2;
        param2 = malloc( (InPar_length + 1)*sizeof(char));
        
        // Save the current parameter string.
        for (h = 0; h < InPar_length; h++)
        {
            param2[h] = line[h];
        }
        param2[InPar_length] = '\0';
        
        // Find match.
        if( (strcmp(param1, param2)) == 0 )
        {
            // Search for equal sign
            for (h = 0; h < bufmax ;h++)
            {
                if(strlen(line)> 0 && line[h]=='=')
                {
                    value = atof(&line[h+1]);
                    free(param1);
                    free(param2);
                    rewind(fin);
                    return value;
                }
            }
        }
        free(param2);
    }
    rewind(fin);
    free(param1);
    return Default;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double ReadMatProps( FILE *fin, char FieldName[], int PhaseID, double Default )
{
    // Some declarations.
    int     bufmax = 1000;
    int     h = 0, ID_value, find = 0, subfind = 0;
    char    line[bufmax], phase_line[bufmax];
    double  value = 0.0;
    
    // Find the current phase ID starting from the beginning of the input file.
    rewind(fin);
    
    // Buffer array to contain the string to compare with
    char *param1;
    asprintf(&param1, "%s", FieldName);
    
    // Loop over all lines in input file.
    while(find == 0)
    {
        // Read new line and check for existence of phases. If not exit.
        fgets ( line, sizeof(line), fin );
        if (feof(fin))
        {
            printf("Warning : No phase ID found! I will exit here.");
            free(param1);
            rewind (fin);
            exit(0);
        }
        
        // Determine the length of the parameter string in the current line.
        int InPar_length = 0;
        
        while(line[InPar_length] != ' ')
        {
            InPar_length = InPar_length + 1;
        }
        
        // Allocate memory to save the current parameter string.
        char *param2;
        param2 = malloc( (InPar_length + 1)*sizeof(char) );
        
        // Save the current parameter string.
        for (h = 0; h < InPar_length; h++)
        {
            param2[h] = line[h];
        }
        param2[InPar_length] = '\0';
        
        // Find match.
        if( (strcmp("ID", param2)) == 0 )
        {
            // Search for equal sign
            for (h = 0; h < bufmax ;h++)
            {
                if(strlen(line)> 0 && line[h] == '=')
                {
                    ID_value = atoi(&line[h+1]);
                    break;
                }
            }
            
            // Read parameters if current phase.
            if (ID_value == PhaseID)
            {
                // Loop over all lines in input file
                while(subfind == 0)
                {
                    // Read new line
                    fgets ( phase_line, sizeof(phase_line), fin );
                    
                    // Determine the length of the parameter string in the current line.
                    InPar_length = 0;
                    while(phase_line[InPar_length] != ' ')
                    {
                        InPar_length = InPar_length + 1;
                    }
                    
                    // Allocate memory to save the current parameter string.
                    char *param3;
                    param3 = malloc( (InPar_length + 1)*sizeof(char));
                    
                    // Save the current parameter string.
                    for (h = 0; h < InPar_length; h++)
                    {
                        param3[h] = phase_line[h];
                    }
                    param3[InPar_length] = '\0';
                    
                    // Break in case the parameter has not been defined for the current phase.
                    if ( strcmp(param3,"ID") == 0 || feof(fin) ) {
                        printf("Warning : Parameter '%s' not found in the setup file, running with default value %.2lf\n", FieldName, Default);
                        rewind (fin);
                        free(param1);
                        free(param2);
                        free(param3);
                        return Default;
                    }
                    
                    // Find match.
                    if( (strcmp(param1, param3)) == 0 )
                    {
                        // Search for equal sign
                        for (h = 0; h < bufmax ;h++)
                        {
                            
                            if(strlen(line)> 0 && phase_line[h]=='=')
                            {
                                value = atof(&phase_line[h+1]);
                                free(param1);
                                free(param2);
                                free(param3);
                                return value;
                            }
                        }
                        subfind = 1;
                    }
                    free(param3);
                }
                free(param1);
            }
        }
        free(param2);
    }
    free(param1);
    return Default;
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
