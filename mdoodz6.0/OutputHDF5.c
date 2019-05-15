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
#include "header_MDOODZ.h"
#include "hdf5.h"
#include "zlib.h"

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

// Scale parameters to physical scale
void ScaleBack(float* FieldF, double scale, int size) {
    int k;
    for (k=0; k<size; k++) {
        FieldF[k] = FieldF[k]*scale;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Scale parameters to physical scale
void ScaleBackD(double* FieldD, double scale, int size) {
    int k;
    for (k=0; k<size; k++) {
        FieldD[k] = FieldD[k]*scale;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Function that cast double arrays to float (used for visualisation purposes only)
void DoubleToFloat(double* FieldD, float* FieldF, int size) {
    int k;
    for (k=0; k<size; k++) {
        FieldF[k] = (float)FieldD[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// The following routines were originally writen by D. A. May, and used in I2VIS //
void create_output_hdf5( const char name[] )
{
	hid_t       file_id;   /* file identifier */
	
	/* Create a new file using default properties. */
	file_id = H5Fcreate( name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
	/* Terminate access to the file. */
	H5Fclose(file_id);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddGroup_to_hdf5( const char filename[], const char group[] )
{
	hid_t       file_id, group_id;  /* identifiers */
	char        *group_name;
	
	asprintf( &group_name, "/%s", group );
	
	/* Open exisiting file */
	file_id = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT);
	
	/* Create group "ParticleGroup" in the root group using absolute name. */
	group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	
	/* Close group. */
	H5Gclose(group_id);
	
	/* Close the file. */
	H5Fclose(file_id);
	
	free( group_name );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddFieldToGroup_generic( int compress, const char filename[], const char group[], const char field[], char d_type, int np, void* data, int dim )
{
	hid_t   file_id, particle_group_id, coord_dataset_id, coord_dataspace_id;  /* identifiers */
	hsize_t length;
	char *dataset_name;
	char *group_name;
	hid_t    plist=0,SET_CREATION_PLIST;
	hsize_t  cdims[2] = {0,0};
	double percentage_chunk;
	double chunk_size;
	int deflation_level, quiet=1;
	
	asprintf( &group_name, "/%s", group );
	asprintf( &dataset_name, "%s/%s", group, field );
	
	length = dim * np;
	
	/* Open exisiting file */
	file_id = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT);
	
	/* Open group "ParticleGroup" */
	particle_group_id = H5Gopen(file_id, group_name,H5P_DEFAULT);
	
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
	percentage_chunk = 5.0;
	chunk_size = percentage_chunk*((double)np)/((double)100.0 );
	if( chunk_size < 1 ) {
		cdims[0] = 1;
	}
	else {
		cdims[0] = (int)chunk_size;
	}
	
	plist  = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, 1, cdims);
	deflation_level = 4;
    H5Pset_deflate( plist, deflation_level);
    
	/* Create the data space for the dataset. */
	coord_dataspace_id = H5Screate_simple( 1, &length, NULL);
	
	if(compress==_TRUE_) {
		SET_CREATION_PLIST = plist;
        if ( quiet == 0 ) {
            printf("*** Compression info *** \n");
            printf("  chunk_size = %f \n", chunk_size );
            printf("  deflation level = %d \n", deflation_level );
        }
	}
	else {
		SET_CREATION_PLIST = H5P_DEFAULT;
	}
	
	if( d_type == 'd' ) {
		/* Create a dataset within "group". */
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_DOUBLE, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
	    
		/* Write the particle dataset. */
		H5Dwrite(coord_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
	}
	
	else if( d_type == 'c' ) {
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_CHAR, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
        H5Dwrite(coord_dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
		
	}
	
	else if( d_type == 'i' ) {
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_STD_I32BE, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
        H5Dwrite(coord_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
		
	}
	else if( d_type == 'f' ) {
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_FLOAT, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
        H5Dwrite(coord_dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
		
	}
	
	else {
		printf("ERROR: Only know how to write doubles (d), ints (i), or chars (c) \n");
		exit(1);
	}
	
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
	/* Close the data space for this particle dataset. */
	//status = H5Sclose(attr_dataspace_id);
    
	/* Close the particle coord dataset. */
	H5Dclose(coord_dataset_id);
    
	/* free(dataset_nameP); */
	/* } */
    
	/* Close datatspace. */
	H5Sclose(coord_dataspace_id);
    
	/* Close list. */
	H5Pclose(plist);
    
	/* Close group. */
	H5Gclose(particle_group_id);
	
	/* Close the file. */
	H5Fclose(file_id);
	
	free(group_name);
	free(dataset_name);
}
// The above routines were originally writen by D. A. May, and used in I2VIS //

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void WriteOutputHDF5( grid *mesh, markers *particles, surface *topo, markers* topo_chain, params model, char *txtout, mat_prop materials, scale scaling ) {
    
    char *name;
    double A[8], E[5];
    int k;
    double *strain, *strain_el, *strain_pl, *strain_pwl, *strain_exp, *strain_lin, *strain_gbs, *X;
    float *Crho_s, *Crho_n, *Ceta_s, *Ceta_n, *CVx, *CVz, *CP, *Csxxd, *Csxz, *Cexxd, *Cexz, *Cstrain, *Cstrain_el, *Cstrain_pl, *Cstrain_pwl, *Cstrain_exp, *Cstrain_lin, *Cstrain_gbs, *CT, *Cd;
    float *Cxg_coord, *Czg_coord, *Cxc_coord, *Czc_coord, *Czvx_coord, *Cxvz_coord;
    float *CeII_el, *CeII_pl, *CeII_pwl, *CeII_exp, *CeII_lin, *CeII_gbs, *CX;
    double *Fxx, *Fxz, *Fzx, *Fzz, *nx, *nz;
    float *CFxx, *CFxz, *CFzx, *CFzz, *Cnx, *Cnz;
    double *T0, *P0, *x0, *z0, *Tmax, *Pmax;
    float *CT0, *CP0, *Cx0, *Cz0, *CTmax, *CPmax;
    
    int    res_fact = 1;
    int    nxviz, nzviz, nxviz_hr, nzviz_hr;
    char  *compo, *compo_hr;
    float *Cxviz, *Czviz, *Cxviz_hr, *Czviz_hr, *Cxtopo, *Cztopo, *Cheight, *Ctopovx, *Ctopovz;
    double *P_total;
    float  *Ccohesion, *Cfriction;

    P_total  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));

    // Build total pressure
    for (k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
        P_total[k] = mesh->p_in[k];
    }

    // ---------------------------------------------------------
    // Genrate phase map with normal resolution
    res_fact = 1;
    nxviz = res_fact*(mesh->Nx-1) + 1;
    nzviz = res_fact*(mesh->Nz-1) + 1;
    double xviz[nxviz], zviz[nzviz];
    // Define visual grid
    xviz[0] = mesh->xg_coord[0];
    zviz[0] = mesh->zg_coord[0];
    
    for (k=1; k<nxviz; k++) {
        xviz[k] = xviz[k-1] + mesh->dx/res_fact;
    }
    for (k=1; k<nzviz; k++) {
        zviz[k] = zviz[k-1] + mesh->dz/res_fact;
    }
    
    compo  = DoodzMalloc( sizeof(char)*(nxviz-1)*(nzviz-1));
    // Closest point interpolation: marker phase --> visual nodes
    Interp_Phase2VizGrid ( *particles, particles->phase, mesh, compo, xviz, zviz, nxviz, nzviz, model, *topo );
    
    // ---------------------------------------------------------
    // Genrate phase map with double resolution
    res_fact = 2;
    nxviz_hr = res_fact*(mesh->Nx-1) + 1;
    nzviz_hr = res_fact*(mesh->Nz-1) + 1;
    double xviz_hr[nxviz_hr], zviz_hr[nzviz_hr];
    // Define visual grid
    xviz_hr[0] = mesh->xg_coord[0];
    zviz_hr[0] = mesh->zg_coord[0];
    
    for (k=1; k<nxviz_hr; k++) {
        xviz_hr[k] = xviz_hr[k-1] + mesh->dx/res_fact;
    }
    for (k=1; k<nzviz_hr; k++) {
        zviz_hr[k] = zviz_hr[k-1] + mesh->dz/res_fact;
    }
    
    compo_hr  = DoodzMalloc( sizeof(char)*(nxviz_hr-1)*(nzviz_hr-1));
    // Closest point interpolation: marker phase --> visual nodes
    Interp_Phase2VizGrid ( *particles, particles->phase, mesh, compo_hr, xviz_hr, zviz_hr, nxviz_hr, nzviz_hr, model, *topo );
    // ---------------------------------------------------------
    // Cast grid arrays
    Crho_s  = DoodzMalloc( sizeof(float)*model.Nx*model.Nz);
    DoubleToFloat( mesh->rho_app_s, Crho_s, model.Nx*model.Nz);
    ScaleBack( Crho_s, scaling.rho, model.Nx*model.Nz );
    
    Crho_n  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->rho_app_n, Crho_n, (model.Nx-1)*(model.Nz-1));
    ScaleBack( Crho_n, scaling.rho, (model.Nx-1)*(model.Nz-1));
    
    Ceta_s  = DoodzMalloc( sizeof(float)*model.Nx*model.Nz);
    DoubleToFloat( mesh->eta_phys_s, Ceta_s, model.Nx*model.Nz);
    ScaleBack( Ceta_s, scaling.eta, model.Nx*model.Nz );
    
    Ceta_n  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eta_phys_n, Ceta_n, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Ceta_n, scaling.eta, (model.Nx-1)*(model.Nz-1) );
    
    CVx  = DoodzMalloc( sizeof(float)*model.Nx*(model.Nz+1));
    DoubleToFloat( mesh->u_in, CVx, model.Nx*(model.Nz+1) );
    ScaleBack( CVx, scaling.V, model.Nx*(model.Nz+1) );
    
    CVz  = DoodzMalloc( sizeof(float)*(model.Nx+1)*model.Nz);
    DoubleToFloat( mesh->v_in, CVz, (model.Nx+1)*model.Nz );
    ScaleBack( CVz, scaling.V, (model.Nx+1)*model.Nz );
    
    CP  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( P_total, CP, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CP, scaling.S, (model.Nx-1)*(model.Nz-1) );
    
    Csxxd  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->sxxd, Csxxd, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Csxxd, scaling.S, (model.Nx-1)*(model.Nz-1) );
    
    Csxz  = DoodzMalloc( sizeof(float)*model.Nx*model.Nz);
    DoubleToFloat( mesh->sxz, Csxz, model.Nx*model.Nz );
    ScaleBack( Csxz, scaling.S, model.Nx*model.Nz );
    
    Cexxd  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->exxd, Cexxd, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Cexxd, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    Cexz  = DoodzMalloc( sizeof(float)*model.Nx*model.Nz);
    DoubleToFloat( mesh->exz, Cexz, model.Nx*model.Nz );
    ScaleBack( Cexz, scaling.E, model.Nx*model.Nz );
    
    CeII_el  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_el, CeII_el, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_el, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    CeII_pl  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_pl, CeII_pl, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_pl, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    CeII_pwl  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_pwl, CeII_pwl, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_pwl, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    CeII_exp  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_exp, CeII_exp, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_exp, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    CeII_lin  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_lin, CeII_lin, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_lin, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    CeII_gbs  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_gbs, CeII_gbs, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_gbs, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    Cd        = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->d, Cd, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Cd, scaling.L, (model.Nx-1)*(model.Nz-1) );
    
    //---------------------------------------------------
    CT  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->T, CT, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CT, scaling.T, (model.Nx-1)*(model.Nz-1) );
    
    //---------------------------------------------------
    Cxg_coord = DoodzMalloc( sizeof(float)*model.Nx);
    DoubleToFloat( mesh->xg_coord, Cxg_coord, model.Nx );
    ScaleBack( Cxg_coord, scaling.L, model.Nx );
    
    Czg_coord = DoodzMalloc( sizeof(float)*model.Nz);
    DoubleToFloat( mesh->zg_coord, Czg_coord, model.Nz );
    ScaleBack( Czg_coord, scaling.L, model.Nz );
    
    Cxc_coord = DoodzMalloc( sizeof(float)*(model.Nx-1));
    DoubleToFloat( mesh->xc_coord, Cxc_coord, model.Nx-1 );
    ScaleBack( Cxc_coord, scaling.L, model.Nx-1 );
    
    Czc_coord = DoodzMalloc( sizeof(float)*(model.Nz-1));
    DoubleToFloat( mesh->zc_coord, Czc_coord, model.Nz-1 );
    ScaleBack( Czc_coord, scaling.L, model.Nz-1 );
    
    Czvx_coord = DoodzMalloc( sizeof(float)*(model.Nz+1));
    DoubleToFloat( mesh->zvx_coord, Czvx_coord, model.Nz+1 );
    ScaleBack( Czvx_coord, scaling.L, model.Nz+1 );
    
    Cxvz_coord = DoodzMalloc( sizeof(float)*(model.Nx+1));
    DoubleToFloat( mesh->xvz_coord, Cxvz_coord, model.Nx+1 );
    ScaleBack( Cxvz_coord, scaling.L, model.Nx+1 );
    
    Cxviz = DoodzMalloc( sizeof(float)*nxviz);
    DoubleToFloat( xviz, Cxviz, nxviz );
    ScaleBack( Cxviz, scaling.L, nxviz );
    
    Czviz = DoodzMalloc( sizeof(float)*nzviz);
    DoubleToFloat( zviz, Czviz, nzviz );
    ScaleBack( Czviz, scaling.L, nzviz );
    
    Cxviz_hr = DoodzMalloc( sizeof(float)*nxviz_hr);
    DoubleToFloat( xviz_hr, Cxviz_hr, nxviz_hr );
    ScaleBack( Cxviz_hr, scaling.L, nxviz_hr );
    
    Czviz_hr = DoodzMalloc( sizeof(float)*nzviz_hr);
    DoubleToFloat( zviz_hr, Czviz_hr, nzviz_hr );
    ScaleBack( Czviz_hr, scaling.L, nzviz_hr );
    
    // Total strain
    strain  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
    Interp_P2C ( *particles,  particles->strain, mesh, strain, mesh->xg_coord, mesh->zg_coord, 1, 0 );
    Cstrain  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain, Cstrain, (model.Nx-1)*(model.Nz-1) );
    
    // Elastic strain
    strain_el  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
    Interp_P2C ( *particles,  particles->strain_el, mesh, strain_el, mesh->xg_coord, mesh->zg_coord, 1, 0 );
    Cstrain_el  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_el, Cstrain_el, (model.Nx-1)*(model.Nz-1) );
    
    // Plastic strain
    strain_pl  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
    Interp_P2C ( *particles,  particles->strain_pl, mesh, strain_pl, mesh->xg_coord, mesh->zg_coord, 1, 0 );
    Cstrain_pl  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_pl, Cstrain_pl, (model.Nx-1)*(model.Nz-1) );
    
    // Power-law strain
    strain_pwl  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
    Interp_P2C ( *particles,  particles->strain_pwl, mesh, strain_pwl, mesh->xg_coord, mesh->zg_coord, 1, 0 );
    Cstrain_pwl  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_pwl, Cstrain_pwl, (model.Nx-1)*(model.Nz-1) );
    
    // Exponential flow strain
    strain_exp  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
    Interp_P2C ( *particles,  particles->strain_exp, mesh, strain_exp, mesh->xg_coord, mesh->zg_coord, 1, 0 );
    Cstrain_exp  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_exp, Cstrain_exp, (model.Nx-1)*(model.Nz-1) );
    
    // Linear flow strain
    strain_lin  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
    Interp_P2C ( *particles,  particles->strain_lin, mesh, strain_lin, mesh->xg_coord, mesh->zg_coord, 1, 0 );
    Cstrain_lin  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_lin, Cstrain_lin, (model.Nx-1)*(model.Nz-1) );
    
    // GBS flow strain
    strain_gbs  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
    Interp_P2C ( *particles,  particles->strain_gbs, mesh, strain_gbs, mesh->xg_coord, mesh->zg_coord, 1, 0 );
    Cstrain_gbs  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_gbs, Cstrain_gbs, (model.Nx-1)*(model.Nz-1) );
    
    if ( model.fstrain == 1 ) {
        // Fxx
        Fxx  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->Fxx, mesh, Fxx, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        CFxx  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Fxx, CFxx, (model.Nx-1)*(model.Nz-1) );
        
        // Fxz
        Fxz  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->Fxz, mesh, Fxz, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        CFxz  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Fxz, CFxz, (model.Nx-1)*(model.Nz-1) );
        
        // Fzx
        Fzx  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->Fzx, mesh, Fzx, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        CFzx  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Fzx, CFzx, (model.Nx-1)*(model.Nz-1) );
        
        // Fzz
        Fzz  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->Fzz, mesh, Fzz, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        CFzz = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Fzz, CFzz, (model.Nx-1)*(model.Nz-1) );
    }
    
    if ( model.aniso == 1 ) {
        
        // nx
        nx  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->nx, mesh, nx, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        Cnx = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( nx, Cnx, (model.Nx-1)*(model.Nz-1) );
        
        // nz
        nz  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->nz, mesh, nz, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        Cnz = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( nz, Cnz, (model.Nx-1)*(model.Nz-1) );
        
    }
    
    if ( model.rec_T_P_x_z == 1 ) {
        // T0
        T0  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->T0, mesh, T0, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        for (k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
            if ( mesh->BCt.type[k] == 30 ) T0[k]   = zeroC/scaling.T;
        }
        CT0  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( T0, CT0, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( CT0, scaling.T, (model.Nx-1)*(model.Nz-1) );
        

        // P0
        P0  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->P0, mesh, P0, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        for (k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
            if ( mesh->BCt.type[k] == 30 ) P0[k]   = 0.0;
        }
        CP0  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( P0, CP0, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( CP0, scaling.S, (model.Nx-1)*(model.Nz-1) );
        
        // Tmax
        Tmax  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->Tmax, mesh, Tmax, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        for (k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
            if ( mesh->BCt.type[k] == 30 ) Tmax[k] = zeroC/scaling.T;
        }
        CTmax  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Tmax, CTmax, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( CTmax, scaling.T, (model.Nx-1)*(model.Nz-1) );
        
        
        // Pmax
        Pmax  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->Pmax, mesh, Pmax, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        for (k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
            if ( mesh->BCt.type[k] == 30 ) Pmax[k] = 0.0;
        }
        CPmax  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Pmax, CPmax, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( CPmax, scaling.S, (model.Nx-1)*(model.Nz-1) );

        
        // x0
        x0  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->x0, mesh, x0, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        Cx0  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( x0, Cx0, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( Cx0, scaling.L, (model.Nx-1)*(model.Nz-1) );
        
        // z0
        z0  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
        Interp_P2C ( *particles,  particles->z0, mesh, z0, mesh->xg_coord, mesh->zg_coord, 1, 0 );
        Cz0 = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( z0, Cz0, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( Cz0, scaling.L, (model.Nx-1)*(model.Nz-1) );

    }
    
    if (model.isPl_soft == 1) {
        Cfriction  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( mesh->phi_n, Cfriction, (model.Nx-1)*(model.Nz-1) );
        Ccohesion  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( mesh->C_n, Ccohesion, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( Ccohesion, scaling.S, (model.Nx-1)*(model.Nz-1) );
    }
    
    // Get X from particles
    X  = DoodzCalloc( sizeof(double),(model.Nx-1)*(model.Nz-1));
    Interp_P2C ( *particles,  particles->X, mesh, X, mesh->xg_coord, mesh->zg_coord, 1, 0 );
    CX  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( X, CX, (model.Nx-1)*(model.Nz-1) );
    
    
    //---------------------------------------------------
    
    // Topography
    if ( model.free_surf == 1 ) {
        
        Cxtopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->x, Cxtopo, topo_chain->Nb_part );
        ScaleBack( Cxtopo, scaling.L, topo_chain->Nb_part );
        
        Cztopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->z, Cztopo, topo_chain->Nb_part );
        ScaleBack( Cztopo, scaling.L, topo_chain->Nb_part );
        
        Cheight = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo->height, Cheight, model.Nx );
        ScaleBack( Cheight, scaling.L, model.Nx );
        
        Ctopovx = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo->vx, Ctopovx, model.Nx );
        ScaleBack( Ctopovx, scaling.V, model.Nx );
        
        Ctopovz = DoodzMalloc( sizeof(float)*(model.Nx+1));
        DoubleToFloat( topo->vz, Ctopovz, model.Nx+1 );
        ScaleBack( Ctopovz, scaling.V, model.Nx+1 );
        
    }
    
    // Generate file name
    asprintf( &name, "%s%05d%s",txtout, model.step, ".gzip.h5");
    create_output_hdf5( name );
    
    // Add groups
    AddGroup_to_hdf5( name, "Model" );
    AddGroup_to_hdf5( name, "Vertices" );
    AddGroup_to_hdf5( name, "Centers" );
    AddGroup_to_hdf5( name, "VxNodes" );
    AddGroup_to_hdf5( name, "VzNodes" );
    AddGroup_to_hdf5( name, "Particles" );
    AddGroup_to_hdf5( name, "VizGrid" );
    AddGroup_to_hdf5( name, "Topo" );
    
    // Model Parameters
    A[0] = (double)(model.time) * scaling.t;
    A[1] = (double)(model.xmax - model.xmin) * scaling.L;
    A[2] = (double)(model.zmax - model.zmin) * scaling.L;
    A[3] = (double)model.Nx;
    A[4] = (double)model.Nz;
    A[5] = (double)model.dx * scaling.L;
    A[6] = (double)model.dz * scaling.L;
    A[7] = (double)model.dt * scaling.t;
    
    E[0] = model.time*scaling.t;
    E[1] = (model.L0-(model.xmax - model.xmin))/model.L0 * 100.0;
    E[2] = mesh->Ut*scaling.S*scaling.L*scaling.L;
    E[3] = mesh->Ue*scaling.S*scaling.L*scaling.L;
    E[4] = mesh->W *scaling.S*scaling.L*scaling.L;
    
    // Parameter array
    AddFieldToGroup_generic( _TRUE_, name, "Model", "Params"   , 'd', 8,  A, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Model", "Energy"   , 'd', 5,  E, 1 );

    // Grid coordinate arrays
    AddFieldToGroup_generic( _TRUE_, name, "Model", "xg_coord" , 'f', model.Nx,    Cxg_coord,  1 );
    AddFieldToGroup_generic( _TRUE_, name, "Model", "zg_coord" , 'f', model.Nz,    Czg_coord,  1 );
    AddFieldToGroup_generic( _TRUE_, name, "Model", "xc_coord" , 'f', model.Nx-1,  Cxc_coord,  1 );
    AddFieldToGroup_generic( _TRUE_, name, "Model", "zc_coord" , 'f', model.Nz-1,  Czc_coord,  1 );
    AddFieldToGroup_generic( _TRUE_, name, "Model", "xvz_coord", 'f', model.Nx+1,  Cxvz_coord, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Model", "zvx_coord", 'f', model.Nz+1,  Czvx_coord, 1 );

    // Visualisation grid
    AddFieldToGroup_generic( _TRUE_, name, "VizGrid", "xviz"    , 'f', nxviz, Cxviz,     1 );
    AddFieldToGroup_generic( _TRUE_, name, "VizGrid", "zviz"    , 'f', nzviz, Czviz,     1 );
    AddFieldToGroup_generic( _TRUE_, name, "VizGrid", "xviz_hr" , 'f', nxviz_hr, Cxviz_hr,  1 );
    AddFieldToGroup_generic( _TRUE_, name, "VizGrid", "zviz_hr" , 'f', nzviz_hr, Czviz_hr,  1 );
    AddFieldToGroup_generic( _TRUE_, name, "VizGrid", "compo"   , 'c', (nxviz-1)*(nzviz-1), compo,    1 );
    AddFieldToGroup_generic( _TRUE_, name, "VizGrid", "compo_hr", 'c', (nxviz_hr-1)*(nzviz_hr-1), compo_hr,    1 );
    
    // Add casted grid fields
    AddFieldToGroup_generic( _TRUE_, name, "Vertices", "rho_s", 'f', model.Nx*model.Nz,         Crho_s, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "rho_n", 'f', (model.Nx-1)*(model.Nz-1), Crho_n, 1 );
    
    AddFieldToGroup_generic( _TRUE_, name, "Vertices", "eta_s", 'f', model.Nx*model.Nz,         Ceta_s, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "eta_n", 'f', (model.Nx-1)*(model.Nz-1), Ceta_n, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "VxNodes" , "Vx"   , 'f', model.Nx*(model.Nz+1),     CVx, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "VzNodes" , "Vz"   , 'f', (model.Nx+1)*model.Nz,     CVz, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "P"    , 'f', (model.Nx-1)*(model.Nz-1), CP, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "sxxd" , 'f', (model.Nx-1)*(model.Nz-1), Csxxd, 1 );
    
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "strain",    'f', (model.Nx-1)*(model.Nz-1), Cstrain, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "strain_el", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_el, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "strain_pl", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_pl, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "strain_pwl", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_pwl, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "strain_exp", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_exp, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "strain_lin", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_lin, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "strain_gbs", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_gbs, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "T",     'f', (model.Nx-1)*(model.Nz-1), CT, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Vertices", "sxz"  , 'f', model.Nx*model.Nz,         Csxz, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "exxd" , 'f', (model.Nx-1)*(model.Nz-1), Cexxd, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Vertices", "exz"  , 'f', model.Nx*model.Nz,         Cexz, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "eII_el" , 'f', (model.Nx-1)*(model.Nz-1), CeII_el, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "eII_pl" , 'f', (model.Nx-1)*(model.Nz-1), CeII_pl, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "eII_pwl" , 'f', (model.Nx-1)*(model.Nz-1), CeII_pwl, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "eII_exp" , 'f', (model.Nx-1)*(model.Nz-1), CeII_exp, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "eII_lin" , 'f', (model.Nx-1)*(model.Nz-1), CeII_lin, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "eII_gbs" , 'f', (model.Nx-1)*(model.Nz-1), CeII_gbs, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "d" , 'f', (model.Nx-1)*(model.Nz-1), Cd, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "X" , 'f', (model.Nx-1)*(model.Nz-1), CX, 1 );
    if ( model.free_surf == 1 ) {
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "height" , 'f', (model.Nx), Cheight, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "vx" , 'f', (model.Nx), Ctopovx, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "vz" , 'f', (model.Nx+1), Ctopovz, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "x" , 'f', topo_chain->Nb_part, Cxtopo, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "z" , 'f', topo_chain->Nb_part, Cztopo, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "phase" , 'i', topo_chain->Nb_part, topo_chain->phase, 1 );
    }
    
    if (model.isPl_soft == 1) {
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "friction", 'f', (model.Nx-1)*(model.Nz-1), Cfriction, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "cohesion", 'f', (model.Nx-1)*(model.Nz-1), Ccohesion, 1 );
    }
    
    if (model.fstrain == 1) {
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "Fxx", 'f', (model.Nx-1)*(model.Nz-1), CFxx, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "Fxz", 'f', (model.Nx-1)*(model.Nz-1), CFxz, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "Fzx", 'f', (model.Nx-1)*(model.Nz-1), CFzx, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "Fzz", 'f', (model.Nx-1)*(model.Nz-1), CFzz, 1 );
    }
    
    if (model.aniso == 1) {
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "nx", 'f', (model.Nx-1)*(model.Nz-1), Cnx, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "nz", 'f', (model.Nx-1)*(model.Nz-1), Cnz, 1 );
    }
    
    if (model.rec_T_P_x_z == 1) {
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "T0", 'f', (model.Nx-1)*(model.Nz-1), CT0, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "P0", 'f', (model.Nx-1)*(model.Nz-1), CP0, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "Tmax", 'f', (model.Nx-1)*(model.Nz-1), CTmax, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "Pmax", 'f', (model.Nx-1)*(model.Nz-1), CPmax, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "x0", 'f', (model.Nx-1)*(model.Nz-1), Cx0, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Centers" , "z0", 'f', (model.Nx-1)*(model.Nz-1), Cz0, 1 );
    }
    
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "BCp" , 'c', (model.Nx-1)*(model.Nz-1), mesh->BCp.type, 1 );

    
    // Freedom
    free( name );
    //------------//
    DoodzFree( Crho_s );
    DoodzFree( Crho_n );
    DoodzFree( Ceta_s );
    DoodzFree( Ceta_n );
    DoodzFree( CVx );
    DoodzFree( CVz );
    DoodzFree( CP );
    DoodzFree( Csxxd );
    DoodzFree( Csxz );
    DoodzFree( Cexxd );
    DoodzFree( Cexz );
    DoodzFree( CeII_el  );
    DoodzFree( CeII_pl  );
    DoodzFree( CeII_pwl );
    DoodzFree( CeII_exp );
    DoodzFree( CeII_lin );
    DoodzFree( CeII_gbs );
    
    DoodzFree( Cxg_coord );
    DoodzFree( Czg_coord );
    DoodzFree( Cxc_coord );
    DoodzFree( Czc_coord );
    DoodzFree( Czvx_coord );
    DoodzFree( Cxvz_coord );
	DoodzFree( Cxviz );
	DoodzFree( Czviz );
    DoodzFree( Cxviz_hr );
    DoodzFree( Czviz_hr );
    
    DoodzFree( strain     );
    DoodzFree( strain_el  );
    DoodzFree( strain_pl  );
    DoodzFree( strain_pwl );
    DoodzFree( strain_exp );
    DoodzFree( strain_lin );
    DoodzFree( strain_gbs );
    DoodzFree( X );
    DoodzFree( Cstrain     );
    DoodzFree( Cstrain_el  );
    DoodzFree( Cstrain_pl  );
    DoodzFree( Cstrain_pwl );
    DoodzFree( Cstrain_exp );
    DoodzFree( Cstrain_lin );
    DoodzFree( Cstrain_gbs );
    DoodzFree( CX );
    DoodzFree( CT );
    DoodzFree( Cd );
    
    if ( model.free_surf == 1 ) {
        DoodzFree( Cxtopo );
        DoodzFree( Cztopo );
        DoodzFree( Cheight );
        DoodzFree( Ctopovx );
        DoodzFree( Ctopovz );
    }
    
    if (model.isPl_soft == 1) {
        DoodzFree( Cfriction );
        DoodzFree( Ccohesion );
    }
    
    DoodzFree( compo    );
    DoodzFree( compo_hr );
    
    if ( model.fstrain == 1 ) {
        DoodzFree( Fxx  );
        DoodzFree( Fxz  );
        DoodzFree( Fzx  );
        DoodzFree( Fzz  );
        DoodzFree( CFxx );
        DoodzFree( CFxz );
        DoodzFree( CFzx );
        DoodzFree( CFzz );
    }
    
     if ( model.aniso == 1 ) {
         DoodzFree( nx  );
         DoodzFree( nz  );
         DoodzFree( Cnx );
         DoodzFree( Cnz );
     }
    
    if ( model.rec_T_P_x_z == 1 ) {
        DoodzFree( T0 );
        DoodzFree( P0 );
        DoodzFree( Tmax );
        DoodzFree( Pmax );
        DoodzFree( x0 );
        DoodzFree( z0 );
        DoodzFree( CT0 );
        DoodzFree( CP0 );
        DoodzFree( CTmax );
        DoodzFree( CPmax );
        DoodzFree( Cx0 );
        DoodzFree( Cz0 );
    }
    
    DoodzFree(P_total);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void WriteOutputHDF5Particles( grid *mesh, markers *particles, surface *topo, markers* topo_chain, surface *topo_ini, markers* topo_chain_ini, params model, char *txtout, mat_prop materials, scale scaling ) {
    
    char *name;
    double A[8];
    char  *part_ph, *part_gen;
    float *part_x, *part_z, *part_Vx, *part_Vz, *part_T, *part_P, *part_sxxd, *part_sxz;
    int d_part, ind=0, Nb_part_viz=particles->Nb_part;
    float *Cxtopo, *Cztopo, *Cvxtopo, *Cvztopo, *Cheight, *Cvxsurf, *Cvzsurf;
    float *Cxtopo_ini, *Cztopo_ini, *Cvxtopo_ini, *Cvztopo_ini, *Cheight_ini, *Cvxsurf_ini, *Cvzsurf_ini;
    int k;
    
    // Only save a give number of particles
    if(particles->Nb_part<Nb_part_viz) {
        Nb_part_viz = particles->Nb_part;
    }
    d_part    = particles->Nb_part / (Nb_part_viz);
    d_part      = particles->Nb_part / (Nb_part_viz);
    part_x      = DoodzMalloc( sizeof(float)*Nb_part_viz );
    part_z      = DoodzMalloc( sizeof(float)*Nb_part_viz );
    part_Vx     = DoodzMalloc( sizeof(float)*Nb_part_viz );
    part_Vz     = DoodzMalloc( sizeof(float)*Nb_part_viz );
    part_ph     = DoodzMalloc( sizeof(char) *Nb_part_viz );
    part_gen    = DoodzMalloc( sizeof(char) *Nb_part_viz );
    part_T      = DoodzMalloc( sizeof(float) *Nb_part_viz );
    part_P      = DoodzMalloc( sizeof(float) *Nb_part_viz );
    part_sxxd   = DoodzMalloc( sizeof(float) *Nb_part_viz );
    part_sxz    = DoodzMalloc( sizeof(float) *Nb_part_viz );
    
    for (k=0; k<Nb_part_viz; k++) {
        part_x[k]     = (float)particles->x[ind];
        part_z[k]     = (float)particles->z[ind];
        part_Vx[k]    = (float)particles->Vx[ind];
        part_Vz[k]    = (float)particles->Vz[ind];
        part_P[k]     = (float)particles->P[ind];
        part_T[k]     = (float)particles->T[ind];
        part_sxxd[k]  = (float)particles->sxxd[ind];
        part_sxz[k]   = (float)particles->sxz[ind];
        part_ph[k]    = (char)particles->phase[ind];
        part_gen[k]   = (char)particles->generation[ind];
        ind += d_part;
    }
    
    ScaleBack( part_x,     scaling.L, Nb_part_viz );
    ScaleBack( part_z,     scaling.L, Nb_part_viz );
    ScaleBack( part_Vx,    scaling.V, Nb_part_viz );
    ScaleBack( part_Vz,    scaling.V, Nb_part_viz );
    ScaleBack( part_T,     scaling.T, Nb_part_viz );
    ScaleBack( part_P,     scaling.S, Nb_part_viz );
    ScaleBack( part_sxxd,  scaling.S, Nb_part_viz );
    ScaleBack( part_sxz,   scaling.S, Nb_part_viz );
     
    // Topography
    if ( model.free_surf == 1 ) {
        
        // Real topo
        Cxtopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->x, Cxtopo, topo_chain->Nb_part );
        ScaleBack( Cxtopo, scaling.L, topo_chain->Nb_part );
        
        Cztopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->z, Cztopo, topo_chain->Nb_part );
        ScaleBack( Cztopo, scaling.L, topo_chain->Nb_part );
        
        Cheight = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo->height, Cheight, model.Nx );
        ScaleBack( Cheight, scaling.L, model.Nx );
        
        Cvxtopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->Vx, Cvxtopo, topo_chain->Nb_part );
        ScaleBack( Cvxtopo, scaling.V, topo_chain->Nb_part );
        
        Cvztopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->Vz, Cvztopo, topo_chain->Nb_part );
        ScaleBack( Cvztopo, scaling.V, topo_chain->Nb_part );
        
        Cvxsurf = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo->vx, Cvxsurf, model.Nx );
        ScaleBack( Cvxsurf, scaling.V, model.Nx );
        
        Cvzsurf = DoodzMalloc( sizeof(float)*(model.Nx+1));
        DoubleToFloat( topo->vz, Cvzsurf, model.Nx+1);
        ScaleBack( Cvzsurf, scaling.V, model.Nx+1 );
        
        // Advected initial topo
        Cxtopo_ini = DoodzMalloc( sizeof(float)*topo_chain_ini->Nb_part);
        DoubleToFloat( topo_chain_ini->x, Cxtopo_ini, topo_chain_ini->Nb_part );
        ScaleBack( Cxtopo_ini, scaling.L, topo_chain_ini->Nb_part );
        
        Cztopo_ini = DoodzMalloc( sizeof(float)*topo_chain_ini->Nb_part);
        DoubleToFloat( topo_chain_ini->z, Cztopo_ini, topo_chain_ini->Nb_part );
        ScaleBack( Cztopo_ini, scaling.L, topo_chain_ini->Nb_part );
        
        Cheight_ini = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo_ini->height, Cheight_ini, model.Nx );
        ScaleBack( Cheight_ini, scaling.L, model.Nx );
        
        Cvxtopo_ini = DoodzMalloc( sizeof(float)*topo_chain_ini->Nb_part);
        DoubleToFloat( topo_chain_ini->Vx, Cvxtopo_ini, topo_chain_ini->Nb_part );
        ScaleBack( Cvxtopo_ini, scaling.V, topo_chain_ini->Nb_part );
        
        Cvztopo_ini = DoodzMalloc( sizeof(float)*topo_chain_ini->Nb_part);
        DoubleToFloat( topo_chain_ini->Vz, Cvztopo_ini, topo_chain_ini->Nb_part );
        ScaleBack( Cvztopo_ini, scaling.V, topo_chain_ini->Nb_part );
        
        Cvxsurf_ini = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo_ini->vx, Cvxsurf_ini, model.Nx );
        ScaleBack( Cvxsurf_ini, scaling.V, model.Nx );
        
        Cvzsurf_ini = DoodzMalloc( sizeof(float)*(model.Nx+1));
        DoubleToFloat( topo_ini->vz, Cvzsurf_ini, model.Nx+1);
        ScaleBack( Cvzsurf_ini, scaling.V, model.Nx+1 );
        
    }
    
    // Generate file name
    asprintf( &name, "%s%05d%s",txtout, model.step, ".gzip.h5");
    create_output_hdf5( name );
    
    // Add groups
    AddGroup_to_hdf5( name, "Model" );
    AddGroup_to_hdf5( name, "Vertices" );
    AddGroup_to_hdf5( name, "Centers" );
    AddGroup_to_hdf5( name, "VxNodes" );
    AddGroup_to_hdf5( name, "VzNodes" );
    AddGroup_to_hdf5( name, "Particles" );
    AddGroup_to_hdf5( name, "VizGrid" );
    AddGroup_to_hdf5( name, "Topo" );
    AddGroup_to_hdf5( name, "Topo_ini" );
    
    // Model Parameters
    A[0] = (double)(model.time) * scaling.t;
    A[1] = (double)(model.xmax - model.xmin) * scaling.L;
    A[2] = (double)(model.zmax - model.zmin) * scaling.L;
    A[3] = (double)model.Nx;
    A[4] = (double)model.Nz;
    A[5] = (double)model.dx * scaling.L;
    A[6] = (double)model.dz * scaling.L;
    A[7] = (double)model.dt * scaling.t;
    
    // Parameter array
    AddFieldToGroup_generic( _TRUE_, name, "Model", "Params"   , 'd', 8,  A, 1 );
    
    // Add casted partial particle fields
    AddFieldToGroup_generic( _TRUE_, name, "Particles", "x"    , 'f', Nb_part_viz, part_x,     1 );
    AddFieldToGroup_generic( _TRUE_, name, "Particles", "z"    , 'f', Nb_part_viz, part_z,     1 );
    AddFieldToGroup_generic( _TRUE_, name, "Particles", "phase", 'c', Nb_part_viz, part_ph,    1 );
    
    AddFieldToGroup_generic( _TRUE_, name, "Particles", "generation", 'c', Nb_part_viz, part_gen,    1 );

    
    AddFieldToGroup_generic( _TRUE_, name, "Particles", "T", 'f', Nb_part_viz, part_T,     1 );
    AddFieldToGroup_generic( _TRUE_, name, "Particles", "P", 'f', Nb_part_viz, part_P,    1 );

    
    AddFieldToGroup_generic( _TRUE_, name, "Particles", "Vx", 'f', Nb_part_viz, part_Vx,    1 );
    AddFieldToGroup_generic( _TRUE_, name, "Particles", "Vz", 'f', Nb_part_viz, part_Vz,    1 );
    
    AddFieldToGroup_generic( _TRUE_, name, "Particles", "sxxd", 'f', Nb_part_viz, part_sxxd,    1 );
    AddFieldToGroup_generic( _TRUE_, name, "Particles", "sxz",  'f', Nb_part_viz, part_sxz,    1 );
    
    if ( model.free_surf == 1 ) {
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "height" , 'f', (model.Nx), Cheight, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "vxsurf" , 'f', (model.Nx), Cvxsurf, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "vzsurf" , 'f', (model.Nx+1), Cvzsurf, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "x" , 'f', topo_chain->Nb_part, Cxtopo, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "z" , 'f', topo_chain->Nb_part, Cztopo, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "vx" , 'f', topo_chain->Nb_part, Cvxtopo, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "vz" , 'f', topo_chain->Nb_part, Cvztopo, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo", "phase" , 'i', topo_chain->Nb_part, topo_chain->phase, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo_ini", "height" , 'f', (model.Nx), Cheight_ini, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo_ini", "vxsurf" , 'f', (model.Nx), Cvxsurf_ini, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo_ini", "vzsurf" , 'f', (model.Nx+1), Cvzsurf_ini, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo_ini", "x" , 'f', topo_chain_ini->Nb_part, Cxtopo_ini, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo_ini", "z" , 'f', topo_chain_ini->Nb_part, Cztopo_ini, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo_ini", "vx" , 'f', topo_chain_ini->Nb_part, Cvxtopo_ini, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo_ini", "vz" , 'f', topo_chain_ini->Nb_part, Cvztopo_ini, 1 );
        AddFieldToGroup_generic( _TRUE_, name, "Topo_ini", "phase" , 'i', topo_chain_ini->Nb_part, topo_chain_ini->phase, 1 );
    }
    
    // Freedom
    free( name );
    //------------//
    DoodzFree( part_x );
    DoodzFree( part_z );
    DoodzFree( part_Vx );
    DoodzFree( part_Vz );
    DoodzFree( part_ph );
    DoodzFree( part_gen );
    DoodzFree( part_T );
    DoodzFree( part_P );
    DoodzFree( part_sxxd );
    DoodzFree( part_sxz  );
    //------------//
    
    if ( model.free_surf == 1 ) {
        DoodzFree( Cxtopo );
        DoodzFree( Cztopo );
        DoodzFree( Cvxtopo );
        DoodzFree( Cvxsurf );
        DoodzFree( Cvzsurf );
        DoodzFree( Cvztopo );
        DoodzFree( Cheight );
        
        DoodzFree( Cxtopo_ini );
        DoodzFree( Cztopo_ini );
        DoodzFree( Cvxtopo_ini );
        DoodzFree( Cvxsurf_ini );
        DoodzFree( Cvzsurf_ini );
        DoodzFree( Cvztopo_ini );
        DoodzFree( Cheight_ini );
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void WriteResiduals( grid Mmesh, params model, Nparams Nmodel, scale scaling ) {
    
    char *name;
    double A[8];
    
    asprintf( &name, "Residuals%05d%s", Nmodel.nit, ".gzip.h5");
    create_output_hdf5( name );
    
    
    // Model Parameters
    A[0] = (double)(model.time) * scaling.t;
    A[1] = (double)(model.xmax - model.xmin) * scaling.L;
    A[2] = (double)(model.zmax - model.zmin) * scaling.L;
    A[3] = (double)model.Nx;
    A[4] = (double)model.Nz;
    A[5] = (double)model.dx * scaling.L;
    A[6] = (double)model.dz * scaling.L;
    A[7] = (double)model.dt * scaling.t;
    
    // Add groups
    AddGroup_to_hdf5( name, "Model" );
    AddGroup_to_hdf5( name, "Vertices" );
    AddGroup_to_hdf5( name, "Centers" );
    AddGroup_to_hdf5( name, "VxNodes" );
    AddGroup_to_hdf5( name, "VzNodes" );
    AddGroup_to_hdf5( name, "Particles" );
    AddGroup_to_hdf5( name, "VizGrid" );
    AddGroup_to_hdf5( name, "Topo" );
    
    // Scaling
    ArrayTimesScalar( Mmesh.eta_phys_n, scaling.eta, (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.rho_app_n,  scaling.rho, (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.eta_phys_s, scaling.eta, (model.Nx)*(model.Nz) );
    ArrayTimesScalar( Mmesh.rho_app_s,  scaling.rho, (model.Nx)*(model.Nz) );
    ArrayTimesScalar( Mmesh.rp,      scaling.E,   (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.ru,      scaling.F,   (model.Nx)*(model.Nz+1) );
    ArrayTimesScalar( Mmesh.rv,      scaling.F,   (model.Nx+1)*(model.Nz) );
    
    // Parameter array
    AddFieldToGroup_generic( _TRUE_, name, "Model", "Params"   , 'd', 8,  A, 1 );
    
    AddFieldToGroup_generic( _TRUE_, name, "VxNodes" , "ru"  , 'd', model.Nx*(model.Nz+1),     Mmesh.ru, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "VzNodes" , "rv"  , 'd', (model.Nx+1)*model.Nz,     Mmesh.rv, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "rp"  , 'd', (model.Nx-1)*(model.Nz-1), Mmesh.rp, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "rho"  , 'd', (model.Nx-1)*(model.Nz-1), Mmesh.rho_app_n, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Vertices" , "rho"  , 'd', (model.Nx)*(model.Nz), Mmesh.rho_app_s, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Centers" , "eta"  , 'd', (model.Nx-1)*(model.Nz-1), Mmesh.eta_phys_n, 1 );
    AddFieldToGroup_generic( _TRUE_, name, "Vertices" , "eta"  , 'd', (model.Nx)*(model.Nz), Mmesh.eta_phys_s, 1 );
    
    // Scaling
    ArrayTimesScalar( Mmesh.eta_phys_n, 1.0/scaling.eta, (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.rho_app_n,  1.0/scaling.rho, (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.eta_phys_s, 1.0/scaling.eta, (model.Nx)*(model.Nz) );
    ArrayTimesScalar( Mmesh.rho_app_s,  1.0/scaling.rho, (model.Nx)*(model.Nz) );
    ArrayTimesScalar( Mmesh.rp,      1.0/scaling.E,   (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.ru,      1.0/scaling.F,   (model.Nx)*(model.Nz+1) );
    ArrayTimesScalar( Mmesh.rv,      1.0/scaling.F,   (model.Nx+1)*(model.Nz) );
    
    free( name );
    
}