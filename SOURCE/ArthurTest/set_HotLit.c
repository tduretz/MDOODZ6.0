#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "head.h"
#include "time.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// SIMPLE
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, Mgrid Mmesh, scale scaling ) {
    
    int k;
    double TopoLevel = -0.0e3/scaling.L; // sets zero initial topography

    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        topo_chain->z[k]     = TopoLevel;
        topo_chain->phase[k] = 0;
    }
    
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {
    
    int np;
    
    // Define dimensions
    double Lx      = (double) (model.xmax - model.xmin) ;
    double Lz      = (double) (model.zmax - model.zmin) ;
    double gsbg    = 2.0e-3/scaling.L;                            // reference grain size
    double HLit    = 98.e3/scaling.L;                           // lithosphere thickness
    double HCrust  = 30.e3/scaling.L;                            // crust thickness
    double Huc     = 2.5e3/scaling.L;                             // Upper crust
    double Hlc     = HCrust - Huc;                            // Lower crust
    
//    double HMtlSup = 40e3/scaling.L;                           // mantle sup
    double Tsurf   = zeroC/scaling.T, Tpart;                   // surface temperature
    double Tmant   = (1330.0+zeroC)/scaling.T;                   // adiabatic mantle temperature
    double rad=model.user0/scaling.L, la= 1.0*rad, sa = 1.0*rad, theta=(0.0)*M_PI/180;             // Dimensions for granite seed
    printf("Rad: \t %.10f\n", rad);
//    double xc = 0.0/scaling.L, zc = -15.0e3/scaling.L;
    double spacing = 2.0e3/scaling.L;
    double zlc     = 14.0e3/scaling.L;
    printf("Model paramters: \t %.4f \t %.4f \t %.4f \t %.4f\n", model.xmax, model.xmin, model.zmax, model.zmin);
    printf("Model extents: \t %.4f \t %.4f\n", Lx, Lz);
    
    //
    //--------------------------------------------------//
    // lire un fichier de T:
    //--------------------------------------------------//
    FILE *read;
    
    int nlayers = ( model.user0 - 1)/2-1, il;
    int s1, s2, nb_pt_prof,iprof;
    double *zprof, *Tprof;
    double Tprofmax, zprofmax;
    double Tco, Toc, distRatio;
    
    s1 = sizeof(int);
    s2 = sizeof(double);
    
    
    // Read the Doudou file
    model.input_file = "ProfCont_30_Tm750.bin";
    if (fopen(model.input_file, "rb")!=NULL){
        read = fopen(model.input_file, "rb");
        //printf("Loading %d particles from file %s...\n", particles->Nb_part, model.input_file );
    }
    else {
        printf("Cannot open file %s, check if the file exists in the current location !\n Exiting", model.input_file);
        exit(1);
    }
    
    fread(&nb_pt_prof, s1, 1, read);
    
    zprof = DoodzMalloc(sizeof(double)*nb_pt_prof);
    Tprof = DoodzMalloc(sizeof(double)*nb_pt_prof);
    
    fread( Tprof,     s2, nb_pt_prof, read);
    fread( zprof,     s2, nb_pt_prof, read);
    
    fclose(read);
    
    printf("nb pt in profile= %d : \n", nb_pt_prof);
    
    MinMaxArray(zprof, 1, nb_pt_prof, "zprof" );
    MinMaxArray(Tprof, 1, nb_pt_prof, "Tprof" );
    
    zprofmax = zprof[nb_pt_prof-1];
    Tprofmax = Tprof[nb_pt_prof-1];
    
    printf("Tmax: %f  /  profmax: %f \n", Tprofmax, zprofmax);
    
    // scale it
    for( iprof=0; iprof<nb_pt_prof; iprof++ ) {
        //printf("T[C]: %f  /  depth[m]: %f \n", Tprof[iprof], zprof[iprof]);
        Tprof[iprof] = (Tprof[iprof]+zeroC)/scaling.T;
        zprof[iprof] = -zprof[iprof]/scaling.L;
    }
    
    zprofmax = zprof[nb_pt_prof-1];
    Tprofmax = Tprof[nb_pt_prof-1];
    
    printf("Tmax: %f  /  profmax: %f \n", Tprofmax, zprofmax);
    
    
    DoodzFree(zprof);
    DoodzFree(Tprof);
    //--------------------------------------------------//
    
// JPoh Addition:===
    //* Perturbation switches *//
    int rand_pert = 0;
    int layering = 0;
    int passive_marker = 1;
    int weak_seed = 1;
    
    //* Creation of crustal seed markers *//
    int seed_number = 15;                               // number of circles
    double zc_1 = -8.0e3/scaling.L;
    double zc_2 = -20.0e3/scaling.L;
    double xc_freq = Lx / seed_number;
    double X, Z, Xn, Zn, X2, Z2, Xn2, Zn2;
    int i;
    
    /* Crustal weak seed */
    double seed_x = -10.0e3/scaling.L;
    double seed_z = -8.5e3/scaling.L;
    double seed_X, seed_Z, seed_Xn, seed_Zn;
    
    // Vectors to store values
    double xc_Layer1[seed_number];
    double xc_Layer2[seed_number];
    double zc_Layer1[seed_number];
    double zc_Layer2[seed_number];
    
    srand((unsigned) time(NULL));
    for ( i = 0; i <= (seed_number); i++) {
        xc_Layer1[i] = (model.xmin + 5.e3/scaling.L) + (i * xc_freq);
        xc_Layer2[i] = (model.xmin + 10.e3/scaling.L) + (i * xc_freq);
            
        zc_Layer1[i] = zc_1;
        zc_Layer2[i] = zc_2;
            
//            printf("%.2f \t %.2f \n", xc_Layer1[i], zc_Layer1[i] );
//            printf("%.2f \t %.2f \n", xc_Layer2[i], zc_Layer2[i] );
    }
    

    //* Creating noise layer at the Moho *//
    double hi_num = 2.0;
    double low_num = 1.0;
    double line_z = HCrust;                            // Change value here to place your noise layer
    double dx_line = Lx/(model.Nx - 1);
    double rand_array[model.Nx];
    double rand_array2[model.Nx];
    double rand_line[model.Nx];
    double m, c, y, y1, y2, x1, x2;
    double m_uc, c_uc, y_uc, y1_uc, y2_uc, x1_uc, x2_uc;
    
    srand((unsigned) time(NULL));
    for ( i = 0; i <= (model.Nx-1) ; i++) {
        // Generate noise layer for crust
        rand_array[i] = ((((double)rand() / (RAND_MAX)) * hi_num - low_num)) * (1.e3/scaling.L);
        rand_array[i] *= 2.;                           // Amplitude of error
        rand_array[i] += line_z;
        
        // Noise layer for 2nd random line
        rand_array2[i] = ((((double)rand() / (RAND_MAX)) * hi_num - low_num)) * (1.e3/scaling.L);
        rand_array2[i] *= 1.;                           // Amplitude of error
        rand_array2[i] += zc_2;
        //rand_array[i] = rand_array[i]/scaling.L;
        rand_line[i] = model.xmin;
        if (i > 0) {
            rand_line[i] = rand_line[i-1] + dx_line;
            //            temp_rand_line[i] = temp_rand_line[i-1] + dx_line;
        }
        //        printf("%.2f \t %.2f \n", rand_array[i], rand_line[i]); // Prints out elements in error array
        //        printf("%.4f \n", temp_rand_array[i]);
    }
    
    //--------------------------------------------------//

    
    // Loop on particles
    for( np=0; np<particles->Nb_part; np++ ) {
        
        // Standart initialisation of particles
        particles->Vx[np]    = -particles->x[np]*model.EpsBG; // set initial particle velocity (unused)
        particles->Vz[np]    =  particles->z[np]*model.EpsBG; // set initial particle velocity (unused)
        particles->phase[np] = 0;                             // same phase number everywhere
        particles->d[np]     = gsbg;                          // same grain size everywhere
        particles->phi[np]   = 0.0;                           // zero porosity everywhere
        
//        X = particles->x[np]-xc;
//        Z = particles->z[np]-zc;
        
        //--------------------------//
        // TEMPERATURE - linear geotherm
//        Tpart = ((Tmant-Tsurf)/HLit)*(-particles->z[np]) + Tsurf;
//        if (Tpart>Tmant) Tpart = Tmant;
//        particles->T[np]       = Tpart;

        Tpart = Tprofmax;
        
        if (particles->z[np] > zprofmax) {
            for( iprof=0; iprof<nb_pt_prof-1; iprof++ ) {
                if(particles->z[np] > zprof[iprof] && particles->z[np]<zprof[iprof-1]){
                    Tpart = Tprof[iprof];
                }
            }
        }
        
        particles->T[np]  = Tpart;

        //--------------------------//
//        // Phases - lithosphere-asthenosphere
//         if (particles->z[np] <-(HCrust))  particles->phase[np] = 1;
//         if (particles->z[np] <-(HLit)  )  particles->phase[np] = 2;
       
         //printf("%2.2e\n", rad*scaling.L);
         //if (pow(particles->z[np]+30e3/scaling.L,2) + pow(particles->x[np]-0.0e3/scaling.L,2) < rad*rad )  particles->phase[np] = 0 ;
        
// POTATO_SARAH: ==
        if (particles->z[np]<-(HLit)) particles->phase[np] = 3;
        if (particles->z[np]>-(HLit)) particles->phase[np] = 2;
        if (particles->z[np]>-(HCrust)) particles->phase[np] = 1;
        
        /* Initiate noise layer between Moho and crust*/
        if (rand_pert == 1) {
            for (i = 0; i <= (model.Nx-1); i++) {
                if (particles->x[np]> rand_line[i] && particles->x[np]<= rand_line[i+1]){
                    //printf("x_part: %.2f; xleft = %.2f; xright = %.2f\n", particles->x[np], rand_line[i], rand_line[i+1]); // Prints out elements in error array
                    y1 = -(rand_array[i]);
                    y2 = -(rand_array[i+1]);
                    x1 = rand_line[i];
                    x2 = rand_line[i+1];
                    
                    m = (y1 - y2) / (x1 - x2);               // Calculate gradient of line
                    c = y1 - (m * x1);                       // Calculation of the y-intercept
                    y = m * particles->x[np] + c;            // Determine where the marker would be on the line
                    
                    if (particles->z[np]>=y){
                        //printf("z_part: %.2f; y = %.2f ; I''m ABOVE \n", particles->z[np], y); // Prints out elements in error array
                        particles->phase[np] = 1;
                    }
                    
                    if (particles->z[np]<y && particles->z[np]>=-(HCrust)){
                        //printf("z_part: %.2f; y = %.2f ; I''m BELOW \n", particles->z[np], y); // Prints out elements in error array
                        particles->phase[np] = 2;
                    }
                }
            }
        }
        
//        // Add distinction between aged craton and juvenile
//        if (particles->z[np]>-(Huc) && particles->x[np]<(0.0e3/scaling.L)) particles->phase[np] = 0;
//        if (particles->z[np]>-(Huc) && particles->x[np]>(0.0e3/scaling.L)) particles->phase[np] = 2;
        
        // Add upper crust layer
        if (particles->z[np]>-(Huc)) particles->phase[np] = 0;
        
//        // Crustal layering:
//        if (layering == 1) {
//            for( il=0; il<= 50; il=il+2 ) {
//                if (particles->z[np]>-(HCrust) && particles->z[np]>=-((il+1)*spacing) && particles->z[np]<-(il*spacing)) particles->phase[np] = 1;
//            }
//        }
        
        // A single passive marker
//        if (pow(particles->z[np]+ (-15.0e3/scaling.L),2) + pow(particles->x[np]- (-10.0e3/scaling.L),2) <= 0.25*(rad*rad) )  particles->phase[np] = 4 ;
//        X = particles->x[np]-(10.e3/scaling.L);
//        Z = particles->z[np]-(15.0e3/scaling.L);
//        Xn = X*cos(theta) + Z*sin(theta);
//        Zn = X*sin(theta) - Z*cos(theta);
//        if ( pow(Xn/la,2) + pow(Zn/sa,2) - 1 < 0 ) particles->phase[np] = 4;
//        if (particles->z[np]>=-(10.0e3/scaling.L) && particles->z[np]<=-(5.0e3/scaling.L)) {
//            if (particles->x[np]>= (-25.0e3/scaling.L) && particles->x[np]<=(25.0e3/scaling.L)) {
//                particles->phase[np] = 4;
//            }
//        }
        
        // Passive markers
        if (passive_marker == 1) {
//            for ( i = 0; i < (seed_number); i ++) {
//            
//                /* First row of markers */
//                X = particles->x[np]-xc_Layer1[i];
//                Z = particles->z[np]-zc_Layer1[i];
//                Xn = X*cos(theta) + Z*sin(theta);
//                Zn = X*sin(theta) - Z*cos(theta);
//                
//                if ( pow(Xn/rad,2.0) + pow(Zn/rad,2.0) - 1.0 < 0.0 ) {
////                    printf("Putting marker at X = %.4g \t Z = %.4g\n", Xn, Zn);
//                    particles->phase[np] = 4;
//                }
//                
////                if ( pow(Xn/rad,2) + pow(Zn/rad,2) - 1 < 0 && particles->x[np]>(0.0e3/scaling.L)) {
//////                    printf("Putting marker at X = %.4g \t Z = %.4g\n", Xn, Zn);
////                    particles->phase[np] = 7;
////                }
//
//                /* Second row of markers */
//                X2 = particles->x[np]-xc_Layer2[i];
//                Z2 = particles->z[np]-zc_Layer2[i];
//                Xn2 = X2*cos(theta) + Z2*sin(theta);
//                Zn2 = X2*sin(theta) - Z2*cos(theta);
//                
//                if ( pow(Xn2/rad,2.0) + pow(Zn2/rad,2.0) - 1.0 < 0.0 ) {
////                    printf("Putting marker at X = %.4g \t Z = %.4g\n", Xn2, Zn2);
//                    particles->phase[np] = 4;
//                }
////                if ( pow(Xn2/rad,2) + pow(Zn2/rad,2) - 1 < 0 && particles->x[np]>(0.0e3/scaling.L)) {
////                    //                    printf("Putting marker at X = %.4g \t Z = %.4g\n", Xn2, Zn2);
////                    particles->phase[np] = 7;
////                }
//            }
            if (particles->phase[np] == 1) {
            for( il=0; il<= 50; il=il+2 ) {
                if (particles->z[np]>-(HCrust) && particles->z[np]>=-((il+1)*spacing) && particles->z[np]<-(il*spacing)) particles->phase[np] = 4;
            }
            }
        }
        
        // Crustal weak seed
        if (weak_seed == 1) {
            seed_X = particles->x[np]-seed_x;
            seed_Z = particles->z[np]-seed_z;
            seed_Xn = seed_X*cos(theta) + seed_Z*sin(theta);
            seed_Zn = seed_X*sin(theta) - seed_Z*cos(theta);
            
            if ( pow(seed_Xn/rad,2.0) + pow(seed_Zn/rad,2.0) - 1.0 < 0.0) {
//                printf("Putting marker at X = %.4g \t Z = %.4g\n", seed_Xn, seed_Zn);
                particles->phase[np] = 5;
            }
        }
// POTATO_SARAH: ==


        
        //--------------------------//
        // DENSITY
        if ( model.eqn_state > 0 ) {
            particles->rho[np] = materials->rho[particles->phase[np]] * (1.0 -  materials->alp[particles->phase[np]] * (Tpart - materials->T0[particles->phase[np]]) );
        }
        else {
            particles->rho[np] = materials->rho[particles->phase[np]];
        }
        
        //--------------------------//
        // SANITY CHECK
        if (particles->phase[np] > model.Nb_phases) {
            printf("Lazy bastard! Fix your particle phase ID! \n");
            exit(144);
        }
        //--------------------------//
        
        //particles->phase[np] = 0;
        
    }
    
//    // Loop on particles
//    for( np=0; np<particles->Nb_part; np++ ) {
//        
//        // Standart initialisation of particles
//        particles->Vx[np]    = -particles->x[np]*model.EpsBG; // set initial particle velocity (unused)
//        particles->Vz[np]    =  particles->z[np]*model.EpsBG; // set initial particle velocity (unused)
//        particles->phase[np] = 2;                             // same phase number everywhere
//        particles->d[np]     = gsbg;                          // same grain size everywhere
//        particles->phi[np]   = 0.0;                           // zero porosity everywhere
//        
//        //--------------------------//
//        // TEMPERATURE - linear geotherm
//        Tpart = Tmant;//((Tmant-Tsurf)/HLit)*(-particles->z[np]) + Tsurf;
//        if (Tpart>Tmant) Tpart = Tmant;
//        if (pow(particles->z[np]+45e3/scaling.L,2) + pow(particles->x[np]-0.0e3/scaling.L,2) < rad*rad )  Tpart = 2*Tpart;
//
//        particles->T[np]       = Tpart;
//        
//
//        
//        //--------------------------//
////        // Phases - lithosphere-asthenosphere
////        if (particles->z[np] <-(HCrust))  particles->phase[np] = 1;
////        if (particles->z[np] <-(HLit)  )  particles->phase[np] = 2;
////        //        printf("%2.2e\n", rad*scaling.L);
////        if (pow(particles->z[np]+45e3/scaling.L,2) + pow(particles->x[np]-0.0e3/scaling.L,2) < rad*rad )  {printf("sada \n");particles->phase[np] = 0;};
//        
//        //--------------------------//
//        // DENSITY
//        if ( model.eqn_state > 0 ) {
//            particles->rho[np] = materials->rho[particles->phase[np]] * (1 -  materials->alp[particles->phase[np]] * (Tpart - materials->T0[particles->phase[np]]) );
//        }
//        else {
//            particles->rho[np] = materials->rho[particles->phase[np]];
//        }
//        
//        //--------------------------//
//        // SANITY CHECK
//        if (particles->phase[np] > model.Nb_phases) {
//            printf("Lazy bastard! Fix your particle phase ID! \n");
//            exit(144);
//        }
//    }

    MinMaxArray(particles->Vx, scaling.V, particles->Nb_part, "Vxp init" );
    MinMaxArray(particles->Vz, scaling.V, particles->Nb_part, "Vzp init" );
    MinMaxArray(particles->T, scaling.T, particles->Nb_part,  "Tp init" );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void MSetFields( Mgrid *mesh, params *model, Mparams *Mmodel, scale scaling, markers* particles, mat_prop *materials ) {
    
    int   kk, k, l, c, c1, np;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
    double dmin, VzBC, width = 1.0 / scaling.L, eta = 1.0e4 / scaling.eta ;
    double Lx, Lz, T1, T2, rate=model->EpsBG,  z_comp=-140.0e3/scaling.L;
    double Vx_r, Vx_l, Vz_b, Vz_t, Vx_tot, Vz_tot;
    
    NX  = mesh->Nx[0];
    NZ  = mesh->Nz[0];
    //double Lxinit = 1400e3/scaling.L, ShortSwitchV0 = 0.40;
    //double Vfix = (3.11576/(1000.0*365.25*24.0*3600.0))/(scaling.L/scaling.t); // [3.11576 == 0.315576 cm/yr]
    
    Lx = sqrt((mesh->xg_coord[0][NX-1]-mesh->xg_coord[0][0])*(mesh->xg_coord[0][NX-1]-mesh->xg_coord[0][0]));
    Lz = sqrt((mesh->zg_coord[0][NZ-1]-mesh->zg_coord[0][0])*(mesh->zg_coord[0][NZ-1]-mesh->zg_coord[0][0]));

    double Vfix = 2.0 * rate * Lx;
    
    double Tchange   = (1329.0+zeroC)/scaling.T;                  // phase 4 => phase 5 @ Tchange;
    
    
    // ---- T-Dependent marker types
    // -------------------- SPECIFIC TO YOANN's SETUP -------------------- //
    
    NX  = mesh->Nx[0];
    NZ  = mesh->Nz[0];
    NCX = NX-1;
    NCZ = NZ-1;
    NXVZ = NX+1;
    NZVX = NZ+1;
    
    X  = malloc (NX*sizeof(double));
    Z  = malloc (NZ*sizeof(double));
    XC = malloc (NCX*sizeof(double));
    ZC = malloc (NCZ*sizeof(double));
    
    for (k=0; k<NX; k++) {
        X[k] = mesh->xg_coord[0][k];
    }
    for (k=0; k<NCX; k++) {
        XC[k] = mesh->xc_coord[0][k];
    }
    for (l=0; l<NZ; l++) {
        Z[l] = mesh->zg_coord[0][l];
    }
    for (l=0; l<NCZ; l++) {
        ZC[l] = mesh->zc_coord[0][l];
    }
    
    // Fix VelocityM
//    Lx = sqrt((mesh->xg_coord[0][NX-1]-mesh->xg_coord[0][0])*(mesh->xg_coord[0][NX-1]-mesh->xg_coord[0][0]));
//    Lz = sqrt((mesh->zg_coord[0][NZ-1]-mesh->zg_coord[0][0])*(mesh->zg_coord[0][NZ-1]-mesh->zg_coord[0][0]));

    Vx_tot = Vfix;
    Vx_r = -0.5*Vx_tot;
    Vx_l = 0.5*Vx_tot;
    
    rate = Vx_tot/Lx;
    
    Vz_tot = -rate*Lz;
    Vz_t = 0.5*rate*Lz;
    Vz_b = Vz_tot;   //-0.5*rate*Lz;
    
    double Vz_tib = -Vx_tot * (0.0 - model->zmin) / (model->xmax - model->xmin);
    
    printf("Vz_b = %2.2e Vz_tib = %2.2e\n", Vz_b, Vz_tib);

    
    /* --------------------------------------------------------------------------------------------------------*/
	/* Set the BCs for Vx on all grid levels                                                                   */
	/* Type  0: Dirichlet point that matches the physical boundary (Vx: left/right, Vz: bottom/top)            */
	/* Type 11: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)       */
	/* Type  2: Neumann point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)         */
	/* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)              */
	/* Type -2: periodic in the x direction (matches the physical boundary)                                    */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
	/* --------------------------------------------------------------------------------------------------------*/
	
	for(kk=0; kk<Mmodel->n_level; kk++) {
		NX  = mesh->Nx[kk];
		NZ  = mesh->Nz[kk];
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;
        
		for (l=0; l<mesh->Nz[kk]+1; l++) {
			for (k=0; k<mesh->Nx[kk]; k++) {
				
				c = k + l*(mesh->Nx[kk]);
                
                if ( mesh->BCu.type[kk][c] != 30 ) {
                    
                    // Internal points:  -1
                    mesh->BCu.type[kk][c] = -1;
                    mesh->BCu.val[kk][c]  =  0;
                    
                    // Matching BC nodes WEST
                    if (k==0 ) {
                        mesh->BCu.type[kk][c] = 0;
                        //mesh->BCu.val[kk][c]  = -mesh->xg_coord[0][k] * model->EpsBG;
                        mesh->BCu.val[kk][c]  = Vx_l;
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx[kk]-1 ) {
                        mesh->BCu.type[kk][c] = 0;
                        //mesh->BCu.val[kk][c]  = -mesh->xg_coord[0][k] * model->EpsBG;
                        mesh->BCu.val[kk][c]  = Vx_r;
                    }
                    
                    // Free slip SOUTH
                    if (l==0  ) {
                        mesh->BCu.type[kk][c] = 13;
                        mesh->BCu.val[kk][c]  =  0.0;
                    }
                    
                    // Free slip NORTH
                    if ( l==mesh->Nz[kk] ) {
                        mesh->BCu.type[kk][c] = 13;
                        mesh->BCu.val[kk][c]  =  0.0;
                    }
                }
                
			}
		}
		
	}
	
	/* --------------------------------------------------------------------------------------------------------*/
	/* Set the BCs for Vz on all grid levels                                                                   */
	/* Type  0: Dirichlet point that matches the physical boundary (Vx: left/right, Vz: bottom/top)            */
	/* Type 11: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)       */
	/* Type  2: Neumann point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)         */
	/* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)              */
	/* Type -2: periodic in the x direction (does not match the physical boundary)                             */
	/* Type-10: useless point (set to zero)                                                                    */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
	/* --------------------------------------------------------------------------------------------------------*/
	
	for(kk=0; kk<Mmodel->n_level; kk++) {
		
		NX  = mesh->Nx[kk];
		NZ  = mesh->Nz[kk];
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;
		
		for (l=0; l<mesh->Nz[kk]; l++) {
			for (k=0; k<mesh->Nx[kk]+1; k++) {
				
				c  = k + l*(mesh->Nx[kk]+1);
                
                if ( mesh->BCv.type[kk][c] != 30 ) {
                    
                    // Internal points:  -1
                    mesh->BCv.type[kk][c] = -1;
                    mesh->BCv.val[kk][c]  =  0;
                    
                    // Matching BC nodes SOUTH
                    if (l==0 ) {
                        mesh->BCv.type[kk][c] = 0;
                        //mesh->BCv.val[kk][c]  = mesh->zg_coord[0][l] * model->EpsBG;
                        mesh->BCv.val[kk][c]  = Vz_tib;
                    }
                    
                    // Matching BC nodes NORTH
                    if (l==mesh->Nz[kk]-1 ) {
                        mesh->BCv.type[kk][c] = 0;
                        //mesh->BCv.val[kk][c]  = mesh->zg_coord[0][l] * model->EpsBG;
                        mesh->BCv.val[kk][c]  = 0.0;
                    }
                    
                    // Non-matching boundary WEST
                    if ( (k==0) ) {
                        mesh->BCv.type[kk][c] =   13;
                        mesh->BCv.val[kk][c]  =   0.0;
                    }
                    
                    // Non-matching boundary EAST
                    if ( (k==mesh->Nx[kk]) ) {
                        mesh->BCv.type[kk][c] =   13;
                        mesh->BCv.val[kk][c]  =   0.0;
                    }
                }
				
			}
		}
		
	}
    
    /* --------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for P on all grid levels                                                                    */
    /* Type  0: Dirichlet within the grid                                                                      */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
    /* Type 31: surface pressure (Dirichlet)                                                                   */
    /* --------------------------------------------------------------------------------------------------------*/
    
    for(kk=0; kk<Mmodel->n_level; kk++) {
        
        NX  = mesh->Nx[kk];
        NZ  = mesh->Nz[kk];
        NCX = NX-1;
        NCZ = NZ-1;
        NXVZ = NX+1;
        NZVX = NZ+1;
        
        for (l=0; l<NCZ; l++) {
            for (k=0; k<NCX; k++) {
                
                c  = k + l*(NCX);
                
            
                
                if (mesh->BCt.type[c] != 30) {
                    
                    
                    
                    // Internal points:  -1
                    mesh->BCp.type[kk][c] = -1;
                    mesh->BCp.val[kk][c]  =  0;
                    
                    if ( (k==0 || k==NCX-1) && l==NCZ-1 ) {
                        mesh->BCp.type[kk][c] =  0;
                        mesh->BCp.val[kk][c]  =  0;
                    }
                }
            }
        }
    }
    
	/* -------------------------------------------------------------------------------------------------------*/
	/* Set the BCs for T on all grid levels                                                                   */
	/* Type  1: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)      */
	/* Type  0: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)             */
    /* Type -2: periodic in the x direction (matches the physical boundary)                                   */
	/* Type -1: not a BC point (tag for inner points)                                                         */
    /* Type 30: not calculated (part of the "air")                                                            */
	/* -------------------------------------------------------------------------------------------------------*/
        
    double TN =  zeroC/scaling.T, TS = 1330.0/scaling.T;
    double TW = 1330/scaling.T, TE = 1330/scaling.T;
    double Tbot, Tleft, Tright;
    
	for(kk=0; kk<1; kk++) {
		NX  = mesh->Nx[kk];
		NZ  = mesh->Nz[kk];
		NCX = NX-1;
		NCZ = NZ-1;
		NXVZ = NX+1;
		NZVX = NZ+1;
		
		for (l=0; l<mesh->Nz[kk]-1; l++) {
			for (k=0; k<mesh->Nx[kk]-1; k++) {
				
				c = k + l*(NCX);
                
                if ( mesh->BCt.type[c] != 30 ) {
                    
                    // WEST
                    if ( k==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typW[l] = 0;
                        mesh->BCt.valW[l] = TW;
                    }
                    
                    // EAST
                    if ( k==NCX-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typE[l] = 0;
                        mesh->BCt.valE[l] = TE;
                    }
                    
                    // SOUTH
                    if ( l==0 ) {
                        mesh->BCt.type[c] = 1;
                        mesh->BCt.typS[k] = 1;
                        mesh->BCt.valS[k] = mesh->T[c];
                    }
                    
                    // NORTH
                    if ( l==NCZ-1 ) {
                        mesh->BCt.type[c] = 1;
                        mesh->BCt.typN[k] = 1;
                        mesh->BCt.valN[k] = TN;//mesh->T[c];
                    }
                    
                    // FREE SURFACE
                    else {
                        if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c+NCX] == 30) {
                            mesh->BCt.type[c] = 1;
                            mesh->BCt.val[c]  = TN;
                        }
                    }
                    
                    
                }
                
			}
		}

    }
    
	free(X);
	free(Z);
	free(XC);
	free(ZC);
	printf("Velocity and pressure were initialised\n");
	printf("Boundary conditions were set up\n");
	
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
