#include "cholmod.h"
#include "cs.h"
//---------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ MACROS DEFINITIONS -------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------//
#define MINV(a,b) ((a)<=(b)? (a):(b))
#define MAXV(a,b) ((a)>=(b)? (a):(b))
#define ABSV(a)   ((a)>=0? (a):-(a))
#define _TRUE_  1
#define _FALSE_ 0
#define DoodzFP double
#define zeroC 273.15
#define Rg    8.314510
#define PI    3.14159265359

//---------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------- STRUCTURE DEFINITIONS -----------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------//

// Stucture scale contains scaling parameters
typedef struct _scale scale;
struct _scale {
	double eta, L, V, T, t, a, E, S, m, rho, F, J, W, Cv, rhoE, k;
};

// Structure surface contains free surface data
typedef struct _surface surface;
struct _surface {
    DoodzFP *a, *b, *height, *vx, *vz, *a0, *b0, *height0;
    int *VertInd;
};

// Stucture scale contains Stokes sparse matrix
typedef struct _OutputSparseMatrix OutputSparseMatrix ;
struct _OutputSparseMatrix  {
    double params[4];
    double *V, *eta_cell, *b;
    int *Ic, *J, *eqn_u, *eqn_v, *eqn_p;
};

// mat_prop contains information related to material phase properties
typedef struct _mat_prop mat_prop;
struct _mat_prop {
	int    Nb_phases;
    DoodzFP R, eta_VP;
    DoodzFP eta0[20], rho[20], mu[20], Cv[20], k[20], Qr[20], C[20], phi[20], Slim[20], n[20], A[20], Ea[20], Va[20], alp[20], bet[20], Qm[20], T0[20], P0[20], drho[20], k_eff[20];
    DoodzFP tpwl[20], Qpwl[20], Vpwl[20], npwl[20], mpwl[20], Apwl[20], apwl[20], fpwl[20], rpwl[20], Fpwl[20], pref_pwl[20];
    DoodzFP texp[20], Qexp[20], Vexp[20], Sexp[20], Eexp[20], Gexp[20], aexp[20], fexp[20], rexp[20], qexp[20], nexp[20];
    DoodzFP tlin[20], Qlin[20], Vlin[20], nlin[20], mlin[20], Alin[20], alin[20], flin[20], rlin[20], Flin[20];
    DoodzFP tgbs[20], Qgbs[20], Vgbs[20], ngbs[20], mgbs[20], Agbs[20], agbs[20], fgbs[20], rgbs[20], Fgbs[20];
    DoodzFP ppzm[20], Kpzm[20], Qpzm[20], Gpzm[20], cpzm[20], Lpzm[20], gs_ref[20];
    int     gs[20], cstv[20], pwlv[20], linv[20], expv[20], gbsv[20], phase_diagram[20], density_model[20];
    DoodzFP C_end[20], phi_end[20], pls_start[20], pls_end[20];
};

// markers is the particles structure
typedef struct _p_markers markers;
struct _p_markers {
	int    Nx_part, Nz_part, Nb_part, Nb_part_max, min_part_cell;
	DoodzFP *x, *z, *Vx, *Vz, *P, *sxxd, *szzd, *sxz, *progress, *rho, *om_p, *T, *d, *phi, *X;
    DoodzFP *strain, *strain_el, *strain_pl, *strain_pwl, *strain_exp, *strain_lin, *strain_gbs;
	int    *phase, *generation;
    markers* marker_chain;
    int    *intag;
    double *rhoUe0;
    double *Fxx, *Fxz, *Fzx, *Fzz, *nx, *nz;
    double *T0, *P0, *x0, *z0, *Tmax, *Pmax;
};

// BC is a boundary condition structure for the mechanical solver
typedef struct _BC BC;
struct _BC {
	char   *type;
	double *val;
};

// BC is a boundary condition structure for the thermal solver
typedef struct _BCT BCT;
struct _BCT {
	char   *type, *typW, *typE, *typS, *typN;
	double *val, *valW, *valE, *valS, *valN;
};

// Contains discrete systems of equation
typedef struct _SparseMat SparseMat;
struct _SparseMat {
	double   *A, *x, *b, *F, *d, *bbc;
	int      *Ic, *J, neq;
    int   *eqn_u, *eqn_v, *eqn_p;
    int   nnz, neq_mom, neq_cont;
};

// params contains the model parameters
typedef struct _params params;
struct _params {
	double  xmin, zmin, xmax, zmax, time, dx, dz, dt, dt0, dt_start, L0;
    double  xmin0, zmin0, xmax0, zmax0;
	double gx, gz;
	int Nx, Nz, Nt, step, Newton;
	int eta_avg, p_avg;
	int ismechanical, isperiodic_x, isinertial, iselastic, isnonnewtonian, isthermal, ispureshear_ale, free_surf, eqn_state, write_markers, write_debug, write_energies;
    double free_surf_stab;
    int dt_constant, moving_front, imp_advection, RK, line_search, thermal_eq, subgrid_diff, adiab_heat, shear_heat, advection, fstrain;
    int isPl_soft, surf_processes, cpc, surf_remesh, loc_iter, therm_pert, surf_ised1, surf_ised2, MantleID, topografix, aniso;
    double EpsBG, user0, user1, user2, user3, user4, user5, user6, user7, user8;
	char *input_file;
    int    Nb_phases;
    int    ncont;
    double Courant, mineta, maxeta;
    // Linear solver
    int decoupled_solve, lsolver, diag_scaling;
    double penalty, abs_tol_div, rel_tol_div, auto_penalty, compressible;
    // Deformation maps
    int nT, nE, nd, def_maps;
    double Pn, Tmin, Tmax, Emin, Emax, dmin, dmax, PrBG;
    // Surface processes
    double surf_diff, surf_sedirate, surf_baselev;
    // Stuff related Pips periodic aggregate deformation
    int cut_noise, rheo_on_cells, DefectCorrectionForm, HsOnly, HomoFields;
    double accu;
    // Initial thermal perturbation
    double therm_pert_x0, therm_pert_z0, therm_pert_dT, therm_pert_rad, cooling_time;
    // For rheological database reasons...
    int    force_act_vol_ast;
    double act_vol_dis_ast, act_vol_dif_ast;
    // Phase diagrams
    int    isPD, num_PD, *PDMnT, *PDMnP;
    double **PDMrho, *PDMTmin, *PDMTmax, *PDMPmin, *PDMPmax;
    // Visualisation
    int rec_T_P_x_z, rm_break;
};

// Nparams contains numerical parameters of the non-linear solver
typedef struct _n_params Nparams;
struct _n_params {
	int    nit, nit_max, stagnated;
    double tol_u, tol_p;
	double resx, resz, resp, rest;
    double resx_f, resz_f, resp_f;
	double vrlx,  prlx, trlx;
};

// grid contains all the fine grid arrays (double *)
typedef struct _grid grid;
struct _grid {
	int    Nx, Nz, NN, NC;
	double dx,dz;
	double *roger_x, *roger_z, *div_u, *u_in, *v_in, *p_in, *sxxd, *szzd, *sxz, *exxd, *ezzd, *exz, *VE_s, *VE_n, *sxxd0, *szzd0, *sxz0, *mu_s, *mu_n, *u_adv, *v_adv, *eta_phys_n, *kx, *kz, *Cv, *Qr, *eta_phys_s, *u_start, *v_start, *p_start;
	int    *iter_smooth;
	int    *nb_part_cell, *nb_part_vert;
	BC     BCu, BCv, BCp;
	BCT    BCt, BCg;
	double *xg_coord, *zg_coord, *xc_coord, *zc_coord, *xvz_coord, *zvx_coord, *xg_coord0, *zg_coord0, *xg_coord_ext, *zg_coord_ext;
	double *eta_s, *eta_n, *rho_s, *rho_n, *rho_app_s, *rho_app_n;
    double *strain_n, *strain_s;
	double *u, *v, *p;
	double *ru, *rv, *rp;
	double *rhs_u, *rhs_v, *rhs_p, *rhs_t;
	double p_scale;
    double *alp, *bet, *p_lith, *dp;
    double *VxVz, *VzVx;
    int    *P2N, *P2C;
    int    *kvx, *lvx, *kvz, *lvz, *kp, *lp, *kn, *ln;
    double **phase_perc_n, **phase_perc_s;
    double *sxxd0_s, *szzd0_s, *sxz0_n, *exxd_s, *ezzd_s, *exz_n, *sxz_n;
    double *rho_app_s0, *rho_app_n0;
    double Ut, Ue, W, *Work, *Uelastic, *Uthermal, *Time, *Short;
    double *T, *dT, *d, *d0, *phi, *X;
    double *eII_el, *eII_pl, *eII_pl_s, *eII_pwl, *eII_exp, *eII_lin, *eII_gbs, *eII_cst, *A2_pwl_n, *eii_n, *eii_s, *tii0_n, *tii0_s;
    double *eII_pwl_s, *A2_pwl_s;
    double *exx_el, *ezz_el, *exz_el, *exx_diss, *ezz_diss, *exz_diss;
    int   *comp_cells;
    
    // To remove
    double *exx_pwl_n, *exz_pwl_n, *exx_pwl_s, *exz_pwl_s, *exx_pl, *exz_pl;
    
    double *cell_min_z, *cell_max_z, *vert_min_z, *vert_max_z;
    double *phi_n, *phi_s, *C_n, *C_s;
    double *rhoUe0;
    double *exz_n_el, *exz_n_diss, *exz_n_pl;
};

// Contains information needed for the direct solver
typedef struct _DirectSolver DirectSolver;
struct _DirectSolver {
    int mtype;           /* Real unsymmetric matrix */
    double res, res0;
    int nrhs;             /* Number of right hand sides */
    int maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    int i, j;
    double ddum;                  /* Double dummy */
    int idum;                 /* Integer dummy */
    char *uplo;
    int n_th;
    int flag;
    int     Analyze;
    cholmod_factor *Lfact;
    cholmod_common c;
};

//---------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------ FUNCTION PROTOTYPES ------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------//

// Miscellaneous functions
void Initialise1DArrayDouble( double*, int, double );

// Input/Output
void ScaleMe( scale* );
void    ReadInputFile( char*, int*, int*, int*, int*, params*, scale*, mat_prop*, markers*, Nparams* );
void    UpdateInputFile( char[], int );
int     ReadInt2( FILE*, char[], int );
double  ReadDou2( FILE*, char[], double );
//float   ReadFlo2( FILE*, char[], float );
double  ReadMatProps( FILE* , char[], int, double );
//double  ReadDou21( FILE*, char[], double );
char*   ReadChar( FILE*, char[], char[] );
char*   ReadPhaseDiagram( FILE*, char FieldName[] );
double* ReadBin( char[], int, int, double );
//void    WriteBin(double*, int, const char[]);

// Allocation
void GridAlloc( grid* , params*  );
void GridFree( grid* , params*  );
void PartAlloc( markers*, params*  );
void PartFree( markers*, params* );

// Initialisation
void InitialiseSolutionFields( grid* , params* );
void GridIndices( grid* );
void SetGridCoordinates( grid*, params*, int, int );

// Memory
void* DoodzMalloc(  size_t );
void* DoodzCalloc( int, size_t );
void* DoodzRealloc( void*, size_t );
void  DoodzFree( void* );
void  AllocMat( SparseMat*, int );
void  FreeMat( SparseMat* );
//void  gridAllocStar( grid*, params*OutputSparseMatrix  );
//void  gridAllocStarStar( grid*, params*OutputSparseMatrix  );
//void  gridFree( grid*, params*OutputSparseMatrix  );
//void  MPartAllocStar( markers*, params* );
//void  MPartReAllocStar( markers );
//void  MPartFree( markers*, params* );
//void  MDefineNumGSiteration( Mparams, grid* );

// General function prototypes
void ArrayPlusArray( double*, double*, int );
void ArrayMinusArray( double*, double*, int );
void ArrayEqualArray( double*, double*, int );
void ArrayEqualArrayI( int*, int*, int );
//void ArrayPlusScalar( double*, double, int );
void ArrayTimesScalar( double*, double, int );
void ArrayTimesScalarArray( double*, double, double*, int );
void ArrayDividedScalarArray( double*, double, double*, int );
void ArrayPlusScalarArray( double*, double, double*, int );
void Initialise2DArray( double*, int, int, double );
void Initialise2DArrayInt( int*, int, int, int );
void MinMaxArray( double*, double, int, char* );
void MinMaxArrayVal( DoodzFP*, int, double*, double* );
void MinMaxArrayTag( DoodzFP*, double, int, char*, char* );
void MinMaxArrayPart( DoodzFP*, double, int, char*, int* ) ;
double SumArray( double*, double, int, char*);
//void MinMaxArrayF( float*, double, int, char*);
void MinMaxArrayI( int*, int, int, char*);
//void SumArrayF( float*, double, int, char*);
//void SumArrayI( int*, double, int, char* );
//void Initialise1DArrayDouble( double*, int, double );
void Initialise1DArrayChar( char*, int, char );
//void Print2DArrayDoubleTag( DoodzFP*, int, int, double, char* );
//void Print2DArrayChar( char*, int, int, double);
void Initialise1DArrayInt( int*, int, int );
void IsNanArray2DFP( DoodzFP*, int );
void IsInfArray2DFP( DoodzFP*, int );
void InterpCentroidsToVerticesDouble( double*, double*, grid*, params*, scale* );
//void InterpVerticesToCentroidsDouble( double*, double*, grid*, params*, scale* );
//
//
//// Grid related function
//void MInitialiseSolutionFields( grid*, params* );
//void MinitMG( grid*, paramsOutputSparseMatrix  );
void SetBCs( grid*, params*, scale, markers*, mat_prop* );
//void MSetRes( grid*, paramsOutputSparseMatrix  );
//void gridFillLevels( grid*, paramsOutputSparseMatrix  );
//void gridSetCoords( grid*, params*OutputSparseMatrix );
//void MSetGrid( grid*, params, int, int, int);
//void eval_anal_Dani( double*, double*, double*, double*, double*, double*, double, double, int, double, double, double );
void ComputeLithostaticPressure( grid*, params*, double, scale, int );

//// Particles
void PutPartInBox( markers*, grid*, params, surface, scale );
void SetParticles( markers*, scale, params, mat_prop* );
void PartInit( markers*, params* );
void Interp_P2U( markers, DoodzFP*, grid*, double*, double*, double*, int, int, int, char* , params* );
void Interp_P2N( markers, DoodzFP*, grid*, double*, double*, double*, int, int, params* );
void Interp_P2C( markers, DoodzFP*, grid*, double*, double*, double*, int, int );
void Interp_Grid2P( markers, DoodzFP*, grid*, double*, double*, double*, int, int, char* );
void Interp_Grid2P_strain( markers, DoodzFP*, grid*, double*, double*, double*, int, int, char* );
//void FreeP2Mesh( grid* );
void Interp_Phase2VizGrid( markers, int*, grid*, char*, double*, double*, int, int, params, surface );
void ParticleInflowCheck ( markers*, grid*, params, surface, int);
//
//// Stokes
//void ResidualCalc2( grid*OutputSparseMatrix *, params, int, double*, double*, double*, double*, double*, double*, int, scale );
//void Gauss_Seidel( grid*OutputSparseMatrix *, params, scale );
//void SaveRHS0( grid* );
void EvaluateRHS( grid*, params, scale, double );
//void DefineResidualScales( Mparams*, grid* );
//void ResidualNorm( gridOutputSparseMatrix  * );
//void StrainStressCalc( grid* );
//void PressureScaling( Mparams, grid*, mat_prop, Eparams );
//
//// Interpolation fine to coarse (restriction)
//void RestrictAll( grid *, int, int );
//void Interp_F2C_U5( grid*, int, double*, double* );
//void Interp_F2C_V5( grid*, int, double*, double* );
//void Interp_F2C_Peta2( grid*, int, double*, double*, double*, double* );
//void Interp_F2C_C4( grid*, int, double*, double*, int );
//void Interp_F2C_N4( grid*, int, double*, double*, int );
//
//// Interpolation coarse to fine (prolongation)
//void ProlongAll( grid*, intOutputSparseMatrix * );
//void Interp_C2F_C( grid*, int, double*, double* );
//void Interp_C2F_U( grid*, int, double*, double* );
//void Interp_C2F_V( grid*, int, double*, double* );
//
//// Visualisation prototypes
void create_output_hdf5( const char[] );
void AddGroup_to_hdf5( const char[], const char[] );
void AddFieldToGroup_generic( int, const char[], const char[], const char[], char, int, void*, int );
//void Myfopen( char*, FILE** );
//void MViz_vtk( grid*, char*OutputSparseMatrix  );
void WriteOutputHDF5( grid*, markers*, surface*, markers*, params, char*, mat_prop, scale );
void WriteOutputHDF5Particles( grid*, markers*, surface*, markers*, surface*, markers* , params, char*, mat_prop, scale );
void WriteResiduals( grid, params, Nparams, scale );
//
//// Viscosity continuation
//void LimitViscosityFields( Mparams, grid*, Eparams );
//void FieldLimiter( grid*, double*, double*, double*, int );
//

//
//// Phase diagrams
//void AllocatePhaseDiagrams( params* );
//void FreePhaseDiagrams( params* );
//
// Output
void MakeBreakpointParticles( markers*, grid*, markers*, markers*, params , surface*, surface*, scale );
void LoadBreakpointParticles( markers*, grid*, markers*, markers*, params*, surface*, surface*, scale );
void DeletePreviousBreakpoint( int, int );
//
// Direct solver
void EvalNumberOfEquations( grid*, SparseMat* );
void SAlloc( SparseMat*, int );
void SFree( SparseMat* );
void BuildStokesOperator( grid*, params, int, double*, double*, double*, SparseMat*, int );
void EvaluateStokesResidual( SparseMat*, Nparams*, grid*, params, scale, int );
//void StokesDirectSolve( grid*OutputSparseMatrix *, params, int, double*, double*, double*, double*, double*, double* );
//void StokesDirectSolveCoarse( grid*OutputSparseMatrix *, params, int, double*, double*, double*, double*, double*, double* );
//void DirectSolverCall( double*, int*, int*, grid*, int );
void SolveStokes( SparseMat*, DirectSolver* );
void SolveStokesDefect( SparseMat*, DirectSolver*, Nparams*, grid*, params*, markers*, markers*, surface*, mat_prop, scale );
void DirectStokes( SparseMat*, DirectSolver*, double* , double* );
void ExtractSolutions( SparseMat*, grid*, params* );
void InitialiseSolutionVector( grid*, SparseMat*, params* );
//
//// Viscoelastoplasticity
//void EffectiveViscosityCalc( markers*, mat_prop, params, scale );
//void EffectiveViscosityCalcGrid(  grid*, mat_prop, params, scaleOutputSparseMatrix  );
//void StrainStressCalcVE( grid* );
void RotateStresses( grid, markers*, params, scale* );
void UpdateParticleStress( grid*, markers*, params*, mat_prop*, scale* );
void ShearModulusGrid( grid*, mat_prop, params, scale );
void CohesionFrictionGrid( grid* , mat_prop, params, scale  );
//
//// Non-Newtonian rheology
//double Viscosity( int, double, double, double, double, double, double, double, double, double, double, double, mat_prop*, params*, scale*, int, double*, double*, double*, double* , double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double* , double, double, double );
//double ViscosityPlast( int, double, double, double, double, double, double, double, double, double, double, double, mat_prop*, params*, scale*, int, double*, double*, double*, double* , double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double* ,double);
void UpdateNonLinearity( grid*, markers*, markers*, surface*, mat_prop, params*, Nparams*, scale, int, double );
double LineSearch( SparseMat*, double*, grid*, params*, Nparams*, markers*, markers*, surface*, mat_prop, scale );
//void StressStrainRateInvariantsMarkers( markers*, grid*, scale );
//void StressStrainRateInvariantsGrid( grid*, scale );
////void NonNewtonianViscosityMarkers( markers*, mat_prop, paramsOutputSparseMatrix , scale, int );
void NonNewtonianViscosityCells( grid*, mat_prop*, params*, Nparams, scale*, int );
void NonNewtonianViscosityGrid( grid*, mat_prop*, params*, Nparams, scale*, int );
//void NonNewtonianViscosityGridPart( grid*, markers*, mat_prop*, params*OutputSparseMatrix , scale*, int );
void StrainRateComponents( grid*, scale, params* );
void GenerateDeformationMaps( grid*, mat_prop*, params*, Nparams, scale*, int );
void UpdateParticleGrainSize( grid*, scale, params, markers*, mat_prop* );
//
// Advection
void DefineInitialTimestep( params*, grid*, markers, mat_prop, scale );
void EvaluateCourantCriterion( double*, double*, params*, scale, grid*, int);
//void EvaluateCourantCriterionParticles( markers, params*, scale);
void RogerGunther( markers*, params, grid, int, scale );
void isout( markers*, params );
void CountPartCell    ( markers*, grid* , params, surface, surface, int, scale );
void CountPartCell_Old( markers*, grid* , params, surface, int, scale );
//void CountPartVertex ( markers*, grid*, params );
void AccumulatedStrain( grid*, scale , params, markers* );
void PureShearALE( params*,  grid*, markers*, scale );
void VelocitiesOnCenters( double*, double*, double*, double*, int, int, scale );
//void VelocitiesOnVertices( double*, double*, double*, double*, int, int, scale );
//void FirstOrderUpwindAdvection( double*, double*, double*, double*, grid*, int, int, params, scale, int );
void VelocitiesToParticles( grid*, markers*, DoodzFP*, DoodzFP*, params, scale );
void DeformationGradient ( grid, scale, params , markers * );
//
//
//// Energy
void UpdateParticleEnergy( grid*, scale, params, markers*, mat_prop* );
//void UpdateParticleVelocity( grid*, scaleOutputSparseMatrix , params, markers* );
void EnergyDirectSolve( grid*, params, double*, double*, double*, double*, markers*, double, int, int, scale, int );
cholmod_factor* FactorEnergyCHOLMOD( cholmod_common*, double*, int*, int*, int, int );
void SolveEnergyCHOLMOD( cholmod_common*, cholmod_factor*, double*, double*, int, int );
//void EnergyTemperatureConvertPart( double*, markers*, mat_prop, params );
void ThermalSteps( grid*, params, double*, double*, double*, double*, markers*, double, scale );
void Energies( grid*, params, scale );
void SetThermalPert( grid*, params, scale );
void UpdateMaxPT ( scale, params, markers* );

// Free surface routines
void BuildInitialTopography( surface*, markers*, params, grid, scale );
void SetTopoChainHorizontalCoords( surface*, markers*, params, grid, scale );
void AllocateMarkerChain( surface*, markers*, params );
void FreeMarkerChain( surface*, markers* );
void CellFlagging( grid*, params, surface, scale );
void ProjectTopography( surface*, markers*, params, grid, scale, double*, int );
//double TopoFun( double, int, surface, scale );
void MarkerChainPolyFit( surface*, markers*, params, grid );
void CleanUpSurfaceParticles( markers*, grid*, surface, scale );
void RemeshMarkerChain( markers*, surface*, params, scale, grid*, int );
void SurfaceDensityCorrection( grid*, params, surface, scale);
void SurfaceVelocity( grid*, params, surface*, markers*, scale );
void UpdateDensity( grid*, markers*, mat_prop*, params*, scale* );
//void AdvectFreeSurf( markers*, params, scale );
//void PhaseGrowth( markers*, params, grid* );
void DiffuseAlongTopography( grid*, params, scale, double*, int, double, double );
void AddPartSed( markers *, mat_prop , markers *, surface *, params , scale , grid *);
void CorrectTopoIni( markers *, mat_prop , markers *, surface *, params , scale , grid *);
//
// Decoupled solver
void KSPStokesDecoupled( SparseMat*,  SparseMat*,  SparseMat*,  SparseMat*, DirectSolver*, double*, double*, double*, params, grid*, scale, SparseMat*, SparseMat*, SparseMat*,  SparseMat*,  SparseMat* );
double LineSearchDecoupled( SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, double*, grid*, params*, Nparams*, markers*, markers*, surface*, mat_prop, scale );
void EvaluateStokesResidualDecoupled( SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, Nparams*, grid*, params, scale, int );
void BuildStokesOperatorDecoupled( grid*, params,int, double*, double*, double*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, int );
void SolveStokesDecoupled( SparseMat*, SparseMat*, SparseMat*,  SparseMat*, SparseMat*, DirectSolver*, params, grid*, scale );
void SolveStokesDefectDecoupled( SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, DirectSolver*, Nparams*, grid*, params*, markers*, markers*, surface*, mat_prop, scale, SparseMat*, SparseMat*, SparseMat* );
void AddCoeff2( int*, double*, int, int, int*, double, int, double, double* );
void MergeParallelMatrix( SparseMat*, double**, int**, int**, grid*, int*, int*, int*, int*, int*, int, char*, int* );
void DirectStokesDecoupled( SparseMat*, SparseMat*, SparseMat*,  SparseMat*, DirectSolver*, double*, double*, double*, params, grid*, scale, SparseMat* );
void DirectStokesDecoupledComp( SparseMat*, SparseMat*, SparseMat*,  SparseMat*, DirectSolver*, double*, double*, double*, params, grid*, scale, SparseMat* );
//void ConvertTo1Based( int*, int );
void DecompressCSRtoTriplets( int, int*, int* );
//void  ArrayEqualScalarArray( DoodzFP*, DoodzFP, DoodzFP*, int );
//void ApplyBC( grid*, params );
//void GridIndices( grid* );
//void ScaleBack          (float*,  double , int );
void ScaleBackD         (double*, double , int );
//void DoubleToFloat      (double*,  float*, int );
//
// Newton
void BuildJacobianOperatorDecoupled( grid*, params,int, double*, double*, double*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, SparseMat*, int );
//
//
// Flow law Database
void ReadDataPowerLaw   ( mat_prop*, params*, int, int, scale* );
void ReadDataLinear     ( mat_prop*, params*, int, int, scale* );
void ReadDataGBS        ( mat_prop*, params*, int, int, scale* );
void ReadDataExponential( mat_prop*, params*, int, int, scale* );
void ReadDataGSE        ( mat_prop*, params*, int, int, scale* );
void AllocatePhaseDiagrams( params* );
void FreePhaseDiagrams( params* );

double DotProduct( DoodzFP*, DoodzFP*, int  );
void BackToSolutionVector( cholmod_dense*, cholmod_dense*, double*, grid*, SparseMat* );
void NormResidualCholmod( double*, double*, cholmod_dense*, cholmod_dense*, int, int, params, scale, int );
void BuildInitialSolutions( double*, double*, grid* );
void copy_cholmod_to_cs_matrix( cholmod_sparse*, cs* );
void copy_cs_to_cholmod_matrix( cholmod_sparse*, cs* );
void copy_vec_to_cholmod_dense( cholmod_dense*, DoodzFP* );
void copy_cholmod_dense_to_cholmod_dense( cholmod_dense*, cholmod_dense* );
void cholmod_dense_plus_cholmod_dense( cholmod_dense*, cholmod_dense* );

void ApplyBC( grid*, params* );
void AssignMarkerProperties (markers*, int, int, params* );

// GLOBAL
//void Interp_P2G( markers, DoodzFP*, grid*, double*, double*, double*, int, int, double, double, int, int, params*, char*  );

void AdvectFreeSurf_BEN( markers*, params, scale );
void BuildInitialTopography_BEN( surface*, markers*, params, grid, scale );
void SetTopoChainHorizontalCoords_BEN( surface*, markers*, params, grid, scale );
void AllocateMarkerChain_BEN( surface*, markers*, params );
void FreeMarkerChain_BEN( surface*, markers* );
void CellFlagging_BEN( grid*, params, surface, scale );
void ProjectTopography_BEN( surface*, markers*, params, grid, scale, double*, int );
//double TopoFun( double, int, surface, scale );
void MarkerChainPolyFit_BEN( surface*, markers*, params, grid );
void CleanUpSurfaceParticles_BEN( markers*, grid*, surface, scale );
void RemeshMarkerChain_BEN( markers*, surface*, params, scale, grid*, int );
void SurfaceDensityCorrection_BEN( grid*, params, surface, scale);
void SurfaceVelocity_BEN( grid*, params, surface*, markers*, scale );
void UpdateDensity_BEN( grid*, markers*, mat_prop*, params*, scale* );
//void AdvectFreeSurf( markers*, params, scale );
//void PhaseGrowth( markers*, params, grid* );
void DiffuseAlongTopography_BEN( grid*, params, scale, double*, int, double, double );
void AddPartSed_BEN( markers *, mat_prop , markers *, surface *, params , scale , grid *);
void CorrectTopoIni_BEN( markers *, mat_prop , markers *, surface *, params , scale , grid *);
void EvaluateRHS_BEN( grid*, params, scale, double );
void UpdateNonLinearity_BEN( grid*, markers*, markers*, surface*, mat_prop, params*, Nparams*, scale, int, double );
void PressureScaling_BEN(  grid *, mat_prop , params  );
void StrainStressCalc_BEN( grid* mesh );
void UpdateParticleVelocity_BEN( grid*, scale, params, markers* );

void BuildStokesOperator_BEN( grid*, params, int, double*, double*, double*, SparseMat*, int );
int EvalNumberOfEquations_BEN( grid* mesh, SparseMat *Stokes );
void ExtractSolutions_BEN( SparseMat*, grid*, params );
void EvaluateStokesResidual_BEN( SparseMat *, Nparams *, grid *, params , scale , int  );
void EvaluateCourantCriterion_BEN( double* , double* , params *, scale , grid*, int  );
void RogerGunther_BEN( markers *, params , grid , int , scale  );

void Interp_P2N_BEN ( markers , DoodzFP* , grid *, double* , double* , double* , int , int, params*  ) ;
void Interp_P2C_BEN ( markers , DoodzFP* , grid *, double* , double* , double* , int , int   ) ;
void Interp_P2U_BEN ( markers , DoodzFP* , grid *, double* , double* , double* , int , int , int , char*  ) ;
void Interp_Grid2P_BEN ( markers , DoodzFP* , grid *, double* , double* , double* , int , int , char * ) ;
void CountPartCell_BEN ( markers* , grid *, params , surface, int, scale  );
void PutPartInBox_BEN( markers *, grid *, params , surface , scale  );
void SetBCs_BEN( grid *, params *, scale , markers* , mat_prop * );

void SetParticles_BEN( markers *, scale , params , mat_prop*   );
void BuildInitialTopography_BEN( surface *, markers *, params , grid , scale  );
void SolveStokes_BEN( SparseMat*, DirectSolver* );

double Grid2P( markers*, double*, double*, double*, int, int , char *, double, double, int );
void RogerGuntherII( markers*, params, grid, int, scale );
void AccumulatedStrainII( grid*, scale, params, markers*, double*, double*, int, int, char * );
void AdvectFreeSurf( markers*, params, scale );

void RotateDirectorVector( grid, markers*, params, scale* );
void UpdateParticlePressure( grid*, scale, params, markers*, mat_prop* );
void DetectCompressibleCells ( grid* , params*  );
void ScaleVelocitiesRHSBack(SparseMat*, double*);


void ExtractDiagonalScale(SparseMat *, SparseMat *, SparseMat *, SparseMat * );

void ScaleMatrix(SparseMat *, SparseMat *, SparseMat *, SparseMat * ) ;
