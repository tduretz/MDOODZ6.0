// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2021  MDOODZ Developper team
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
#include "cholmod.h"
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

double ViscosityTest( int phase, double G, double T, double P, double d, double phi, double X0, double Exx, double Ezz, double Exz, double Gxx, double Gzz, double Gxz, double Txx0, double Tzz0, double Txz0, mat_prop* materials, params *model, scale *scaling, double *txxn, double *tzzn, double *txzn, double* etaVE, double* VEcoeff, double* Eii_el, double* Eii_pl, double* Eii_pwl, double* Eii_exp , double* Eii_lin, double* Eii_gbs, double* Eii_cst, double* Exx_el, double* Ezz_el, double* Exz_el, double* Exx_diss, double* Ezz_diss, double* Exz_diss, double *d1, double strain_acc, double dil, double fric, double C, double *detadexx, double *detadezz, double *detadexz, double *detadp, double P0,  double *X1, double *OverS, double *ddivpdexx, double *ddivpdezz, double *ddivpdexz, double *ddivpdp, double *Pcorr, double *drhodp, double *rho, double beta, double div, double *div_el, double *div_pl, double *div_r ) {
    
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
    double Eii_vis=0.0, Eii= 0.0, Gii = 0.0, f_ani = 0.0;
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
    double K = 1.0/beta, dQdP=0.0;

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
    double retro = 1.0;
    double ddivrdpc = 0.0, divr = 0.0, drhodX = 0.0;
    double dCgbsdP = 0.0, dClindP = 0.0, dCpwldP = 0.0, Jii, alpha, rho_ref, drho_ref_dP;
    
    alpha  = materials->alp[phase];
    
    if (model->diffuse_X==0) constant_mix  = 0;
    
    //------------------------------------------------------------------------//
    
    // Initialise strain rate invariants to 0
    *Eii_exp = 0.0; *Eii_lin = 0.0; *Eii_pl = 0.0; *Eii_pwl = 0.0; *Eii_el = 0.0, *Eii_gbs=0, *Eii_cst=0.0;
    *txxn=0.0; *tzzn=0.0; *txzn=0.0; *etaVE=0.0; *VEcoeff=0.0, *Exx_el=0.0, *Ezz_el=0.0, *Exz_el=0.0, *Exx_diss=0.0, *Ezz_diss=0.0, *Exz_diss=0.0, *d1=0.0;
    *detadexx=0.0;  *detadezz=0.0;  *detadexz=0.0;  *detadp=0.0;
    *ddivpdexx=0.0; *ddivpdezz=0.0; *ddivpdexz=0.0; *ddivpdp=0.0;
    *X1=0.0; *drhodp = 0.0; *rho = 0.0;
    *OverS = 0.0; *Pcorr = 0.0; *div_el = 0.0; *div_pl = 0.0; *div_r = 0.0;
    
    //------------------------------------------------------------------------//
    
    // Invariants
    double Tyy, Eyy  = -( Exx  + Ezz  ), Gyy  = -( Gxx  + Gzz  ), Tyy0 = -( Txx0 + Tzz0 ); // definition of deviatoric tensor
    Eii   = sqrt(1.0/2.0*(Exx*Exx + Ezz*Ezz + Eyy*Eyy) + Exz*Exz);
    Gii   = sqrt(1.0/2.0*(Gxx*Gxx + Gzz*Gzz + Gyy*Gyy) + Gxz*Gxz);
    f_ani = Gii/Eii;
    if (Eii*scaling->E<1e-30) Eii=1e-30/scaling->E;
    
    // P corr will be corrected if plasticity feedbacks on pressure (dilation)
    *Pcorr = P;
    
    //------------------------------------------------------------------------//

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
        C_pwl   = pow(2.0*f_ani*B_pwl, -n_pwl);
        dCpwldP = -C_pwl*Va_pwl/R/T;
    }
    if ( diffusion == 1 ) {
        if (m_lin>0.0 && d<1e-13/scaling->L){
            printf("Cannot run with grain size dependent viscosity if grain size is set to 0 --> d = %2.2e!!!\n", d*scaling->L);
            exit(1);
        };
        B_lin = F_lin * pow(A_lin,-1.0/n_lin) * exp( (Ea_lin + P*Va_lin)/R/n_lin/T ) * pow(f_lin, -r_lin/n_lin) * exp(-a_lin*phi/n_lin); // * pow(d, m_lin/n_lin) !!!!!!!!!!!!!!!!!!!!!!!!
        C_lin = pow(2.0*f_ani*B_lin, -n_lin);
        dClindP = -C_lin*Va_lin/R/T;
    }
    if ( gbs == 1 ) {
        B_gbs = F_gbs * pow(A_gbs,-1.0/n_gbs) * exp( (Ea_gbs + P*Va_gbs)/R/n_gbs/T ) * pow(d, m_gbs/n_gbs) * pow(f_gbs, -r_gbs/n_gbs) * exp(-a_gbs*phi/n_gbs);
        C_gbs = pow(2.0*f_ani*B_gbs, -n_gbs);
        dCgbsdP = -C_gbs*Va_gbs/R/T;
    }
    if ( peierls   == 1 ) {
//        ST                           = Ea_exp/R/T * pow((1.0-gamma),(q-1.0)) * q*gamma;
        ST                           = Ea_exp/R/T * 2.0*gamma*(1-gamma);
//        if ( (int)t_exp == 0) F_exp  = 1.0;
//        if ( (int)t_exp == 1) F_exp  = 1.0/6.0*pow(2.0,1.0/(ST+n_exp)) * pow(3.0,(ST+n_exp-1.0)/2.0/(ST+n_exp));
//        if ( (int)t_exp == 2) F_exp  = 1.0/4.0*pow(2,1.0/(ST+n_exp));
        // old
//        B_exp                   = F_exp * pow(E_exp*exp(-Ea_exp/R/T*pow(1.0-gamma,2.0)), -1.0/(ST+n_exp)) * pow(gamma*S_exp, ST/(ST+n_exp));
//        C_exp                   = pow(2.0*f_ani*B_exp, -(ST+n_exp));
        // committed in April 2020 - should work without aniso (and no correction)
//        C_exp = E_exp *exp(-Ea_exp/R/T * pow(1.0-gamma,2.0)) * pow(gamma*S_exp,-ST);  // ajouter Fexp
//        B_exp = 0.5*pow(C_exp, -1./(n_exp+ST) );
        // new
        double Arr_exp = exp(-Ea_exp/R/T*pow(1.0-gamma,2.0));
        F_exp = pow( pow(2.0,1.0-ST-n_exp) / pow(sqrt(3.0), ST+n_exp+1.0), 1.0/(ST+n_exp));
        B_exp = F_exp * ( pow(gamma*S_exp, ST/(ST+n_exp)) / pow( E_exp*Arr_exp, 1.0/(ST+n_exp)) );
        C_exp = pow(2.0*f_ani*B_exp, -(ST+n_exp));
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
            
            if ( X<X0 && NoReturn == 1 ) {
                retro  = 0.0;
                X      = X0;
            }
        }
        // paprameters of end-members
        ndis1  = materials->npwl[phase];               ndis2  = materials->npwl[materials->reac_phase[phase]];
        Adis1  = materials->Apwl[phase];               Adis2  = materials->Apwl[materials->reac_phase[phase]];
        Qdis1  = materials->Qpwl[phase];               Qdis2  = materials->Qpwl[materials->reac_phase[phase]];
        rho1   = materials->rho[phase];                rho2   = materials->rho[materials->reac_phase[phase]];

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

        // Proper choices of corrections factors
        if ( (int)t_pwl == 0 ) {
            F_pwl = 1.0;
        }
        if ( (int)t_pwl == 1 ) {
            F_pwl = 1.0/6.0*pow(2.0,1.0/n_pwl) * pow(3.0,(n_pwl-1.0)/2.0/n_pwl);
        }
        if ( (int)t_pwl == 2 ) {
            F_pwl = 1.0/4.0*pow(2,1.0/n_pwl);
        }
        
        // Override power-law flow law parameters
        B_pwl  = pre_factor * F_pwl * pow(A_pwl,-1.0/n_pwl) * exp( (Ea_pwl)/R/n_pwl/T );
        C_pwl  = pow(2.0*B_pwl, -n_pwl);
    }
    
    // set pointer value
    *X1 = X;
    
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
//    eta_ve                  = 0.5*(eta_up+eta_lo);
    eta_ve = eta_up;

    // Local iterations
    for (it=0; it<nitmax; it++) {
        
        // Function evaluation at current effective viscosity
        Tii = 2.0 * f_ani * eta_ve * Eii;
        if ( constant    == 1 ) *Eii_cst = Tii/f_ani/2.0/eta_cst;
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
    if ( it>10 ) printf("Warnung: more that 10 local iterations, there might be a problem...\n");
    
    // Recalculate stress components
    Tii                  = 2.0*eta_ve*f_ani*Eii;
    
    //------------------------------------------------------------------------//
    
    // Check yield stress
    F_trial = Tii - Tyield;
    
    double Tyield_trial = Tyield, Tiic;
    double F_corr1=F_trial, F_corr2=F_trial, Tc, Pc_chk, Tc_chk;
    
    // Select appropriate dilation angle for tensile domain, see SimpleYields.m
    if (tens == 1) {
        eta_vp   = eta_vp0 * pow(Eii, 1.0/n_vp - 1);
        Pc_chk    = -(C*cos_fric) / (- Tii/P + sin_fric);
        Tc_chk    = Tii/P*Pc_chk;
        sin_dil = (P*eta_ve + P*eta_vp - Pc_chk*eta_ve - Pc_chk*eta_vp)/(K*dt*(C*cos_fric + Pc_chk*sin_fric - Tii));
    }
    
    if (F_trial > 1e-17) {
        
        // Initial guess - eta_vp = 0
        is_pl    = 1;
        eta_vp   = eta_vp0 * pow(Eii, 1.0/n_vp - 1);
        gdot     = F_trial / ( eta_ve + eta_vp + K*dt*sin_fric*sin_dil);
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
            Tiic    = Tii - eta_ve*gdot;
            F_trial = Tiic - Tyield;
            
            // Residual check
            res = fabs(F_trial);
            if (noisy>0 ) printf("%02d Visco-Plastic iterations It., tens = %d F = %2.2e Frel = %2.2e --- n_vp = %2.2e, eta_vp = %2.2e\n", it, tens, res, res/F_trial0, n_vp, eta_vp*scaling->eta);
            if ( res < tol || res/F_trial0 < tol ) break;
            dFdgdot  = - eta_ve - eta_vp/n_vp - K*dt*sin_fric*sin_dil;
            gdot    -= F_trial / dFdgdot;
            
        }
        if ( it==nitmax-1 && (res > tol || res/F_trial0 > tol)  ) { printf("Visco-Plastic iterations failed!\n"); exit(0);}
        
        eta_vep = Tiic / (2.0*Eii);
    }
    
    
    // ----------------- Reaction volume changes
    
    if (ProgressiveReaction == 1) {
        rho_ref      = (1.0-X)*rho1 + X*rho2;
        *rho         = rho_ref * exp(1.0/K * Pc - alpha * T);
    }
    else {
        rho_ref      = rho1;
        *rho         = rho_ref * exp(1.0/K * Pc - alpha * T );
    }

    if (is_pl == 0) {
        (*etaVE)    = eta_ve;
        (*div_pl)   = 0.0;
    }
    else {
        (*etaVE)    = eta_vep;
        (*div_pl)   = divp;
    }
    
    /*----------------------------------------------------*/
    /*----------------------------------------------------*/
    /*----------------------------------------------------*/
    
    // Viscosity limiter
    if( *etaVE > maxEta ) {
        *etaVE = maxEta;
    }
    
    if( *etaVE < minEta ) {
        *etaVE = minEta;
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void PhaseRheologyLoop( int centroid, double sign, double denom, double Exx, double Ezz, double Exz, double P, double ani, double d0, double d1, int c0, double** vol,
               double* G, double* T, double* P0, double* gs0, double* phi0, double* X0, double* txx0, double* tzz0, double* txz0, double* beta, double* div,
               double* strain, double* dil, double* fric, double* C,
               params* model, mat_prop* materials, scale* scaling,
               double* detadE, double* ddivpdE, double* eta ) {
    
    int cond;
    double gxz   = 2.0*Exz;
    double Gxx   = Exx*(1.0 - ani*d0) + Ezz*ani*d0 + gxz*ani*d1;
    double Gzz   = Ezz*(1.0 - ani*d0) + Exx*ani*d0 - gxz*ani*d1;
    double Gxz   = Exx*ani*d1 - Exx*ani*d1 + gxz*(ani*(d0 - 0.5) + 0.5);  // Gxz = Exz if isotropic
    
    double txx1, tzz1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, dnew, detadexx, detadezz, detadexz, detadp, Xreac, OverS, ddivpdexx, ddivpdezz, ddivpdexz, ddivpdp, Pcorr, drhodp, rho, div_el, div_pl, div_r;

    // Loop on phases
    for (int p=0; p<model->Nb_phases; p++) {
                        
        cond =  fabs(vol[p][c0])>1.0e-13;

        if ( cond == 1 ) {
            ViscosityTest( p, G[c0], T[c0], P, gs0[c0], phi0[c0], X0[c0], Exx, Ezz, Exz, Gxx, Gzz, Gxz, txx0[c0], tzz0[c0], txz0[c0], materials, model, scaling, &txx1, &tzz1, &txz1, &etaVE, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &exx_el, &ezz_el, &exz_el, &exx_diss, &ezz_diss, &exz_diss, &dnew, strain[c0], dil[c0], fric[c0], C[c0], &detadexx, &detadezz, &detadexz, &detadp, P0[c0], &Xreac, &OverS, &ddivpdexx, &ddivpdezz, &ddivpdexz, &ddivpdp, &Pcorr, &drhodp, &rho, beta[c0], div[c0], &div_el, &div_pl, &div_r );
        }
        
        if ( model->eta_avg == 0) { // ARITHMETIC AVERAGE
            if ( cond == 1 ) detadE[c0]      += sign*vol[p][c0] * etaVE/(denom);
        }
        
        if ( model->eta_avg == 1) { // HARMONIC AVERAGE
            if ( cond == 1 ) detadE[c0]      += sign*vol[p][c0] * etaVE/(denom) / pow(etaVE,2.0);
        }
        
        if ( model->eta_avg == 2) { // GEOMETRIC AVERAGE
            if ( cond == 1 ) detadE[c0]      += sign*vol[p][c0] * etaVE/(denom) / etaVE;
        }
        // OTHERs
        if (centroid == 1) {
            if ( cond == 1 ) ddivpdE[c0]     += sign*vol[p][c0] * div_pl/(denom);
        }
    }
    
    // ----------------------------------------------------------------------------//
    if ( model->eta_avg == 1) { // HARMONIC AVERAGE
        detadE[c0] *= pow(eta[c0],2.0);
    }
    
    if ( model->eta_avg == 2) { // GEOMETRIC AVERAGE
        detadE[c0] *= eta[c0];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ViscosityDerivatives( grid *mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling ) {

    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1, k1, cond;
    double eta, txx1, tzz1, txz1, Pn, Tn, etaVE, VEcoeff=0.0, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, div_el, div_pl, div_r;
    double exx_pwl, exz_pwl, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss;
    int average = model->eta_avg, UnsplitDiffReac = model->UnsplitDiffReac;
    double detadexx, detadezz, detadexz, detadp;
    double Xreac;
    double OverS;
    double ddivpdexx, ddivpdezz, ddivpdexz, ddivpdp, Pcorr, drhodp, rho;
    double Exx, Ezz, Exz, gxz, Gxx, Gzz, Gxz, el, etae, ani, d0, d1, nx, nz;
    double Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det;
    double tol = 1e-5;

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
#pragma omp parallel for shared( mesh  ) private( cond, k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, detadexx, detadezz, detadexz, detadp, Xreac, OverS, ddivpdexx, ddivpdezz, ddivpdexz, ddivpdp, Pcorr, drhodp, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, gxz, Gxx, Gzz, Gxz, el, etae, ani, d0, d1, nx, nz, Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det ) firstprivate( UnsplitDiffReac, materials, scaling, average, model, Ncx, Ncz, tol )
    for ( k1=0; k1<Ncx*Ncz; k1++ ) {

        k      = mesh->kp[k1];
        l      = mesh->lp[k1];
        c0     = k  + l*(Ncx);

        // Initialise arrays to 0
        mesh->detadexx_n[c0]       = 0.0;
        mesh->ddivpdexx_n[c0]      = 0.0;
        mesh->detadezz_n[c0]       = 0.0;
        mesh->ddivpdezz_n[c0]      = 0.0;
        mesh->detadgxz_n[c0]       = 0.0;
        mesh->ddivpdgxz_n[c0]      = 0.0;
        mesh->detadp_n[c0]         = 0.0;
        mesh->ddivpdp_n[c0]        = 0.0;

        // Loop on grid nodes
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31 ) {
            
            //----------------------------------------------------------//
            if ( model->iselastic==1 ) etae      = model->dt*mesh->mu_n[c0];
            else                       etae      = 1.0; // set to arbitrary value to avoid division by 0.0
            //----------------------------------------------------------//
            if ( model->aniso == 0 ) {
                ani = 0.0; d0   = 0.0; d1   = 0.0;
            }
            else {
                // Director
                nx = mesh->nx0_n[c0];
                nz = mesh->nz0_n[c0];
                // See Anisotropy_v2.ipynb
                if ( model->aniso_fstrain  == 0 ) ani = 1.0 - 1.0 / mesh->aniso_factor_n[c0];
                if ( model->aniso_fstrain  == 1 ) ani = 1.0 - 1.0 / mesh->FS_AR_n[c0];
                d0   =  2.0*pow(nx, 2.0)*pow(nz, 2.0);
                d1   = nx*nz*(-pow(nx, 2.0) + pow(nz, 2.0));
            }
            //----------------------------------------------------------//
//            Exx = mesh->exxd[c0]  + mesh->sxxd0[c0] /etae/2.0;
//            Ezz = mesh->ezzd[c0]  + mesh->szzd0[c0] /etae/2.0;
//            Exz = mesh->exz_n[c0] + mesh->sxz0_n[c0]/etae/2.0;
//            gxz = 2.0*Exz;
            
            Da11  = 2.0 - 2.0*ani*d0;
            Da12  = 2.0*ani*d0;
            Da13  = 2.0*ani*d1;
            Da22  = 2.0 - 2.0*ani*d0;
            Da23  =-2.0*ani*d1;
            Da33  = 1.0  + 2.0*ani*(d0 - 0.5);
            a11   = Da33 * Da22 - pow(Da23,2);
            a12   = Da13 * Da23 - Da33 * Da12;
            a13   = Da12 * Da23 - Da13 * Da22;
            a22   = Da33 * Da11 - pow(Da13,2);
            a23   = Da12 * Da13 - Da11 * Da23;
            a33   = Da11 * Da22 - pow(Da12,2);
            det   = (Da11 * a11) + (Da12 * a12) + (Da13 * a13);
            iDa11 = a11/det; iDa12 = a12/det; iDa13 = a13/det;
            iDa22 = a22/det; iDa23 = a23/det;
            iDa33 = a33/det;
            
            double Exx_ref   = mesh->exxd[c0]      + (iDa11*mesh->sxxd0[c0] + iDa12*mesh->szzd0[c0] + iDa13*mesh->sxz0_n[c0])/etae;
            double Ezz_ref   = mesh->ezzd[c0]      + (iDa12*mesh->sxxd0[c0] + iDa22*mesh->szzd0[c0] + iDa23*mesh->sxz0_n[c0])/etae;
            double Exz_ref   = mesh->exz_n[c0]     + (iDa13*mesh->sxxd0[c0] + iDa23*mesh->szzd0[c0] + iDa33*mesh->sxz0_n[c0])/2.0/etae;
            double P_ref     = mesh->p_in[c0];
            double gxz_ref   = 2.0*Exz_ref;
            double pert_xx   = tol*fabs(Exx_ref) + tol/1e3;
            double pert_zz   = tol*fabs(Ezz_ref) + tol/1e3;
            double pert_xz   = tol*fabs(Exz_ref) + tol/1e3;
            double pert_p    = tol*fabs(P_ref)   + tol/1e3;
            
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
            
            // 1) Positive perturbation in Exx
            PhaseRheologyLoop(  1, 1.0, 2.0*pert_xx, Exx_ref+pert_xx, Ezz_ref, Exz_ref, P_ref, ani, d0, d1, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n,
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadexx_n, mesh->ddivpdexx_n, mesh->eta_n );
            
            // 2) Negative perturbation in Exx
            PhaseRheologyLoop( 1, -1.0, 2.0*pert_xx, Exx_ref-pert_xx, Ezz_ref, Exz_ref, P_ref, ani, d0, d1, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n,
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadexx_n, mesh->ddivpdexx_n, mesh->eta_n );
            
            //----------------------------------------------------------------------------------------------------------------------------------------------------//

            // 1) Positive perturbation in Ezz
            PhaseRheologyLoop( 1, 1.0, 2.0*pert_zz, Exx_ref, Ezz_ref+pert_zz, Exz_ref, P_ref, ani, d0, d1, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n,
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadezz_n, mesh->ddivpdezz_n, mesh->eta_n );

            // 2) Negative perturbation in Ezz
            PhaseRheologyLoop( 1, -1.0, 2.0*pert_zz, Exx_ref, Ezz_ref-pert_zz, Exz_ref, P_ref, ani, d0, d1, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n,
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadezz_n, mesh->ddivpdezz_n, mesh->eta_n );

            //----------------------------------------------------------------------------------------------------------------------------------------------------//

            // 1) Positive perturbation in Exz ---- NOTE THE FACTOR 2 due to Gxz = 2*Exz
            PhaseRheologyLoop(  1, 1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref+pert_xz, P_ref, ani, d0, d1, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n,
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadgxz_n, mesh->ddivpdgxz_n, mesh->eta_n );

            // 2) Negative perturbation in Exz ---- NOTE THE FACTOR 2 due to Gxz = 2*Exz
            PhaseRheologyLoop( 1, -1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref-pert_xz, P_ref, ani, d0, d1, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n,
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadgxz_n, mesh->ddivpdgxz_n, mesh->eta_n );

            //----------------------------------------------------------------------------------------------------------------------------------------------------//

            // 1) Positive perturbation in P
            PhaseRheologyLoop(  1, 1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref+pert_p, ani, d0, d1, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n,
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadp_n, mesh->ddivpdp_n, mesh->eta_n );

            // 2) Negative perturbation in P
            PhaseRheologyLoop( 1, -1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref-pert_p, ani, d0, d1, c0, mesh->phase_perc_n,
                              mesh->mu_n, mesh->T, mesh->p0_n,
                              mesh->d0_n, mesh->phi0_n, mesh->X0_n,
                              mesh->sxxd0, mesh->szzd0, mesh->sxz0_n,
                              mesh->bet_n, mesh->div_u,
                              mesh->strain_n, mesh->dil_n, mesh->fric_n, mesh->C_n,
                              model, materials, scaling,
                              mesh->detadp_n, mesh->ddivpdp_n, mesh->eta_n );
      
        }
    }

    #pragma omp parallel for shared( mesh ) private( cond, k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, etaVE, VEcoeff, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, exx_el, ezz_el, exz_el, exx_diss, ezz_diss, exz_diss, detadexx, detadezz, detadexz, detadp, Xreac, OverS, ddivpdexx, ddivpdezz, ddivpdexz, ddivpdp, Pcorr, drhodp, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, gxz, Gxx, Gzz, Gxz, el, etae, ani, d0, d1, nx, nz, Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det ) firstprivate( UnsplitDiffReac, materials, scaling, average, model, Nx, Nz, tol )
    for ( k1=0; k1<Nx*Nz; k1++ ) {

        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;

        // Initialise arrays to 0
        mesh->detadexx_s[c1]       = 0.0;
        mesh->detadezz_s[c1]       = 0.0;
        mesh->detadgxz_s[c1]       = 0.0;
        mesh->detadp_s[c1]         = 0.0;

        if ( mesh->BCg.type[c1] != 30 ) {

            //----------------------------------------------------------//
            if ( model->iselastic==1   ) etae      = model->dt*mesh->mu_s[c1];
            else           etae      = 1.0; // set to arbitrary value to avoid division by 0.0
            //----------------------------------------------------------//
            if ( model->aniso == 0 ) {
                ani = 0.0; d0   = 0.0; d1   = 0.0;
            }
            else {
                // Director
                nx = mesh->nx0_s[c1];
                nz = mesh->nz0_s[c1];
                // See Anisotropy_v2.ipynb
                if ( model->aniso_fstrain  == 0 ) ani = 1.0 - 1.0 / mesh->aniso_factor_s[c1];
                if ( model->aniso_fstrain  == 1 ) ani = 1.0 - 1.0 / mesh->FS_AR_s[c1];
                d0   =  2.0*pow(nx, 2.0)*pow(nz, 2.0);
                d1   = nx*nz*(-pow(nx, 2.0) + pow(nz, 2.0));
            }
            //----------------------------------------------------------//
            //            Exx = mesh->exxd_s[c1] + mesh->sxxd0_s[c1]/etae/2.0;
            //            Ezz = mesh->ezzd_s[c1] + mesh->szzd0_s[c1]/etae/2.0;
            //            Exz = mesh->exz[c1]    + mesh->sxz0[c1]   /etae/2.0;
            //            gxz = 2.0*Exz;

            Da11  = 2.0 - 2.0*ani*d0;
            Da12  = 2.0*ani*d0;
            Da13  = 2.0*ani*d1;
            Da22  = 2.0 - 2.0*ani*d0;
            Da23  =-2.0*ani*d1;
            Da33  = 1.0  + 2.0*ani*(d0 - 0.5);
            a11   = Da33 * Da22 - pow(Da23,2);
            a12   = Da13 * Da23 - Da33 * Da12;
            a13   = Da12 * Da23 - Da13 * Da22;
            a22   = Da33 * Da11 - pow(Da13,2);
            a23   = Da12 * Da13 - Da11 * Da23;
            a33   = Da11 * Da22 - pow(Da12,2);
            det   = (Da11 * a11) + (Da12 * a12) + (Da13 * a13);
            iDa11 = a11/det; iDa12 = a12/det; iDa13 = a13/det;
            iDa22 = a22/det; iDa23 = a23/det;
            iDa33 = a33/det;
            double Exx_ref = mesh->exxd_s[c1]  + (iDa11*mesh->sxxd0_s[c1] + iDa12*mesh->szzd0_s[c1] + iDa13*mesh->sxz0[c1])/etae;
            double Ezz_ref = mesh->ezzd_s[c1]  + (iDa12*mesh->sxxd0_s[c1] + iDa22*mesh->szzd0_s[c1] + iDa23*mesh->sxz0[c1])/etae;
            double Exz_ref = mesh->exz[c1]     + (iDa13*mesh->sxxd0_s[c1] + iDa23*mesh->szzd0_s[c1] + iDa33*mesh->sxz0[c1])/2.0/etae;
            double P_ref     = mesh->P_s[c1];
            double gxz_ref   = 2.0*Exz_ref;
            double pert_xx   = tol*fabs(Exx_ref) + tol/1e3;
            double pert_zz   = tol*fabs(Ezz_ref) + tol/1e3;
            double pert_xz   = tol*fabs(Exz_ref) + tol/1e3;
            double pert_p    = tol*fabs(P_ref)   + tol/1e3;

            //----------------------------------------------------------------------------------------------------------------------------------------------------//

            // 1) Positive perturbation in Exx
            PhaseRheologyLoop(  0, 1.0, 2.0*pert_xx, Exx_ref+pert_xx, Ezz_ref, Exz_ref, P_ref, ani, d0, d1, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadexx_s, NULL, mesh->eta_s );

            // 2) Negative perturbation in Exx
            PhaseRheologyLoop( 0, -1.0, 2.0*pert_xx, Exx_ref-pert_xx, Ezz_ref, Exz_ref, P_ref, ani, d0, d1, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadexx_s, NULL, mesh->eta_s );
            
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
            
            // 1) Positive perturbation in Ezz
            PhaseRheologyLoop(  0, 1.0, 2.0*pert_zz, Exx_ref, Ezz_ref+pert_zz, Exz_ref, P_ref, ani, d0, d1, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadezz_s, NULL, mesh->eta_s );
            
            // 2) Negative perturbation in Ezz
            PhaseRheologyLoop( 0, -1.0, 2.0*pert_zz, Exx_ref, Ezz_ref-pert_zz, Exz_ref, P_ref, ani, d0, d1, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadezz_s, NULL, mesh->eta_s );
            
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
            
            // 1) Positive perturbation in Exz
            PhaseRheologyLoop(  0, 1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref+pert_xz, P_ref, ani, d0, d1, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadgxz_s, NULL, mesh->eta_s );
        
            // 2) Negative perturbation in Exz
            PhaseRheologyLoop( 0, -1.0, 4.0*pert_xz, Exx_ref, Ezz_ref, Exz_ref-pert_xz, P_ref, ani, d0, d1, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadgxz_s, NULL, mesh->eta_s );
            
            //----------------------------------------------------------------------------------------------------------------------------------------------------//
            // 1) Positive perturbation in P
            PhaseRheologyLoop(  0, 1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref+pert_p, ani, d0, d1, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadp_s, NULL, mesh->eta_s );
            
            // 2) Negative perturbation in Ezz
            PhaseRheologyLoop( 0, -1.0, 2.0*pert_p, Exx_ref, Ezz_ref, Exz_ref, P_ref-pert_p, ani, d0, d1, c1, mesh->phase_perc_s,
                              mesh->mu_s, mesh->T_s, mesh->p0_s,
                              mesh->d0_s, mesh->phi0_s, mesh->X0_s,
                              mesh->sxxd0_s, mesh->szzd0_s, mesh->sxz0,
                              mesh->bet_s, mesh->div_u_s,
                              mesh->strain_s, mesh->dil_s, mesh->fric_s, mesh->C_s,
                              model, materials, scaling,
                              mesh->detadp_s, NULL, mesh->eta_s );


        }

    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
