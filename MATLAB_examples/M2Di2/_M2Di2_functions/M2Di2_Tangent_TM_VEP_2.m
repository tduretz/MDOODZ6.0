function [ Dv ] = M2Di2_Tangent_TM_VEP_2( etac, etav, litmax, littol, noisy, Exxc, Eyyc, Exyc, Txxco, Tyyco, Txyco, Tc, Exxv, Eyyv, Exyv, Txxvo, Tyyvo, Txyvo, Tv, T0, mc, mv, etaec, etaev, Ptc, Ptv, phic, phiv, Cc, Cv, Eiic2, Eiiv2, Exxc_el, Eyyc_el, Exyc_el, Exxv_el, Eyyv_el, Exyv_el  )
% Engineering convention
Gxyc    = 2*Exyc;
Gxyv    = 2*Exyv;
Gxyc_el = 2*Exyc_el;
Gxyv_el = 2*Exyv_el;
%% Derivative of viscosity versus solution fields
% Perturbation - centroids
toljac = 1e-5;
Eiic = sqrt(Eiic2);
dExx = toljac*Eiic;
dEyy = toljac*Eiic;
dExy = toljac*Eiic;
dP   = toljac*Ptc;
dT   = toljac*Tc;
[ etac_xx, ~, ~, ~, Exxc_el_xx, Eyyc_el_xx, Exyc_el_xx ] = M2Di2_LocalIteration_VEP_TM_2( litmax, littol, noisy, Exxc+dExx, Eyyc, Exyc, Txxco, Tyyco, Txyco, Tc, T0, mc, etaec, Ptc, phic, Cc );
[ etac_yy, ~, ~, ~, Exxc_el_yy, Eyyc_el_yy, Exyc_el_yy ] = M2Di2_LocalIteration_VEP_TM_2( litmax, littol, noisy, Exxc,      Eyyc+dEyy, Exyc, Txxco, Tyyco, Txyco, Tc, T0, mc, etaec, Ptc, phic, Cc );
[ etac_xy, ~, ~, ~, Exxc_el_xy, Eyyc_el_xy, Exyc_el_xy ] = M2Di2_LocalIteration_VEP_TM_2( litmax, littol, noisy, Exxc,      Eyyc, Exyc+dExy, Txxco, Tyyco, Txyco, Tc, T0, mc, etaec, Ptc, phic, Cc );
[ etac_p , ~, ~, ~, Exxc_el_p , Eyyc_el_p , Exyc_el_p  ] = M2Di2_LocalIteration_VEP_TM_2( litmax, littol, noisy, Exxc,      Eyyc, Exyc, Txxco, Tyyco, Txyco, Tc, T0, mc, etaec, Ptc+dP, phic, Cc );
[ etac_T , ~, ~, ~, Exxc_el_T , Eyyc_el_T , Exyc_el_T  ] = M2Di2_LocalIteration_VEP_TM_2( litmax, littol, noisy, Exxc,      Eyyc, Exyc, Txxco, Tyyco, Txyco, Tc+dT, T0, mc, etaec, Ptc, phic, Cc );
% Viscosity derivatives
detadexxc = (etac_xx - etac) ./ dExx;
detadeyyc = (etac_yy - etac) ./ dEyy;
detadgxyc = 1/2*(etac_xy - etac) ./ dExy;
detadpc   = (etac_p  - etac) ./ dP;
detadTc   = (etac_T  - etac) ./ dT;
% Elastic strain rate derivatives
dexxeldexxc = (Exxc_el_xx - Exxc_el) ./ dExx;
deyyeldexxc = (Eyyc_el_xx - Eyyc_el) ./ dExx;
dgxyeldexxc = (Exyc_el_xx - Exyc_el) ./ dExx;
dexxeldeyyc = (Exxc_el_yy - Exxc_el) ./ dEyy;
deyyeldeyyc = (Eyyc_el_yy - Eyyc_el) ./ dEyy;
dgxyeldeyyc = (Exyc_el_yy - Exyc_el) ./ dEyy;
dexxeldgxyc = 1/2*(Exxc_el_xy - Exxc_el) ./ dExy;
deyyeldgxyc = 1/2*(Eyyc_el_xy - Eyyc_el) ./ dExy;
dgxyeldgxyc = 1/2*(Exyc_el_xy - Exyc_el) ./ dExy;
dexxeldTc = (Exxc_el_T - Exxc_el) ./ dT;
deyyeldTc = (Eyyc_el_T - Eyyc_el) ./ dT;
dgxyeldTc = (Exyc_el_T - Exyc_el) ./ dT;
% Perturbation - vertices
Eiiv = sqrt(Eiiv2);
dExx = toljac*Eiiv;
dEyy = toljac*Eiiv;
dExy = toljac*Eiiv;
dP   = toljac*Ptv;
dT   = toljac*Tv;
[ etav_xx, ~, ~, ~, Exxv_el_xx, Eyyv_el_xx, Exyv_el_xx ] = M2Di2_LocalIteration_VEP_TM_2( litmax, littol, noisy, Exxv+dExx, Eyyv, Exyv, Txxvo, Tyyvo, Txyvo, Tv, T0, mv, etaev, Ptv, phiv, Cv );
[ etav_yy, ~, ~, ~, Exxv_el_yy, Eyyv_el_yy, Exyv_el_yy ] = M2Di2_LocalIteration_VEP_TM_2( litmax, littol, noisy, Exxv, Eyyv+dEyy, Exyv, Txxvo, Tyyvo, Txyvo, Tv, T0, mv, etaev, Ptv, phiv, Cv );
[ etav_xy, ~, ~, ~, Exxv_el_xy, Eyyv_el_xy, Exyv_el_xy ] = M2Di2_LocalIteration_VEP_TM_2( litmax, littol, noisy, Exxv, Eyyv, Exyv+dExy, Txxvo, Tyyvo, Txyvo, Tv, T0, mv, etaev, Ptv, phiv, Cv );
[ etav_p , ~, ~, ~, Exxv_el_p, Eyyv_el_p, Exyv_el_p ] = M2Di2_LocalIteration_VEP_TM_2( litmax, littol, noisy, Exxv, Eyyv, Exyv, Txxvo, Tyyvo, Txyvo, Tv, T0, mv, etaev, Ptv+dP, phiv, Cv );
[ etav_T , ~, ~, ~, Exxv_el_T, Eyyv_el_T, Exyv_el_T ] = M2Di2_LocalIteration_VEP_TM_2( litmax, littol, noisy, Exxv, Eyyv, Exyv, Txxvo, Tyyvo, Txyvo, Tv+dT, T0, mv, etaev, Ptv, phiv, Cv );
% Viscosity derivatives
detadexxv = (etav_xx - etav) ./ dExx;
detadeyyv = (etav_yy - etav) ./ dEyy;
detadgxyv = 1/2*(etav_xy - etav) ./ dExy;
detadpv   = (etav_p  - etav) ./ dP;
detadTv   = (etav_T  - etav) ./ dT;
% Elastic strain rate derivatives
dexxeldexxv = (Exxv_el_xx - Exxv_el) ./ dExx;
deyyeldexxv = (Eyyv_el_xx - Eyyv_el) ./ dExx;
dgxyeldexxv = (Exyv_el_xx - Exyv_el) ./ dExx;
dexxeldeyyv = (Exxv_el_yy - Exxv_el) ./ dEyy;
deyyeldeyyv = (Eyyv_el_yy - Eyyv_el) ./ dEyy;
dgxyeldeyyv = (Exyv_el_yy - Exyv_el) ./ dEyy;
dexxeldgxyv = 1/2*(Exxv_el_xy - Exxv_el) ./ dExy;
deyyeldgxyv = 1/2*(Eyyv_el_xy - Eyyv_el) ./ dExy;
dgxyeldgxyv = 1/2*(Exyv_el_xy - Exyv_el) ./ dExy;
dexxeldTv   = (Exxv_el_T - Exxv_el) ./ dT;
deyyeldTv   = (Eyyv_el_T - Eyyv_el) ./ dT;
dgxyeldTv   = (Exyv_el_T - Exyv_el) ./ dT;
%% Rheological coefficients for the Jacobian
Exxc1   = Exxc + Txxco./2./etaec;
Eyyc1   = Eyyc + Tyyco./2./etaec;
Gxyc1   = Gxyc + Txyco./1./etaec;
Exxv1   = Exxv + Txxvo./2./etaev;
Eyyv1   = Eyyv + Tyyvo./2./etaev;
Gxyv1   = Gxyv + Txyvo./1./etaev;
% Centroids
Dv.D11c = 2*etac +   2*detadexxc.*Exxc1;       Dv.D12c =            2*detadeyyc.*Exxc1;       Dv.D13c =           2*detadgxyc.*Exxc1;       Dv.D14c = 2*detadTc.*Exxc1;
Dv.D21c =            2*detadexxc.*Eyyc1;       Dv.D22c = 2*etac +   2*detadeyyc.*Eyyc1;       Dv.D23c =           2*detadgxyc.*Eyyc1;       Dv.D24c = 2*detadTc.*Eyyc1;
Dv.D31c =            1*detadexxc.*Gxyc1;       Dv.D32c =            1*detadeyyc.*Gxyc1;       Dv.D33c = 1*etac +  1*detadgxyc.*Gxyc1;       Dv.D34v = 1*detadTc.*Gxyc1;
Dv.D41c = 2.*detadexxc.*Exxc1.*(Exxc-Exxc_el) + 2.*etac.*(Exxc-Exxc_el) + 2.*etac.*Exxc1.*(1-dexxeldexxc) + 2.*detadexxc.*Eyyc1.*(Eyyc-Eyyc_el) - 2.*etac.*Eyyc1.*deyyeldexxc + detadexxc.*Gxyc1.*(Gxyc-Gxyc_el)-etac.*Gxyc1.*dgxyeldexxc;
Dv.D42c = 2.*detadeyyc.*Exxc1.*(Exxc-Exxc_el) - 2.*etac.*Exxc1.*dexxeldeyyc + 2.*detadeyyc.*Eyyc1.*(Eyyc-Eyyc_el) + 2.*etac.*(Eyyc-Eyyc_el) + 2.*etac.*Eyyc1.*(1-deyyeldeyyc) + detadeyyc.*Gxyc1.*(Gxyc-Gxyc_el)-etac.*Gxyc1.*dgxyeldeyyc;
Dv.D43c = 2.*detadgxyc.*Exxc1.*(Exxc-Exxc_el) - 2.*etac.*Exxc1.*dexxeldgxyc + 2.*detadgxyc.*Eyyc1.*(Eyyc-Eyyc_el) - 2.*etac.*Eyyc1.*deyyeldgxyc + detadgxyc.*Gxyc1.*(Gxyc-Gxyc_el) + etac.*(Gxyc-Gxyc_el) + etac.*Gxyc1.*(1-dgxyeldgxyc);
Dv.D44c = 2.*detadTc.*Exxc1.*(Exxc-Exxc_el) - 2.*etac.*Exxc1.*dexxeldTc + 2.*detadTc.*Eyyc1.*(Eyyc-Eyyc_el) - 2.*etac.*Eyyc1.*deyyeldTc + detadTc.*Gxyc1.*(Gxyc-Gxyc_el) - etac.*Gxyc1.*dgxyeldTc;
% Vertices
Dv.D11v = 2*etav +   2*detadexxv.*Exxv1;       Dv.D12v =            2*detadeyyv.*Exxv1;       Dv.D13v =            2*detadgxyv.*Exxv1;       Dv.D14v = 2*detadTv.*Exxv1;
Dv.D21v =            2*detadexxv.*Eyyv1;       Dv.D22v = 2*etav +   2*detadeyyv.*Eyyv1;       Dv.D23v =            2*detadgxyv.*Eyyv1;       Dv.D24v = 2*detadTv.*Eyyv1;
Dv.D31v =            1*detadexxv.*Gxyv1;       Dv.D32v =            1*detadeyyv.*Gxyv1;       Dv.D33v = 1*etav +   1*detadgxyv.*Gxyv1;       Dv.D34v = 1*detadTv.*Gxyv1;
Dv.D41v = 2.*detadexxv.*Exxv1.*(Exxv-Exxv_el) + 2.*etav.*(Exxv-Exxv_el) + 2.*etav.*Exxv1.*(1-dexxeldexxv) + 2.*detadexxv.*Eyyv1.*(Eyyv-Eyyv_el) - 2.*etav.*Eyyv1.*deyyeldexxv + detadexxv.*Gxyv1.*(Gxyv-Gxyv_el)-etav.*Gxyv1.*dgxyeldexxv;
Dv.D42v = 2.*detadeyyv.*Exxv1.*(Exxv-Exxv_el) - 2.*etav.*Exxv1.*dexxeldeyyv + 2.*detadeyyv.*Eyyv1.*(Eyyv-Eyyv_el) + 2.*etav.*(Eyyv-Eyyv_el) + 2.*etav.*Eyyv1.*(1-deyyeldeyyv) + detadeyyv.*Gxyv1.*(Gxyv-Gxyv_el)-etav.*Gxyv1.*dgxyeldeyyv;
Dv.D43v = 2.*detadgxyv.*Exxv1.*(Exxv-Exxv_el) - 2.*etav.*Exxv1.*dexxeldgxyv + 2.*detadgxyv.*Eyyv1.*(Eyyv-Eyyv_el) - 2.*etav.*Eyyv1.*deyyeldgxyv + detadgxyv.*Gxyv1.*(Gxyv-Gxyv_el) + etav.*(Gxyv-Gxyv_el) + etav.*Gxyv1.*(1-dgxyeldgxyv);
Dv.D44v = 2.*detadTv.*Exxv1.*(Exxv-Exxv_el) - 2.*etav.*Exxv1.*dexxeldTv + 2.*detadTv.*Eyyv1.*(Eyyv-Eyyv_el) - 2.*etav.*Eyyv1.*deyyeldTv + detadTv.*Gxyv1.*(Gxyv-Gxyv_el) - etav.*Gxyv1.*dgxyeldTv;

end

