function [ eta_0, Txx, Tyy, Txy, exx_el, eyy_el, exy_el ] = M2Di2_LocalIteration_VEP_TM_2( litmax, tol, noisy, Exx, Eyy, Exy, Txxo, Tyyo, Txyo, T, T0, m_pwl, eta_el, P, phi, C )

% Strain rate components and invariant
exx    = Exx;
eyy    = Eyy;
exy    = Exy;
eII    = sqrt(1/2*(exx.^2 + eyy.^2) + exy.^2);

B_pwl  = exp( -T.*(1./(1 + T./T0)) );

% Viscoelastic strain rate components and invariant
Exx    = Exx + Txxo./2./eta_el;
Eyy    = Eyy + Tyyo./2./eta_el;
Exy    = Exy + Txyo./2./eta_el;
Eii    = sqrt(1/2*(Exx.^2 + Eyy.^2) + Exy.^2);
%eta_pl  = (Cref(h) + P0ref*sin(phiref(h))) / (2*Eii);

% Use n instead of m
n_pwl   = 1./(2*m_pwl+1);
C_pwl   = (2.*B_pwl).^(-n_pwl);

% Isolated viscosities
eta_pwl = B_pwl.*eII.^(1./n_pwl-1);
eta_pl  = (C + P.*sin(phi) ) ./ (2*Eii);

% Define viscosity bounds
eta_up = eta_pwl;
eta_up = min(eta_up, eta_el);
eta_lo = 1./( 1./eta_pwl  + 1./eta_el );

% Start at midpoint (initial guess)
eta_0     = 0.5*(eta_up+eta_lo);

if noisy>1, fprintf('Starting local iterations:\n'); end

% Function evaluation at plastic effective viscosity
Tii       = 2*eta_pl.*Eii;
Eii_pwl_0 = C_pwl.*Tii.^(n_pwl);
Eii_vis_0 = Eii_pwl_0;
r_eta_0   = Eii - Tii./(2*eta_el) - Eii_vis_0;

pl = r_eta_0 >= 0;
eta_0(pl==1) = eta_pl(pl==1);

% figure(9), clf, imagesc(pl),colorbar, drawnow

% Local iterations
for lit=1:litmax
    
    % Function evaluation at current effective viscosity
    Tii(pl==0)       = 2*eta_0(pl==0).*Eii(pl==0);
    Eii_pwl_0(pl==0) = C_pwl(pl==0).*Tii(pl==0).^(n_pwl(pl==0));  
    Eii_vis_0(pl==0) = Eii_pwl_0(pl==0);
    r_eta_0(pl==0)   = Eii(pl==0) - Tii(pl==0)./(2*eta_el(pl==0)) - Eii_vis_0(pl==0);
    
    % Residual check
    
    res = max(max(abs(r_eta_0(pl==0))));
    if lit == 1; res0 = res; end
    if noisy>1, fprintf('It. %02d, r = %2.2e\n', lit, res/res0); end
    if res/res0 < tol
        break;
    end
    
    % Exact derivative
    drdeta = -2*Eii./(2*eta_el) - C_pwl .* Tii.^n_pwl.*n_pwl./eta_0;
    
    % Update viscosity
    eta_0(pl==0)  = eta_0(pl==0) - r_eta_0(pl==0) ./ drdeta(pl==0);
    
end

% Update individual deviatoric stress components
Txx = 2*eta_0.*Exx;
Tyy = 2*eta_0.*Eyy;
Txy = 2*eta_0.*Exy;
% Eii_pwl2 = Eii_pwl_0.^2;
% eta_pwl = B_pwl.*Eii_pwl_0.^(1./n_pwl-1);
% eta_pwl1 = B_pwl.*Eii_pwl_0.^(2*m_pwl);
% eta_pwl2 = (2.^(-2.*m_pwl).*B_pwl.*Tii.^(2.*m_pwl)) .^(1./(2.*m_pwl+1));
% diffeta = eta_0 - (1./eta_el + 1./eta_pwl).^(-1);
% 
% % Strain rates
exx_el  = (Txx - Txxo) ./ (2*eta_el);
eyy_el  = (Tyy - Tyyo) ./ (2*eta_el);
exy_el  = (Txy - Txyo) ./ (2*eta_el);
% exx_pwl = Eii_pwl_0.*exx./eII;
% eyy_pwl = Eii_pwl_0.*eyy./eII;
% exy_pwl = Eii_pwl_0.*exy./eII;
% exx_net = exx-exx_el-exx_pwl;
% eyy_net = eyy-eyy_el-eyy_pwl;
% exy_net = exy-exy_el-exy_pwl;
% 
% exx_pwl1 = exx-exx_el;
% eyy_pwl1 = eyy-eyy_el;
% exy_pwl1 = exy-exy_el;
% 
% eII_pwl1  = sqrt(1/2*(exx_pwl1.^2 + eyy_pwl1.^2) + exy_pwl1.^2);
% 
% Eii2 = Eii.^2;
% dexx_pwldexx     = (1  - eta_0./eta_el);
% dexx_pwldexx     = 2.*Eii_pwl_0./eII - 4.*exx.^2.*Eii_pwl_0./eII.^(3/2);
% deta_pwldexx_pwl = exx_pwl1.*m_pwl.*B_pwl.*Eii_pwl2.^(m_pwl-1);
% deta_pwldexx     = dexx_pwldexx .* deta_pwldexx_pwl;
% 
% deta_pwldexx    = -2*eta_pwl .*m_pwl.*exx_pwl1.*eta_el.*(-eta_el + eta_0) ./ (eta_el.^2.*exx_pwl1.^2+eta_el.^2.*eyy_pwl1.^2+2.*eta_el.^2.*exy_pwl1.^2+ 2.*eta_pwl.* m_pwl.*exx_pwl1.*exx.*eta_el + eta_pwl .*m_pwl.*exx_pwl1.*Txxo);
% 
% deta_0deta_pwl   = eta_0.^2 .* eta_pwl.^(-2);
% 
% 
% Eii_pwl2         = Eii_pwl_0.^2;
% deta_pwldexx_pwl1 = eta_pwl.*m_pwl.*exx_pwl1 ./ Eii_pwl2;
% detadExx         = -2.*deta_pwldexx_pwl1.*(-eta_el+eta_0).*eta_el./(2.*1./(eta_0).^2.*eta_pwl.^2.*eta_el.^2+2.*deta_pwldexx_pwl1.*exx.*eta_el+deta_pwldexx_pwl1.*Txxo);
% 
% detadExx0        = -2.*eta_pwl.*m_pwl.*exx_pwl1.*eta_el.*(-eta_el+eta_0)./(eta_el.^2.*exx_pwl1.^2+eta_el.^2.*eyy_pwl1.^2+2.*eta_el.^2.*exy_pwl1.^2+2.*eta_pwl.*eta_el.*exx_pwl1.^2+2.*eta_pwl.*eta_el.*eyy_pwl1.^2+4.*eta_pwl.*eta_el.*exy_pwl1.^2+eta_pwl.^2.*exx_pwl1.^2+eta_pwl.^2.*eyy_pwl1.^2+2.*eta_pwl.^2.*exy_pwl1.^2+2.*eta_pwl.*m_pwl.*exx_pwl1.*exx.*eta_el+eta_pwl.*m_pwl.*exx_pwl1.*Txxo);
% 
% detadExx1         = eta_0.^2 .* (deta_pwldexx./eta_pwl.^2);
% 
% detadExx2         = (deta_pwldexx.*eta_el.^2) ./ (eta_pwl + eta_el).^2;
% 
% detadExx3         = deta_0deta_pwl .* deta_pwldexx;
% 
%  detadExx = 0;%  (m_pwl.*(exx_pwl1)) .* eta_0.^2 ./ ( eta_pwl.*Eii_pwl2);
%  detadEyy = 0;%  (m_pwl.*(eyy-eyy_el)) .* eta_0.^2 ./ ( eta_pwl.*Eii_pwl2);
%  detadExy = 0;%2*(m_pwl.*(exy-exy_el)) .* eta_0.^2 ./ ( eta_pwl.*Eii_pwl2);


% figure(89), clf
% subplot(311)
% imagesc(eta_pwl), colorbar
% subplot(312)
% imagesc(eta_pwl1), colorbar
% subplot(313)
% imagesc(eta_pwl2), colorbar
% drawnow

% figure(89), clf
% subplot(211)
% imagesc(eta_0), colorbar
% subplot(212)
% imagesc(diffeta), colorbar
% drawnow

% figure(90), clf
% subplot(311)
% imagesc(exx_net), colorbar
% subplot(312)
% imagesc(eyy_net), colorbar
% subplot(313)
% imagesc(eII_pwl1-Eii_pwl_0), colorbar
% drawnow

end