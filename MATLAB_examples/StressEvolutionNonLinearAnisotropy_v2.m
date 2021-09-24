clear 
% close all
clc

nonlinear = 1; % 0: power-law exponent is 1.0 - 1: larger than 1.0 
two       = 1; % 1: uses engineering strain - 2 is standard form

% Constants
R   = 8.314;

% Conditions 
exx = 1e-14;
eyy = -exx;
ezz = 0;    % incompressible
exy = exx/2;
eii = sqrt(1/2*(exx^2+eyy^2+ezz^2) + exy^2);
T   = 773;

% Numerics
dt     = 1e11;
nt     = 2;
litmax = 10;
noisy  = 3;
tol    = 1e-12;

% Mat props
aniso_fact = 100;
angle      = 65*pi/180;
G          = 2e10;
Q          = 276e3;
if nonlinear == 1
    A      = 3.2e-20;
    n      = 3.0;
else
    A      = 2*(2.237e8)^(-1)*exp(Q/R/T)*eii; % same flow stress as non-linear model
    n      = 1.0;
end

% Precompute
B    = A^(-1/n)*exp(Q/n/R/T);
C    = (2*B)^-n;

% Isolated viscosity
eta_pwl  = B*eii^(1/n-1);
eta_el   = G*dt;

% Anisotropy
ani      = (1 - 1 / aniso_fact);
nx       = cos(angle); ny = sin(angle);
Director = [nx, ny];
a0       = 2*Director(1)^2*Director(2)^2;
a1       = Director(1)*Director(2) * (-Director(1)^2 + Director(2)^2);
C_ANI    = [-a0 a0 two*a1; a0 -a0 -two*a1; a1 -a1 -two*(1/2+a0)];
C_ISO    = [1 0 0; 0 1 0; 0 0 two*1/2];
Dani     = (2* C_ISO + 2 * C_ANI * ani); 

% Isotropic matrix
Diso = [2.0 0 0;...
0 2 0;...
0 0 two*1];

% Initial stresses
Txx = 0; Tyy = 0; Tzz = 0; Txy = 0; time = 0;
Txx_iso = 0; Tyy_iso = 0; Tzz_iso = 0; Txy_iso = 0;
T = zeros(nt,1); T_iso = zeros(nt,1); t = zeros(nt,1);

for it=1:nt
    % Old guys
    Txx0 = Txx; Tyy0 = Tyy; Tzz0 = Tzz; Txy0 = Txy;

    % New effective strain rate
    Exx  = exx + Txx0 / 2 / eta_el;
    Eyy  = eyy + Tyy0 / 2 / eta_el;
    Ezz  = ezz + Tzz0 / 2 / eta_el;
    Exy  = exy + Txy0 / 2 / eta_el;
    Gxy  = 2*Exy;
    Eii  = sqrt(1/2*(Exx^2+Eyy^2+Ezz^2) + Exy^2);

    % Anisotropic viscosity matrix
    Dv = eta_pwl * Dani;

    % Anisotropic elastic matrix
    De = G * Dani * dt;

    % Quasi harmonic mean
    Dve = inv( inv(Dv) + inv(De) );
    
    E = [Exx; Eyy; (1/two)*Gxy ];
    Tau = Dve*E;
    Txx = Tau(1);
    Tyy = Tau(2);
    Txy = Tau(3);
    Tzz = 0;
    
    % so far
    Tii  = sqrt(1/2*(Txx^2+Tyy^2+Tzz^2) + Txy^2);
 
    %%%%%%%%%%%%%%%%% CURRENT MDOODZ6.0
    
    % No local iteration
    eta_ve = (1/eta_pwl+1/eta_el)^(-1);
    
    % New deviatoric stress
    Txx0 = Txx_iso; Tyy0 = Tyy_iso; Tzz0 = Tzz_iso; Txy0 = Txy_iso;
    Exx  = exx + Txx0 / 2 / eta_el;
    Eyy  = eyy + Tyy0 / 2 / eta_el;
    Ezz  = ezz + Tzz0 / 2 / eta_el;
    Exy  = exy + Txy0 / 2 / eta_el;
    Gxy  = 2*Exy;
    Txx_iso = 2*eta_ve*Exx;
    Tyy_iso = 2*eta_ve*Eyy;
    Tzz_iso = 2*eta_ve*Ezz;
    Txy_iso = 2*eta_ve*Exy;
    Tii_iso = 2*eta_ve*Eii;
    
    % Sanity check
    Tii_iso  = sqrt(1/2*(Txx_iso^2+Tyy_iso^2+Tzz_iso^2) + Txy_iso^2);
    norm(Tii_iso-Tii)
    %%%%%%%%%%%%%%%%% CURRENT MDOODZ6.0
    
    time  = time + dt;
    T_iso(it) = Tii_iso;
    T(it) = Tii;
    t(it) = time;
    
end

% Visualisation step #1
figure(1), clf
hold on
plot(t, 2*eta_pwl*eii*ones(size(T)), '-.k')
plot(t, T, '+k')
plot(t, T_iso, '-r')

%%%%%%%%%%%%%%%%%%%%%% NOW WITH LOCAL ITERATIONS

% Initial stresses
Txx = 0; Tyy = 0; Tzz = 0; Txy = 0; time = 0;
Txx_iso = 0; Tyy_iso = 0; Tzz_iso = 0; Txy_iso = 0;
T = zeros(nt,1); T_iso = zeros(nt,1); t = zeros(nt,1);

for it=1:nt
    % Old guys
    Txx0 = Txx; Tyy0 = Tyy; Tzz0 = Tzz; Txy0 = Txy;
    Tau0 = [Txx0; Tyy0; Txy0];
    e    = [exx; eyy; (1/two)*(2*exy) ];
    
    % Anisotropic elastic matrix
    De = G * Dani * dt;
      
    % New effective strain rate
    E = e + inv(De) * Tau0;
    
    % Only works for isotropic 
%     Exx  = exx + Txx0 / 2 / eta_el;
%     Eyy  = eyy + Tyy0 / 2 / eta_el;
%     Ezz  = ezz + Tzz0 / 2 / eta_el;
%     Exy  = exy + Txy0 / 2 / eta_el;
    Gxy  = 2*Exy;
    Eii  = sqrt(1/2*(Exx^2+Eyy^2+Ezz^2) + Exy^2);
%     E = [Exx; Eyy; (1/two)*Gxy ];
    
    figure(2), clf, hold on
    for iter =1:30
        
        % Stress vector
        Tau = [Txx; Tyy; Txy];
        
        % Anisotropic viscosity matrix
        Dv = eta_pwl * Dani;
                 
%         % Quasi harmonic mean
%         Dve = inv( inv(Dv) + inv(De) );
        
        % Superb simplification, if same viscous and elastic anisotropy:
        Dve = eta_ve .* Dani;
        
        % Additive strain decomposition
        fr    = E - ( inv(Dv) + inv(De) )* Tau;
        fr   = e - inv(Dv)* Tau - inv(De)*(Tau-Tau0);
        plot(iter, log10(norm(fr)), 'ok')

        Tau = Dve*E;
        Txx = Tau(1);
        Tyy = Tau(2);
        Txy = Tau(3);
        Tzz = 0;
        
        % so far
        Tii       = sqrt(1/2*(Txx^2+Tyy^2+Tzz^2) + 1/aniso_fact^2*Txy^2);
        Eii_pwl_0 = C.*Tii.^(n);
        eta_pwl   = B* Eii_pwl_0^(1/n-1);
        eta_ve    = (1/eta_pwl+1/eta_el)^(-1);
        
    end
    
    % Check decomposition for components and invariants
   
    
    ee = inv(De)*(Tau-Tau0);
    ev = inv(Dv)*Tau;
    e - ee - ev
    
    exxv = ev(1); eyyv = ev(2); ezzv = -exxv-eyyv; exyv = (two/2)*ev(3);
    eiiv = sqrt(1/2*(exxv^2+eyyv^2+ezzv^2) + 1/aniso_fact^2*exyv^2);
    exxe = ee(1); eyye = ee(2); ezze = -exxe-eyye; exye = (two/2)*ee(3);
    eii = sqrt(1/2*(exxe^2+eyye^2+ezze^2) + exye^2);
    
%     exxe = (Txx-Txx0)/2/eta_e;
%     exxv = Eii_pwl_0*Txx/Tii;
%     eyye = (Tyy-Tyy0)/2/eta_e;
%     eyyv = Eii_pwl_0*Tyy/Tii;
%     exye = (Txy-Txy0)/2/eta_e;
%     exyv = Eii_pwl_0*Txy/Tii;
%     ezzv = -exxv-eyyv;
%     ezze = -exxe-eyye;
%     eiiv = sqrt(1/2*(exxv^2+eyyv^2+ezzv^2) + exyv^2);
%     eiie = sqrt(1/2*(exxe^2+eyye^2+ezze^2) + exye^2);
%    
%     exxe1 = eiie*Txx/Tii;
%     
%     Eii - Tii./(2*eta_e) - Eii_pwl_0
    Eii_pwl_0  - eiiv
%     exx - exxe - exxv
%     eyy - eyye - eyyv
%     exy - exye - exyv
%     eii - eiie - eiiv
%     (eii - eiie) - (Eii - Tii./(2*eta_e))
%     exxe - exxe1
    
    T(it) = Tii;
    
    %%%%%%%%%%%%%%%%% CURRENT MDOODZ6.0
    Txx0 = Txx_iso; Tyy0 = Tyy_iso; Tzz0 = Tzz_iso; Txy0 = Txy_iso;
    Exx  = exx + Txx0 / 2 / eta_el;
    Eyy  = eyy + Tyy0 / 2 / eta_el;
    Ezz  = ezz + Tzz0 / 2 / eta_el;
    Exy  = exy + Txy0 / 2 / eta_el;
    Gxy  = 2*Exy;
    Eii  = sqrt(1/2*(Exx^2+Eyy^2+Ezz^2) + Exy^2);
    eta_pwl   = B* eii^(1/n-1);
    
    % with local iteration
    eta_ve = (1/eta_pwl+1/eta_el)^(-1);
    
    % Local iterations
    for lit=1:litmax
        
        % Function evaluation at current effective viscosity
        Tii       = 2*eta_ve.*Eii;
        Eii_pwl_0 = C.*Tii.^(n);
        r_eta_0   = Eii - Tii./(2*eta_el) - Eii_pwl_0;
        
        % Residual check
        res = max(max(abs(r_eta_0)));
        if lit == 1; res0 = res; end
        if noisy>2, fprintf('It. %02d, r abs. = %2.2e r rel. = %2.2e\n', lit, res, res/res0); end
        if res/res0 < tol || res/eii<tol
            break;
        end
        
        % Exact derivative
        drdeta = -2*Eii./(2*eta_el) - C .* Tii.^n.*n./eta_ve;
        
        % Update viscosity
        eta_ve  = eta_ve - r_eta_0 ./ drdeta;
        
    end

    % New deviatoric stress
    Txx_iso = 2*eta_ve*Exx;
    Tyy_iso = 2*eta_ve*Eyy;
    Tzz_iso = 2*eta_ve*Ezz;
    Txy_iso = 2*eta_ve*Exy;
    Tii_iso = 2*eta_ve*Eii;
    %%%%%%%%%%%%%%%%% CURRENT MDOODZ6.0
    
    % Sanity check
    Tii_iso  = sqrt(1/2*(Txx_iso^2+Tyy_iso^2+Tzz_iso^2) + Txy_iso^2);
    norm(Tii_iso-Tii);
    %%%%%%%%%%%%%%%%%
    
    time  = time + dt;
    T_iso(it) = Tii_iso;
    t(it) = time;
    
end

% Visualisation step #2
figure(1)
plot(t, T, '-dr')
plot(t, T_iso, '-g')
ylabel('$\tau_\mathrm{II}$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
legend('Flow stress', 'No iterations', 'Local iterations', 'location', 'southeast')
set(gca, 'fontsize', 18)