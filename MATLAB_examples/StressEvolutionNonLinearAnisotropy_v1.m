clear 
close all

nonlinear = 0; % 0: power-law exponent is 1.0 - 1: larger than 1.0 

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
nt     = 30;
litmax = 10;
noisy  = 3;
tol    = 1e-12;

% Mat props
aniso_fact = 1;
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
ani   = (1 - 1 / aniso_fact);
deta  = eta_pwl*ani;
nx    = cos(angle); ny = cos(angle);
d0   = 2.0*power(nx, 2.0)*power(ny, 2.0);
d1   = nx*ny*(-power(nx, 2.0) + power(ny, 2.0));

Dv = [2.0*eta_pwl -2.0*deta*d0 2.0*deta*d0;...
2.0*deta*d0 2.0*eta_pwl - 2.0*deta*d0 -2.0*deta*d1;...
2.0*deta*d1 -2.0*deta*d1 eta_pwl+2.0*deta*(d0 - 0.5);];
Dv_inv = inv(Dv);

Dani = [2.0 -2.0*ani*d0 2.0*ani*d0;...
2.0*ani*d0 2.0 - 2.0*ani*d0 -2.0*ani*d1;...
2.0*ani*d1 -2.0*ani*d1 1+2.0*ani*(d0 - 0.5);];
Dv2 = eta_pwl * Dani;


Dv_inv2 = 1/eta_pwl * inv(Dani); % Inverse of anisotropic viscosity matrix

Dve_inv = Dv_inv;

IGdt = [1/2/eta_el 0 0; 0 1/2/eta_el 0; 0 0 1/1/eta_el];

Dve_inv = Dv_inv + IGdt ;

% Dve_inv(1,1) = Dve_inv(1,1) + 1/2/eta_el;
% Dve_inv(2,2) = Dve_inv(2,2) + 1/2/eta_el;
% Dve_inv(3,3) = Dve_inv(3,3) + 1/1/eta_el;  % Inverse of anisotropic visco-elastic matrix

Dve = inv(Dve_inv)

% % inverse of matrix sum: inv(A + B) =  inv(A) - 1/(1+g) * inv(A)*B*inv(A)
g    = trace( IGdt*Dv ); %Dv = inv(A) et IGdt  = B
Dve2 = Dv - 1/(1+g) * Dv*IGdt*Dv


I = eye(3,3);

Dve2 = Dv - Dv*I*inv(inv(IGdt) + I*Dv*I)*I*Dv
Dve2 = Dv - Dv*inv(inv(IGdt) + Dv)*Dv

Dve2 = Dv - Dv*I*inv(I + IGdt*I*Dv*I) * IGdt*I*Dv



% Initial stresses
Txx = 0; Tyy = 0; Tzz = 0; Txy = 0; time = 0;
T = zeros(nt,1); t = zeros(nt,1);

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
    
    Eps = [Exx; Eyy; Gxy ];
%     Tau = Dve_inv \ Eps;
    Tau = Dve2 * Eps;
    Txx = Tau(1);
    Tyy = Tau(2);
    Txy = Tau(3);
    Tzz = 0;
    
    % so far
    Tii  = sqrt(1/2*(Txx^2+Tyy^2+Tzz^2) + Txy^2);
 
    % No local iteration
    eta_ve = (1/eta_pwl+1/eta_el)^(-1);
    
%     % New deviatoric stress
%     Txx1 = 2*eta_ve*Exx;
%     Tyy1 = 2*eta_ve*Eyy;
%     Tzz1 = 2*eta_ve*Ezz;
%     Txy1 = 2*eta_ve*Exy;
%     Tii1 = 2*eta_ve*Eii;
%     
%     % Sanity check
%     Tii_chk  = sqrt(1/2*(Txx^2+Tyy^2+Tzz^2) + Txy^2);
%     norm(Tii_chk-Tii)
%     
%     Txx1 - Txx
%     Tyy1 - Tyy
%     Txy1 - Txy
    
    
    time  = time + dt;
    T(it) = Tii;
    t(it) = time;
    
end

% Visualisation step #1
figure(1), clf
hold on
plot(t, 2*eta_pwl*eii*ones(size(T)), '-.k')
plot(t, T, '+k')

% %%%%%%%%%%%%%%%%%%%%%% NOW WITH LOCAL ITERATIONS
% 
% % Initial stresses
% Txx = 0; Tyy = 0; Tzz = 0; Txy = 0; time = 0;
% T = zeros(nt,1); t = zeros(nt,1);
% 
% for it=1:nt
%     % Old guys
%     Txx0 = Txx; Tyy0 = Tyy; Tzz0 = Tzz; Txy0 = Txy;
%     
%     % New effective strain rate
%     Exx  = exx + Txx0 / 2 / eta_el;
%     Eyy  = eyy + Tyy0 / 2 / eta_el;
%     Ezz  = ezz + Tzz0 / 2 / eta_el;
%     Exy  = exy + Txy0 / 2 / eta_el;
%     Eii  = sqrt(1/2*(Exx^2+Eyy^2+Ezz^2) + Exy^2);
%     
%     % Trial viscosity
%     eta_ve = (1/eta_pwl+1/eta_el)^(-1);
%     
%     % Local iterations
%     for lit=1:litmax
%         
%         % Function evaluation at current effective viscosity
%         Tii       = 2*eta_ve.*Eii;
%         Eii_pwl_0 = C.*Tii.^(n);
%         r_eta_0   = Eii - Tii./(2*eta_el) - Eii_pwl_0;
%         
%         % Residual check
%         res = max(max(abs(r_eta_0)));
%         if lit == 1; res0 = res; end
%         if noisy>2, fprintf('It. %02d, r abs. = %2.2e r rel. = %2.2e\n', lit, res, res/res0); end
%         if res/res0 < tol || res/eii<tol
%             break;
%         end
%         
%         % Exact derivative
%         drdeta = -2*Eii./(2*eta_el) - C .* Tii.^n.*n./eta_ve;
%         
%         % Update viscosity
%         eta_ve  = eta_ve - r_eta_0 ./ drdeta;
%         
%     end
%     
%     % New deviatoric stress
%     Txx = 2*eta_ve*Exx;
%     Tyy = 2*eta_ve*Eyy;
%     Tzz = 2*eta_ve*Ezz;
%     Txy = 2*eta_ve*Exy;
%     Tii = 2*eta_ve*Eii;
%     
%     % Sanity check
%     Tii_chk  = sqrt(1/2*(Txx^2+Tyy^2+Tzz^2) + Txy^2);
%     norm(Tii_chk-Tii)
%     
%     time  = time + dt;
%     T(it) = Tii;
%     t(it) = time;
%     
% end
% 
% % Visualisation step #2
% plot(t, T, '-r')
% ylabel('$\tau_\mathrm{II}$', 'interpreter', 'latex')
% xlabel('$t$', 'interpreter', 'latex')
% legend('Flow stress', 'No iterations', 'Local iterations', 'location', 'southeast')
% set(gca, 'fontsize', 18)