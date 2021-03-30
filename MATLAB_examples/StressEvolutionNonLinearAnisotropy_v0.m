clear 
close all

nonlinear = 1; % 0: power-law exponent is 1.0 - 1: larger than 1.0 

% Constants
R   = 8.314;

% Conditions 
exx = 1e-14;
eyy = -exx;
ezz = 0;    % incompressible
exy = 0;
eii = sqrt(1/2*(exx^2+eyy^2+ezz^2) + exy^2);
T   = 773;

% Numerics
dt     = 1e11;
nt     = 50;
litmax = 10;
noisy  = 3;
tol    = 1e-12;

% Mat props
G   = 2e10;
Q   = 276e3;
if nonlinear == 1
    A   = 3.2e-20;
    n   = 3.0;
else
    A   = 2*(2.237e8)^(-1)*exp(Q/R/T)*eii; % same flow stress as non-linear model
    n   = 1.0;
end

% Precompute
B    = A^(-1/n)*exp(Q/n/R/T);
C    = (2*B)^-n;

% Isolated viscosity
eta_pwl  = B*eii^(1/n-1);
eta_el   = G*dt;

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
    Eii  = sqrt(1/2*(Exx^2+Eyy^2+Ezz^2) + Exy^2);
    
    % No local iteration
    eta_ve = (1/eta_pwl+1/eta_el)^(-1);
    
    % New deviatoric stress
    Txx = 2*eta_ve*Exx;
    Tyy = 2*eta_ve*Eyy;
    Tzz = 2*eta_ve*Ezz;
    Txy = 2*eta_ve*Exy;
    Tii = 2*eta_ve*Eii;
    
    % Sanity check
    Tii_chk  = sqrt(1/2*(Txx^2+Tyy^2+Tzz^2) + Txy^2);
    norm(Tii_chk-Tii)
    
    time  = time + dt;
    T(it) = Tii;
    t(it) = time;
    
end

% Visualisation step #1
figure(1), clf
hold on
plot(t, 2*eta_pwl*eii*ones(size(T)), '-.k')
plot(t, T, '+k')

%%%%%%%%%%%%%%%%%%%%%% NOW WITH LOCAL ITERATIONS

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
    Eii  = sqrt(1/2*(Exx^2+Eyy^2+Ezz^2) + Exy^2);
    
    % Trial viscosity
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
    Txx = 2*eta_ve*Exx;
    Tyy = 2*eta_ve*Eyy;
    Tzz = 2*eta_ve*Ezz;
    Txy = 2*eta_ve*Exy;
    Tii = 2*eta_ve*Eii;
    
    % Sanity check
    Tii_chk  = sqrt(1/2*(Txx^2+Tyy^2+Tzz^2) + Txy^2);
    norm(Tii_chk-Tii)
    
    time  = time + dt;
    T(it) = Tii;
    t(it) = time;
    
end

% Visualisation step #2
plot(t, T, '-r')
ylabel('$\tau_\mathrm{II}$', 'interpreter', 'latex')
xlabel('$t$', 'interpreter', 'latex')
legend('Flow stress', 'No iterations', 'Local iterations', 'location', 'southeast')
set(gca, 'fontsize', 18)