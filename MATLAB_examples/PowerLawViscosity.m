clear

% Gas constant
R   = 8.314;

% Ambient  conditions
expEii = -18:0.1:-11;
Eii = 10.^(expEii);
T   = 500+273;
P   = 250e8;
d   = 2e-3;
phi = 0.0;

% lithosphere
eta0 = 1e9;

% case 1
tpwl = 0;
npwl = 10;
mpwl = 0.0;
rpwl = 0.0;
Qpwl = 0;
Vpwl = 0.0e-6;
Apwl = eta0^(-npwl);
fpwl = 0.0;
apwl = 0;

% asthenosphere
% eta0 = 1e10;
% 
% % case 1
% tpwl = 0;
% npwl = 3;
% mpwl = 0.0;
% rpwl = 0.0;
% Qpwl = 0;
% Vpwl = 0.0e-6;
% Apwl = eta0^(-npwl);
% fpwl = 0.0;
% apwl = 0;

Eapwl = Qpwl;
Vapwl = Vpwl;

Fpwl    = 1.0;
A1pwl   = Fpwl * power(Apwl,-1/npwl) * exp( (Eapwl + P*Vapwl)/R/npwl/T ) * power(d, mpwl/npwl) * power(fpwl, -rpwl/npwl) * exp(-apwl*phi/npwl);
eta_pwl = A1pwl * power( Eii, 1.0/npwl - 1.0 );

figure(1), clf
plot(log10(Eii), log10(eta_pwl), 'd')
xlabel('Eii', 'interpreter', 'latex')
ylabel('Viscosity - $\eta_{\textrm{pwl}}$', 'interpreter', 'latex')
drawnow

figure(2), clf
plot(Eii, (eta_pwl).*2.*Eii/1e3, 'd')
xlabel('Eii', 'interpreter', 'latex')
ylabel('Stress $', 'interpreter', 'latex')
drawnow



