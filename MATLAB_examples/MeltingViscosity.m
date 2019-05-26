clear

% Gas constant
R   = 8.314;

% Ambient  conditions
Eii = 1e-15;
T   = 500+273;
P   = 250e8;
d   = 2e-3;
phi = 0:0.001:0.5;

% Westerly granite
tpwl = 1;
npwl = 3.3;
mpwl = 0.0;
rpwl = 0.0;
Qpwl = 186.5e3;
Vpwl = 0.0e-6;
Apwl = 3.1623e-26;
fpwl = 0.0;
apwl = 80.0;

Eapwl = Qpwl;
Vapwl = Vpwl;

Fpwl    = 1.0/6.0*power(2.0,1.0/npwl) * power(3.0,(npwl-1.0)/2.0/npwl);
A1pwl   = Fpwl * power(Apwl,-1/npwl) * exp( (Eapwl + P*Vapwl)/R/npwl/T ) * power(d, mpwl/npwl) * power(fpwl, -rpwl/npwl) * exp(-apwl*phi/npwl);
eta_pwl = A1pwl * power( Eii, 1.0/npwl - 1.0 );

figure(1), clf
plot(phi, (eta_pwl)*2*Eii/1e3, 'd')
xlabel('Melt content - $\phi$', 'interpreter', 'latex')
ylabel('Viscosity - $\eta_{\textrm{pwl}}$', 'interpreter', 'latex')
drawnow


