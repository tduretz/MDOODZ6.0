clear
close all

R   = 8.314510; % [J.mol-1.K-1]
d   = 1e-3;

factor2 = 1; % 1 = as in rheology papers, 2 = as in code

% dislocation creep
ndis = 4; % Shelton 1981
mdis = 1;
Adis = 10^0.9*(1e-6)^(ndis);
Edis = 431e3;

% diffusion creep
ndif = 1; % Rybacki et Dresen 2000
mdif = 3;
Adif = 10^(5.1)*(1e-6)^(ndif)*(1e-6)^(mdif);
Edif = 268e3;

% Peierls creep
npei = 2; % Azuma et al., 2014
Apei = 10^(3.48)*(1e-6)^(npei);
Spei = 9831e6;
Epei = 431e3;
p    = 1;
q    = 2;
gam  = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Smin =   5;
Smax =  10;
Tmin = 500;
Tmax = 1500;
nS   = 401;
nT   = 401;
dS   = (Smax-Smin)/(nS-1);
dT   = (Tmax-Tmin)/(nT-1);
S    =  Smin:dS:Smax;
S    = 10.^S;
T    =  Tmin:dT:Tmax;
[T2d,S2d] = meshgrid(T,S);

E_dif = Adif.*d.^(-mdif).*S2d.^ndif.*exp(-Edif./R./T2d);
E_dis = Adis*S2d.^ndis.*exp(-Edis./R./T2d);
E_pei = Apei*S2d.^npei.*exp(-Epei./R./T2d .* (1 - (S2d./Spei).^p).^q );
E_pei(S2d<100e6) = 0; %% NEED TO CUT OF PEIERLS
St        = Epei./R./T2d*(1-gam)^(q-1)*q*gam;
E_pei_reg = Apei.*S2d.^npei.*exp(-Epei./R./T2d .*(1-gam)^q) .* (S2d./gam./Spei).^St;

E     = max(max(E_dif,E_dis), E_pei);
E_reg = max(max(E_dif,E_dis), E_pei_reg);

dom   = zeros(size(T2d));
dom(E_pei>E_dif & E_pei>E_dis) = 3;
dom(E_dis>E_dif & E_dis>E_pei) = 2;
dom(E_dif>E_dis & E_dif>E_pei) = 1;

dom_reg   = zeros(size(T2d));
dom_reg(E_pei_reg>E_dif & E_pei_reg>E_dis) = 3;
dom_reg(E_dis>E_dif & E_dis>E_pei_reg) = 2;
dom_reg(E_dif>E_dis & E_dif>E_pei_reg) = 1;

figure(1), clf
subplot(121)
hold on
imagesc(T, log10(S/1e6), dom), shading interp, axis xy
contour(T2d, log10(S2d/1e6), E, [1e-15 1e-15], 'w')
caxis([1 3])
xlabel('T [K]')
ylabel('log_{10} S [MPa]')
title('Dry plagioclase (Azuma et al., 2016)')
axis square

subplot(122)
hold on
imagesc(T, log10(S/1e6), dom_reg), shading interp, axis xy
contour(T2d, log10(S2d/1e6), E_reg, [1e-15 1e-15], 'w')
contour(T2d, log10(S2d/1e6), E_reg, [1e-12 1e-12], 'k', 'LineWidth', 2)
plot((550+273)* ones(size(S)), log10(S2d/1e6), '-k' , 'LineWidth', 2 )
caxis([1 3])
xlabel('T [K]')
ylabel('log_{10} S [MPa]')
title('Dry plagioclase regularized Peierls')
axis square

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Emin = -30;
Emax = -1;
nE   = 401;
dE   = (Emax-Emin)/(nS-1);
E    =  Emin:dE:Emax;
E    = 10.^E;

[T2d,E2d] = meshgrid(T,E);

St      = Epei./R./T2d*(1-gam)^(q-1)*q*gam;

% Updated in MDOODZ6.0 
Cdif = Adif.*d.^(-mdif).*exp(-Edif./R./T2d);
Cdis = Adis.*exp(-Edis./R./T2d);
Cpei = Apei.*exp(-Epei./R./T2d.*power(1.0-gam,2.0)) .* (gam*Spei).^(-St);

Bdif = 1/factor2*Cdif.^(-1./ndif);
Bdis = 1/factor2*Cdis.^(-1./ndis);
Bpei = 1/factor2*Cpei.^(-1./(npei+St));

S_dif     = Bdif .* E2d.^(1/ndif-1) .* factor2.*E2d;
S_dis     = Bdis .* E2d.^(1/ndis-1) .* factor2.*E2d;
S_pei_reg = Bpei .* E2d.^(1./(npei+St)-1) .* factor2.*E2d;

S_reg   = min(S_dis,min(S_dis,S_pei_reg));
dom_reg = zeros(size(T2d));
dom_reg(S_pei_reg<S_dif & S_pei_reg<S_dis) = 3;
dom_reg(S_dis<S_dif & S_dis<S_pei_reg) = 2;
dom_reg(S_dif<S_dis & S_dif<S_pei_reg) = 1;

S_target = 10^(3.05+6);

figure(2), clf
hold on
imagesc(T-273, log10(E), dom_reg), shading interp, axis xy
% imagesc(T-273, log10(E), S_pei_reg), shading interp, axis xy

plot((550)* ones(size(E)), log10(E), '-k' , 'LineWidth', 2 )
plot(T-273, log10(1e-12)*ones(size(T)), '-k' , 'LineWidth', 2 )
contour(T2d-273, log10(E2d), S_reg, [S_target S_target], 'w')
caxis([1 3])
% contour(T2d, log10(S2d/1e6), E_reg, [1e-15 1e-15], 'w')
xlabel('T [C]')
ylabel('log_{10} E [1/s]')
title('Dry plagioclase regularized Peierls')
axis square

