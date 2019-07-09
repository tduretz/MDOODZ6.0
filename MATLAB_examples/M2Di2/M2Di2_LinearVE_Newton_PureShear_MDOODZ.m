% =========================================================================
% Nonlinear power-law Stokes flow with analytical Newton solver in pure shear.

% Copyright (C) 2017  Ludovic Raess, Thibault Duretz, Yury podladchikov

% This file is part of M2Di.

% M2Di is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% M2Di is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with M2Di.  If not, see <http://www.gnu.org/licenses/>.
% =========================================================================

function M2Di2_PowerLawVE_Newton_PureShear
addpath('./_M2Di2_functions/')
%% Switches
SuiteSparse = 0;                                                             % 0=NO SuiteSparse  1=SuiteSparse
noisy       = 2;                                                             % 0=NO Output       1=Std output     2=Very noisy  3=Very noisy + kspgcr output
Newton      = 1;                                                             % 0=Picard          1=Newton
mat_pert    = 1;                                                             % Thermal or material perturbation
%% Physics
xmin        =-0.5;                                                           % min x
xmax        = 0.5;                                                           % max x
ymin        =-0.5;                                                           % min y
ymax        = 0.5;                                                           % max y
Rg          = 1;                                                             % Gas constant
r           = 0.05;                                                          % Inclusion radius
Ebg         = 1;                                                             % Background strain rate
Tbg         = 1;                                                             % Background temperature
% Material 1 - MATRIX
ndis(1)     = 1;                                                            % Power Law exponent
Qdis(1)     = 20;                                                            % Activation energy
Adis(1)     = 1^ndis(1)*exp(Qdis(1)/Rg/Tbg);                                 % Pre-exponent (Power law model)
mu_r(1)     = 1*exp(-Qdis(1)/ndis(1)/Rg/Tbg);                                % Reference viscosity (Carreau model)
mu_0(1)     = 1e3*mu_r(1);                                                   % Viscosity for zero strain rate (Carreau model)
mu_i(1)     = 1e-6;                                                          % Viscosity for infinite strain rate (Carreau model)
ksi(1)      = (mu_r(1)/(mu_0(1)-mu_i(1)))^(1/(1/ndis(1)-1));
G(1)        = 1;
% Material 2 - INCLUSION
ndis(2)     = 1;
Qdis(2)     = 1;
Adis(2)     = (1/1e1)^ndis(2)*exp(Qdis(2)/Rg/Tbg);
mu_r(2)     = 1e1*exp(-Qdis(2)/ndis(2)/Rg/Tbg);
mu_0(2)     = 1e0*mu_r(2);
mu_i(2)     = 1e0*mu_r(2);
ksi(2)      = (mu_r(2)/(mu_0(2)-mu_i(2)))^(1/(1/ndis(2)-1));
G(2)        = 1;
%% Numerics
nx          = 101;                                                           % Number of cells in x
ny          = nx;                                                            % Number of cells in y
nt          = 20;
dt          = 1/2;
% Non-linear solver settings
nmax_glob   = 200;                                                           % Max. number of non-linear iterations
tol_glob    = 1e-8;                                                          % Global nonlinear tolerance
eps_kspgcr  = 1e-4;                                                          % KSPGCR solver tolerance
LineSearch  = 1;                                                             % Line search algorithm
alpha       = 1.0;                                                           % Default step (in case LineSearch = 0)
% Linear solver settings
lsolver     = 1;                                                             % 0=Backslash 1=Powell-Hestenes
gamma       = 1e-5;                                                          % Numerical compressibility
nPH         = 30;                                                            % Max. number of linear iterations
tol_linu    = tol_glob/5000;                                                 % Velocity tolerance
tol_linp    = tol_glob/1000;                                                 % Divergence tolerance
%% Dimensionalisation
h     = 1; % Material phase used for scaling
% Characterictic values
muc   = Adis(h)^(-1/ndis(h))*Ebg^(1/ndis(h)-1)*exp(Qdis(h)/ndis(h)/Rg/Tbg);
Tc    = Qdis(h)/Rg;                                                          % Kelvins
Lc    = (xmax-xmin);                                                         % Meters
tc    = 1/Ebg;                                                               % Seconds
tauc  = muc*(1/tc);                                                          % Pascals
mc    = tauc*Lc*tc^2;                                                        % Kilograms
Jc    = mc*Lc^2/tc^2;                                                        % Joules
% Scaling
Adis  = Adis./(tauc.^-ndis.*1./tc);
Qdis  = Qdis./Jc;
Rg    = Rg/(Jc/Tc);
Ebg   = Ebg/(1/tc);
Tbg   = Tbg/Tc;
r     = r/Lc;
mu_0  = mu_0/muc;
mu_i  = mu_i/muc;
xmin  = xmin/Lc;  xmax = xmax/Lc;
ymin  = ymin/Lc;  ymax = ymax/Lc;
%% Preprocessing
tic
Lx = xmax-xmin;     dx = Lx/nx;                                              % Grid spacing x
Ly = ymax-ymin;     dy = Ly/ny;                                              % Grid spacing y
xn = xmin:dx:xmax;  xc = xmin+dx/2:dx:xmax-dx/2;                             % Nodal coordinates x
yn = ymin:dy:ymax;  yc = ymin+dy/2:dy:ymax-dy/2;                             % Nodal coordinates y
[xn2,  yn2] = ndgrid(xn,yn);
[xc2,  yc2] = ndgrid(xc,yc);
[xvx2,yvx2] = ndgrid(xn,yc);
[xvy2,yvy2] = ndgrid(xc,yn);
rxvec  = zeros(nmax_glob,1); cpu = zeros(9,1);
%% Initial arrays
Vx     = -Ebg.*xvx2;                                                         % Initial solutions for velocity in x
Vy     =  Ebg.*yvy2;                                                         % Initial solutions for velocity in y
Pt     =    0.* xc2;                                                         % Initial solutions for pressure
phc    =     ones(nx  ,ny  );                                                % Material phase on centroids
phv    =     ones(nx+1,ny+1);                                                % Material phase on vertices
Tec    = Tbg*ones(nx  ,ny  );                                                % Temperature on centroids
Txxc   =    zeros(nx  ,ny  );
Tyyc   =    zeros(nx  ,ny  );
Txyv   =    zeros(nx+1,ny+1);
stress = zeros(nt,1);
etaec  = G(phc).*dt.*ones(nx  ,ny  );
etaev  = G(phv).*dt.*ones(nx+1,ny+1);
if mat_pert==0, Tec((xc2+0*Lx/2).^2+(yc2-0*Ly/2).^2<r^2) = Tbg + 0.1*Tbg;    % Thermal perturbation
else            phc((xc2+0*Lx/2).^2+(yc2-0*Ly/2).^2<r^2) = 2;                % Material perturbation
                phv((xn2+0*Lx/2).^2+(yn2-0*Ly/2).^2<r^2) = 2;                % Material perturbation
end
%% Numbering Pt and Vx,Vy
NumVx  = reshape(1:(nx+1)*ny,nx+1,ny  );
NumVy  = reshape(1:nx*(ny+1),nx  ,ny+1); NumVyG = NumVy + max(NumVx(:));
NumPt  = reshape(1:nx*ny    ,nx  ,ny  ); NumPtG = NumPt + max(NumVyG(:));
cpu(1)=toc;
%% Define BC's ---------------- NEW STUFF
% Free slip / no slip setting for x momentum
BC.nsxS    =  zeros(size(Vx)); %BC.nsxS( :         , 1  ) = 1;
BC.nsxN    =  zeros(size(Vx)); %BC.nsxN( :         , end) = 1;
BC.nsxW    =  zeros(size(Vx)); %BC.nsxW([1     2  ], :  ) = 1;
BC.nsxE    =  zeros(size(Vx)); %BC.nsxE([end-1 end], :  ) = 1;
BC.fsxS    =  zeros(size(Vx)); BC.fsxS( :         , 1  ) = 1;  % Free slip S
BC.fsxN    =  zeros(size(Vx)); BC.fsxN( :         , end) = 1;  % Free slip N
BC.fsxW    =  zeros(size(Vx)); BC.fsxW([1     2  ], :  ) = 1;  % Free slip W
BC.fsxE    =  zeros(size(Vx)); BC.fsxE([end-1 end], :  ) = 1;  % Free slip E
% Free slip / no slip setting for y momentum
BC.nsyW    =  zeros(size(Vy)); %BC.nsyW( 1 , :         ) = 1;
BC.nsyE    =  zeros(size(Vy)); %BC.nsyE(end, :         ) = 1;
BC.nsyS    =  zeros(size(Vy)); %BC.nsyS( : ,[1     2]  ) = 1;
BC.nsyN    =  zeros(size(Vy)); %BC.nsyN( : ,[end-1 end]) = 1;
BC.fsyW    =  zeros(size(Vy)); BC.fsyW( 1 , :         ) = 1;  % Free slip W
BC.fsyE    =  zeros(size(Vy)); BC.fsyE(end, :         ) = 1;  % Free slip E
BC.fsyS    =  zeros(size(Vy)); BC.fsyS( : ,[1     2  ]) = 1;  % Free slip S
BC.fsyN    =  zeros(size(Vy)); BC.fsyN( : ,[end-1 end]) = 1;  % Free slip N
% Free surface
BC.frxN    =  zeros(size(Vx)); %BC.frxN(2:end-1,end) = 1;
BC.fryN    =  zeros(size(Vy)); %BC.fryN( :     ,end) = 1;
%% Boundary values - Dirichlets
Vx_W   = -Ebg.*xvx2(1  ,:)'.*ones(1 ,ny  )';                                 % BC value Vx West
Vx_E   = -Ebg.*xvx2(end,:)'.*ones(1 ,ny  )';                                 % BC value Vx East
Vy_S   =  Ebg.*yvy2(:,  1) .*ones(nx,1   );                                  % BC value Vy South
Vy_N   =  Ebg.*yvy2(:,end) .*ones(nx,1   );                                  % BC value Vy North
Vx_S = zeros(nx+1,1);  Vx_N = zeros(nx+1,1);
Vy_W = zeros(1,ny+1)'; Vy_E = zeros(1,ny+1)';
%% Construct BC structure ---------------- NEW STUFF
BC.Ux_W = Vx_W; BC.Ux_E = Vx_E; BC.Ux_S = Vx_S; BC.Ux_N = Vx_N;
BC.Uy_S = Vy_S; BC.Uy_N = Vy_N; BC.Uy_W = Vy_W; BC.Uy_E = Vy_E;
cpu(2)=toc;
%% Assemble pressure gradient and divergence operator
tic
[ grad, div, BcD ] = M2Di2_AssembleDivGrad( BC, dx, dy, nx, ny,NumVx, NumVyG, NumPt, SuiteSparse );
%% Assemble Block PP
iPt   = NumPt;  I = iPt(:)';  J = I;                                     % Eq. index center (pressure diagonal)
V     = gamma*ones(size(Pt(:)));                                         % Center coeff.
if SuiteSparse==1, PP = sparse2(I,J,V);                                  % Matrix assembly
else               PP = sparse (I,J,V); end
PPI   = spdiags(1./diag(PP),0,PP);                                       % Trivial inverse of diagonal PP matrix
cpu(3)=toc;
%% Timesteps
for it=1:nt
%     Vx     = -Ebg.*xvx2;                                                         % Initial solutions for velocity in x
%     Vy     =  Ebg.*yvy2;                                                         % Initial solutions for velocity in y
%     Pt     =    0.* xc2;                                                         % Initial solutions for pressure
    Txxco  = Txxc;
    Tyyco  = Tyyc;
    Txyvo  = Txyv;
    %% Non linear iterations
    iter=0; resnlu=2*tol_glob; resnlu0=2*tol_glob;
    du = [0*Vx(:); 0*Vy(:)];         dp = zeros(max(NumPt(:)),1);
    % Constant 2D field (temperature perturbation) - Interpolate from centers to nodes and extrapolates boundaries
    Tei = 0.25*(Tec(1:end-1,1:end-1) + Tec(2:end,1:end-1) + Tec(1:end-1,2:end) + Tec(2:end,2:end));
    Tev = zeros(nx+1,ny+1);  Tev(2:end-1,2:end-1)=Tei; Tev([1,end],:)=Tev([2,end-1],:); Tev(:,[1,end])=Tev(:,[2,end-1]);
    while( resnlu/resnlu0 > tol_glob && iter < nmax_glob )
        iter = iter+1; if noisy>=1, fprintf('\n *** Nonlin iter %d *** \n', iter ); end
        tic
        % Initial guess or iterative solution for strain increments
        divV        = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
        Exxc        = diff(Vx,1,1)/dx - 1/3*divV;
        Eyyc        = diff(Vy,1,2)/dy - 1/3*divV;
        Vx_exp      = [2*BC.nsxS(:,1).*BC.Ux_S-BC.nsxS(:,1).*Vx(:,1) + BC.fsxS(:,1).*Vx(:,1), Vx, 2*BC.nsxN(:,end).*BC.Ux_N-BC.nsxN(:,end).*Vx(:,end) + BC.fsxN(:,end).*Vx(:,end)];
        Vy_exp      = [2*BC.nsyW(1,:).*BC.Uy_W'-BC.nsyW(1,:).*Vy(1,:) + BC.fsyW(1,:).*Vy(1,:); Vy; 2*BC.nsyE(end,:).*BC.Uy_E'-BC.nsyE(end,:).*Vy(end,:) + BC.fsyE(end,:).*Vy(end,:)];
        dVxdy       = diff(Vx_exp,1,2)/dy;
        dVydx       = diff(Vy_exp,1,1)/dx;
        Exyv        = 0.5*( dVxdy + dVydx );
        % Extrapolate trial strain components
        Exyc      = 0.25*(Exyv (1:end-1,1:end-1) + Exyv (2:end,1:end-1) + Exyv (1:end-1,2:end) + Exyv (2:end,2:end));
        Txyco     = 0.25*(Txyvo(1:end-1,1:end-1) + Txyvo(2:end,1:end-1) + Txyvo(1:end-1,2:end) + Txyvo(2:end,2:end));
        [ Exxv  ] = M2Di2_centroids2vertices( Exxc );
        [ Eyyv  ] = M2Di2_centroids2vertices( Eyyc );
        [ Txxvo ] = M2Di2_centroids2vertices( Txxco );
        [ Tyyvo ] = M2Di2_centroids2vertices( Tyyco );
        % Engineering convention
        Gxyc = 2*Exyc;
        Gxyv = 2*Exyv;
        % Invariants
        Eiic2    = 1/2*(Exxc.^2 + Eyyc.^2) + Exyc.^2;
        Eiiv2    = 1/2*(Exxv.^2 + Eyyv.^2) + Exyv.^2;
        % Viscosity
        mc   = 1/2*( 1./ndis(phc) - 1);
        mv   = 1/2*( 1./ndis(phv) - 1);
        Bc   = Adis(phc).^(-1./ndis(phc)) .* exp(Qdis(phc)/Rg./Tec./ndis(phc));
        Bv   = Adis(phv).^(-1./ndis(phv)) .* exp(Qdis(phv)/Rg./Tev./ndis(phv));
        e1 = Bc.*Eiic2.^mc;
        etac = (1./(Bc.*Eiic2.^mc) + 1./etaec).^(-1);
        etav = (1./(Bv.*Eiiv2.^mv) + 1./etaev).^(-1);
        
        min(e1(:))*muc
        max(e1(:))*muc
        
        min(etac(:))*muc
        max(etac(:))*muc
        %% Rheological coefficients for the Picard operator
        % Centroids
        D.D11c = 2*etac; D.D12c = 0*etac; D.D13c = 0*etac;
        D.D21c = 0*etac; D.D22c = 2*etac; D.D23c = 0*etac;
        D.D31c = 0*etac; D.D32c = 0*etac; D.D33c = 1*etac;
        % Vertices
        D.D11v = 2*etav; D.D12v = 0*etav; D.D13v = 0*etav;
        D.D21v = 0*etav; D.D22v = 2*etav; D.D23v = 0*etav;
        D.D31v = 0*etav; D.D32v = 0*etav; D.D33v = 1*etav;
        %% Assemble Picard operator
        [ K, BcK ] = M2Di2_GeneralAnisotropicVVAssembly( D, BC, NumVx, NumVy, NumVyG, nx, ny, dx, dy, 1, Vx, Vy, SuiteSparse );
        %% Rheological coefficients for the Jacobian
        detadexxc =     Bc.*Eiic2.^(mc-1).*mc.*Exxc.*etaec.^2 ./ (etaec + Bc.*Eiic2.^mc).^2; detadexxv =     Bv.*Eiiv2.^(mv-1).*mv.*Exxv.*etaev.^2 ./ (etaev + Bv.*Eiiv2.^mv).^2;
        detadeyyc =     Bc.*Eiic2.^(mc-1).*mc.*Eyyc.*etaec.^2 ./ (etaec + Bc.*Eiic2.^mc).^2; detadeyyv =     Bv.*Eiiv2.^(mv-1).*mv.*Eyyv.*etaev.^2 ./ (etaev + Bv.*Eiiv2.^mv).^2;
        detadgxyc = 1/2*Bc.*Eiic2.^(mc-1).*mc.*Gxyc.*etaec.^2 ./ (etaec + Bc.*Eiic2.^mc).^2; detadgxyv = 1/2*Bv.*Eiiv2.^(mv-1).*mv.*Gxyv.*etaev.^2 ./ (etaev + Bv.*Eiiv2.^mv).^2;
        % Centroids
        D.D11c = 2*etac +   2*detadexxc.*Exxc + 1/1*detadexxc.*Txxco./etaec; D.D12c =            2*detadeyyc.*Exxc + 1/1*detadeyyc.*Txxco./etaec; D.D13c =            2*detadgxyc.*Exxc + 1/1*detadgxyc.*Txxco./etaec;
        D.D21c =            2*detadexxc.*Eyyc + 1/1*detadexxc.*Tyyco./etaec; D.D22c = 2*etac +   2*detadeyyc.*Eyyc + 1/1*detadeyyc.*Tyyco./etaec; D.D23c =            2*detadgxyc.*Eyyc + 1/1*detadgxyc.*Tyyco./etaec;
        D.D31c =            1*detadexxc.*Gxyc + 1/1*detadexxc.*Txyco./etaec; D.D32c =            1*detadeyyc.*Gxyc + 1/1*detadeyyc.*Txyco./etaec; D.D33c = 1*etac +   1*detadgxyc.*Gxyc + 1/1*detadgxyc.*Txyco./etaec;
        % Vertices
        D.D11v = 2*etav +   2*detadexxv.*Exxv + 1/1*detadexxv.*Txxvo./etaev; D.D12v =            2*detadeyyv.*Exxv + 1/1*detadeyyv.*Txxvo./etaev; D.D13v =            2*detadgxyv.*Exxv + 1/1*detadgxyv.*Txxvo./etaev;
        D.D21v =            2*detadexxv.*Eyyv + 1/1*detadexxv.*Tyyvo./etaev; D.D22v = 2*etav +   2*detadeyyv.*Eyyv + 1/1*detadeyyv.*Tyyvo./etaev; D.D23v =            2*detadgxyv.*Eyyv + 1/1*detadgxyv.*Tyyvo./etaev;
        D.D31v =            1*detadexxv.*Gxyv + 1/1*detadexxv.*Txyvo./etaev; D.D32v =            1*detadeyyv.*Gxyv + 1/1*detadeyyv.*Txyvo./etaev; D.D33v = 1*etav +   1*detadgxyv.*Gxyv + 1/1*detadgxyv.*Txyvo./etaev;
        %% Assemble Jacobian
        [ J,  ~  ] = M2Di2_GeneralAnisotropicVVAssembly( D, BC, NumVx, NumVy, NumVyG, nx, ny, dx, dy, 1, Vx, Vy, SuiteSparse );
        %% Elasticity
        Fex = zeros(size(Vx));
        Fex(2:end-1,:) =       diff(etac./etaec.*Txxco,1,1)./dx;
        Fex            = Fex + diff(etav./etaev.*Txyvo,1,2)./dy;
        Fey = zeros(size(Vy));
        Fey(:,2:end-1) =       diff(etac./etaec.*Tyyco,1,2)./dy;
        Fey            = Fey + diff(etav./etaev.*Txyvo,1,1)./dx;
        % Evaluate non-Linear residual using matrix-vector products
        bu = BcK + [Fex(:); Fey(:)]; bp = BcD; u = [Vx(:); Vy(:)]; p = Pt(:);  % Assemble velocity/pressure vectors
        fu = -(   K*u + grad*p - bu );  resnlu = norm(fu)/length(fu);          % Velocity non-linear residuals
        fp = -( div*u          - bp );  resnlp = norm(fp)/length(fp);          % Pressure non-linear residuals
        if iter==1, resnlu0=resnlu; resnlp0=resnlp; end
        if noisy>=1,fprintf('Chk: NonLin res. ||res.u||=%2.4e, ||res.u||/||res.u0||=%2.4e\n', resnlu, resnlu/resnlu0 );
            fprintf('Chk: NonLin res. ||res.p||=%2.4e, ||res.p||/||res.p0||=%2.4e\n', resnlp, resnlp/resnlp0 ); rxvec(iter) = resnlu/resnlu0; end
        if (resnlu/resnlu0 < tol_glob), break; end
        cpu(4)=cpu(4)+toc;
        %% Linear solver - obtain velocity and pressure corrections for current non linear iteration
        if lsolver==0 % Backslash coupled solver ------------------------------
            tic
            Ms  = [ J   , grad ; ...
                div , 0*PP ];                                                % Assemble entire Jacobian
            f   = [ fu  ; fp   ];                                                % Assemble entire residual vector
            cpu(5)=cpu(5)+toc;
            tic
            dX  = Ms\f;                                                          % Call direct solver
            du  = dX(1:max(NumVyG(:)));                                          % Extract velocity correction
            dp  = dX(NumPtG(:));                                                 % Extract pressure correction
            f1  = Ms*dX - f;                                                     % Compute entire linear residual
            fu1 = f1(1:max(NumVyG(:)));                                          % Extract velocity linear residual
            fp1 = f1(NumPtG(:));                                                 % Extract pressure linear residual
            cpu(6)=cpu(6)+toc;
        elseif lsolver==1 % Powell-Hestenes INCOMPRESSIBLE - Decoupled/segregated solve --------
            if Newton==0, J = K; end
            tic
            Kt  = K - grad*(PPI*div);                                            % Velocity Schur complement of Picard operator (Kt)
            Jt  = J - grad*(PPI*div);                                            % Velocity Schur complement of Jacobian operator
            [Kc,e,s] = chol(Kt,'lower','vector');                                % Choleski factorization of Kt
            cpu(5)=cpu(5)+toc;
            tic
            % Powell-Hestenes iterations
            fu0 = fu;                                                                    % Save linear norm 0
            for itPH=1:nPH
                fut  = fu - grad*dp - grad*PPI*fp;                                       % Iterative right hand side
                [du,norm_r,its] = M2Di2_kspgcr_m(Jt,fut,du,Kc,s,eps_kspgcr,noisy,SuiteSparse); % Apply inverse of Schur complement
                dp   = dp + PPI*(fp -  div*du);                                          % Pressure corrctions
                fu1  = fu -   J*du  - grad*dp;                                           % Compute linear velocity residual
                fp1  = fp - div*du;                                                      % Compute linear pressure residual
                if noisy>1, fprintf('--- iteration %d --- \n',itPH);
                    fprintf('  Res. |u| = %2.2e \n',norm(fu1)/length(fu1));
                    fprintf('  Res. |p| = %2.2e \n',norm(fp1)/length(fp1));
                    fprintf('  KSP GCR: its=%1.4d, Resid=%1.4e \n',its,norm_r); end
                if ((norm(fu1)/length(fu1)) < tol_linu) && ((norm(fp1)/length(fp1)) < tol_linp), break; end
                if ((norm(fu1)/length(fu1)) > (norm(fu0)/length(fu1)) && norm(fu1)/length(fu1) < tol_glob), fprintf(' > Linear residuals do no converge further:\n'); break; end
                fu0 = fu1;
            end
            cpu(6)=cpu(6)+toc;
        end
        tic
        if noisy>=1, fprintf(' - Linear res. ||res.u||=%2.2e\n', norm(fu1)/length(fu1) );
                     fprintf(' - Linear res. ||res.p||=%2.2e\n', norm(fp1)/length(fp1) ); end
        % Call line search - globalization procedure
        if LineSearch==1, [alpha,~] = LineSearch_Direct(BC,nx,ny,dx,dy,u,p,du,dp,Tec,Tev,phc,phv,Rg,Adis,Qdis,ndis,mu_i,mu_0,ksi,0,Newton,0,1/3,noisy,0,resnlu0,resnlp0,resnlu,resnlp,etaec,etaev,Txxco,Tyyco,Txyvo); end
        u    = u + alpha*du;                                                     % Velocity update from non-linear iterations
        p    = p + alpha*dp;                                                     % Pressure update from non-linear iterations
        Vx   = reshape(u(NumVx(:)) ,[nx+1,ny  ]);                                % Velocity in x (2D array)
        Vy   = reshape(u(NumVyG(:)),[nx  ,ny+1]);                                % Velocity in y (2D array)
        Pt   = reshape(p(NumPt(:)) ,[nx  ,ny  ]);                                % Pressure      (2D array)
        cpu(7)=cpu(7)+toc;
    end
    Txxc      = 2*etac.*Exxc + etac./etaec.*Txxco;
    Tyyc      = 2*etac.*Eyyc + etac./etaec.*Tyyco;
    Txyv      = 2*etav.*Exyv + etav./etaev.*Txyvo;
    Txyc      = 0.25*(Txyv (1:end-1,1:end-1) + Txyv (2:end,1:end-1) + Txyv (1:end-1,2:end) + Txyv (2:end,2:end));
    Tiic      = sqrt(1/2*Txxc.^2 + 1/2*Tyyc.^2 + 1/1*Txyc.^2);
    stress(it)= norm(Tiic);
    tic
    %% Post-processing
    figure(4),clf,colormap('jet'),set(gcf,'Color','white')
    plot((1:it)*dt, (stress(1:it)/rxvec(1)),'bx-')
    xlabel('Time'),ylabel('Stress'),drawnow
    figure(1),clf,colormap('jet'),set(gcf,'Color','white')
    plot(1:iter, log10(rxvec(1:iter)/rxvec(1)),'rx-')
    xlabel('Iterations'),ylabel('Residuals'),drawnow
    figure(2),clf,colormap('jet'),set(gcf,'Color','white')
    imagesc(flipud(log10(sqrt(Eiic2)'*(1/tc)))),colorbar,axis image,title(['max Eii = ', num2str(max(sqrt(Eiic2(:))*(1/tc)), '%2.4e')])
    xlabel('x','FontSize',15), ylabel('y','FontSize',15)
    title(['$\dot{\epsilon}_{II}$', ' min = ', num2str((min(sqrt(Eiic2(:))*(1/tc))), '%2.2e'  ), ' s$^{-1}$', ' max = ', num2str((max(sqrt(Eiic2(:))*(1/tc))), '%2.2e'  ), ' s$^{-1}$' ], 'interpreter', 'latex'),set(gca,'FontSize',15)
    figure(5),clf,colormap('jet'),set(gcf,'Color','white')
    imagesc(flipud(((Tiic)'*(tauc)))),colorbar,axis image,title(['max Tii = ', num2str(max((Tiic(:))*(tauc)), '%2.4e')])
    xlabel('x','FontSize',15), ylabel('y','FontSize',15)
    title(['${\tau}_{II}$', ' min = ', num2str((min((Tiic(:))*(tauc))), '%2.2e'  ), ' Pa', ' max = ', num2str((max(sqrt(Tiic(:))*(tauc))), '%2.2e'  ), ' Pa' ], 'interpreter', 'latex'),set(gca,'FontSize',15)
    figure(3),clf,colormap('jet'),set(gcf,'Color','white')
    imagesc(flipud(log10(etac'*muc))),colorbar,axis image,title(['max eta = ', num2str(max(etac(:)*muc), '%2.4e')])
    xlabel('x','FontSize',15),ylabel('y','FontSize',15)
    title(['$\eta$', ' min = ', num2str((min(etac(:)*muc)), '%2.2e'  ), ' s$^{-1}$', ' max = ', num2str((max(etac(:)*muc)), '%2.2e'  ), ' s$^{-1}$' ], 'interpreter', 'latex'),set(gca, 'FontSize', 15)
    cpu(8)=toc; cpu(9)=sum(cpu(1:7));
    display([' Time preprocess   = ', num2str(cpu(1))]);
    display([' Time BC           = ', num2str(cpu(2))]);
    display([' Time cst block    = ', num2str(cpu(3))]);
    display([' Time assemble     = ', num2str(cpu(4))]);
    display([' Time precondition = ', num2str(cpu(5))]);
    display([' Time solve        = ', num2str(cpu(6))]);
    display([' => WALLTIME      == ', num2str(cpu(9))]);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    FUNCTIONS USED IN MAIN CODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [alpha, LSsuccess] = LineSearch_Direct(BC,nx,ny,dx,dy,u,p,du,dp,Tec,Tev,phc,phv,Rg,Adis,Qdis,ndis,mu_i,mu_0,ksi,KRO,Newton,comp,pls,noisy,INV,resnlu0,resnlp0,resnlu,resnlp,etaec,etaev,Txxco,Tyyco,Txyvo)
% Line search explicit algorithm
nsteps = 11;                                   % number of steps
amin   = 0.0;                                  % minimum step
if Newton==1, amax = 1.0; else amax = 2.0; end % maximum step
dalpha = (amax-amin)/(nsteps-1);
alphav = amin:dalpha:amax;
nVx    = (nx+1)*ny; nRMe   = zeros(nsteps,1); nRCTe  = zeros(nsteps,1);
u0 = u; p0 = p;
% Compute non linear residual for different steps
for ils=1:nsteps;
    % Primitive variables with correction step
    u       = u0 + alphav(ils).*du;
    p       = p0 + alphav(ils).*dp;
    % Evaluate non-Linear residual in matrix-free form
    Pt      = reshape(p           ,[nx  ,ny  ]);
    Vx      = reshape(u(1:nVx)    ,[nx+1,ny  ]);
    Vy      = reshape(u(nVx+1:end),[nx  ,ny+1]);
    % Initial guess or iterative solution for strain increments
    divV        = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
    Exxc        = diff(Vx,1,1)/dx - 1/3*divV;
    Eyyc        = diff(Vy,1,2)/dy - 1/3*divV;
    Vx_exp      = [2*BC.nsxS(:,1).*BC.Ux_S-BC.nsxS(:,1).*Vx(:,1) + BC.fsxS(:,1).*Vx(:,1), Vx, 2*BC.nsxN(:,end).*BC.Ux_N-BC.nsxN(:,end).*Vx(:,end) + BC.fsxN(:,end).*Vx(:,end)];
    Vy_exp      = [2*BC.nsyW(1,:).*BC.Uy_W'-BC.nsyW(1,:).*Vy(1,:) + BC.fsyW(1,:).*Vy(1,:); Vy; 2*BC.nsyE(end,:).*BC.Uy_E'-BC.nsyE(end,:).*Vy(end,:) + BC.fsyE(end,:).*Vy(end,:)];
    dVxdy       = diff(Vx_exp,1,2)/dy;
    dVydx       = diff(Vy_exp,1,1)/dx;
    Exyv        = 0.5*( dVxdy + dVydx );
    % Extrapolate trial strain components
    Exyc     = 0.25*(Exyv(1:end-1,1:end-1) + Exyv(2:end,1:end-1) + Exyv(1:end-1,2:end) + Exyv(2:end,2:end));
    [ Exxv ] = M2Di2_centroids2vertices( Exxc );
    [ Eyyv ] = M2Di2_centroids2vertices( Eyyc );
    % Engineering convention
    Gxyc = 2*Exyc;
    Gxyv = 2*Exyv;
    % Invariants
    Eiic2    = 1/2*(Exxc.^2 + Eyyc.^2) + Exyc.^2;
    Eiiv2    = 1/2*(Exxv.^2 + Eyyv.^2) + Exyv.^2;
    % Viscosity
    mc   = 1/2*( 1./ndis(phc) - 1);
    mv   = 1/2*( 1./ndis(phv) - 1);
    Bc   = Adis(phc).^(-1./ndis(phc)) .* exp(Qdis(phc)/Rg./Tec./ndis(phc));
    Bv   = Adis(phv).^(-1./ndis(phv)) .* exp(Qdis(phv)/Rg./Tev./ndis(phv));
    etac = (1./(Bc.*Eiic2.^mc) + 1./etaec).^(-1);
    etav = (1./(Bv.*Eiiv2.^mv) + 1./etaev).^(-1);
    % Momentum residual
    tau_xx     = 2*etac.*Exxc + etac./etaec.*Txxco;
    tau_yy     = 2*etac.*Eyyc + etac./etaec.*Tyyco;
    tau_xy     = 2*etav.*Exyv + etav./etaev.*Txyvo;
    Res_x      = diff(-Pt + tau_xx,1,1)/dx + diff(tau_xy(2:end-1,:),1,2)/dy;
    Res_y      = diff(-Pt + tau_yy,1,2)/dy + diff(tau_xy(:,2:end-1),1,1)/dx;
    RMe        = [Res_x(:) ; Res_y(:)];
    nRMe(ils)  = norm(RMe)/((nx+1)*ny + (ny+1)*nx);
    % Continuity residual
    RCTe       = -(diff(Vx,1,1)/dx+diff(Vy,1,2)/dy);
    nRCTe(ils) = norm(RCTe(:))/(nx*ny);
    
    if ils == 1
        fprintf(' LS: it. = %02d\n', ils );
        fprintf(' LS: NonLin res. ||res.u||=%2.4e, ||res.u||/||res.u0||=%2.4e\n', nRMe (ils), nRMe (ils)/resnlu0 );
        fprintf(' LS: NonLin res. ||res.p||=%2.4e, ||res.p||/||res.p0||=%2.4e\n', nRCTe(ils), nRCTe(ils)/resnlp0 );
    end
end
% Find optimal step (i.e. yielding to lowest residuals)
[~,ibestM] = min(nRMe);
if ibestM==1, LSsuccess=0; alpha=0; fprintf('No descent found - continuing ...')
    nRMe(1) = 2*max(nRMe); [~,ibestM] = min(nRMe);
    LSsuccess=1; alpha = alphav(ibestM);
else          LSsuccess=1; alpha = alphav(ibestM); end
if noisy>=1
    fprintf(' LS: Selected alpha = %2.2f\n', alpha );
    fprintf(' LS: NonLin res. ||res.u||=%2.4e, ||res.u||/||res.u0||=%2.4e\n', nRMe (ibestM), nRMe (ibestM)/resnlu0 );
    fprintf(' LS: NonLin res. ||res.p||=%2.4e, ||res.p||/||res.p0||=%2.4e\n', nRCTe(ibestM), nRCTe(ibestM)/resnlp0 );
end
end