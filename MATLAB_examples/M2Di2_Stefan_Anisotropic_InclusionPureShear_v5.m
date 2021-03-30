% =========================================================================
% Anisotropic Stokes flow in pure shear.

% Copyright (C) 2019 Thibault Duretz, Stefan Schmalholz

% This file is part of M2Di2.

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

function M2Di2_Stefan_Anisotropic_InclusionPureShear_v5
addpath('./_M2Di2_functions/')
addpath('./M2Di2/_M2Di2_functions/')
aniso             = 1;
director_angle    = 35;
anisotropy_factor = 10;
%% Physics                                                               
xmin    = -0.5;
xmax    =  0.5;
ymin    = -0.5;
ymax    =  0.5;
mu0     = 10;                                                               % Background viscosity
mu_i    = 1;                                                             % Inclusion viscosity
rad     = 0.2337;                                                               % Inclusion radius
%% Numerics
nx          = 20;                                                         % Grid points in x
ny          = nx;                                                          % Grid points in y
gamma       = 1e5;                                                         % Numerical compressibility                                                            % Max number of linear iterations
SuiteSparse = 0;
lsolver     = 0;
nPH         = 20;
eps_kspgcr  = 5e-4;                                                        % KSPGCR solver tolerance
tol_glob    = 1e-10;                                                       % Global nonlinear tolerance
tol_linu    = tol_glob/5000;                                               % Velocity tolerance  
tol_linp    = tol_glob/10000;
noisy       = 2;
chol_ksp    = 0;   % uses LU (UMFPACK) in Powell-Hestenes iterations instead of KSP + Cholesky (CHOLMOD) as preconditioner  
%% Preprocessing
tic
dx = (xmax-xmin)/nx;                                                                 % cell size in x
dy = (ymax-ymin)/ny;                                                                 % cell size in y
xv = xmin:dx:xmax;           yv = ymin:dy:ymax;           [xv2, yv2] = ndgrid(xv, yv); % cell coord. grid
xc = xmin+dx/2:dx:xmax-dx/2; yc = ymin+dy/2:dy:ymax-dy/2; [xc2, yc2] = ndgrid(xc, yc); % cell coord. grid
[xvx2,yvx2] = ndgrid(xv,yc);
[xvy2,yvy2] = ndgrid(xc,yv);
%% Initial conditions
Ebg    = -1;
Vx     =   Ebg.*xvx2;                                                         % Initial solutions for velocity in x
Vy     =  -Ebg.*yvy2;                                                         % Initial solutions for velocity in y
Pt     =    0.* xc2;
muc                               =  mu0*ones(nx  ,ny  );                  % viscosity at cell center
muv                               =  mu0*ones(nx+1,ny+1);                  % viscosity at vertex
muc(sqrt(xc2.^2 + yc2.^2) < rad ) = mu_i;                                  % define inclusion
muv(sqrt(xv2.^2 + yv2.^2) < rad ) = mu_i;
% muc(sqrt((xc2-xmin).^2 + (yc2-ymin).^2) < rad ) = mu_i;                                  % define inclusion
% muv(sqrt((xv2-xmin).^2 + (yv2-ymin).^2) < rad ) = mu_i;
% muc(sqrt((xc2-xmin).^2 + (yc2+ymin).^2) < rad ) = mu_i;                                  % define inclusion
% muv(sqrt((xv2-xmin).^2 + (yv2+ymin).^2) < rad ) = mu_i;
% muc(sqrt((xc2+xmin).^2 + (yc2-ymin).^2) < rad ) = mu_i;                                  % define inclusion
% muv(sqrt((xv2+xmin).^2 + (yv2-ymin).^2) < rad ) = mu_i;
% muc(sqrt((xc2+xmin).^2 + (yc2+ymin).^2) < rad ) = mu_i;                                  % define inclusion
% muv(sqrt((xv2+xmin).^2 + (yv2+ymin).^2) < rad ) = mu_i;
delta_ani_c                       = anisotropy_factor*ones(size(muc)); % eta_normal - eta_shear
% delta_ani_c(sqrt(xc2.^2 + yc2.^2) > rad ) = dani;
% delta_ani_c([1 end],:) = 1; delta_ani_c(:,[1 end]) = 1; % make boundaries isotropic - to be improved
delta_ani_v                       = anisotropy_factor*ones(size(muv)); % eta_normal - eta_shear
% delta_ani_v(sqrt(xv2.^2 + yv2.^2) > rad ) = dani; 
% delta_ani_v([1 end],:) = 1; delta_ani_v(:,[1 end]) = 1; % make boundaries isotropic - to be improved
n1c                               = cosd(director_angle).*ones(size(muc));
n2c                               = sind(director_angle).*ones(size(muc));
n1v                               = cosd(director_angle).*ones(size(muv));
n2v                               = sind(director_angle).*ones(size(muv));
normc  =  sqrt(n1c.^2 + n2c.^2);  
normv  =  sqrt(n1v.^2 + n2v.^2);
n1c    = n1c./normc; n2c    = n2c./normc;
n1v    = n1v./normv; n2v    = n2v./normv;
%% Numbering Pt and Vx,Vy
NumVx  = reshape(1:(nx+1)*ny,nx+1,ny  );
NumVy  = reshape(1:nx*(ny+1),nx  ,ny+1); NumVyG = NumVy + max(NumVx(:));    % G stands for Gobal numbering
NumPt  = reshape(1:nx*ny    ,nx  ,ny  ); NumPtG = NumPt + max(NumVyG(:));   % G stands for Gobal numbering
cpu(1)=toc;
%% Boundary values - Dirichlets
Vx_W   =  xmin*Ebg.*ones(1 ,ny  )';                           % BC value Vx West
Vx_E   =  xmax*Ebg.*ones(1 ,ny  )';                           % BC value Vx East
Vy_S   = -ymin*Ebg.*ones(nx,1   );                            % BC value Vy South
Vy_N   = -ymax*Ebg.*ones(nx,1   );                            % BC value Vy South
Vx_S = yv(1)*ones(nx+1,1);  Vx_N =  yv(end)*ones(nx+1,1);  
Vy_W =     zeros(1,ny+1)';  Vy_E =     zeros(1,ny+1)';
%% Define BC's  ---------------- NEW STUFF
% Free slip / no slip setting for x momentum
BC.nsxS    =  zeros(size(Vx)); %BC.nsxS( :         , 1  ) = 1;
BC.nsxN    =  zeros(size(Vx)); %BC.nsxN( :         , end) = 1;
BC.nsxW    =  zeros(size(Vx)); %BC.nsxW([1     2  ], :  ) = 1;
BC.nsxE    =  zeros(size(Vx)); %BC.nsxE([end-1 end], :  ) = 1;
BC.fsxS    =  zeros(size(Vx)); BC.fsxS( :         , 1  ) = 1;  % Free slip S
BC.fsxN    =  zeros(size(Vx)); BC.fsxN( :         , end) = 1;  % Free slip N
BC.fsxW    =  zeros(size(Vx)); BC.fsxW([1     2  ], :  ) = 1;  % Free slip W
BC.fsxE    =  zeros(size(Vx)); BC.fsxE([end-1 end], :  ) = 1;  % Free slip E
BC.NeuW    =  zeros(size(Vx));
BC.NeuE    =  zeros(size(Vx));
BC.Dirx    =  zeros(size(Vx)); BC.Dirx([1 end], :) = 1;
% Free slip / no slip setting for y momentum
BC.nsyW    =  zeros(size(Vy)); %BC.nsyW( 1 , :         ) = 1;
BC.nsyE    =  zeros(size(Vy)); %BC.nsyE(end, :         ) = 1;
BC.nsyS    =  zeros(size(Vy)); %BC.nsyS( : ,[1     2]  ) = 1;
BC.nsyN    =  zeros(size(Vy)); %BC.nsyN( : ,[end-1 end]) = 1;
BC.fsyW    =  zeros(size(Vy)); BC.fsyW( 1 , :         ) = 1;  % Free slip W
BC.fsyE    =  zeros(size(Vy)); BC.fsyE(end, :         ) = 1;  % Free slip E
BC.fsyS    =  zeros(size(Vy)); BC.fsyS( : ,[1     2  ]) = 1;  % Free slip S
BC.fsyN    =  zeros(size(Vy)); BC.fsyN( : ,[end-1 end]) = 1;  % Free slip N
BC.NeuS    =  zeros(size(Vy));
BC.NeuN    =  zeros(size(Vy));
BC.Diry    =  zeros(size(Vy)); BC.Diry( :,[1 end]) = 1;
% Free surface
BC.frxN    =  zeros(size(Vx)); %BC.frxN(2:end-1,end) = 1;
BC.fryN    =  zeros(size(Vy)); %BC.fryN( :     ,end) = 1;
%% Construct BC structure  ---------------- NEW STUFF
BC.Ux_W = Vx_W; BC.Ux_E = Vx_E; BC.Ux_S = Vx_S; BC.Ux_N = Vx_N; 
BC.Uy_S = Vy_S; BC.Uy_N = Vy_N; BC.Uy_W = Vy_W; BC.Uy_E = Vy_E;
%% Assemble pressure gradient and divergence operator
tic
[ grad, div, BcD ] = M2Di2_AssembleDivGrad( BC, dx, dy, nx, ny, NumVx, NumVyG, NumPt, SuiteSparse );
cpu(3)=toc;
%% PP block
iPt    = NumPt;  I = iPt(:)';  J = I;                                       % Eq. index center (pressure diagonal)
V      = ones(nx*ny,1)./gamma;                                              % Center coeff.
if SuiteSparse==1, PP = sparse2(I,J,V); else                                % Matrix assembly
                   PP =  sparse(I,J,V); end
PPI   = spdiags(1./diag(PP),0,PP);
%% Rheological coefficients  ---------------- NEW STUFF
% Centroids - pressure and normal points
detac  = aniso*(muc - muc./delta_ani_c);
nc     = sqrt(n1c.^2 + n2c.^2);
n1c    = n1c./nc;
n2c    = n2c./nc;
d0     = 2*n1c.^2.0.*n2c.^2.0;
d1     = n1c.*n2c.*(-n1c.^2 + n2c.^2);
D.D11c = -2.0*detac.*d0 + 2.0*muc;
D.D12c =  2.0*detac.*d0;
D.D13c =  2.0*detac.*d1;
D.D14c =  0.0.*d0;
D.D21c =  2.0*detac.*d0;
D.D22c = -2.0*detac.*d0 + 2.0*muc;
D.D23c = -2.0*detac.*d1;
D.D24c =  0.0.*d0;
% Vertices - nodal points - shear stress points
detav  = aniso*(muv - muv./delta_ani_v);
nv     = sqrt(n1v.^2 + n2v.^2);
n1v    = n1v./nv;
n2v    = n2v./nv;
d0     = 2*n1v.^2.0.*n2v.^2.0;
d1     = n1v.*n2v.*(-n1v.^2 + n2v.^2);
D.D31v =  2.0*detav.*d1;
D.D32v = -2.0*detav.*d1;
D.D33v =  2.0*detav.*(d0 - 0.5) + muv ;
D.D34v =  0.0*d0;
% D.D31v([1 end], :) = 0; D.D31v(:, [1 end]) = 0;
% D.D32v([1 end], :) = 0; D.D32v(:, [1 end]) = 0;
% D.D33v([1 end], :) = 0; D.D33v(:, [1 end]) = 0;
%% Assemble Picard operator  ---------------- NEW STUFF
% ACHTUNG: Matrix symmetrisation is turned off
% [ K, BcK ] = M2Di2_GeneralAnisotropicVVAssembly( D, BC, NumVx, NumVy, NumVyG, nx, ny, dx, dy, 0, Vx, Vy, SuiteSparse );
ps = 1/3;
% Residual
divV        = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
Exxc        = diff(Vx,1,1)/dx - ps*divV;
Eyyc        = diff(Vy,1,2)/dy - ps*divV;
Vx_exp      = [2*BC.nsxS(:,1).*BC.Ux_S-BC.nsxS(:,1).*Vx(:,1) + BC.fsxS(:,1).*Vx(:,1), Vx, 2*BC.nsxN(:,end).*BC.Ux_N-BC.nsxN(:,end).*Vx(:,end) + BC.fsxN(:,end).*Vx(:,end)];
Vy_exp      = [2*BC.nsyW(1,:).*BC.Uy_W'-BC.nsyW(1,:).*Vy(1,:) + BC.fsyW(1,:).*Vy(1,:); Vy; 2*BC.nsyE(end,:).*BC.Uy_E'-BC.nsyE(end,:).*Vy(end,:) + BC.fsyE(end,:).*Vy(end,:)];
dVxdy       = diff(Vx_exp,1,2)/dy;
dVydx       = diff(Vy_exp,1,1)/dx;
Exyv        = 0.5*( dVxdy + dVydx );
% Extrapolate trial strain components
Exyc     = 0.25*(Exyv(1:end-1,1:end-1) + Exyv(2:end,1:end-1) + Exyv(1:end-1,2:end) + Exyv(2:end,2:end));
[ Exxv ] = M2Di2_centroids2vertices( Exxc );
[ Eyyv ] = M2Di2_centroids2vertices( Eyyc );
% Momentum residual
tau_xx      = D.D11c.*Exxc + D.D12c.*Eyyc + D.D13c.*2.*Exyc;
tau_yy      = D.D21c.*Exxc + D.D22c.*Eyyc + D.D23c.*2.*Exyc;
tau_xy      = D.D31v.*Exxv + D.D32v.*Eyyv + D.D33v.*2.*Exyv;
Res_x       = zeros(size(Vx)); Res_y       = zeros(size(Vy));
Res_x(2:end-1, :)  = diff(-Pt + tau_xx,1,1)/dx + diff(tau_xy(2:end-1,:),1,2)/dy;
Res_y(:, 2:end-1)  = diff(-Pt + tau_yy,1,2)/dy + diff(tau_xy(:,2:end-1),1,1)/dx;
% Assembly
Fu = [Res_x(:); Res_y(:)];

Fp = divV(:);
[ K ] = M2Di2_GeneralAnisotropicVVAssembly_v3( D, BC, NumVx, NumVy, NumVyG, nx, ny, dx, dy, 0, Vx, Vy, ps, SuiteSparse );
cpu(4)=toc; display(['Time Build Blocks = ', num2str(cpu(4))])

save('AnisotropyMatrixM2Di','K', 'grad', 'div', 'PP', 'Fu', 'Fp', 'dx', 'dy', 'NumVx', 'NumVyG')

%% Linear solver
tic
if lsolver == 0       % One that sucks
    PPd      = 0*PP;  % set one pressure constrain 
    PPd(1,1) = gamma; % set one pressure constrain 
    BcD(1)   = 0;     % set one pressure constrain 
    div(1,:) = 0;     % set one pressure constrain 
    grad(:,1)= 0;     % set one pressure constrain 
    M   = [ K   , grad ; ...
            div ,  PPd   ];
    b   = [ Fu;   Fp  ];
    up  = M\b;
elseif lsolver == 1   % Killer solver
    J   = K;
    fu  = BcK;  fp  = BcD;
    du  = 0*fu; dp  = 0*fp;
    tic
    Js  = 1/2*(J'+ J);                                                     % Symmetrisation of Jacobian
    Jt  = J  - grad*(PPI*div);                                            % Velocity Schur complement of Jacobian operator
    Jts = Js - grad*(PPI*div);                                            % Velocity Schur complement of Jacobian operator
    if chol_ksp == 1
        if SuiteSparse == 0
            [Jcs,e,s] = chol(Jts,'lower','vector');                               % Choleski factorization of Jts
        else
            s     = amd2(Jts);
            Jts   = cs_transpose(Jts);
            Jts   = cs_symperm(Jts,s);
            Jts   = cs_transpose(Jts);
            Jcs   = lchol(Jts);
        end
    else
        if( exist('umfpac2k')==3 ) [L1, U1, P1, Q1] = umfpack2(Jt); end
        if( exist('umfpack' )==3 ) [L1, U1, P1, Q1] = umfpack (Jt); end
        perm = 1:size(Jt,1);
    end
    norm_r = 0; its = 0;
    % Powell-Hestenes iterations
    fu0 = fu;                                                                    % Save linear norm 0
    fp0 = fp;  
    for itPH=1:nPH
        fut  = fu - grad*dp - grad*PPI*fp;                                       % Iterative right hand side
        if chol_ksp==1, 
            [du,norm_r,its] = M2Di2_kspgcr_m(Jt,fut,du,Jcs,s,eps_kspgcr,noisy,SuiteSparse);  % Apply inverse of Schur complement
        else
            du   = cs_usolve( U1, cs_lsolve(L1, P1*fut(perm))  );
        	du(perm)    = Q1*du;
        end
        dp   = dp + PPI*(fp -  div*du);                                          % Pressure corrctions
        fu1  = fu -   J*du  - grad*dp;                                           % Compute linear velocity residual
        fp1  = fp - div*du;                                                      % Compute linear pressure residual
        if noisy>1, fprintf('--- iteration %d --- \n',itPH);
            fprintf('  Res. |u| = %2.2e \n',norm(fu1)/length(fu1));
            fprintf('  Res. |p| = %2.2e \n',norm(fp1)/length(fp1));
            fprintf('  KSP GCR: its=%1.4d, Resid=%1.4e \n',its,norm_r); end
        if ((norm(fu1)/length(fu1)) < tol_linu) && ((norm(fp1)/length(fp1)) < tol_linp), break; end
        if ((norm(fp1)/length(fp1)) > (norm(fp0)/length(fp1)) && norm(fp1)/length(fp1) < tol_glob), fprintf(' > Linear residuals do no converge further:\n'); break; end
%         if ((norm(fu1)/length(fu1)) > (norm(fu0)/length(fu1)) && norm(fu1)/length(fu1) < tol_glob), fprintf(' > Linear residuals do no converge further:\n'); break; end
        fu0 = fu1;
        fp0 = fp1;
    end
    up  = [du;dp];
end
cpu(5)=toc;
%% Post-process
Pt  = Pt + reshape(up(NumPtG(:)),[nx  ,ny  ]); Pt = Pt - mean(Pt(:));
Vx  = Vx + reshape(up(NumVx(:)) ,[nx+1,ny  ]);
Vy  = Vy + reshape(up(NumVyG(:)),[nx  ,ny+1]);
u   = up([NumVx(:);NumVyG(:)]);
% Check residuals % errors
fuv = Fu -   K*u - grad*Pt(:);
fpt = Fp - div*u;
error_mom = norm(fuv)/length(fuv); display([' Error momentum        = ', num2str(error_mom)]);
error_div = norm(fpt)/length(fpt); display([' Error divergence      = ', num2str(error_div)]);
% DivV and Mom X & Y Check
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
Eii      = sqrt( 1/2*Exxc.^2 + Eyyc.^2 + Exyc.^2 );
% Stresses
Txx      = D.D11c.*Exxc + D.D12c.*Eyyc + D.D13c.*2.*Exyc;
Tyy      = D.D21c.*Exxc + D.D22c.*Eyyc + D.D23c.*2.*Exyc;
Txy      = D.D31v.*Exxv + D.D32v.*Eyyv + D.D33v.*2.*Exyv;
% Residuals
Res_x       = diff(-Pt + Txx,1,1)/dx + diff(Txy(2:end-1,:),1,2)/dy;
Res_y       = diff(-Pt + Tyy,1,2)/dy + diff(Txy(:,2:end-1),1,1)/dx;
% Viz
figure(1),clf,colormap('jet'),set(gcf,'Color','white')
subplot(221),imagesc(xv,yc,Vx' ),axis image, colorbar, title('Vx'), set(gca,'ydir','normal')
subplot(222),imagesc(xc,yv,Vy' ),axis image, colorbar, title('Vy'), set(gca,'ydir','normal')
subplot(223),imagesc(xc,yc,Pt' ),axis image, colorbar, title('P'), set(gca,'ydir','normal'), %caxis([-2 2])
subplot(224),imagesc(xc,yc,Eii'),axis image, colorbar, title('Eii'), set(gca,'ydir','normal')
%
figure(2),clf,colormap('jet'),set(gcf,'Color','white')
subplot(311),imagesc(xc,yc,Txx' ),axis image, colorbar, title(['mean. Txx = ', num2str(mean(Txx(:))) ] ), set(gca,'ydir','normal')
subplot(312),imagesc(xc,yc,Tyy' ),axis image, colorbar, title(['mean. Tyy = ', num2str(mean(Tyy(:))) ] ), set(gca,'ydir','normal')
subplot(313),imagesc(xv,yv,Txy' ),axis image, colorbar, title(['mean. Txy = ', num2str(mean(Txy(:))) ] ), set(gca,'ydir','normal')
hold on
hv = quiver(mean(xc2(:)), mean(yc2(:)), mean(n1c(:)), mean(n2c(:)), 'k', 'Linewidth', 3);
set(hv, 'MaxHeadSize',2.0, 'AutoScaleFactor', 2, 'ShowArrowHead', 'off');
hv = quiver(mean(xc2(:)), mean(yc2(:)), -mean(n1c(:)), -mean(n2c(:)), 'k', 'Linewidth', 3);
set(hv, 'MaxHeadSize',2.0, 'AutoScaleFactor', 2, 'ShowArrowHead', 'off');
axis([xmin xmax ymin ymax])
%
figure(3),clf,colormap('jet'),set(gcf,'Color','white')
subplot(311),imagesc(xc,yc,divV'),colorbar,axis image, ylabel('Div'  ), set(gca,'ydir','normal'),title('Residuals')
subplot(312),imagesc(xv,yc,Res_x'    ),colorbar,axis image, ylabel('Res X'), set(gca,'ydir','normal')
subplot(313),imagesc(xc,yv,Res_y'    ),colorbar,axis image, ylabel('Res Y'), set(gca,'ydir','normal'),drawnow
cpu(11)=toc; display(['WALL TIME = ', num2str(sum(cpu(1:10))),' seconds (without postprocessing)']);
display(['SOLVE TIME = ', num2str((cpu(5))),' seconds']);
end