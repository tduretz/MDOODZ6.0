function [alpha, LSsuccess] = LineSearch_Direct_TM(BC,nx,ny,dx,dy,u,p,T,du,dp,dT,Tc,Tv,phc,phv,npow,T0,Newton,comp,pls,noisy,INV,resnlu0,resnlp0,resnlT0,resnlu,resnlp,resnlT,Tco,N1,N2,dt,etaec,etaev,Txxco,Tyyco,Txyco,Txxvo,Tyyvo,Txyvo,phic,phiv,Cc,Cv,lit,litmax,littol)
% Line search explicit algorithm
nsteps = 11;                                   % number of steps
amin   = 0.0;                                  % minimum step
if Newton==1, amax = 1.0; else amax = 2.0; end % maximum step
dalpha = (amax-amin)/(nsteps-1);
alphav = amin:dalpha:amax;
nVx    = (nx+1)*ny; nRMe   = zeros(nsteps,1); nRCTe  = zeros(nsteps,1); nRTe  = zeros(nsteps,1);
u0 = u; p0 = p; Tini = T;
% Compute non linear residual for different steps
for ils=1:nsteps;
    % Primitive variables with correction step
    u       = u0 + alphav(ils).*du;
    p       = p0 + alphav(ils).*dp;
    T       = Tini + alphav(ils).*dT;
    % Evaluate non-Linear residual in matrix-free form
    Ptc     = reshape(p           ,[nx  ,ny  ]);
    Vx      = reshape(u(1:nVx)    ,[nx+1,ny  ]);
    Vy      = reshape(u(nVx+1:end),[nx  ,ny+1]);
    Tc      = reshape(T           ,[nx  ,ny  ]);
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
    Txyco     = 0.25*(Txyvo(1:end-1,1:end-1) + Txyvo(2:end,1:end-1) + Txyvo(1:end-1,2:end) + Txyvo(2:end,2:end));
    [ Exxv ] = M2Di2_centroids2vertices( Exxc );
    [ Eyyv ] = M2Di2_centroids2vertices( Eyyc );
    [ Tv   ] = M2Di2_centroids2vertices( Tc   );
    [ Txxvo ] = M2Di2_centroids2vertices( Txxco );
    [ Tyyvo ] = M2Di2_centroids2vertices( Tyyco );
    [  Ptv  ] = M2Di2_centroids2vertices(  Ptc );
    % Engineering convention
    Gxyc = 2*Exyc;
    Gxyv = 2*Exyv;
    % Invariants
    Eiic2    = 1/2*(Exxc.^2 + Eyyc.^2) + Exyc.^2;
    Eiiv2    = 1/2*(Exxv.^2 + Eyyv.^2) + Exyv.^2;
    % Viscosity
    mc   = 1/2*( 1./npow - 1) * ones(nx  , ny  );
    mv   = 1/2*( 1./npow - 1) * ones(nx+1, ny+1);
    [ etac, Txxc, Tyyc, Txyc, detadexxc, detadeyyc, detadexyc ] = M2Di2_LocalIteration_VEP_TM( litmax, littol, noisy, Exxc, Eyyc, Exyc, Txxco, Tyyco, Txyco, Tc, T0, mc, etaec, Ptc, phic, Cc );
    [ etav, Txxv, Tyyv, Txyv, detadexxv, detadeyyv, detadexyv ] = M2Di2_LocalIteration_VEP_TM( litmax, littol, noisy, Exxv, Eyyv, Exyv, Txxvo, Tyyvo, Txyvo, Tv, T0, mv, etaev, Ptv, phiv, Cv );
    % Momentum residual
    Res_x      = diff(-Ptc + Txxc,1,1)/dx + diff(Txyv(2:end-1,:),1,2)/dy;
    Res_y      = diff(-Ptc + Tyyc,1,2)/dy + diff(Txyv(:,2:end-1),1,1)/dx;
    RMe        = [Res_x(:) ; Res_y(:)];
    nRMe(ils)  = norm(RMe)/((nx+1)*ny + (ny+1)*nx);
    % Continuity residual
    RCTe       = -(diff(Vx,1,1)/dx+diff(Vy,1,2)/dy);
    nRCTe(ils) = norm(RCTe(:))/(nx*ny);
    % Thermal residual
    Hs        = Txxc.*Exxc + Tyyc.*Eyyc + 2*Txyc.*Exyc;
    qx        = -diff(Tc,1,1)/dx; qxa = zeros(nx+1,ny); qxa(2:end-1,:) = qx;
    qy        = -diff(Tc,1,2)/dy; qya = zeros(nx,ny+1); qya(:,2:end-1) = qy;
    dTdt      = -(diff(qxa,1,1)/dx + diff(qya,1,2)/dy) + N2*Hs;
    RTe       = N1*(Tc(:)-Tco(:))./dt - dTdt(:); 
    nRTe(ils) = norm(RTe)/length(RTe);
    if ils == 1
        fprintf(' LS: it. = %02d\n', ils );
        fprintf(' LS: NonLin res. ||res.u||=%2.4e, ||res.u||/||res.u0||=%2.4e\n', nRMe (ils), nRMe (ils)/resnlu0 );
        fprintf(' LS: NonLin res. ||res.p||=%2.4e, ||res.p||/||res.p0||=%2.4e\n', nRCTe(ils), nRCTe(ils)/resnlp0 );
        fprintf(' LS: NonLin res. ||res.T||=%2.4e, ||res.T||/||res.T0||=%2.4e\n', nRTe (ils), nRTe (ils)/resnlT0 );
    end
end
% Find optimal step (i.e. yielding to lowest residuals)
[~,ibestM] = min(nRMe + nRTe);
if ibestM==1, LSsuccess=0; alpha=0; fprintf('No descent found - continuing ...')
    nRMe(1) = 2*max(nRMe); [~,ibestM] = min(nRMe);
    LSsuccess=1; alpha = alphav(ibestM);
else          LSsuccess=1; alpha = alphav(ibestM); end

if noisy>=1
    fprintf(' LS: Selected alpha = %2.2f\n', alpha );
    fprintf(' LS: NonLin res. ||res.u||=%2.4e, ||res.u||/||res.u0||=%2.4e\n', nRMe (ibestM), nRMe (ibestM)/resnlu0 );
    fprintf(' LS: NonLin res. ||res.p||=%2.4e, ||res.p||/||res.p0||=%2.4e\n', nRCTe(ibestM), nRCTe(ibestM)/resnlp0 );
    fprintf(' LS: NonLin res. ||res.T||=%2.4e, ||res.T||/||res.T0||=%2.4e\n', nRTe(ibestM), nRTe(ibestM)/resnlT0 );
end
end