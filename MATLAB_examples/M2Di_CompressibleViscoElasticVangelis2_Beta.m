close all
SuiteSparse = 1;
noisy       = 1;
no_slip     = 0;
solver      = 2;
% M2Di Linear Stokes: Thibault Duretz, Ludovic Raess, Yuri Podladchikov  - Unil 2017
%% Physics
Lx      = 1;                                                                % Box width
Ly      = 1;                                                                % Box height
mus0    = 1e-1;                                                             % Background viscosity
mus_i   = 1e0;                                                              % Inclusion viscosity
rad     = Ly/15;                                                            % Inclusion radius
B0      = 1e5;                                                              % Background bulk modulus
B_i     = 1e5;                                                              % Inclusion bulk modulus
G0      = 1e10; 
dt      = 0.1*mus_i*B0;
VE      = 1;
Ebg     = 1;
%% Numerics
nx      = 50;                                                              % Grid points in x
ny      = nx;                                                              % Grid points in y
nt      = 50;
time    = 0;
tolL    = 1e-11;
MaxIterL= 100;
%% Preprocessing
tic
dx = Lx/nx;                                                                 % cell size in x
dy = Ly/ny;                                                                 % cell size in y
xv = -Lx/2:dx:Lx/2;           yv = -Ly/2:dy:Ly/2;           [XV2, YV2] = ndgrid(xv, yv); % cell coord. grid
xc = -Lx/2+dx/2:dx:Lx/2-dx/2; yc = -Ly/2+dy/2:dy:Ly/2-dy/2; [XC2, YC2] = ndgrid(xc, yc); % cell coord. grid
%% Initial conditions
Txx                               =      zeros(nx  ,ny  );
Tyy                               =      zeros(nx  ,ny  );
Txy                               =      zeros(nx+1,ny+1);
Tzz                               =      zeros(nx  ,ny  );
divU                              =      zeros(nx  ,ny  );  
Pt                                =      zeros(nx  ,ny  );                  % Pressure
Vx                                =      zeros(nx+1,ny  );                  % x velocity
Vy                                =      zeros(nx  ,ny+1);                  % y velocity
musc                              =  mus0*ones(nx  ,ny  );                  % viscosity at cell center
musv                              =  mus0*ones(nx+1,ny+1);
Gc                                =    G0*ones(nx  ,ny  );                  % viscosity at cell center
Gv                                =    G0*ones(nx+1,ny+1);% viscosity at vertex
Bc                                 =   B0*ones(nx  ,ny  );
% musc(sqrt((XC2-0*Lx/2).^2 + YC2.^2) < rad) = mus_i;                                  % define inclusion
% musv(sqrt((XV2-0*Lx/2).^2 + YV2.^2) < rad) = mus_i;
% Kc  ( sqrt(XC2        .^2 + YC2.^2) < rad) = K_i;

% Vangelis layer
musc( abs(YC2)<rad ) = mus_i;                                  % define inclusion
musv( abs(YV2)<rad ) = mus_i;
Bc  ( abs(YC2)<rad ) = B_i;

%% Numbering Pt and Vx,Vy
NumVx  = reshape(1:(nx+1)*ny,nx+1,ny  );
NumVy  = reshape(1:nx*(ny+1),nx  ,ny+1); NumVyG = NumVy + max(NumVx(:));    % G stands for Gobal numbering
NumPt  = reshape(1:nx*ny    ,nx  ,ny  ); NumPtG = NumPt + max(NumVyG(:));   % G stands for Gobal numbering
cpu(1)=toc;
%% Boundary Conditions on velocities [W E S N]
tic
ibcVxW = NumVx(1,:)';         ibcVxE = NumVx(end,:)';                       % Indexes of BC nodes for Vx East/West
ibcVxS = NumVx(2:end-1,1);    ibcVxN = NumVx(2:end-1,end);                  % Indexes of BC nodes for Vx North/South
ibcVyS = NumVyG(:,1);         ibcVyN = NumVyG(:,end);                       % Indexes of BC nodes for Vy North/South
ibcVyW = NumVyG(1,2:end-1)';  ibcVyE = NumVyG(end,2:end-1)';                % Indexes of BC nodes for Vy East/West
ibc    = [ ibcVxW; ibcVxE; ibcVyS; ibcVyN ];                                % Group all indexes
ibcNC  = [ ibcVxS; ibcVxN; ibcVyW; ibcVyE ];                                % Non Confroming to the physical boundary
Vx_W   = Ebg*( Lx/2)*ones(size(ibcVxW));
Vx_E   = Ebg*(-Lx/2)*ones(size(ibcVxE));
Vx_S   = 0*ones(size(ibcVxS));
Vx_N   = 0*ones(size(ibcVxN));
Vy_S   = Ebg*(-Lx/2)*ones(size(ibcVyS));
Vy_N   = Ebg*( Lx/2)*ones(size(ibcVyN));
Vy_W   = 0*ones(size(ibcVyW));
Vy_E   = 0*ones(size(ibcVyE));
fprintf('VxW = %2.2e\n', Vx_W)

vBc    = [ Vx_W;          Vx_E;          Vy_S;          Vy_N         ];     % Group all values
vBcNC  = [ Vx_S; Vx_N; Vy_W; Vy_E];     % Non Confroming to the physical boundary
cpu(3)=toc;
%% grad and div blocs
tic
iVxC   = NumVx;                                                             % dP/dx
iPtW   =  ones(size(iVxC));    iPtW(2:end-0,:) = NumPt;
iPtE   =  ones(size(iVxC));    iPtE(1:end-1,:) = NumPt;
cPtW   = -ones(size(iVxC))/dx; cPtW([1 end],:) = 0;
cPtE   =  ones(size(iVxC))/dx; cPtE([1 end],:) = 0;
Idx    = [ iVxC(:); iVxC(:) ]';
Jdx    = [ iPtW(:); iPtE(:) ]';
Vdx    = [ cPtW(:); cPtE(:) ]';
iVyC   = NumVyG;                                                            % dP/dy
iPtS   =  ones(size(iVyC));    iPtS(:,2:end-0) = NumPt;
iPtN   =  ones(size(iVyC));    iPtN(:,1:end-1) = NumPt;
cPtS   = -ones(size(iVyC))/dy; cPtS(:,[1 end]) = 0;
cPtN   =  ones(size(iVyC))/dy; cPtN(:,[1 end]) = 0;
Idy    = [ iVyC(:); iVyC(:) ]';
Jdy    = [ iPtS(:); iPtN(:) ]';
Vdy    = [ cPtS(:); cPtN(:) ]';
%% PP block
iPt    = NumPt;  I = iPt(:)';  J = I;                                       % Eq. index center (pressure diagonal)
V      = ones(nx*ny,1).*Bc(:)./dt;
if SuiteSparse==1, PP = sparse2(I,J,V); else                                % Matrix assembly
    PP =  sparse(I,J,V); end
for istep=1:nt
    Pto  = Pt; Txxo = Txx; Tyyo = Tyy; Txyo = Txy; Tzzo = Tzz;
    if VE==1, etac =  1./(1./musc + 1./(Gc*dt) ); else etac = musc; end
    if VE==1, etav =  1./(1./musv + 1./(Gv*dt) ); else etav = musv; eta_axy = etav(2:end-1,2:end-1);   end
    if VE==1, VEc  = 1./(Gc*dt) .* etac; else VEc = 0*etac; end
    if VE==1, VEv  = 1./(Gv*dt) .* etav; else VEv = 0*etav;end
    %% Block UU
    mus_W  = zeros(size(Vx));  mus_W(2:end  , :     ) = etac;
    mus_E  = zeros(size(Vx));  mus_E(1:end-1, :     ) = etac;
    mus_S  = zeros(size(Vx));  mus_S(2:end-1,2:end  ) = etav(2:end-1,2:end-1);
    mus_N  = zeros(size(Vx));  mus_N(2:end-1,1:end-1) = etav(2:end-1,2:end-1);
    iVx    = NumVx;                                                              % Eq. index for Vx (C,W,E,S,N)
    iVxW    = ones(size(Vx));   iVxW(2:end  , :     ) =  NumVx(1:end-1, :     );
    iVxS    = ones(size(Vx));   iVxS( :     ,2:end  ) =  NumVx( :     ,1:end-1);
    cVxC   = 4/3*(mus_W+mus_E)/dx/dx + (mus_S+mus_N)/dy/dy;                      % Center coeff.
    scVx   = max(cVxC(:));                               cVxC([1,end],:) = scVx; % Scaling factor for Vx Dirichlet values
    cVxW   = -4/3*mus_W/dx/dx;  cVxW([1,end],:) = 0;     cVxW  ([1 end],:) = 0;  cVxW( 2   , : ) = 0; % West coeff.
    cVxS   = -mus_S/dy/dy;  cVxS([1,end],:) = 0;         cVxS  ([1 end],:) = 0;  cVxS( :   , 1 ) = 0; % South coeff.
    Iuu    = [  iVx(:);  iVx(:);  iVx(:) ]';                                     % Triplets [I,J,V]
    Juu    = [  iVx(:); iVxW(:); iVxS(:) ]';
    Vuu    = [ cVxC(:); cVxW(:); cVxS(:) ]';
    %% Block VV
    mus_W = zeros(size(Vy));  mus_W(2:end  ,2:end-1) = etav(2:end-1,2:end-1);
    mus_E = zeros(size(Vy));  mus_E(1:end-1,2:end-1) = etav(2:end-1,2:end-1);
    mus_S = zeros(size(Vy));  mus_S( :     ,2:end  ) = etac;
    mus_N = zeros(size(Vy));  mus_N( :     ,1:end-1) = etac;
    iVy    = NumVyG;                                                             % Eq. index for Vy (C,W,E,S,N)
    iVyW    = ones(size(Vy));   iVyW(2:end  , :     ) = NumVyG(1:end-1, :     );
    iVyS    = ones(size(Vy));   iVyS( :     ,2:end  ) = NumVyG( :     ,1:end-1);
    cVyC   = (mus_W+mus_E)/dx/dx + 4/3*(mus_S+mus_N)/dy/dy;                      % Center coeff.
    scVy   = max(cVyC(:));                              cVyC(:,[1,end]) = scVy;  % Scaling factor for Vy Dirichlet values
    cVyW   = -mus_W/dx/dx;                              cVyW  (:,[1 end]) = 0;  cVyW( 1 , :   ) = 0; % West coeff.
    cVyS   = -4/3*mus_S/dy/dy;  cVyS(:,[1,end]) = 0;    cVyS  (:,[1 end]) = 0;  cVyS( : , 2   ) = 0; % South coeff.
    Ivv    = [  iVy(:);  iVy(:);  iVy(:) ]';                                    % Triplets [I,J,V]
    Jvv    = [  iVy(:); iVyW(:); iVyS(:) ]';
    Vvv    = [ cVyC(:); cVyW(:); cVyS(:) ]';
    %% Block VU
    mus_W  = zeros(size(Vy));  mus_W( :, :     ) = etav(1:end-1, :);            % Viscosities (W,E,S,N)
    mus_E  = zeros(size(Vy));  mus_E( :, :     ) = etav(2:end  , :);
    mus_S  = zeros(size(Vy));  mus_S( :,2:end  ) = etac;
    mus_N  = zeros(size(Vy));  mus_N( :,1:end-1) = etac;
    iVy    = NumVyG( :    ,2:end-1);                                            % Eq. index for VyC
    iVxSW  = NumVx(1:end-1,1:end-1);                                            % Eq. index for Vx (SW,SE,NW,NE)
    iVxSE  = NumVx(2:end  ,1:end-1);
    iVxNW  = NumVx(1:end-1,2:end  );
    iVxNE  = NumVx(2:end  ,2:end  );
    cVxSW  = (-mus_W(:,2:end-1) + 2/3*mus_S(:,2:end-1))/(dx*dy);  cVxSW(1  ,:) = 0; % Coeff. for Vx (SW,SE,NW,NE)
    cVxSE  = ( mus_E(:,2:end-1) - 2/3*mus_S(:,2:end-1))/(dx*dy);  cVxSE(end,:) = 0;
    cVxNW  = ( mus_W(:,2:end-1) - 2/3*mus_N(:,2:end-1))/(dx*dy);  cVxNW(1  ,:) = 0;
    cVxNE  = (-mus_E(:,2:end-1) + 2/3*mus_N(:,2:end-1))/(dx*dy);  cVxNE(end,:) = 0;
    Ivu    = [   iVy(:);   iVy(:);   iVy(:);   iVy(:) ]';                       % Triplets [I,J,V]
    Jvu    = [ iVxSW(:); iVxSE(:); iVxNW(:); iVxNE(:) ]';
    Vvu    = [ cVxSW(:); cVxSE(:); cVxNW(:); cVxNE(:) ];
    cpu(4)=toc; display(['Time Build Blocks = ', num2str(cpu(4))]);
    %% Assemble Blocs
    tic
    if SuiteSparse==1, K    = sparse2( [Iuu(:); Ivv(:); Ivu(:)], [Juu(:); Jvv(:); Jvu(:)], [Vuu(:); Vvv(:); Vvu(:)], (nx+1)*ny+(ny+1)*nx, (nx+1)*ny+(ny+1)*nx );
        grad = sparse2( [Idx(:); Idy(:)], [Jdx(:); Jdy(:)], [Vdx(:); Vdy(:)], (nx+1)*ny+(ny+1)*nx, nx*ny );                                      else
        K    =  sparse( [Iuu(:); Ivv(:); Ivu(:)], [Juu(:); Jvv(:); Jvu(:)], [Vuu(:); Vvv(:); Vvu(:)], (nx+1)*ny+(ny+1)*nx, (nx+1)*ny+(ny+1)*nx );
        grad =  sparse( [Idx(:); Idy(:)], [Jdx(:); Jdy(:)], [Vdx(:); Vdy(:)], (nx+1)*ny+(ny+1)*nx, nx*ny );                                      end
    div  = -grad';
    cpu(5)=toc; display(['Time Assemble     = ', num2str(cpu(5))]);
    % Elastic force - x momentum
    %     exxpW = zeros(size(Vx));  exxpW(2:end  , :     ) = exx_plc;
    %     exxpE = zeros(size(Vx));  exxpE(1:end-1, :     ) = exx_plc;
    %     exypS = zeros(size(Vx));  exypS(2:end-1,1:end-1) = exy_plv(2:end-1,2:end-1);
    %     exypN = zeros(size(Vx));  exypN(2:end-1,2:end  ) = exy_plv(2:end-1,2:end-1);
    Fxx   = VEc.*Txxo;% - 2*etac.*exx_plc ;
    Fxy   = VEv.*Txyo;% - 2*etav.*exy_plv ;
    Fexi  = diff(Fxx,1,1)/dx + diff(Fxy(2:end-1,:),1,2)/dy;
    Fex   = zeros(size(Vx)); Fex(2:end-1,:) = Fexi;
    % Elastic force - y momentum
    %     exypW = zeros(size(Vy));  exypW(1:end-1,2:end-1)  = exy_plv(2:end-1,2:end-1);
    %     exypE = zeros(size(Vy));  exypE(2:end  ,2:end-1)  = exy_plv(2:end-1,2:end-1);
    %     eyypS = zeros(size(Vy));  eyypS( :     ,2:end  )  = eyy_plc;
    %     eyypN = zeros(size(Vy));  eyypN( :     ,1:end-1)  = eyy_plc;
    Fyy   = VEc.*Tyyo;% - 2*etac.*eyy_plc ;
    Fxy   = VEv.*Txyo;% - 2*etav.*exy_plv ;
    Feyi  = diff(Fyy,1,2)/dy + diff(Fxy(:,2:end-1),1,1)/dx;
    Fey   = zeros(size(Vy)); Fey(:,2:end-1) = Feyi;
    %% BC's on K and DivV
    tic
    BcK                  = zeros(size(K,1),1);
    BcK(NumVx(:))        = BcK(NumVx(:))        + Fex(:);
    BcK(NumVyG(:))       = BcK(NumVyG(:))       + Fey(:);
    BcK(NumVx(2    ,:))  = BcK(NumVx(2    ,:))  + 4/3*etac(1  ,:  )'/dx/dx.*Vx_W;
    BcK(NumVx(end-1,:))  = BcK(NumVx(end-1,:))  + 4/3*etac(end,:  )'/dx/dx.*Vx_E;
    BcK(NumVyG(:, 2   )) = BcK(NumVyG(:, 2   )) + 4/3*etac(:  ,1  ) /dy/dy.*Vy_S;
    BcK(NumVyG(:,end-1)) = BcK(NumVyG(:,end-1)) + 4/3*etac(:  ,end) /dy/dy.*Vy_N;
    BcK(ibcVxN) = BcK(ibcVxN) +  (  2/3*etac(1:end-1,end) - etav(2:end-1,end) )/dx/dy.*Vy_N(1:end-1) ;   % VyNW
    BcK(ibcVxN) = BcK(ibcVxN) +  ( -2/3*etac(2:end  ,end) + etav(2:end-1,end) )/dx/dy.*Vy_N(2:end  ) ;   % VyNE
    BcK(ibcVxS) = BcK(ibcVxS) +  ( -2/3*etac(1:end-1,1  ) + etav(2:end-1, 1 ) )/dx/dy.*Vy_S(1:end-1) ;   % VySW
    BcK(ibcVxS) = BcK(ibcVxS) +  (  2/3*etac(2:end  ,1  ) - etav(2:end-1, 1 ) )/dx/dy.*Vy_S(2:end  ) ;   % VySE
    BcK(ibcVyE) = BcK(ibcVyE) + ((  2/3*etac(end,1:end-1) - etav(end,2:end-1) )/dx/dy.*Vx_E(1:end-1)')'; % VxSE
    BcK(ibcVyE) = BcK(ibcVyE) + (( -2/3*etac(end,2:end  ) + etav(end,2:end-1) )/dx/dy.*Vx_E(2:end  )')'; % VxNE
    BcK(ibcVyW) = BcK(ibcVyW) + (( -2/3*etac(1  ,1:end-1) + etav(1  ,2:end-1) )/dx/dy.*Vx_W(1:end-1)')'; % VxSW
    BcK(ibcVyW) = BcK(ibcVyW) + ((  2/3*etac(1  ,2:end  ) - etav(1  ,2:end-1) )/dx/dy.*Vx_W(2:end  )')'; % VxNW
    % Conforming Dirichlets
    BcK([ibcVxW; ibcVxE; ibcVyS; ibcVyN]) = [Vx_W*scVx; Vx_E*scVx; Vy_S*scVy; Vy_N*scVy];
    if no_slip == 1
        % Non-conforming Dirichlets
        cNC        = [eta_axy(:,1)/dy/dy; eta_axy(:,end)/dy/dy; eta_axy(1,:)'/dx/dx; eta_axy(end,:)'/dx/dx];
        d0         = spdiags(K,0);
        d0(ibcNC)  = d0(ibcNC)  + 2*cNC;
        BcK(ibcNC) = BcK(ibcNC) + 2*cNC.*vBcNC;
        K          = spdiags(d0,0,K);
    end
    % BC on div
    BcD                 = zeros(size(div,1),1);
    BcD                 = BcD + Pto(:).*Bc(:)./dt;
    BcD(NumPt(1  ,:  )) = BcD(NumPt(1  ,:  )) + 1/dx*Vx_W;
    BcD(NumPt(end,:  )) = BcD(NumPt(end,:  )) - 1/dx*Vx_E;
    BcD(NumPt( : ,1  )) = BcD(NumPt(:  ,1  )) + 1/dy*Vy_S;
    BcD(NumPt( : ,end)) = BcD(NumPt(:  ,end)) - 1/dy*Vy_N;
    cpu(6)=toc; display(['Time BC           = ', num2str(cpu(6))]);
    %% Prepare solve
    tic
    K = K + K' - diag(diag(K));                                                 % Build full from tri-lower (for checking)
    cpu(7)=toc;
    tic
    PPI   = spdiags(1./diag(PP   ), 0, PP   );
    Kt    = K - grad*(PPI*div);                                          % PPI*DivV = PP\DivV
    [Kchol,e,s] = chol(Kt,'lower','vector');                                   % Cholesky factorization
    cpu(8)=toc; display(['Time CHOLESKY     = ', num2str(cpu(8))]);
    tic
    cpu(9)=toc;
    %% Solve
    tic
    p   =   Pt(:); u   =  [Vx(:) ; Vy(:)];
    if solver == 1 % Penalty
        Rhs = BcK  - grad*(PPI*BcD + p);                                         % Powell-Hestenes
        if SuiteSparse==1, u(s) = cs_ltsolve(Kchol,cs_lsolve(Kchol,Rhs(s))); else % Powell-Hestenes
            u(s) = Kchol'\(Kchol\Rhs(s));                     end  % Powell-Hestenes Matlab
        p   = p + PPI*( BcD - div*u - 0*PP*p);
        fu  = BcK  -   K*u - grad*p;
        fp  = BcD  - div*u -   PP*p; %   Compressible continuity residual
        if noisy==1,
            fprintf('   Res. |fu| = %2.2e \n',norm(fu)/length(fu));
            fprintf('   Res. |fp| = %2.2e \n',norm(fp)/length(fp));
        end
        cpu(10)=toc; display(['Time Backsubs     = ', num2str(cpu(10))]);
    else
        % Pre-conditionning (~Jacobi)
        u = 0*u; p = 0*p;  du = 0*u;
        Fu            = BcK - grad*p -   K*u;
        Fp            = BcD -   PP*p - div*u;
        PPi           = spdiags(1./diag(PP),0,size(PP,1),size(PP,2));
        AscP          = K - grad*(PPi*div);
        [AsccP,pA,sA] = chol(AscP,'lower','vector');                      % Cholesky factors
        resu=1; resp=resu; iR=0;
        
        fprintf('min u = %2.6e max u = %2.6e\n', min(u), max(u))
        fprintf('min p = %2.6e max p = %2.6e\n', min(p), max(p))
        fprintf('min bu = %2.6e max bu = %2.6e\n', min(BcK), max(BcK))
        fprintf('min bp = %2.6e max bp = %2.6e\n', min(BcD), max(BcD))
        fprintf('min fu = %2.6e max fu = %2.6e\n', min(Fu), max(Fu))
        fprintf('min fp = %2.6e max fp = %2.6e\n', min(Fp), max(Fp))
        
        while (resu > tolL || resp > tolL ) && iR<MaxIterL
            iR       = iR+1;
            ru       = -(        K*u +   grad*p  - Fu ); resu = norm(ru)/length(ru);
            rp       = -(      div*u +     PP*p  - Fp ); resp = norm(rp)/length(rp);
            rusc     = ru - grad*(PPi*rp);
            du(sA)   = cs_ltsolve(AsccP,cs_lsolve(AsccP,rusc(sA)));
            dp       = PPi*(rp-div*du);
            u        = u + du; % Updates
            p        = p + dp;
        end%iR
        if iR < MaxIterL %success of linear solve
            fprintf('\n >> Linear solve [OK] in %02d iterations - %3d s\n', iR, toc);
            %                 failed=0;
            success = 1;
        else %fail of linear solve
            fprintf(' -------------------------------------------\n')
            fprintf(' >> Linear solver [FAIL] after %02d iterations\n', iR);
            failed=failed+1; restart=1; iter=0; nfail=nfail+1;
            success = 0;
            fprintf(' -------------------------------------------\n')
        end
        fprintf(' Linear:   Momentum resid = %2.6e\n', norm(ru)/length(Fu) )
        fprintf('         Continuity resid = %2.6e\n', norm(rp)/length(Fp) )
    end
    tic
    XPH = [u ; p];
    %% Post-processing
    up   = XPH;
    Pt   = reshape(up(NumPtG(:)),[nx  ,ny  ]);
    Vx   = reshape(up(NumVx(:)) ,[nx+1,ny  ]);
    Vy   = reshape(up(NumVyG(:)),[nx  ,ny+1]);
    u    = up([NumVx(:);NumVyG(:)]);
    divV = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
    exx  = diff(Vx,1,1)/dx;
    eyy  = diff(Vy,1,2)/dy;
    Vx_e = [ Vx(:,1), Vx, Vx(:,end) ];                                       % expand array using BC's - Free slip
    Vy_e = [ Vy(1,:); Vy; Vy(end,:) ];                                       % expand array using BC's - Free slip
    exy  = 0.5*( diff(Vx_e,1,2)/dy + diff(Vy_e,1,1)/dx );
    Txy  = 2*etav.*exy + VEv.*Txyo;
    Sxx  = -Pto + dt*1./Bc.*(exx+eyy) + 2*etac.*( (2/3)*exx-(1/3)*eyy) + VEc.*Txxo;
    Syy  = -Pto + dt*1./Bc.*(exx+eyy) + 2*etac.*(-(1/3)*exx+(2/3)*eyy) + VEc.*Tyyo;
    Szz  = -Pto + dt*1./Bc.*(exx+eyy) + 2*etac.*(-(1/3)*exx-(1/3)*eyy) + VEc.*Tzzo;
    
    Sxx  = -Pto + VEc.*Txxo + (dt*1./Bc+4/3*etac).*exx + (dt*1./Bc-2/3*etac).*eyy;
    Syy  = -Pto + VEc.*Tyyo + (dt*1./Bc-2/3*etac).*exx + (dt*1./Bc+4/3*etac).*eyy;
    Szz  = -Pto + VEc.*Tzzo + (dt*1./Bc-2/3*etac).*exx + (dt*1./Bc-2/3*etac).*eyy;
   
    Pt1  = -1/3*(Sxx+Syy+Szz);
    
    Txx  = Pt + Sxx; Tyy  = Pt + Syy; Tzz  = Pt + Szz;
    divU = divU + divV*dt;
    time = time + dt;
    % Check residuals % errors
    fuv = BcK -   K*u - grad*p;
    fpt = BcD - div*u -   PP*p;
    error_mom = norm(fuv)/length(fuv); display([' Error momentum        = ', num2str(error_mom)]);
    error_div = norm(fpt)/length(fpt); display([' Error divergence      = ', num2str(error_div)]);
    figure(1),clf,colormap('jet'),set(gcf,'Color','white')
    subplot(221),imagesc(xc,yc,flipud(Pt'  )),axis image,colorbar,caxis([min(Pt(:)) max(Pt(:))]), title('Pt')
    subplot(223),imagesc(xv,yc,flipud(Vx'  )),axis image,colorbar, title('Vx')
    subplot(224),imagesc(xc,yv,flipud(Vy'  )),axis image,colorbar, title('Vy')
    subplot(222),imagesc(xc,yc,flipud(log10(abs(divV')))),axis image,colorbar, title('div V')
    figure(2),clf,colormap('jet'),set(gcf,'Color','white')
    imagesc(xc,yc,flipud(Pt1'-Pt'  )),axis image,colorbar
    
    figure(3), set(gcf,'Color','white'), hold on,
    tmax = (mus_i*B_i);
    dPa  = 1 - exp(-3/4*time/tmax);
    dTxxa= -( 1 - 1/2*exp(-3/4*time/tmax));
    dmus = musc(1,fix(ny/2)) - musc(1,1);
    ffs  = abs( 2*exx(1,fix(ny/2))* dmus);
    dTxx = Txx(1,fix(ny/2)) - Txx(1,1);
    dP   = Pt(1,fix(ny/2)) - Pt(1,1);
    
    subplot(211),hold on
    plot(time/(mus_i*B0), dTxx/ffs,'bo', time/(mus_i*B0), dTxxa,'xr')
    xlabel('$t$','interpreter','latex'); ylabel('$d\tau{xx}$','interpreter','latex'); l1=legend('$\tau{xx}$','$\tau{xx}$ anal.'); set(l1,'interpreter','latex', 'box', 'off')
    axis([0 5 -1 -0.5])
    subplot (212),hold on
    plot(time/(mus_i*B0), dP/ffs,'bo',time/(mus_i*B0), dPa, 'xr') 
    axis([0 5 0 1])

    xlabel('$t$','interpreter','latex'); ylabel('$dP$','interpreter','latex'); l2=legend('$dP$', '$dP$ anal.'); set(l2,'interpreter','latex','location','southeast', 'box', 'off')
    drawnow
    fprintf('Strain rate = %2.2e, flow stress = %2.2e\n\n', Ebg, Ebg*mus0*2)
  
end