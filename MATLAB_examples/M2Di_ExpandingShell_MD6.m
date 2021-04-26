clear
SuiteSparse = 1;
noisy       = 1;
no_slip     = 0;
solver      = 2;
MD6         = 1;
rho1 = 2700;
rho2 = 3200;
% M2Di Linear Stokes: Thibault Duretz, Ludovic Raess, Yuri Podladchikov  - Unil 2017
%% Physics
Lx      = 2;                                                                % Box width
Ly      = 2;                                                                % Box height
mus0    = 1e20;                                                             % Background viscosity
mus_i   = 1e20;                                                              % Inclusion viscosity
rad     = 0.25;                                                            % Inclusion radius
B0      = 1/4e10;                                                              % Background bulk modulus
B_i     = 1/4e10;                                                              % Inclusion bulk modulus
G0      = 3e10; 
G_i     = 3e10; 
dt      = 1%0.1*mus_i*B0;
VE      = 1;
Ebg     = 0;
%% Dimensionalisation
% Characterictic units
muc   = G0*dt;
Tc    = 1;                                                                 % Kelvins
Lc    = Lx;                                                         % Meters
tc    = dt;     
% derived units
tauc  = muc*(1/tc);                                                          % Pascals
mc    = tauc*Lc*tc^2;                                                        % Kilograms
% Scaling
Ebg    = Ebg/(1/tc);
dt     = dt/tc;
Lx     = Lx/Lc;
Ly     = Ly / Lc;
mus0    = mus0/(tauc*tc);                                                             % Background viscosity
mus_i   = mus_i/(tauc*tc);                                                              % Inclusion viscosity
rad     = rad/Lc;                                                            % Inclusion radius
B0      = B0*tauc;                                                              % Background bulk modulus
B_i     = B_i*tauc;                                                              % Inclusion bulk modulus
G0      = G0/tauc; 
G_i     = G_i/tauc; 
rho1    = rho1/(mc/Lc^3);
rho2    = rho2/(mc/Lc^3);
%% Numerics
nx      = 200;                                                              % Grid points in x
ny      = nx;                                                              % Grid points in y
nt      = 1;
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
Gc(sqrt((XC2-0*Lx/2).^2 + YC2.^2) < rad) = mus_i;                                  % define inclusion
Gv(sqrt((XV2-0*Lx/2).^2 + YV2.^2) < rad) = mus_i;
Gc(sqrt((XC2-0*Lx/2).^2 + YC2.^2) < rad) = G_i;                                  % define inclusion
Gv(sqrt((XV2-0*Lx/2).^2 + YV2.^2) < rad) = G_i;
Bi( sqrt(XC2        .^2 + YC2.^2) < rad) = B_i;


f    = -15.65e-2              % percentage vol change
c    = sqrt(f + 1) - 1; 
divr                               = zeros(nx  ,ny  );
divr( sqrt(XC2.^2 + YC2.^2) < rad) = 2*c;

rho0                               = rho1*ones(nx  ,ny  );
X                                  = zeros(nx  ,ny  );
X( sqrt(XC2.^2 + YC2.^2) < rad)    = 1;
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
for istep=1:nt
    Pto  = Pt; Txxo = Txx; Tyyo = Tyy; Txyo = Txy; Tzzo = Tzz;
    if VE==1, etac =  1./(1./musc + 1./(Gc*dt) ); else etac = musc; end
    if VE==1, etav =  1./(1./musv + 1./(Gv*dt) ); else etav = musv; eta_axy = etav(2:end-1,2:end-1);   end
    if VE==1, VEc  = 1./(Gc*dt) .* etac; else VEc = 0*etac; end
    if VE==1, VEv  = 1./(Gv*dt) .* etav; else VEv = 0*etav; end
    for iter=1:2
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
    if MD6 == 0
        BcD = BcD + Pto(:).*Bc(:)./dt;
        BcD = BcD + divr(:);
        %% PP block
        iPt    = NumPt;  I = iPt(:)';  J = I;                                       % Eq. index center (pressure diagonal)
        V      = ones(nx*ny,1).*Bc(:)./dt;
        if SuiteSparse==1, PP = sparse2(I,J,V); 
        else                                % Matrix assembly
            PP =  sparse(I,J,V); 
        end
        
    else
        BcD  = BcD + log(rho0(:))./dt;
        
        BcD  = BcD + X(:).*(log(rho1) - log(rho2))./dt;

        rho_ref      = (1.0-X).*rho1 + X.*rho2;
        drhodX       = rho2 - rho1;
        dXdP         = 0;
        drho_ref_dP  = drhodX .* 0;
        rho          = rho_ref .* exp(Bc .* Pt);
        drhodp       = (rho).*Bc + exp(Bc .* Pt) .* drho_ref_dP;
        %% PP block
        iPt    = NumPt;  I = iPt(:)';  J = I;                                       % Eq. index center (pressure diagonal)
        V      = ones(nx*ny,1).*drhodp(:)./(rho(:).*dt);
        if SuiteSparse==1, PP = sparse2(I,J,V); 
        else                                % Matrix assembly
            PP =  sparse(I,J,V); 
        end
    end
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
            fprintf(' Linear:   Momentum resid = %2.6e\n', norm(ru)/length(Fu) )
            fprintf('         Continuity resid = %2.6e\n', norm(rp)/length(Fp) )
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
    XPH = [u ; p];
    up   = XPH;

    Pt   = reshape(up(NumPtG(:)),[nx  ,ny  ]);
    Vx   = reshape(up(NumVx(:)) ,[nx+1,ny  ]);
    Vy   = reshape(up(NumVyG(:)),[nx  ,ny+1]);
    end
    tic

    %% Post-processing


    u    = up([NumVx(:);NumVyG(:)]);
    divV = diff(Vx,1,1)/dx + diff(Vy,1,2)/dy;
    exx  = diff(Vx,1,1)/dx;
    eyy  = diff(Vy,1,2)/dy;
    Vx_e = [ Vx(:,1), Vx, Vx(:,end) ];                                       % expand array using BC's - Free slip
    Vy_e = [ Vy(1,:); Vy; Vy(end,:) ];                                       % expand array using BC's - Free slip
    exy  = 0.5*( diff(Vx_e,1,2)/dy + diff(Vy_e,1,1)/dx );
    Exx = exx - 1/3*divV;
    Txx = 2*etac.*Exx;
    Sxx = -Pt + Txx;

    divU = divU + divV*dt;
    time = time + dt;
    % Check residuals % errors
    fuv = BcK -   K*u - grad*p;
    fpt = BcD - div*u -   PP*p;
    error_mom = norm(fuv)/length(fuv); display([' Error momentum        = ', num2str(error_mom)]);
    error_div = norm(fpt)/length(fpt); display([' Error divergence      = ', num2str(error_div)]);
    
    Fx = zeros(nx+1,ny);
    Fx(2:end-1,:) = diff(Sxx,1,1)/dx + diff(Txy(2:end-1,:),1,2)/dy;

    
    %%%%%%%%% ANALYTICS
    rmin = 1e-10;
    rmax = 1.0;
    nr   = 1000;
    dr   = (rmax-rmin)/(nr-1);
    r    = rmin:dr:rmax;

    r1 = 0.25;
    r2 = 1.0;
    K1 = 4e10;
    G1 = 3e10;
    K2 = 4e10;
    G2 = 3e10;

    ph       = ones(nr,1)'; 
    ph(r>=r1) = 2;

    C1 = c;

    A1 = -3.0*C1.*K1.*(r1.^2 - r2.^2)./(-G1.*r1.^2 + G2.*r1.^2 + 3.0*G2.*r2.^2 - 3.0*K1.*r1.^2 + 3.0*K2.*r1.^2 + r2.^2.*(G1 + 3.0*K1));
    A2 = -3.0*C1.*K1.*r1.^2./(-G1.*r1.^2 + G2.*r1.^2 + 3.0*G2.*r2.^2 - 3.0*K1.*r1.^2 + 3.0*K2.*r1.^2 + r2.^2.*(G1 + 3.0*K1));
    B2 = -6.0*C1.*K1.*r1.^2.*r2.^2./(-G1.*r1.^2 + G2.*r1.^2 + 3.0*G2.*r2.^2 - 3.0*K1.*r1.^2 + 3.0*K2.*r1.^2 + r2.^2.*(G1 + 3.0*K1));
    B1 = 0;
    
% SxxBC = 4.617e8;
% A1 = 1.5*(-B1.*G2.*(G2 + 3.0*K2).*(r1.^2 - r2.^2) + SxxBC.*r1.^2.*r2.^2.*(4.0*G2 + 3.0*K2) - (B1.*G1 - 2.0*C1.*K1.*r1.^2).*(3.0*G2.*r1.^2 + G2.*r2.^2 + 3.0*K2.*r2.^2))./(r1.^2.*(3.0*G2.*r1.^2.*(G1 + 3.0*K1) - 3.0*G2.*r1.^2.*(G2 + 3.0*K2) + 3.0*G2.*r2.^2.*(G2 + 3.0*K2) + r2.^2.*(G1 + 3.0*K1).*(G2 + 3.0*K2)));
% A2 = 1.5*(-B1.*G2.*(G1 + 3.0*K1) - 3.0*G2.*(B1.*G1 - 2.0*C1.*K1.*r1.^2) + SxxBC.*r2.^2.*(G1 + 3.0*G2 + 3.0*K1))./(3.0*G2.*r1.^2.*(G1 + 3.0*K1) - 3.0*G2.*r1.^2.*(G2 + 3.0*K2) + 3.0*G2.*r2.^2.*(G2 + 3.0*K2) + r2.^2.*(G1 + 3.0*K1).*(G2 + 3.0*K2));
% B2 = r2.^2.*(B1.*(G1 + 3.0*K1).*(G2 + 3.0*K2) + 3.0*SxxBC.*r1.^2.*(G1 - G2 + 3.0*K1 - 3.0*K2) + 3.0*(G2 + 3.0*K2).*(B1.*G1 - 2.0*C1.*K1.*r1.^2))./(3.0*G2.*r1.^2.*(G1 + 3.0*K1) - 3.0*G2.*r1.^2.*(G2 + 3.0*K2) + 3.0*G2.*r2.^2.*(G2 + 3.0*K2) + r2.^2.*(G1 + 3.0*K1).*(G2 + 3.0*K2));
    Atab = [A1 A2];
    Btab = [B1  B2];
    Ktab = [K1 K2];
    Gtab = [G1 G2];
    Ctab = [ c  0];
    A    = Atab(ph);
    B    = Btab(ph);
    K    = Ktab(ph);
    G    = Gtab(ph);
    C    = Ctab(ph);

    ur1  = r;
    ur2  = -1/2./r;
    ur   = A.*ur1 + B.*ur2;

    Err = A + 0.5*B./r.^2;
    Epp = ur./r;
    trr = (K+4/3*G).*Err + (K-2/3*G).*Epp - 2*K.*C;    %2 pour cylindrique
    tpp = (K+4/3*G).*Epp + (K-2/3*G).*Err - 2*K.*C;
    tzz = (K-2/3*G).*Epp + (K-2/3*G).*Err - 2*K.*C;
    div  = Err + Epp;
    p    = -1/3*(trr+tpp+tzz);
    
    trr(end)
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
        
    ixv1 = fix((nx+1)/2);
    ixc1 = fix((nx  )/2);
    iyc1 = fix((nx  )/2);
   
    figure(1), clf
    subplot(311), hold on
    plot(r, ur)
    plot(xv(ixv1:end)*Lc, Vx(ixv1:end,iyc1)*dt*Lc, '.')
    title('Ux')
    
    subplot(312), hold on
    plot(r, p-p(1))
    plot(xc(ixc1:end)*Lc, Pt(ixc1:end,iyc1)*tauc - Pt(ixc1,iyc1)*tauc, '.')
    title('P')
    
    subplot(313), hold on
    plot(r, trr-trr(1))
    plot(xc(ixc1:end)*Lc, Sxx(ixc1:end,iyc1)*tauc - Sxx(ixc1,iyc1)*tauc, '.')
    title('Sxx')
    
%     figure(1),clf,colormap('jet'),set(gcf,'Color','white')
%     subplot(221),imagesc(xc,yc,flipud(Pt'  )),axis image,colorbar,caxis([min(Pt(:)) max(Pt(:))]), title('Pt')
%     subplot(223),imagesc(xv,yc,flipud(Vx'  )),axis image,colorbar, title('Vx')
%     subplot(224),imagesc(xc,yv,flipud(Vy'  )),axis image,colorbar, title('Vy')
%     subplot(222),imagesc(xc,yc,flipud(log10(abs(divV')))),axis image,colorbar, title('div V')
%     figure(2),clf,colormap('jet'),set(gcf,'Color','white')
%     imagesc(xc,yc,flipud(Pt1'-Pt'  )),axis image,colorbar
    
%     figure(3), set(gcf,'Color','white'), hold on,
%     tmax = (mus_i*B_i);
%     dPa  = 1 - exp(-3/4*time/tmax);
%     dTxxa= -( 1 - 1/2*exp(-3/4*time/tmax));
%     dmus = musc(1,fix(ny/2)) - musc(1,1);
%     ffs  = abs( 2*exx(1,fix(ny/2))* dmus);
%     dTxx = Txx(1,fix(ny/2)) - Txx(1,1);
%     dP   = Pt(1,fix(ny/2)) - Pt(1,1);
    
%     subplot(211),hold on
%     plot(time/(mus_i*B0), dTxx/ffs,'bo', time/(mus_i*B0), dTxxa,'xr')
%     xlabel('$t$','interpreter','latex'); ylabel('$d\tau{xx}$','interpreter','latex'); l1=legend('$\tau{xx}$','$\tau{xx}$ anal.'); set(l1,'interpreter','latex', 'box', 'off')
%     axis([0 5 -1 -0.5])
%     subplot (212),hold on
%     plot(time/(mus_i*B0), dP/ffs,'bo',time/(mus_i*B0), dPa, 'xr') 
%     axis([0 5 0 1])

%     xlabel('$t$','interpreter','latex'); ylabel('$dP$','interpreter','latex'); l2=legend('$dP$', '$dP$ anal.'); set(l2,'interpreter','latex','location','southeast', 'box', 'off')
%     drawnow
%     fprintf('Strain rate = %2.2e, flow stress = %2.2e\n\n', Ebg, Ebg*mus0*2)
  
end