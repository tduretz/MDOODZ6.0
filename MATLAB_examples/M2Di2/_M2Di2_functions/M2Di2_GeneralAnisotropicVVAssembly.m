function [ K, BcK ] = M2Di_v2_GeneralAnisotropicVVAssembly( D, BC, NumUx, NumUy, NumUyG, nx, ny, dx, dy, symmetry, Ux, Uy, SuiteSparse )
% RHS
BcK = zeros(max(NumUyG(:)),1);
% Normal stress
frxN = BC.frxN; fryN = BC.fryN;
free_surf = 0;
if (sum(frxN(:,end))>0), free_surf = 1; end
% Boundary condition flags [W E S N]
fsxW = BC.fsxW; fsxE = BC.fsxE;
fsxS = BC.fsxS; fsxN = BC.fsxN;
nsxW = BC.nsxW; nsxE = BC.nsxE;
nsxS = BC.nsxS; nsxN = BC.nsxN;
% Boundary condition flags [W E S N]
fsyW = BC.fsyW; fsyE = BC.fsyE;
fsyS = BC.fsyS; fsyN = BC.fsyN;
nsyW = BC.nsyW; nsyE = BC.nsyE;
nsyS = BC.nsyS; nsyN = BC.nsyN;
% Stencil weights for normal strains interpolation
wSW   = ones(nx+1,ny+1); wSW(end,:) = 2; wSW(:,end) = 2; wSW(end,end) = 4; wSW(1  ,:) = 0; wSW(:,  1) = 0;
wSE   = ones(nx+1,ny+1); wSE(  1,:) = 2; wSE(:,end) = 2; wSE(  1,end) = 4; wSE(end,:) = 0; wSE(:,  1) = 0;
wNW   = ones(nx+1,ny+1); wNW(end,:) = 2; wNW(:,  1) = 2; wNW(end,  1) = 4; wNW(1  ,:) = 0; wNW(:,end) = 0;
wNE   = ones(nx+1,ny+1); wNE(  1,:) = 2; wNE(:,  1) = 2; wNE(  1,  1) = 4; wNE(end,:) = 0; wNE(:,end) = 0;
%% Block UU --- v2
% Stencil weights for normal strains interpolation
wS_SW = wSW(:,1:end-1); wS_SE = wSE(:,1:end-1); wS_W  = wNW(:,1:end-1); wS_E  = wNE(:,1:end-1); 
wN_W  = wSW(:,2:end-0); wN_E  = wSE(:,2:end-0); wN_NW = wNW(:,2:end-0); wN_NE = wNE(:,2:end-0);
% Elastic or elasto-plastic operators West-East
D11W    = zeros(size(Ux));  D11W( 2:end-1, : ) = D.D11c(1:end-1,:);
D11E    = zeros(size(Ux));  D11E( 2:end-1, : ) = D.D11c(2:end,:);
D12W    = zeros(size(Ux));  D12W( 2:end-1, : ) = D.D12c(1:end-1,:);
D12E    = zeros(size(Ux));  D12E( 2:end-1, : ) = D.D12c(2:end,:);
D13W    = zeros(size(Ux));  D13W( 2:end-1, : ) = D.D13c(1:end-1,:);
D13E    = zeros(size(Ux));  D13E( 2:end-1, : ) = D.D13c(2:end,:);
% Elastic or elasto-plastic operators South-North
D31S    = zeros(size(Ux));  D31S( 2:end-1 , :) = D.D31v(2:end-1, 1:end-1);  %D31S(:,1) = 0;
D31N    = zeros(size(Ux));  D31N( 2:end-1 , :) = D.D31v(2:end-1, 2:end  );  %D31N(:,end) = 0;
D32S    = zeros(size(Ux));  D32S( 2:end-1 , :) = D.D32v(2:end-1, 1:end-1);  %D32S(:,1) = 0;
D32N    = zeros(size(Ux));  D32N( 2:end-1 , :) = D.D32v(2:end-1, 2:end  );  %D32N(:,end) = 0;
D33S    = zeros(size(Ux));  D33S( 2:end-1 , :) = D.D33v(2:end-1, 1:end-1);  %D33S(:,1) = 0;
D33N    = zeros(size(Ux));  D33N( 2:end-1 , :) = D.D33v(2:end-1, 2:end  );  %D33N(:,end) = 0;
% x - mommentum
c1UxC = [-0.1e1 ./ dx .* (-0.2e1 ./ 0.3e1 .* D11E ./ dx + D12E ./ dx ./ 0.3e1 + D13E .* (((1 - fsxN) .* (-1 - nsxN) ./ dy) ./ 0.4e1 + ((1 - fsxS) .* (1 + nsxS) ./ dy) ./ 0.4e1) - 0.2e1 ./ 0.3e1 .* D11W ./ dx + D12W ./ dx ./ 0.3e1 - D13W .* (((1 - fsxN) .* (-1 - nsxN) ./ dy) ./ 0.4e1 + ((1 - fsxS) .* (1 + nsxS) ./ dy) ./ 0.4e1)) - 0.1e1 ./ dy .* (D31N .* (-wN_E ./ dx ./ 0.6e1 + wN_W ./ dx ./ 0.6e1) + D32N .* (wN_E ./ dx ./ 0.12e2 - wN_W ./ dx ./ 0.12e2) + (D33N .* (1 - fsxN) .* (-1 - nsxN) ./ dy) - D31S .* (-wS_E ./ dx ./ 0.6e1 + wS_W ./ dx ./ 0.6e1) - D32S .* (wS_E ./ dx ./ 0.12e2 - wS_W ./ dx ./ 0.12e2) - (D33S .* (1 - fsxS) .* (1 + nsxS) ./ dy))];
c1UxW = [-0.1e1 ./ dx .* (0.2e1 ./ 0.3e1 .* D11W ./ dx - D12W ./ dx ./ 0.3e1 - D13W .* (((1 - fsxN) .* (-1 - nsxN) ./ dy) ./ 0.4e1 + ((1 - fsxS) .* (1 + nsxS) ./ dy) ./ 0.4e1)) - 0.1e1 ./ dy .* (-D31N .* wN_W ./ dx ./ 0.6e1 + D32N .* wN_W ./ dx ./ 0.12e2 + D31S .* wS_W ./ dx ./ 0.6e1 - D32S .* wS_W ./ dx ./ 0.12e2)];
c1UxE = [-0.1e1 ./ dx .* (0.2e1 ./ 0.3e1 .* D11E ./ dx - D12E ./ dx ./ 0.3e1 + D13E .* (((1 - fsxN) .* (-1 - nsxN) ./ dy) ./ 0.4e1 + ((1 - fsxS) .* (1 + nsxS) ./ dy) ./ 0.4e1)) - 0.1e1 ./ dy .* (D31N .* wN_E ./ dx ./ 0.6e1 - D32N .* wN_E ./ dx ./ 0.12e2 - D31S .* wS_E ./ dx ./ 0.6e1 + D32S .* wS_E ./ dx ./ 0.12e2)];
c1UxS = [-0.1e1 ./ dx .* ((D13E .* (1 - fsxS) .* (-1 + nsxS) ./ dy) ./ 0.4e1 - (D13W .* (1 - fsxS) .* (-1 + nsxS) ./ dy) ./ 0.4e1) - 0.1e1 ./ dy .* (-D31S .* (wS_SW ./ dx ./ 0.6e1 - wS_SE ./ dx ./ 0.6e1) - D32S .* (-wS_SW ./ dx ./ 0.12e2 + wS_SE ./ dx ./ 0.12e2) - (D33S .* (1 - fsxS) .* (-1 + nsxS) ./ dy))];
c1UxN = [-0.1e1 ./ dx .* ((D13E .* (1 - fsxN) .* (1 - nsxN) ./ dy) ./ 0.4e1 - (D13W .* (1 - fsxN) .* (1 - nsxN) ./ dy) ./ 0.4e1) - 0.1e1 ./ dy .* (D31N .* (wN_NW .* (1 - frxN) ./ dx ./ 0.6e1 - wN_NE .* (1 - frxN) ./ dx ./ 0.6e1) + D32N .* (-wN_NW .* (1 - frxN) ./ dx ./ 0.12e2 + wN_NE .* (1 - frxN) ./ dx ./ 0.12e2) + (D33N .* (1 - fsxN) .* (1 - nsxN) ./ dy))];
c1UySW = [-0.1e1 ./ dx .* (-D13E ./ dx ./ 0.4e1 - D11W ./ dy ./ 0.3e1 + 0.2e1 ./ 0.3e1 .* D12W ./ dy - D13W .* (-0.1e1 ./ dx ./ 0.4e1 + (1 - fsxW) .* (1 + nsxW) ./ dx ./ 0.4e1)) - 0.1e1 ./ dy .* (D31N .* wN_W ./ dy ./ 0.12e2 - D32N .* wN_W ./ dy ./ 0.6e1 - D31S .* (wS_W ./ dy ./ 0.12e2 - wS_SW ./ dy ./ 0.12e2) - D32S .* (-wS_W ./ dy ./ 0.6e1 + wS_SW ./ dy ./ 0.6e1) + D33S ./ dx)];
c1UySE = [-0.1e1 ./ dx .* (D11E ./ dy ./ 0.3e1 - 0.2e1 ./ 0.3e1 .* D12E ./ dy + D13E .* (0.1e1 ./ dx ./ 0.4e1 + (1 - fsxE) .* (-1 - nsxE) ./ dx ./ 0.4e1) - D13W ./ dx ./ 0.4e1) - 0.1e1 ./ dy .* (D31N .* wN_E ./ dy ./ 0.12e2 - D32N .* wN_E ./ dy ./ 0.6e1 - D31S .* (wS_E ./ dy ./ 0.12e2 - wS_SE ./ dy ./ 0.12e2) - D32S .* (-wS_E ./ dy ./ 0.6e1 + wS_SE ./ dy ./ 0.6e1) - D33S ./ dx)];
c1UyNW = [-0.1e1 ./ dx .* (-D13E ./ dx ./ 0.4e1 + D11W ./ dy ./ 0.3e1 - 0.2e1 ./ 0.3e1 .* D12W ./ dy - D13W .* (-0.1e1 ./ dx ./ 0.4e1 + (1 - fsxW) .* (1 + nsxW) ./ dx ./ 0.4e1)) - 0.1e1 ./ dy .* (D31N .* (-wN_W ./ dy ./ 0.12e2 + wN_NW .* (1 - frxN) ./ dy ./ 0.12e2) + D32N .* (wN_W ./ dy ./ 0.6e1 - wN_NW .* (1 - frxN) ./ dy ./ 0.6e1) - D33N ./ dx + D31S .* wS_W ./ dy ./ 0.12e2 - D32S .* wS_W ./ dy ./ 0.6e1)];
c1UyNE = [-0.1e1 ./ dx .* (-D11E ./ dy ./ 0.3e1 + 0.2e1 ./ 0.3e1 .* D12E ./ dy + D13E .* (0.1e1 ./ dx ./ 0.4e1 + (1 - fsxE) .* (-1 - nsxE) ./ dx ./ 0.4e1) - D13W ./ dx ./ 0.4e1) - 0.1e1 ./ dy .* (D31N .* (-wN_E ./ dy ./ 0.12e2 + wN_NE .* (1 - frxN) ./ dy ./ 0.12e2) + D32N .* (wN_E ./ dy ./ 0.6e1 - wN_NE .* (1 - frxN) ./ dy ./ 0.6e1) + D33N ./ dx + D31S .* wS_E ./ dy ./ 0.12e2 - D32S .* wS_E ./ dy ./ 0.6e1)];
c1UxSW = [(1 ./ dx .* D13W .* (1 - fsxS) .* (-1 + nsxS) ./ dy) ./ 0.4e1 - 0.1e1 ./ dy .* ((D31S .* wS_SW ./ dx) ./ 0.6e1 - (D32S .* wS_SW ./ dx) ./ 0.12e2)];
c1UxSE = [-(1 ./ dx .* D13E .* (1 - fsxS) .* (-1 + nsxS) ./ dy) ./ 0.4e1 - 0.1e1 ./ dy .* (-(D31S .* wS_SE ./ dx) ./ 0.6e1 + (D32S .* wS_SE ./ dx) ./ 0.12e2)];
c1UxNW = [(1 ./ dx .* D13W .* (1 - fsxN) .* (1 - nsxN) ./ dy) ./ 0.4e1 - 0.1e1 ./ dy .* (-(D31N .* wN_NW .* (1 - frxN) ./ dx) ./ 0.6e1 + (D32N .* wN_NW .* (1 - frxN) ./ dx) ./ 0.12e2)];
c1UxNE = [-(1 ./ dx .* D13E .* (1 - fsxN) .* (1 - nsxN) ./ dy) ./ 0.4e1 - 0.1e1 ./ dy .* ((D31N .* wN_NE .* (1 - frxN) ./ dx) ./ 0.6e1 - (D32N .* wN_NE .* (1 - frxN) ./ dx) ./ 0.12e2)];
c1UySSW = [-0.1e1 ./ dy .* (-D31S .* wS_SW ./ dy ./ 0.12e2 + D32S .* wS_SW ./ dy ./ 0.6e1)];
c1UySSE = [-0.1e1 ./ dy .* (-D31S .* wS_SE ./ dy ./ 0.12e2 + D32S .* wS_SE ./ dy ./ 0.6e1)];
c1UySWW = [(1 ./ dx .^ 2 .* D13W .* (1 - fsxW) .* (-1 + nsxW)) ./ 0.4e1];
c1UySEE = [-(1 ./ dx .^ 2 .* D13E .* (1 - fsxE) .* (1 - nsxE)) ./ 0.4e1];
c1UyNWW = [(1 ./ dx .^ 2 .* D13W .* (1 - fsxW) .* (-1 + nsxW)) ./ 0.4e1];
c1UyNEE = [-(1 ./ dx .^ 2 .* D13E .* (1 - fsxE) .* (1 - nsxE)) ./ 0.4e1];
c1UyNNW = [-0.1e1 ./ dy .* (-(D31N .* wN_NW .* (1 - frxN) ./ dy) ./ 0.12e2 + (D32N .* wN_NW .* (1 - frxN) ./ dy) ./ 0.6e1)];
c1UyNNE = [-0.1e1 ./ dy .* (-(D31N .* wN_NE .* (1 - frxN) ./ dy) ./ 0.12e2 + (D32N .* wN_NE .* (1 - frxN) ./ dy) ./ 0.6e1)];
if symmetry == 1
    % Contribution from normal stresses (x-momentum)
    BcK(NumUx( 2   , : ))  = BcK(NumUx(2,:))      - c1UxW( 2    , : )'.*BC.Ux_W;
    BcK(NumUx(end-1, : ))  = BcK(NumUx(end-1,:))  - c1UxE( end-1, : )'.*BC.Ux_E;
    % Contribution from shear stresses (x-momentum)
    BcK(NumUx(2:end-1,end))      = BcK(NumUx(2:end-1,end))  -  c1UyNW (2:end-1,end).*BC.Uy_N(1:end-1)   ; % UyNW
    BcK(NumUx(2:end-1,end))      = BcK(NumUx(2:end-1,end))  -  c1UyNE (2:end-1,end).*BC.Uy_N(2:end  )   ; % UyNE
    BcK(NumUx(2:end-1, 1 ))      = BcK(NumUx(2:end-1, 1 ))  -  c1UySW (2:end-1,1  ).*BC.Uy_S(1:end-1)   ; % UySW
    BcK(NumUx(2:end-1, 1 ))      = BcK(NumUx(2:end-1, 1 ))  -  c1UySE (2:end-1,1  ).*BC.Uy_S(2:end  )   ; % UySE
end
% Non conforming Dirichlets BCs
nsxN               = nsxN.*0;
nsxS               = nsxS.*0;
c1UxS1 = [-0.1e1 ./ dx .* ((D13E .* (1 - fsxS) .* (-1 + nsxS) ./ dy) ./ 0.4e1 - (D13W .* (1 - fsxS) .* (-1 + nsxS) ./ dy) ./ 0.4e1) - 0.1e1 ./ dy .* (-D31S .* (wS_SW ./ dx ./ 0.6e1 - wS_SE ./ dx ./ 0.6e1) - D32S .* (-wS_SW ./ dx ./ 0.12e2 + wS_SE ./ dx ./ 0.12e2) - (D33S .* (1 - fsxS) .* (-1 + nsxS) ./ dy))];
c1UxN1 = [-0.1e1 ./ dx .* ((D13E .* (1 - fsxN) .* (1 - nsxN) ./ dy) ./ 0.4e1 - (D13W .* (1 - fsxN) .* (1 - nsxN) ./ dy) ./ 0.4e1) - 0.1e1 ./ dy .* (D31N .* (wN_NW .* (1 - frxN) ./ dx ./ 0.6e1 - wN_NE .* (1 - frxN) ./ dx ./ 0.6e1) + D32N .* (-wN_NW .* (1 - frxN) ./ dx ./ 0.12e2 + wN_NE .* (1 - frxN) ./ dx ./ 0.12e2) + (D33N .* (1 - fsxN) .* (1 - nsxN) ./ dy))];
if sum(BC.nsxN(:))>0, BcK(NumUx(:,end)) = BcK(NumUx(:,end)) - 2*c1UxN1(NumUx(:,end)).*BC.Ux_N; end
if sum(BC.nsxS(:))>0, BcK(NumUx(:,  1)) = BcK(NumUx(:,  1)) - 2*c1UxS1(NumUx(:,  1)).*BC.Ux_S; end
% Set boundary conditions
sc1Ux = max(c1UxC(:));
c1UxC  ([1 end],:) = sc1Ux;
c1UxW  ([1 end],:) = 0;
c1UxE  ([1 end],:) = 0;
c1UxS  ([1 end],:) = 0; c1UxS ( :   , 1 ) = 0;
c1UxN  ([1 end],:) = 0; c1UxN ( :   ,end) = 0;
c1UxSW ([1 end],:) = 0; c1UxSW( :   , 1 ) = 0;
c1UxSE ([1 end],:) = 0; c1UxSE( :   , 1 ) = 0;
c1UxNW ([1 end],:) = 0; c1UxNW( :   ,end) = 0;
c1UxNE ([1 end],:) = 0; c1UxNE( :   ,end) = 0;
c1UySW ([1 end],:) = 0;
c1UySE ([1 end],:) = 0;
c1UyNW ([1 end],:) = 0;
c1UyNE ([1 end],:) = 0;
c1UySWW([1 end],:) = 0; c1UySWW([1 2],:) = 0;
c1UySEE([1 end],:) = 0; c1UySEE([end end-1],:) = 0;
c1UyNWW([1 end],:) = 0; c1UyNWW([1 2],:) = 0;
c1UyNEE([1 end],:) = 0; c1UyNEE([end end-1 ],:) = 0;
c1UySSW([1 end],:) = 0; c1UySSW(:,1  ) = 0;
c1UySSE([1 end],:) = 0; c1UySSE(:,1  ) = 0;
c1UyNNW([1 end],:) = 0; c1UyNNW(:,end) = 0;
c1UyNNE([1 end],:) = 0; c1UyNNE(:,end) = 0;
% Symmetry
if symmetry == 1
    c1UxW  (    2,:) = 0;
    c1UxE  (end-1,:) = 0;
    c1UySW (  :,  1) = 0;
    c1UySE (  :,  1) = 0;
    c1UyNW (  :,end) = 0;
    c1UyNE (  :,end) = 0;
end
% Equation indexes
iUx     = NumUx(:);
iUxW    = ones(size(Ux));   iUxW(2:end-1, :     ) =  NumUx(1:end-2, :     );
iUxE    = ones(size(Ux));   iUxE(2:end-1, :     ) =  NumUx(3:end  , :     );
iUxS    = ones(size(Ux));   iUxS(2:end-1,2:end  ) =  NumUx(2:end-1,1:end-1);
iUxN    = ones(size(Ux));   iUxN(2:end-1,1:end-1) =  NumUx(2:end-1,2:end  );
iUxSW   = ones(size(Ux));  iUxSW(2:end-1,2:end  ) =  NumUx(1:end-2,1:end-1);
iUxSE   = ones(size(Ux));  iUxSE(2:end-1,2:end  ) =  NumUx(3:end-0,1:end-1);
iUxNW   = ones(size(Ux));  iUxNW(2:end-1,1:end-1) =  NumUx(1:end-2,2:end  );
iUxNE   = ones(size(Ux));  iUxNE(2:end-1,1:end-1) =  NumUx(3:end-0,2:end  );
iUySW   = ones(size(Ux));  iUySW(2:end-1, :     ) = NumUyG(1:end-1,1:end-1);
iUySE   = ones(size(Ux));  iUySE(2:end-1, :     ) = NumUyG(2:end-0,1:end-1);
iUyNW   = ones(size(Ux));  iUyNW(2:end-1, :     ) = NumUyG(1:end-1,2:end  );
iUyNE   = ones(size(Ux));  iUyNE(2:end-1, :     ) = NumUyG(2:end-0,2:end  );
iUySWW  = ones(size(Ux)); iUySWW(3:end-1, :     ) = NumUyG(1:end-2,1:end-1);
iUySEE  = ones(size(Ux)); iUySEE(2:end-2, :     ) = NumUyG(3:end  ,1:end-1);
iUyNWW  = ones(size(Ux)); iUyNWW(3:end-1, :     ) = NumUyG(1:end-2,2:end  );
iUyNEE  = ones(size(Ux)); iUyNEE(2:end-2, :     ) = NumUyG(3:end  ,2:end  );
iUySSW  = ones(size(Ux)); iUySSW(2:end-1, 2:end ) = NumUyG(1:end-1,1:end-2);
iUySSE  = ones(size(Ux)); iUySSE(2:end-1, 2:end ) = NumUyG(2:end-0,1:end-2);
iUyNNW  = ones(size(Ux)); iUyNNW(2:end-1 ,1:end-1) = NumUyG(1:end-1,3:end  );
iUyNNE  = ones(size(Ux)); iUyNNE(2:end-1,1:end-1) = NumUyG(2:end-0,3:end  );
IuuJ_2    = [   iUx(:);   iUx(:);   iUx(:);   iUx(:);   iUx(:);    iUx(:);    iUx(:);    iUx(:);    iUx(:);    iUx(:);    iUx(:);    iUx(:);    iUx(:);     iUx(:);     iUx(:);     iUx(:);     iUx(:);     iUx(:);     iUx(:);     iUx(:);     iUx(:) ]';
JuuJ_2    = [   iUx(:);  iUxW(:);  iUxS(:);  iUxE(:);  iUxN(:);  iUxSW(:);  iUxSE(:);  iUxNW(:);  iUxNE(:);  iUySW(:);  iUySE(:);  iUyNW(:);  iUyNE(:);  iUySWW(:);  iUySEE(:);  iUyNWW(:);  iUyNEE(:);  iUySSW(:);  iUySSE(:);  iUyNNW(:);  iUyNNE(:) ]';
VuuJ_2    = [ c1UxC(:); c1UxW(:); c1UxS(:); c1UxE(:); c1UxN(:); c1UxSW(:); c1UxSE(:); c1UxNW(:); c1UxNE(:); c1UySW(:); c1UySE(:); c1UyNW(:); c1UyNE(:); c1UySWW(:); c1UySEE(:); c1UyNWW(:); c1UyNEE(:); c1UySSW(:); c1UySSE(:); c1UyNNW(:); c1UyNNE(:) ]';
%% Block VV --- v2
% Elastic or elasto-plastic operators South-North
if free_surf == 0
    D21S    = zeros(size(Uy));  D21S( : ,2:end-1) = D.D21c(:,1:end-1);% D21S(:,1) = 0;
    D21N    = zeros(size(Uy));  D21N( : ,2:end-1) = D.D21c(:,2:end-0); %D21N(:,end) = 0;
    D22S    = zeros(size(Uy));  D22S( : ,2:end-1) = D.D22c(:,1:end-1);% D22S(:,1) = 0;
    D22N    = zeros(size(Uy));  D22N( : ,2:end-1) = D.D22c(:,2:end-0); %D21N(:,1) = 0;
    D23S    = zeros(size(Uy));  D23S( : ,2:end-1) = D.D23c(:,1:end-1);% D23S(:,1) = 0;
    D23N    = zeros(size(Uy));  D23N( : ,2:end-1) = D.D23c(:,2:end-0);% D23N(:,end) = 0
else
    D21S    = zeros(size(Uy));  D21S( : ,2:end-0) = D.D21c(:,1:end-0);% D21S(:,1) = 0;
    D21N    = zeros(size(Uy));  D21N( : ,2:end-1) = D.D21c(:,2:end-0); %D21N(:,end) = 0;
    D22S    = zeros(size(Uy));  D22S( : ,2:end-0) = D.D22c(:,1:end-0);% D22S(:,1) = 0;
    D22N    = zeros(size(Uy));  D22N( : ,2:end-1) = D.D22c(:,2:end-0); %D21N(:,1) = 0;
    D23S    = zeros(size(Uy));  D23S( : ,2:end-0) = D.D23c(:,1:end-0);% D23S(:,1) = 0;
    D23N    = zeros(size(Uy));  D23N( : ,2:end-1) = D.D23c(:,2:end-0);% D23N(:,end) = 0;    
end
% Elastic or elasto-plastic operators West-East
D31W    = zeros(size(Uy));  D31W( : ,2:end-1 ) = D.D31v(1:end-1, 2:end-1); %D31W(1,:) = 0;
D31E    = zeros(size(Uy));  D31E( : ,2:end-1 ) = D.D31v(2:end  , 2:end-1); %D31E(end,:) = 0;
D32W    = zeros(size(Uy));  D32W( : ,2:end-1 ) = D.D32v(1:end-1, 2:end-1); %D32W(1,:) = 0;
D32E    = zeros(size(Uy));  D32E( : ,2:end-1 ) = D.D32v(2:end  , 2:end-1); %D32E(end,:) = 0;
D33W    = zeros(size(Uy));  D33W( : ,2:end-1 ) = D.D33v(1:end-1, 2:end-1); %D33W(1,:) = 0;
D33E    = zeros(size(Uy));  D33E( : ,2:end-1 ) = D.D33v(2:end  , 2:end-1); %D33E(end,:) = 0;
% Stencil weights for normal strains interpolation
wW_SW = wSW(1:end-1,:); wW_NW = wNW(1:end-1,:); wW_S  = wSE(1:end-1,:); wW_N  = wNW(1:end-1,:); 
wE_S  = wSW(2:end-0,:); wE_N  = wSW(2:end-0,:); wE_SE = wSE(2:end-0,:); wE_NE = wNE(2:end-0,:);
% y - momentum
c2UyC = [-0.1e1 ./ dy .* ((D21N .* (1 - fryN) ./ dy) ./ 0.3e1 - 0.2e1 ./ 0.3e1 .* D22N .* (1 - fryN) ./ dy + D23N .* (((1 - fsyW) .* (1 + nsyW) ./ dx) ./ 0.4e1 + ((1 - fsyE) .* (-1 - nsyE) ./ dx) ./ 0.4e1) + (D21S ./ dy) ./ 0.3e1 - 0.2e1 ./ 0.3e1 .* D22S ./ dy - D23S .* (((1 - fsyW) .* (1 + nsyW) ./ dx) ./ 0.4e1 + ((1 - fsyE) .* (-1 - nsyE) ./ dx) ./ 0.4e1)) - 0.1e1 ./ dx .* (D31E .* (-(wE_S ./ dy) ./ 0.12e2 + (wE_N .* (1 - fryN) ./ dy) ./ 0.12e2) + D32E .* ((wE_S ./ dy) ./ 0.6e1 - (wE_N .* (1 - fryN) ./ dy) ./ 0.6e1) + (D33E .* (1 - fsyE) .* (-1 - nsyE) ./ dx) - D31W .* (-(wW_S ./ dy) ./ 0.12e2 + (wW_N .* (1 - fryN) ./ dy) ./ 0.12e2) - D32W .* ((wW_S ./ dy) ./ 0.6e1 - (wW_N .* (1 - fryN) ./ dy) ./ 0.6e1) - (D33W .* (1 - fsyW) .* (1 + nsyW) ./ dx))];
c2UyW = [-0.1e1 ./ dy .* ((D23N .* (1 - fsyW) .* (-1 + nsyW) ./ dx) ./ 0.4e1 - (D23S .* (1 - fsyW) .* (-1 + nsyW) ./ dx) ./ 0.4e1) - 0.1e1 ./ dx .* (-D31W .* (wW_NW .* (1 - fryN) ./ dy ./ 0.12e2 - wW_SW ./ dy ./ 0.12e2) - D32W .* (-wW_NW .* (1 - fryN) ./ dy ./ 0.6e1 + wW_SW ./ dy ./ 0.6e1) - (D33W .* (1 - fsyW) .* (-1 + nsyW) ./ dx))];
c2UyE = [-0.1e1 ./ dy .* ((D23N .* (1 - fsyE) .* (1 - nsyE) ./ dx) ./ 0.4e1 - (D23S .* (1 - fsyE) .* (1 - nsyE) ./ dx) ./ 0.4e1) - 0.1e1 ./ dx .* (D31E .* (wE_NE .* (1 - fryN) ./ dy ./ 0.12e2 - wE_SE ./ dy ./ 0.12e2) + D32E .* (-wE_NE .* (1 - fryN) ./ dy ./ 0.6e1 + wE_SE ./ dy ./ 0.6e1) + (D33E .* (1 - fsyE) .* (1 - nsyE) ./ dx))];
c2UyS = [-0.1e1 ./ dy .* (-D21S ./ dy ./ 0.3e1 + 0.2e1 ./ 0.3e1 .* D22S ./ dy - D23S .* (((1 - fsyW) .* (1 + nsyW) ./ dx) ./ 0.4e1 + ((1 - fsyE) .* (-1 - nsyE) ./ dx) ./ 0.4e1)) - 0.1e1 ./ dx .* (D31E .* wE_S ./ dy ./ 0.12e2 - D32E .* wE_S ./ dy ./ 0.6e1 - D31W .* wW_S ./ dy ./ 0.12e2 + D32W .* wW_S ./ dy ./ 0.6e1)];
c2UyN = [-0.1e1 ./ dy .* (-(D21N .* (1 - fryN) ./ dy) ./ 0.3e1 + 0.2e1 ./ 0.3e1 .* D22N .* (1 - fryN) ./ dy + D23N .* (((1 - fsyW) .* (1 + nsyW) ./ dx) ./ 0.4e1 + ((1 - fsyE) .* (-1 - nsyE) ./ dx) ./ 0.4e1)) - 0.1e1 ./ dx .* (-(D31E .* wE_N .* (1 - fryN) ./ dy) ./ 0.12e2 + (D32E .* wE_N .* (1 - fryN) ./ dy) ./ 0.6e1 + (D31W .* wW_N .* (1 - fryN) ./ dy) ./ 0.12e2 - (D32W .* wW_N .* (1 - fryN) ./ dy) ./ 0.6e1)];
c2UxSW = [-0.1e1 ./ dy .* (-D23N ./ dy ./ 0.4e1 + 0.2e1 ./ 0.3e1 .* D21S ./ dx - D22S ./ dx ./ 0.3e1 - D23S .* (-0.1e1 ./ dy ./ 0.4e1 + (1 - fsyS) .* (1 + nsyS) ./ dy ./ 0.4e1)) - 0.1e1 ./ dx .* (-D31E .* wE_S ./ dx ./ 0.6e1 + D32E .* wE_S ./ dx ./ 0.12e2 - D31W .* (-wW_S ./ dx ./ 0.6e1 + wW_SW ./ dx ./ 0.6e1) - D32W .* (wW_S ./ dx ./ 0.12e2 - wW_SW ./ dx ./ 0.12e2) + D33W ./ dy)];
c2UxSE = [-0.1e1 ./ dy .* (-D23N ./ dy ./ 0.4e1 - 0.2e1 ./ 0.3e1 .* D21S ./ dx + D22S ./ dx ./ 0.3e1 - D23S .* (-0.1e1 ./ dy ./ 0.4e1 + (1 - fsyS) .* (1 + nsyS) ./ dy ./ 0.4e1)) - 0.1e1 ./ dx .* (D31E .* (wE_S ./ dx ./ 0.6e1 - wE_SE ./ dx ./ 0.6e1) + D32E .* (-wE_S ./ dx ./ 0.12e2 + wE_SE ./ dx ./ 0.12e2) - D33E ./ dy - D31W .* wW_S ./ dx ./ 0.6e1 + D32W .* wW_S ./ dx ./ 0.12e2)];
c2UxNW = [-0.1e1 ./ dy .* (-0.2e1 ./ 0.3e1 .* D21N .* (1 - fryN) ./ dx + (D22N .* (1 - fryN) ./ dx) ./ 0.3e1 + D23N .* (0.1e1 ./ dy ./ 0.4e1 + (1 - fsyN) .* (-1 - nsyN) ./ dy ./ 0.4e1) - D23S ./ dy ./ 0.4e1) - 0.1e1 ./ dx .* (-(D31E .* wE_N .* (1 - fryN) ./ dx) ./ 0.6e1 + (D32E .* wE_N .* (1 - fryN) ./ dx) ./ 0.12e2 - D31W .* (-(wW_N .* (1 - fryN) ./ dx) ./ 0.6e1 + (wW_NW .* (1 - fryN) ./ dx) ./ 0.6e1) - D32W .* ((wW_N .* (1 - fryN) ./ dx) ./ 0.12e2 - (wW_NW .* (1 - fryN) ./ dx) ./ 0.12e2) - D33W ./ dy)];
c2UxNE = [-0.1e1 ./ dy .* (0.2e1 ./ 0.3e1 .* D21N .* (1 - fryN) ./ dx - (D22N .* (1 - fryN) ./ dx) ./ 0.3e1 + D23N .* (0.1e1 ./ dy ./ 0.4e1 + (1 - fsyN) .* (-1 - nsyN) ./ dy ./ 0.4e1) - D23S ./ dy ./ 0.4e1) - 0.1e1 ./ dx .* (D31E .* ((wE_N .* (1 - fryN) ./ dx) ./ 0.6e1 - (wE_NE .* (1 - fryN) ./ dx) ./ 0.6e1) + D32E .* (-(wE_N .* (1 - fryN) ./ dx) ./ 0.12e2 + (wE_NE .* (1 - fryN) ./ dx) ./ 0.12e2) + D33E ./ dy - (D31W .* wW_N .* (1 - fryN) ./ dx) ./ 0.6e1 + (D32W .* wW_N .* (1 - fryN) ./ dx) ./ 0.12e2)];
c2UySW = [(1 ./ dy .* D23S .* (1 - fsyW) .* (-1 + nsyW) ./ dx) ./ 0.4e1 - 0.1e1 ./ dx .* (-(D31W .* wW_SW ./ dy) ./ 0.12e2 + (D32W .* wW_SW ./ dy) ./ 0.6e1)];
c2UySE = [(1 ./ dy .* D23S .* (1 - fsyE) .* (1 - nsyE) ./ dx) ./ 0.4e1 - 0.1e1 ./ dx .* ((D31E .* wE_SE ./ dy) ./ 0.12e2 - (D32E .* wE_SE ./ dy) ./ 0.6e1)];
c2UyNW = [-(1 ./ dy .* D23N .* (1 - fsyW) .* (-1 + nsyW) ./ dx) ./ 0.4e1 - 0.1e1 ./ dx .* ((D31W .* wW_NW .* (1 - fryN) ./ dy) ./ 0.12e2 - (D32W .* wW_NW .* (1 - fryN) ./ dy) ./ 0.6e1)];
c2UyNE = [-(1 ./ dy .* D23N .* (1 - fsyE) .* (1 - nsyE) ./ dx) ./ 0.4e1 - 0.1e1 ./ dx .* (-(D31E .* wE_NE .* (1 - fryN) ./ dy) ./ 0.12e2 + (D32E .* wE_NE .* (1 - fryN) ./ dy) ./ 0.6e1)];
c2UxSSW = [(1 ./ dy .^ 2 .* D23S .* (1 - fsyS) .* (-1 + nsyS)) ./ 0.4e1];
c2UxSSE = [(1 ./ dy .^ 2 .* D23S .* (1 - fsyS) .* (-1 + nsyS)) ./ 0.4e1];
c2UxSWW = [-0.1e1 ./ dx .* (D31W .* wW_SW ./ dx ./ 0.6e1 - D32W .* wW_SW ./ dx ./ 0.12e2)];
c2UxSEE = [-0.1e1 ./ dx .* (D31E .* wE_SE ./ dx ./ 0.6e1 - D32E .* wE_SE ./ dx ./ 0.12e2)];
c2UxNWW = [-0.1e1 ./ dx .* ((D31W .* wW_NW .* (1 - fryN) ./ dx) ./ 0.6e1 - (D32W .* wW_NW .* (1 - fryN) ./ dx) ./ 0.12e2)];
c2UxNEE = [-0.1e1 ./ dx .* ((D31E .* wE_NE .* (1 - fryN) ./ dx) ./ 0.6e1 - (D32E .* wE_NE .* (1 - fryN) ./ dx) ./ 0.12e2)];
c2UxNNW = [-(1 ./ dy .^ 2 .* D23N .* (1 - fsyN) .* (1 - nsyN)) ./ 0.4e1];
c2UxNNE = [-(1 ./ dy .^ 2 .* D23N .* (1 - fsyN) .* (1 - nsyN)) ./ 0.4e1];
if symmetry == 1
    % Contribution from normal stresses
    BcK(NumUyG( : , 2   )) = BcK(NumUyG(:,2))     - c2UyS( : , 2   ).*BC.Uy_S;
    BcK(NumUyG( : ,end-1)) = BcK(NumUyG(:,end-1)) - c2UyN( : ,end-1).*BC.Uy_N;
    % Contribution from shear stresses (y-momentum)
    BcK(NumUyG(end,2:end-1))     = BcK(NumUyG(end,2:end-1)) - (c2UxSE (end,2:end-1).*BC.Ux_E(1:end-1)')'; % UxSE
    BcK(NumUyG(end,2:end-1))     = BcK(NumUyG(end,2:end-1)) - (c2UxNE (end,2:end-1).*BC.Ux_E(2:end  )')'; % UxNE
    BcK(NumUyG( 1 ,2:end-1))     = BcK(NumUyG( 1 ,2:end-1)) - (c2UxSW (1  ,2:end-1).*BC.Ux_W(1:end-1)')'; % UxSW
    BcK(NumUyG( 1 ,2:end-1))     = BcK(NumUyG( 1 ,2:end-1)) - (c2UxNW (1  ,2:end-1).*BC.Ux_W(2:end  )')'; % UxNW
end
% Non conforming Dirichlets BCs
nsyW               = 0.*nsyW;
nsyE               = 0.*nsyE;
c2UyW1 = [-0.1e1 ./ dy .* ((D23N .* (1 - fsyW) .* (-1 + nsyW) ./ dx) ./ 0.4e1 - (D23S .* (1 - fsyW) .* (-1 + nsyW) ./ dx) ./ 0.4e1) - 0.1e1 ./ dx .* (-D31W .* (wW_NW .* (1 - fryN) ./ dy ./ 0.12e2 - wW_SW ./ dy ./ 0.12e2) - D32W .* (-wW_NW .* (1 - fryN) ./ dy ./ 0.6e1 + wW_SW ./ dy ./ 0.6e1) - (D33W .* (1 - fsyW) .* (-1 + nsyW) ./ dx))];
c2UyE1 = [-0.1e1 ./ dy .* ((D23N .* (1 - fsyE) .* (1 - nsyE) ./ dx) ./ 0.4e1 - (D23S .* (1 - fsyE) .* (1 - nsyE) ./ dx) ./ 0.4e1) - 0.1e1 ./ dx .* (D31E .* (wE_NE .* (1 - fryN) ./ dy ./ 0.12e2 - wE_SE ./ dy ./ 0.12e2) + D32E .* (-wE_NE .* (1 - fryN) ./ dy ./ 0.6e1 + wE_SE ./ dy ./ 0.6e1) + (D33E .* (1 - fsyE) .* (1 - nsyE) ./ dx))];
if sum(BC.nsyE(:))>0, BcK(NumUyG(end, :)) = BcK(NumUyG(end, :)) - 2*c2UyE1(NumUy(end, :))'.*BC.Uy_E; end
if sum(BC.nsyW(:))>0, BcK(NumUyG(  1, :)) = BcK(NumUyG(  1, :)) - 2*c2UyW1(NumUy(  1, :))'.*BC.Uy_W; end
% Set boundary conditions and symmetrise
if free_surf == 0
    sc2Uy = max(c2UyC(:));
    c2UyC  (:,[1 end]) = sc2Uy;
    c2UyW  (:,[1 end]) = 0; c2UyW  (1,  :) = 0;
    c2UyE  (:,[1 end]) = 0; c2UyE  (end,:) = 0;
    c2UyS  (:,[1 end]) = 0;
    c2UyN  (:,[1 end]) = 0;
    c2UySW (:,[1 end]) = 0; c2UySW(1  , :) = 0;
    c2UySE (:,[1 end]) = 0; c2UySE(end, :) = 0;
    c2UyNW (:,[1 end]) = 0; c2UyNW(1  , :) = 0;
    c2UyNE (:,[1 end]) = 0; c2UyNE(end, :) = 0;
    c2UxSW (:,[1 end]) = 0;
    c2UxSE (:,[1 end]) = 0;
    c2UxNW (:,[1 end]) = 0;
    c2UxNE (:,[1 end]) = 0;
    c2UxSSW(:,[1 end]) = 0; c2UxSSW(:,[1 2]   ) = 0;
    c2UxSSE(:,[1 end]) = 0; c2UxSSE(:,[1 2]   ) = 0;
    c2UxNNW(:,[1 end]) = 0; c2UxNNW(:,[end end-1]) = 0;
    c2UxNNE(:,[1 end]) = 0; c2UxNNE(:,[end end-1]) = 0;
    c2UxSWW(:,[1 end]) = 0; c2UxSWW(1  ,:)   = 0;
    c2UxNWW(:,[1 end]) = 0; c2UxNWW(1  ,:)   = 0;
    c2UxSEE(:,[1 end]) = 0; c2UxSEE(end,:)   = 0;
    c2UxNEE(:,[1 end]) = 0; c2UxNEE(end,:)   = 0;
else
    sc2Uy = max(c2UyC(:));
    c2UyC  (:,[1    ]) = sc2Uy;
    c2UyW  (:,[1    ]) = 0; c2UyW  (1,  :) = 0;
    c2UyE  (:,[1    ]) = 0; c2UyE  (end,:) = 0;
    c2UyS  (:,[1    ]) = 0;
    c2UyN  (:,[1 end]) = 0;
    c2UySW (:,[1    ]) = 0; c2UySW(1  , :) = 0;
    c2UySE (:,[1    ]) = 0; c2UySE(end, :) = 0;
    c2UyNW (:,[1 end]) = 0; c2UyNW(1  , :) = 0;
    c2UyNE (:,[1 end]) = 0; c2UyNE(end, :) = 0;
    c2UxSW (:,[1    ]) = 0;
    c2UxSE (:,[1    ]) = 0;
    c2UxNW (:,[1 end]) = 0;
    c2UxNE (:,[1 end]) = 0;
    c2UxSSW(:,[1    ]) = 0; c2UxSSW(:,[1 2]   ) = 0;
    c2UxSSE(:,[1    ]) = 0; c2UxSSE(:,[1 2]   ) = 0;
    c2UxNNW(:,[1 end]) = 0; c2UxNNW(:,[end]) = 0;
    c2UxNNE(:,[1 end]) = 0; c2UxNNE(:,[end ]) = 0;
    c2UxSWW(:,[1    ]) = 0; c2UxSWW(1  ,:)   = 0;
    c2UxNWW(:,[1 end]) = 0; c2UxNWW(1  ,:)   = 0;
    c2UxSEE(:,[1    ]) = 0; c2UxSEE(end,:)   = 0;
    c2UxNEE(:,[1 end]) = 0; c2UxNEE(end,:)   = 0;
end
% Symmetry
if symmetry == 1
    c2UyS( : , 2   ) = 0;
    c2UyN( : ,end-1) = 0;
    c2UxSW (1  ,  :) = 0;
    c2UxSE (end,  :) = 0;
    c2UxNW (1  ,  :) = 0;
    c2UxNE (end,  :) = 0;
end
% Equation indexes
iUy     = NumUyG(:);
iUyW    = ones(size(Uy));   iUyW(2:end  ,2:end-1) = NumUyG(1:end-1,2:end-1);
iUyE    = ones(size(Uy));   iUyE(1:end-1,2:end-1) = NumUyG(2:end  ,2:end-1);
iUyS    = ones(size(Uy));   iUyS( :     ,2:end-1) = NumUyG( :     ,1:end-2);
iUyN    = ones(size(Uy));   iUyN( :     ,2:end-1) = NumUyG( :     ,3:end  );
iUySW   = ones(size(Uy));  iUySW(2:end  ,2:end-1) = NumUyG(1:end-1,1:end-2);
iUySE   = ones(size(Uy));  iUySE(1:end-1,2:end-1) = NumUyG(2:end  ,1:end-2);
iUyNW   = ones(size(Uy));  iUyNW(2:end  ,2:end-1) = NumUyG(1:end-1,3:end-0);
iUyNE   = ones(size(Uy));  iUyNE(1:end-1,2:end-1) = NumUyG(2:end  ,3:end-0);
iUxSW   = ones(size(Uy));  iUxSW( :     ,2:end-1) = NumUx (1:end-1,1:end-1);
iUxNW   = ones(size(Uy));  iUxNW( :     ,2:end-1) = NumUx (1:end-1,2:end-0);
iUxSE   = ones(size(Uy));  iUxSE( :     ,2:end-1) = NumUx (2:end  ,1:end-1);
iUxNE   = ones(size(Uy));  iUxNE( :     ,2:end-1) = NumUx (2:end  ,2:end-0);
iUxSSW  = ones(size(Uy)); iUxSSW( :     ,3:end-1) = NumUx (1:end-1,1:end-2);
iUxSSE  = ones(size(Uy)); iUxSSE( :     ,3:end-1) = NumUx (2:end  ,1:end-2);
iUxNNW  = ones(size(Uy)); iUxNNW( :     ,2:end-2) = NumUx (1:end-1,3:end  );
iUxNNE  = ones(size(Uy)); iUxNNE( :     ,2:end-2) = NumUx (2:end  ,3:end  );
iUxSWW  = ones(size(Uy)); iUxSWW(2:end  ,2:end-1) = NumUx (1:end-2,1:end-1);
iUxNWW  = ones(size(Uy)); iUxNWW(2:end  ,2:end-1) = NumUx (1:end-2,2:end-0);
iUxSEE  = ones(size(Uy)); iUxSEE(1:end-1,2:end-1) = NumUx (3:end  ,1:end-1);
iUxNEE  = ones(size(Uy)); iUxNEE(1:end-1,2:end-1) = NumUx (3:end  ,2:end-0);
%% Assembly of the sparse matrix
% Triplets
IvvJ_2   = [   iUy(:);   iUy(:);   iUy(:);   iUy(:);   iUy(:);    iUy(:);    iUy(:);    iUy(:);    iUy(:);    iUy(:);    iUy(:);    iUy(:);    iUy(:);     iUy(:);     iUy(:);     iUy(:);     iUy(:);     iUy(:);     iUy(:);     iUy(:);     iUy(:) ]';
JvvJ_2    = [   iUy(:);  iUyW(:);  iUyS(:);  iUyE(:);  iUyN(:);  iUySW(:);  iUySE(:);  iUyNW(:);  iUyNE(:);  iUxSW(:);  iUxSE(:);  iUxNW(:);  iUxNE(:);  iUxSSW(:);  iUxSSE(:);  iUxNNW(:);  iUxNNE(:);  iUxSWW(:);  iUxNWW(:);  iUxSEE(:);  iUxNEE(:) ]';
VvvJ_2    = [ c2UyC(:); c2UyW(:); c2UyS(:); c2UyE(:); c2UyN(:); c2UySW(:); c2UySE(:); c2UyNW(:); c2UyNE(:); c2UxSW(:); c2UxSE(:); c2UxNW(:); c2UxNE(:); c2UxSSW(:); c2UxSSE(:); c2UxNNW(:); c2UxNNE(:); c2UxSWW(:); c2UxNWW(:); c2UxSEE(:); c2UxNEE(:) ]';
if SuiteSparse==1, K  = sparse2([IuuJ_2(:); IvvJ_2(:);], [JuuJ_2(:); JvvJ_2(:);], [VuuJ_2(:); VvvJ_2(:);], (nx+1)*ny+(ny+1)*nx, (nx+1)*ny+(ny+1)*nx );
else               K  = sparse ([IuuJ_2(:); IvvJ_2(:);], [JuuJ_2(:); JvvJ_2(:);], [VuuJ_2(:); VvvJ_2(:);], (nx+1)*ny+(ny+1)*nx, (nx+1)*ny+(ny+1)*nx );
end
% RHS
% Indices of conforming BC nodes
ibcUxW = NumUx ( 1     ,:)';      ibcUxE = NumUx (end, : )';
ibcUyS = NumUyG( :     ,1) ;      ibcUyN = NumUyG( : ,end) ; % global Numbering for Uy !
if free_surf == 0
    ibc      = [ ibcUxW; ibcUxE; ibcUyS; ibcUyN ];
    BcK(ibc) = [BC.Ux_W*sc1Ux; BC.Ux_E*sc1Ux; BC.Uy_S*sc2Uy; BC.Uy_N*sc2Uy];
else
    ibc    = [ ibcUxW; ibcUxE; ibcUyS ];
    BcK(ibc) = [BC.Ux_W*sc1Ux; BC.Ux_E*sc1Ux; BC.Uy_S*sc2Uy];
end
end

