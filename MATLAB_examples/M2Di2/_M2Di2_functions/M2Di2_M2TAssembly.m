function M2T = M2TAssembly( BC, D, dx, dy, nx, ny, ck, N2, NumVx, NumVyG, NumTe, Vx, Vy, SuiteSparse )
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
%% ------------------------ M2T x ------------------------ %
% Stencil weights for normal strains interpolation
wS_SW = wSW(:,1:end-1); wS_SE = wSE(:,1:end-1); wS_W  = wNW(:,1:end-1); wS_E  = wNE(:,1:end-1);
wN_W  = wSW(:,2:end-0); wN_E  = wSE(:,2:end-0); wN_NW = wNW(:,2:end-0); wN_NE = wNE(:,2:end-0);
% Elastic or elasto-plastic operators West-East
D14W    = zeros(size(Vx));  D14W( 2:end-1, : ) = D.D14c(1:end-1,:);
D14E    = zeros(size(Vx));  D14E( 2:end-1, : ) = D.D14c(2:end,:);
% Elastic or elasto-plastic operators South-North
D34S    = zeros(size(Vx));  D34S( 2:end-1 , :) = D.D34v(2:end-1, 1:end-1);  %D31S(:,1) = 0;
D34N    = zeros(size(Vx));  D34N( 2:end-1 , :) = D.D34v(2:end-1, 2:end  );  %D31N(:,end) = 0;
% Coefficicients
c1TSW = ((1 - fsxS) .* D34S .* wS_SW ./ dy) ./ 0.4e1;
c1TSE = ((1 - fsxS) .* D34S .* wS_SE ./ dy) ./ 0.4e1;
c1TW = D14W ./ dx - (((1 - fsxN) .* D34N .* wN_W) ./ 0.4e1 - ((1 - fsxS) .* D34S .* wS_W) ./ 0.4e1) ./ dy;
c1TE = -D14E ./ dx - (((1 - fsxN) .* D34N .* wN_E) ./ 0.4e1 - ((1 - fsxS) .* D34S .* wS_E) ./ 0.4e1) ./ dy;
c1TNW = -((1 - fsxN) .* D34N .* wN_NW ./ dy) ./ 0.4e1;
c1TNE = -((1 - fsxN) .* D34N .* wN_NE ./ dy) ./ 0.4e1;
c1TSW ([1 end],:) = 0;
c1TSE ([1 end],:) = 0;
c1TW  ([1 end],:) = 0; %c1UxS ( :   , 1 ) = 0;
c1TE  ([1 end],:) = 0; %c1UxN ( :   ,end) = 0;
c1TNW ([1 end],:) = 0; %c1UxSW( :   , 1 ) = 0;
c1TNE ([1 end],:) = 0; %c1UxSE( :   , 1 ) = 0;
iVx   = NumVx(:);
iTW   = ones(size(Vx)); iTW(2:end  ,:) = NumTe;
iTE   = ones(size(Vx)); iTE(1:end-1,:) = NumTe;
iTSW  = ones(size(Vx)); iTSW(2:end  , 2:end  ) = NumTe(:,1:end-1);
iTSE  = ones(size(Vx)); iTSE(1:end-1, 2:end  ) = NumTe(:,1:end-1);
iTNW  = ones(size(Vx)); iTNW(2:end  , 1:end-1) = NumTe(:,2:end);
iTNE  = ones(size(Vx)); iTNE(1:end-1, 1:end-1) = NumTe(:,2:end);
IuTJ  = [   iVx(:);   iVx(:);   iVx(:);   iVx(:);   iVx(:);    iVx(:) ];
JuTJ  = [   iTSW(:);  iTSE(:);  iTW(:);   iTE(:);   iTNW(:);   iTNE(:)];
VuTJ  = [  c1TSW(:); c1TSE(:); c1TW(:);  c1TE(:);  c1TNW(:);  c1TNE(:)];
%% ------------------------ M2T y ------------------------ %
% Stencil weights for normal strains interpolation
wW_SW = wSW(1:end-1,:); wW_NW = wNW(1:end-1,:); wW_S  = wSE(1:end-1,:); wW_N  = wNW(1:end-1,:);
wE_S  = wSW(2:end-0,:); wE_N  = wSW(2:end-0,:); wE_SE = wSE(2:end-0,:); wE_NE = wNE(2:end-0,:);
% Elastic or elasto-plastic operators South-North
D24S    = zeros(size(Vy));  D24S( : ,2:end-1) = D.D24c(:,1:end-1);% D21S(:,1) = 0;
D24N    = zeros(size(Vy));  D24N( : ,2:end-1) = D.D24c(:,2:end-0); %D21N(:,end) = 0;
% Elastic or elasto-plastic operators West-East
D34W    = zeros(size(Vy));  D34W( : ,2:end-1 ) = D.D34v(1:end-1, 2:end-1); %D31W(1,:) = 0;
D34E    = zeros(size(Vy));  D34E( : ,2:end-1 ) = D.D34v(2:end  , 2:end-1); %D31E(end,:) = 0;
% Coefficicients
c2TSW = ((1 - fsyW) .* D34W .* wW_SW ./ dx) ./ 0.4e1;
c2TS = D24S ./ dy - (((1 - fsyE) .* D34E .* wE_S) ./ 0.4e1 - ((1 - fsyW) .* D34W .* wW_S) ./ 0.4e1) ./ dx;
c2TSE = -((1 - fsyE) .* D34E .* wE_SE ./ dx) ./ 0.4e1;
c2TNW = ((1 - fsyW) .* D34W .* wW_NW ./ dx) ./ 0.4e1;
c2TN = -D24N ./ dy - (((1 - fsyE) .* D34E .* wE_N) ./ 0.4e1 - ((1 - fsyW) .* D34W .* wW_N) ./ 0.4e1) ./ dx;
c2TNE = -((1 - fsyE) .* D34E .* wE_NE ./ dx) ./ 0.4e1;
c2TSW (:,[1 end]) = 0;
c2TSE (:,[1 end]) = 0;
c2TS  (:,[1 end]) = 0; %c1UxS ( :   , 1 ) = 0;
c2TN  (:,[1 end]) = 0; %c1UxN ( :   ,end) = 0;
c2TNW (:,[1 end]) = 0; %c1UxSW( :   , 1 ) = 0;
c2TNE (:,[1 end]) = 0; %c1UxSE( :   , 1 ) = 0;
% Equation indexes
iVy   = NumVyG(:);
iTSW  = ones(size(Vy)); iTSW(2:end  ,2:end  ) = NumTe(1:end-1,: );
iTSE  = ones(size(Vy)); iTSE(1:end-1,2:end  ) = NumTe(2:end,: );
iTS   = ones(size(Vy)); iTS (:,2:end  ) = NumTe;
iTN   = ones(size(Vy)); iTN (:,1:end-1) = NumTe;
iTNW  = ones(size(Vy)); iTNW(2:end  ,1:end-1) = NumTe(1:end-1,: );
iTNE  = ones(size(Vy)); iTNE(1:end-1,1:end-1) = NumTe(2:end,: );
%% Assembly of the sparse matrix
IvTJ  =  [   iVy(:);   iVy(:);   iVy(:);   iVy(:);   iVy(:);     iVy(:)];
JvTJ  =  [  iTSW(:);  iTSE(:);   iTS(:);   iTN(:);  iTNW(:);    iTNE(:)];
VvTJ  =  [ c2TSW(:); c2TSE(:);  c2TS(:);  c2TN(:); c2TNW(:);   c2TNE(:)];
if SuiteSparse == 1, M2T   = sparse2( [IuTJ(:); IvTJ(:)], [JuTJ(:); JvTJ(:)], [VuTJ(:); VvTJ(:)]);
else                 M2T   = sparse ( [IuTJ(:); IvTJ(:)], [JuTJ(:); JvTJ(:)], [VuTJ(:); VvTJ(:)]);
end
end