function T2M = T2MAssembly( BC, D, dx, dy, nx, ny, ck, N2, NumVx, NumVyG, NumTe, Tc, SuiteSparse )
D41c = D.D41c;
D42c = D.D42c;
D43c = D.D43c;
% BC flags
fs = zeros(nx+1,ny+1);
fs(:,1) = BC.fsxS(:,1);  fs(:,end) = BC.fsxN(:,end);
fs(1,:) = BC.fsyW(1,:);  fs(end,:) = BC.fsyE(end,:);
ns = zeros(nx+1,ny+1);
ns(:,1) = BC.nsxS(:,1);  ns(:,end) = BC.nsxN(:,end);
ns(1,:) = BC.nsyW(1,:);  ns(end,:) = BC.nsyE(end,:);
fsSW = fs(1:end-1,1:end-1);
fsSE = fs(2:end-0,1:end-1);
fsNW = fs(1:end-1,2:end-0);
fsNE = fs(2:end-0,2:end-0);
nsSW = ns(1:end-1,1:end-1);
nsSE = ns(2:end-0,1:end-1);
nsNW = ns(1:end-1,2:end-0);
nsNE = ns(2:end-0,2:end-0);
% Coefficients
c3VxW  = -N2 .* (1 - ck) .* (-0.2e1 ./ 0.3e1 ./ dx .* D41c + 0.1e1 ./ dx .* D42c ./ 0.3e1 + (((1 - fsNW) .* (-1 - nsNW) ./ dy) ./ 0.4e1 + ((1 - fsSW) .* (1 + nsSW) ./ dy) ./ 0.4e1) .* D43c);
c3VxE  = -N2 .* (1 - ck) .* (0.2e1 ./ 0.3e1 ./ dx .* D41c - 0.1e1 ./ dx .* D42c ./ 0.3e1 + (((1 - fsNE) .* (-1 - nsNE) ./ dy) ./ 0.4e1 + ((1 - fsSE) .* (1 + nsSE) ./ dy) ./ 0.4e1) .* D43c);
c3VyS  = -N2 .* (1 - ck) .* (0.1e1 ./ dy .* D41c ./ 0.3e1 - 0.2e1 ./ 0.3e1 ./ dy .* D42c + (((1 - fsSW) .* (1 + nsSW) ./ dx) ./ 0.4e1 + ((1 - fsSE) .* (-1 - nsSE) ./ dx) ./ 0.4e1) .* D43c);
c3VyN  = -N2 .* (1 - ck) .* (-0.1e1 ./ dy .* D41c ./ 0.3e1 + 0.2e1 ./ 0.3e1 ./ dy .* D42c + (((1 - fsNW) .* (1 + nsNW) ./ dx) ./ 0.4e1 + ((1 - fsNE) .* (-1 - nsNE) ./ dx) ./ 0.4e1) .* D43c);
c3VxSW = -(N2 .* (1 - ck) .* (1 - fsSW) .* (-1 + nsSW) ./ dy .* D43c) ./ 0.4e1;
c3VxSE = -(N2 .* (1 - ck) .* (1 - fsSE) .* (-1 + nsSE) ./ dy .* D43c) ./ 0.4e1;
c3VxNW = -(N2 .* (1 - ck) .* (1 - fsNW) .* (1 - nsNW) ./ dy .* D43c) ./ 0.4e1;
c3VxNE = -(N2 .* (1 - ck) .* (1 - fsNE) .* (1 - nsNE) ./ dy .* D43c) ./ 0.4e1;
c3VySW = -(N2 .* (1 - ck) .* (1 - fsSW) .* (-1 + nsSW) ./ dx .* D43c) ./ 0.4e1;
c3VySE = -(N2 .* (1 - ck) .* (1 - fsSE) .* (1 - nsSE) ./ dx .* D43c) ./ 0.4e1;
c3VyNW = -(N2 .* (1 - ck) .* (1 - fsNW) .* (-1 + nsNW) ./ dx .* D43c) ./ 0.4e1;
c3VyNE = -(N2 .* (1 - ck) .* (1 - fsNE) .* (1 - nsNE) ./ dx .* D43c) ./ 0.4e1;
% Equation indexes
iTe   = NumTe(:);
iVxW  = NumVx(1:end-1,:); iVyS  = NumVyG(:,1:end-1);
iVxE  = NumVx(2:end  ,:); iVyN  = NumVyG(:,2:end,:);
iVxSW = ones(size(Tc));  iVxSW(:,2:end  ) = NumVx (1:end-1, 1:end-1);
iVxSE = ones(size(Tc));  iVxSE(:,2:end  ) = NumVx (2:end  , 1:end-1);
iVxNW = ones(size(Tc));  iVxNW(:,1:end-1) = NumVx (1:end-1, 2:end  );
iVxNE = ones(size(Tc));  iVxNE(:,1:end-1) = NumVx (2:end  , 2:end  );
iVySW = ones(size(Tc));  iVySW(2:end  ,:) = NumVyG(1:end-1, 1:end-1);
iVySE = ones(size(Tc));  iVySE(1:end-1,:) = NumVyG(2:end  , 1:end-1);
iVyNW = ones(size(Tc));  iVyNW(2:end  ,:) = NumVyG(1:end-1, 2:end  );
iVyNE = ones(size(Tc));  iVyNE(1:end-1,:) = NumVyG(2:end  , 2:end  );
%% Assembly of the sparse matrix
It2m   = [   iTe(:);  iTe(:);  iTe(:);  iTe(:);   iTe(:);   iTe(:);   iTe(:);   iTe(:);   iTe(:);   iTe(:);   iTe(:);   iTe(:)  ]';
Jt2m   = [  iVxW(:); iVxE(:); iVyS(:); iVyN(:); iVxSW(:); iVxSE(:); iVxNW(:); iVxNE(:); iVySW(:); iVySE(:); iVyNW(:); iVyNE(:)  ]';
Vt2m   = [ c3VxW(:); c3VxE(:); c3VyS(:); c3VyN(:); c3VxSW(:); c3VxSE(:); c3VxNW(:); c3VxNE(:); c3VySW(:); c3VySE(:); c3VyNW(:); c3VyNE(:)  ]';
if SuiteSparse == 1, T2M    = sparse2( It2m, Jt2m, Vt2m );
else                 T2M    = sparse ( It2m, Jt2m, Vt2m );
end
end