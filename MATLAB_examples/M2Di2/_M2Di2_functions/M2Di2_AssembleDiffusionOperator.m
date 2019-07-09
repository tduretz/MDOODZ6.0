function [TM, BcTe] = M2Di2_AssembleDiffusionOperator( BC, dt, dx, dy, ck, nx, ny, kt, NumTe, Tc, Toc, N1, N2, SuiteSparse)
d_dt  = ones(size(Tc))*N1/dt;
k_W   = zeros(size(Tc)); k_W(2:end  ,:) = (kt(1:end-1,:)   + kt(2:end,:))/2;
k_E   = zeros(size(Tc)); k_E(1:end-1,:) = (kt(1:end-1,:)   + kt(2:end,:))/2;
k_S   = zeros(size(Tc)); k_S(:,2:end  ) = (kt(:,1:end-1)   + kt(:,2:end))/2;
k_N   = zeros(size(Tc)); k_N(:,1:end-1) = (kt(:,1:end-1)   + kt(:,2:end))/2;
iTe   = NumTe;
iTeW  = ones(size(Tc)); iTeW(2:end  ,:) = NumTe(1:end-1, :     );
iTeE  = ones(size(Tc)); iTeE(1:end-1,:) = NumTe(2:end  , :     );
iTeS  = ones(size(Tc)); iTeS(:,2:end  ) = NumTe( :     ,1:end-1);
iTeN  = ones(size(Tc)); iTeN(:,1:end-1) = NumTe( :     ,2:end  );
cTeC  = d_dt + (1-ck)*(k_W+k_E)/dx/dx + (1-ck)*(k_S+k_N)/dy/dy;
cTeW  = -(1-ck)*k_W/dx/dx;
cTeS  = -(1-ck)*k_S/dy/dy;
cTeE  = -(1-ck)*k_E/dx/dx;
cTeN  = -(1-ck)*k_N/dy/dy;
I     = [  iTe(:);  iTe(:);  iTe(:);   ]';
J     = [  iTe(:);  iTeW(:); iTeS(:);  ]';
V     = [ cTeC(:);  cTeW(:); cTeS(:);  ]';
if SuiteSparse == 1, TM    = sparse2(I,     J,     V,        nx*ny,nx*ny);
else                 TM    = sparse (I,     J,     V,        nx*ny,nx*ny);
end
TM    = TM  + TM'  - diag(diag(TM));
% BC's on TM
BcTe                 = zeros(size(TM,1),1);
BcTe(NumTe(:))       = d_dt(:).*Toc(:);
% Crank-Nicolson stuff
% BcTe(NumTe(:))       = BcTe(NumTe(:)) + ck*diffo(:);
% BcTe(NumTe(:))       = BcTe(NumTe(:)) + ck*N2*ShearHeat*Hso(:);
% % RHS form of shear heating --- unsused
% Hs1            = 0*(4*etac.*Eiic2);
% BcTe(NumTe(:)) = BcTe(NumTe(:)) + (1-ck)*N2*Hs1(:);
end