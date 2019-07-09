function [ grad, div, BcD ] = M2Di2_AssembleDivGrad( BC, dx, dy, nx, ny,NumVx, NumVyG, NumPt, SuiteSparse )
iVxC   = NumVx;                                                              % dP/dx
iPtW   =  ones(size(iVxC));    iPtW(1:end-1,:) = NumPt;
iPtE   =  ones(size(iVxC));    iPtE(2:end  ,:) = NumPt;
cPtW   =  ones(size(iVxC))/dx; cPtW([1 end],:) = 0;
cPtE   = -ones(size(iVxC))/dx; cPtE([1 end],:) = 0;
Idx    = [ iVxC(:); iVxC(:) ]';
Jdx    = [ iPtW(:); iPtE(:) ]';
Vdx    = [ cPtW(:); cPtE(:) ]';
iVyC   = NumVyG;                                                             % dP/dy
iPtS   =  ones(size(iVyC));    iPtS(:,1:end-1) = NumPt;
iPtN   =  ones(size(iVyC));    iPtN(:,2:end  ) = NumPt;
cPtS   =  ones(size(iVyC))/dy; cPtS(:,[1 end]) = 0;
cPtN   = -ones(size(iVyC))/dy; cPtN(:,[1 end]) = 0;
Idy    = [ iVyC(:); iVyC(:) ]';
Jdy    = [ iPtS(:); iPtN(:) ]';
Vdy    = [ cPtS(:); cPtN(:) ]';
%% Assemble grad and divV
if SuiteSparse==1, grad = sparse2( [Idx(:); Idy(:)], [Jdx(:); Jdy(:)], [Vdx(:); Vdy(:)], (nx+1)*ny+(ny+1)*nx, nx*ny );
else               grad = sparse ( [Idx(:); Idy(:)], [Jdx(:); Jdy(:)], [Vdx(:); Vdy(:)], (nx+1)*ny+(ny+1)*nx, nx*ny ); end
div    = -grad';
%% Build BC's for divergence
BcD = zeros(size(div,1),1);
BcD(NumPt( 1 , : )) = BcD(NumPt(1  ,:  )) + 1/dx*BC.Ux_W;
BcD(NumPt(end, : )) = BcD(NumPt(end,:  )) - 1/dx*BC.Ux_E;
BcD(NumPt(:  , 1 )) = BcD(NumPt(:  ,1  )) + 1/dy*BC.Uy_S;
BcD(NumPt(:  ,end)) = BcD(NumPt(:  ,end)) - 1/dy*BC.Uy_N;
end