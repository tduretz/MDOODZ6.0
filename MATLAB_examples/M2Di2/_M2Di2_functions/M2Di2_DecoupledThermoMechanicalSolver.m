function [x,y,z,rx,ry,rz, failed,restart,iter ] = DecoupledThermoMechanicalSolver(tol_linu,tol_linp,eps_kspgcr,nPH,noisy,SuiteSparse, BcK,BcD,BcT, KM,grad,div,PP,T2M,TM, JM, M2TJ,T2MJ,TMJ, fu,fp,ft, failed)
% Newton correction arrays
x  = zeros(size(BcK));  y = zeros(size(BcD));  z = zeros(size(BcT)); % prepare decoupled solve
% Linear solver correction array
xx = zeros(size(BcK)); yy = zeros(size(BcD)); zz = zeros(size(BcT));
% Picard operator
A=KM;     B=grad;  C=0*PP; rx=BcK;
D=div;    E=PP;    F=0*PP; ry=BcD;
G=T2M;    H=0*PP;  I=TM;   rz=BcT;
% Newton operator
AJ=JM;     BJ=grad;  CJ=M2TJ;
DJ=div;    EJ=PP;    FJ=0*PP;
GJ=T2MJ;   HJ=0*PP;  IJ=TM+TMJ;
Fx = fu;%-( A*Sx + B*Sy - Rx );        resnlX = norm(Fx)/length(Fx) % nonlinear residuals
Fy = fp;%-( D*Sx + E*Sy - Ry );        resnlY = norm(Fy)/length(Fy)
Fz = ft;%-( G*Sx + I*Sz - Rz );        resnlZ = norm(Fz)/length(Fz)
tic
% Pre-conditionning (~Jacobi)
EJi           = spdiags(1./diag(EJ),0,size(EJ,1),size(EJ,2));  % trivial inverse
Asc           = AJ-BJ*(EJi*DJ);
AscP          =  A- B*(EJi*D );
[Icc  ,pI,sI] = chol(IJ  ,'lower','vector');                     % Cholesky factors
[AsccP,pA,sA] = chol(AscP,'lower','vector');                     % Cholesky factors
resx=1; resy=resx; resz=resx; iR=0;
while (resx > tol_linu || resy > tol_linp || resz > tol_linu) && iR < nPH
    iR = iR+1;
    rx   = -( AJ*x +   BJ*y + CJ*z - Fx ); resx = norm(rx)/length(rx);
    ry   = -( DJ*x + 0*EJ*y + FJ*z - Fy ); resy = norm(ry)/length(ry);
    rz   = -( GJ*x +   HJ*y + IJ*z - Fz ); resz = norm(rz)/length(rz);
    rxsc = rx - BJ*(EJi*ry) - CJ*zz;
    [xx,norm_r,its] = M2Di2_kspgcr_m(Asc,rxsc,xx,AsccP,sA,eps_kspgcr,noisy,SuiteSparse);
    yy       = EJi*(ry-DJ*xx);
    rzsc     = rz - GJ*xx;
    zz(sI,1) = cs_ltsolve(Icc, cs_lsolve(Icc, rzsc(sI)));
    x   = x + xx; % Updates
    y   = y + yy;
    z   = z + zz;
end%iR
if iR < nPH %success of linear solve
    fprintf('\n >> Linear solve [OK] in %02d iterations - %3d s\n', iR, toc);
    restart=0;
else %fail of linear solve
    fprintf(' -------------------------------------------\n')
    fprintf(' >> Linear solver [FAIL] after %02d iterations\n', iR);
    failed=failed+1; restart=1; %nfail=nfail+1;
    fprintf(' -------------------------------------------\n')
end
% fprintf(' Linear:   Momentum resid = %2.6e\n'  , norm(rx)/length(Fx) )
% fprintf('         Continuity resid = %2.6e\n'  , norm(ry)/length(Fy) )
% fprintf('            Thermal resid = %2.6e\n\n', norm(rz)/length(Fz) )
end