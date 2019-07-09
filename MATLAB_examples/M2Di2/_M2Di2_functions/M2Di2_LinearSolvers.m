function [du,dp,dT,fu1,fp1,ft1, failed,restart ] = M2Di2_LinearSolvers(lsolver,Newton,tol_linu,tol_linp,eps_kspgcr,nPH,noisy,SuiteSparse, BcK,BcD,BcT, KM,grad,div,PP,PPI,T2M,TM, JM, M2TJ,T2MJ,TMJ, fu,fp,ft, failed, NumVyG,NumPtG,NumTeG)
%% Linear solver - obtain velocity/pressure/temperature corrections for current non-linear iteration
restart = 0;
if lsolver == -2 % Decoupled thermo-mechanical solver V-P-T ------------------------------
    [du,dp,dT,fu1,fp1,ft1, failed,restart ] = M2Di2_DecoupledThermoMechanicalSolver(tol_linu,tol_linp,eps_kspgcr,nPH,noisy,SuiteSparse, BcK,BcD,BcT, KM,grad,div,PP,T2M,TM, JM, M2TJ,T2MJ,TMJ, fu,fp,ft, failed);
elseif lsolver == -1 % Backslash coupled solver V-P-T ------------------------------
    if Newton == 1, % Assemble entire Jacobian
        Ms  = [   JM,  grad,   M2TJ; ...
                 div,  0*PP,   0*PP;
                T2MJ,  0*PP, TM+TMJ;];
    else % Assemble Picard operator
        Ms  = [   KM , grad, 0*grad; ...
                  div, 0*PP,   0*PP;
                0*T2M, 0*PP,     TM;];
    end
    f   = [ fu; fp; ft ];                                                % Assemble entire residual vector
    dX  = Ms\f;                                                          % Call direct solver
    du  = dX(1:max(NumVyG(:)));                                          % Extract velocity correction
    dp  = dX(NumPtG(:));                                                 % Extract pressure correction
    dT  = dX(NumTeG(:));
    f1  = Ms*dX - f;                                                     % Compute entire linear residual
    fu1 = f1(1:max(NumVyG(:)));                                          % Extract velocity linear residual
    fp1 = f1(NumPtG(:));                                                 % Extract pressure linear residual
    ft1 = f1(NumTeG(:));
elseif lsolver == 0 % Backslash coupled solver V-P ------------------------------
    Ms  = [  JM, grad ; ...
            div, 0*PP ];                                                 % Assemble entire Jacobian
    f   = [ fu  ; fp   ];                                                % Assemble entire residual vector
    dX  = Ms\f;                                                          % Call direct solver
    du  = dX(1:max(NumVyG(:)));                                          % Extract velocity correction
    dp  = dX(NumPtG(:));                                                 % Extract pressure correction
    dT  = 0*dp;
    f1  = Ms*dX - f;                                                     % Compute entire linear residual
    fu1 = f1(1:max(NumVyG(:)));                                          % Extract velocity linear residual
    fp1 = f1(NumPtG(:));                                                 % Extract pressure linear residual
    ft1 = 0*fp1;
elseif lsolver == 1 % Powell-Hestenes solver V-P ------------------------------    
    if Newton==0, JM = KM; end
    Kt  = KM - grad*(PPI*div);                                            % Velocity Schur complement of Picard operator (Kt)
    Jt  = JM - grad*(PPI*div);                                            % Velocity Schur complement of Jacobian operator
    [Kc,e,s] = chol(Kt,'lower','vector');                                 % Choleski factorization of Kt
    % Powell-Hestenes iterations
    fu0 = fu;                                                                    % Save linear norm 0
    dp  = zeros(size(BcD));
    du  = zeros(size(BcK));
    for itPH=1:nPH
        fut  = fu - grad*dp - grad*PPI*fp;                                       % Iterative right hand side
        [du,norm_r,its] = M2Di2_kspgcr_m(Jt,fut,du,Kc,s,eps_kspgcr,noisy,SuiteSparse); % Apply inverse of Schur complement
        dp   = dp + PPI*(fp -  div*du);                                          % Pressure corrctions
        fu1  = fu - JM*du  - grad*dp;                                            % Compute linear velocity residual
        fp1  = fp - div*du;                                                      % Compute linear pressure residual
        if noisy>1, fprintf('--- iteration %d --- \n',itPH);
            fprintf('  Res. |u| = %2.2e \n',norm(fu1)/length(fu1));
            fprintf('  Res. |p| = %2.2e \n',norm(fp1)/length(fp1));
            fprintf('  KSP GCR: its=%1.4d, Resid=%1.4e \n',its,norm_r); end
        if ((norm(fu1)/length(fu1)) < tol_linu) && ((norm(fp1)/length(fp1)) < tol_linp), break; end
        if ((norm(fu1)/length(fu1)) > (norm(fu0)/length(fu1)) && norm(fu1)/length(fu1) < tol_linu*5e3), fprintf(' > Linear residuals do no converge further:\n'); break; end
        fu0 = fu1;
    end
    dT  = 0*dp;
    ft1 = 0*fp1;
end
end