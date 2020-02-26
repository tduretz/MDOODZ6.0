clear

step = 5;

% Data from M2Di
path = '/Users/tduretz/Dropbox/Thibault_Sam_Ldx/M2Di_v2/M2Di_v2/progress2/';
Mc = load([path 'MatlabStokes_step' num2str(step,'%02d'), '_iter00']);

% Data from MDOODZ:
path = '/Users/tduretz/REPO_GIT/MDOODZ6.0_PreHolidays/SOURCE/';
filename = [path 'Stokes_01cpu_step', num2str(step,'%02d'), '_iter00.gzip.h5'];
[ MA, MB, MC, MD, rSt, rPt, fSt, fPt ] = ExtractDecoupledMatrixFromMDoodz( filename );

% Matrices
M1 = [Mc.K Mc.grad; Mc.div 0*Mc.PP ];
M2 = [MA MB; MC 0*MD ] ;

% Multiply M2Di matrix entries by cell volume
M1 = M1 * (Mc.dx*Mc.dy);

% Kill dirichlet row/cols from the M2Di code
west  = Mc.NumVx(1,:);
east  = Mc.NumVx(end,:);
south = Mc.NumVyG(:,1);
north = Mc.NumVyG(:,end);
all = [west(:); east(:); south(:); north(:)];
M1(all,:) = [];
M1(:,all) = [];

% RHSs
% R1 = Mc.L;
% R2 = [Md.rSt; Md.rPt];

diff = abs(M1-M2);
diff(diff<1e-4) = 0;

figure(1)
subplot(131), spy(M1)
subplot(132), spy(M2)
subplot(133), spy(diff)

fprintf('Matlab')
M1(1,:)
fprintf('M2Di2')
M2(1,:)

fprintf('Matlab')
M1(699,:)
fprintf('M2Di2')
M2(699,:)

fprintf('Matlab')
M1(2600,:)
fprintf('M2Di2')
M2(2600,:)

% norm(R1-R2)

% Mc = load('CoupledMatrix_BENALLAL');
% Md = load('CoupledMatrix_BENALLAL2');
% 
% % Matrices
% M1 = Mc.M;
% M2 = Md.M;
% 
% % RHSs
% R1 = Mc.L;
% R2 = Md.L;
% 
% figure(1)
% subplot(131), spy(M1)
% subplot(132), spy(M2)
% 
% 
% Mdiff = M1-M2;
% Mdiff = 0;
% subplot(133), spy(Mdiff)
% 
% norm(R1-R2)