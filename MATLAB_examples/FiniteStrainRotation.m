function main
clear
close all

% Deformation gradient --- simple shear
gamma = 1;
Fxx   = 1;
Fxz   = gamma;
Fzx   = 0;
Fzz   = 1;

% Rotation matrix
theta = atan(gamma/2);
R     = [cos(theta) sin(theta); -sin(theta) cos(theta)];

% Director vector: initial orientation
Ndx = 0;  Ndz = 1;
N   = [Ndx; Ndz];
a0N = atan(Ndz/Ndx); % initial angle of vector

% Stress tensor: initial
sxx = 0;  sxz = 1;
szx = 1;  szz = 0;
S   = [sxx sxz; szx szz];
[v,e] = eig(S);
a0S = atan(v(2,2)/v(1,2)); % initial angle of vector

% Rotate vector
N   = R*N;
Ndx = N(1);
Ndz = N(2);

% Rotate stress vector: this is a passive rotation of the stress tensor
S   = R*S*R';
sxx = S(1,1);  sxz = S(1,2);
szx = S(2,1);  szz = S(2,2);

% Pre-process
F    = [Fxx Fxz; Fzx Fzz];
CG   = F'*F;
CGxx = CG(1,1); CGxz = CG(1,2);
CGzx = CG(2,1); CGzz = CG(2,2);

%% Strain tensor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Thibault's version

% Stretch tensor = sqrt(CG)
[ U0xx,U0xz,U0zx,U0zz]  = Sqrt_2x2( CGxx, CGxz, CGzx, CGzz );
[e1,e2,Uxx,Uxz,Uzx,Uzz] = Eigs_2x2( U0xx, U0xz, U0zx, U0zz ); % strecth
[wxx,wxz,wzx,wzz]       =  Inv_2x2( U0xx, U0xz, U0zx, U0zz );
[rxx,rxz,rzx,rzz]       =  AxB_2x2( Fxx, Fxz, Fzx, Fzz, wxx, wxz, wzx, wzz );
[Uxx,Uxz,Uzx,Uzz]       =  AxB_2x2( rxx, rxz, rzx, rzz, Uxx, Uxz, Uzx, Uzz );
n                       = sqrt(Uxx.^2 + Uzx.^2);
v11                     = Uxx./n * e1;
v12                     = Uzx./n * e1;
n                       = sqrt(Uxz.^2 + Uzz.^2);
v21                     = Uxz./n * e2;
v22                     = Uzz./n * e2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Richard's version
CG      = [CGxx(1), CGxz(1); CGzx(1), CGzz(1) ];
F       = [ Fxx(1),  Fxz(1);  Fzx(1),  Fzz(1) ];
U0       = sqrtm(CG);
[UV,UE] = eigs(U0);
R       = F*inv(U0); % compute rotation
Urot    = R*UV;      % rotate
rv1     = Urot(:,2); % no more minus / no more flipud
rv2     = Urot(:,1);
rv1     = rv1 / sqrt( rv1(1)^2 + rv1(2)^2 ) * e2; % Normalise and scale vectors of Spitz
rv2     = rv2 / sqrt( rv2(1)^2 + rv2(2)^2 ) * e1; % Normalise and scale vectors of Spitz

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Marcel
angle  = 0.5*( atan( (2*( Fxx.*Fzx + Fxz.*Fzz)) ./ (Fxx.^2 + Fxz.^2 - Fzx.^2 - Fzz.^2 ) ) );
xell   = [ e1/4.*cos(angle) -e1/4.*cos(angle) e2/4.*cos(pi/2+angle) -e2/4.*cos(pi/2+angle)];
zell   = [ e1/4.*sin(angle) -e1/4.*sin(angle) e2/4.*sin(pi/2+angle) -e2/4.*sin(pi/2+angle)];
x1     = [0 + xell(1,1) 0 + xell(1,2)];
z1     = [0 + zell(1,1) 0 + zell(1,2)];
x2     = [0 + xell(1,3) 0 + xell(1,4)];
z2     = [0 + zell(1,3) 0 + zell(1,4)];

%% Stress tensor
[L1,L2,s11,s12,s21,s22] = Eigs_2x2( sxx, sxz, szx, szz );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Differences in finite strain axis lengths
de         = e1 - e2;

% Check angles of finite strain vectors for simple shear
ang_ana1   = atan(2/gamma)/2            * 180/pi;
ang_mar    = angle(1,1)                 * 180/pi;
ang_ric    = atan(rv2(2)/rv2(1))        * 180/pi;
ang_tib    = atan(v12(1)/v11(1))        * 180/pi;

% Check angles for vector angle
ang_ana2   = (a0N  - atan(gamma/2)       ) * 180/pi;
ang_dir    = (       atan(Ndz(1)/Ndx(1)) ) * 180/pi;

% Check angles for sigma 1 vector angle
ang_ana3   = (a0S  - atan(gamma/2)      ) * 180/pi;
ang_s1     =-(      atan(s12(1)/s11(1)) ) * 180/pi;

fprintf('\nFinite strain - stretch:\n')
fprintf('e1 - e2:  %2.2f  error: %2.2e\n', de, abs(de - gamma))

fprintf('\nFinite strain - principal strain angle:\n')
fprintf('Marcel:   %2.2f  error: %2.2e\n', ang_mar, abs(ang_mar - ang_ana1))
fprintf('Richard:  %2.2f  error: %2.2e\n', ang_mar, abs(ang_ric - ang_ana1))
fprintf('Thibault: %2.2f  error: %2.2e\n', ang_mar, abs(ang_tib - ang_ana1))

fprintf('\nDirector vector  angle:\n')
fprintf('Angle:   %2.2f  error: %2.2e\n', ang_dir, abs(ang_dir - ang_ana2))

fprintf('\nStress - principal stress angle:\n')
fprintf('Angle:   %2.2f  error: %2.2e\n', ang_s1 , abs(ang_s1 - ang_ana3))
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1), hold on

% Plot Richard ellipse
quiver( 0, 0,  rv1(1), rv1(2), 0.25,'b','Autoscale','off',  'LineWidth', 10, 'ShowArrowHead','off' )
quiver( 0, 0,  rv2(1), rv2(2), 0.25,'b','Autoscale','off',  'LineWidth', 10, 'ShowArrowHead','off' )
quiver( 0, 0, -rv1(1), -rv1(2), 0.25,'b','Autoscale','off',  'LineWidth', 10, 'ShowArrowHead','off' )
quiver( 0, 0, -rv2(1), -rv2(2), 0.25,'b','Autoscale','off',  'LineWidth', 10, 'ShowArrowHead','off' )

% Plot Marcel ellipse
plot(x1,z1, '-r', 'LineWidth', 5 )
plot(x2,z2, '-r', 'LineWidth', 5)

% Plot Thibault ellipse
quiver( 0, 0,      v11,      v12,  0.25, 'g', 'Autoscale', 'on',  'LineWidth', 1, 'ShowArrowHead','off')
quiver( 0, 0,     -v11,     -v12,  0.25, 'g', 'Autoscale', 'on',  'LineWidth', 1, 'ShowArrowHead','off')
quiver( 0, 0,      v21,      v22,  0.25, 'g', 'Autoscale', 'on',  'LineWidth', 1, 'ShowArrowHead','off')
quiver( 0, 0,     -v21,     -v22,  0.25, 'g', 'Autoscale', 'on',  'LineWidth', 1, 'ShowArrowHead','off')

% plot vector
quiver( 0, 0,      Ndx,      Ndz,  0.25, 'k', 'Autoscale', 'on',  'LineWidth', 1, 'ShowArrowHead','off')
quiver( 0, 0,     -Ndx,     -Ndz,  0.25, 'k', 'Autoscale', 'on',  'LineWidth', 1, 'ShowArrowHead','off')

% plot stress ellipse
quiver( 0, 0,      s11,      s12,  0.25, 'm', 'Autoscale', 'on',  'LineWidth', 1, 'ShowArrowHead','off')
quiver( 0, 0,     -s11,     -s12,  0.25, 'm', 'Autoscale', 'on',  'LineWidth', 1, 'ShowArrowHead','off')
quiver( 0, 0,      s21,      s22,  0.25, 'm', 'Autoscale', 'on',  'LineWidth', 1, 'ShowArrowHead','off')
quiver( 0, 0,     -s21,     -s22,  0.25, 'm', 'Autoscale', 'on',  'LineWidth', 1, 'ShowArrowHead','off')

axis equal xy, xlabel('x [-]'), ylabel('z [-]')

end

function  [aixx,aizx,aixz,aizz] = Inv_2x2( axx, axz, azx, azz )
Det  = axx(:).*azz(:) - axz(:).*azx(:);
aixx = 1./Det .*azz;
aixz =-1./Det .*axz;
aizx =-1./Det .*azx;
aizz = 1./Det .*axx;
end

function [c11,c12,c21,c22] =  AxB_2x2( a11, a12, a21, a22, b11,b12,b21,b22 )
c11 = a11.*b11 + a12.*b21;
c12 = a11.*b12 + a12.*b22;
c21 = a21.*b11 + a22.*b21;
c22 = a21.*b12 + a22.*b22;
end

function [b11,b12,b21,b22] =  Sqrt_2x2( a11, a12, a21, a22 )
Tr  = a11(:) + a22(:);
Det = a11(:).*a22(:) - a12(:).*a21(:);
s   = sqrt(Det);
t   = sqrt(Tr + 2*s);
b11 = 1/t*(a11 + s);
b12 = 1/t*a12;
b21 = 1/t*a21;
b22 = 1/t*(a22 + s);
end

function  [L1,L2,s11,s21,s12,s22] = Eigs_2x2( axx, axz, azx, azz )
% Eigenvalues
Tr   = axx(:) + azz(:);
Det  = axx(:).*azz(:) - axz(:).*azx(:);
L1   = Tr./2 + sqrt( Tr.^2/4 - Det );
L2   = Tr./2 - sqrt( Tr.^2/4 - Det );

% Extract principal strain axis 2 -- finite strain
a22 = azz(:) - L2;
a21 = azx(:);
alp = (a22./a21);    % should be correct for a normal basis vector
s11 = ones(size(a21)) .* 1./sqrt(alp.^2+1) .* L2;
s12 = alp             .* s11;
s11 = reshape(s11,1,1)'; s11 = s11(1:1:size(s11,1), 1:1:size(s11,2));
s12 = reshape(s12,1,1)'; s12 = s12(1:1:size(s12,1), 1:1:size(s12,2));

% Extract principal strain axis 2 -- finite strain
a22 = azz(:) - L1;
a21 = azx(:);
alp = (a22./a21);    % should be correct for a normal basis vector
s21 = ones(size(a21)) .* 1./sqrt(alp.^2+1) .* L1;
s22 = alp             .* s21;
s21 = reshape(s21,1,1)'; s21 = s21(1:1:size(s21,1), 1:1:size(s21,2));
s22 = reshape(s22,1,1)'; s22 = s22(1:1:size(s22,1), 1:1:size(s22,2));
end