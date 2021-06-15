clear

style = 1; % 1: sinus, 2:Gaussian

% Domain
xmin = -2;
xmax =  2;
ymin = -1;
ymax =  1;

% Discretisation
ncx = 400;
ncy = 200;
dx  = (xmax-xmin)/ncx;
dy  = (ymax-ymin)/ncy;
xc  = xmin+dx/2:dx:xmax-dx/2;
yc  = ymin+dy/2:dy:ymax-dy/2;
[xc2,yc2] = ndgrid(xc,yc);

% Layer parameters
h   = 0.1;
tet = 35*pi/180;
Lx  = xmin-xmax;

% Perturbation
a   = 1*h/2;
sig = Lx/20;

% Coordinates of upper layer limit
X1 = (yc2+h*cos(tet)).*sin(tet) + (xc2-h*sin(tet)).*cos(tet);
Y1 = (yc2+h*cos(tet)).*cos(tet) - (xc2-h*sin(tet)).*sin(tet);

% Coordinates of lower layer limit
X2 = (yc2-h*cos(tet)).*sin(tet) + (xc2+h*sin(tet)).*cos(tet);
Y2 = (yc2-h*cos(tet)).*cos(tet) - (xc2+h*sin(tet)).*sin(tet);

if style==1
    % Sinus
    f1 = Y1 - a*sin(6*pi*(X1)/Lx);
    f2 = Y2 - a*sin(6*pi*(X2)/Lx);
else
    % Gaussian
    f1 = Y1 - a*exp(-(X1/sig).^2);
    f2 = Y2 + a*exp(-(X2/sig).^2);
end

% Phase
phase  = zeros(size(xc2));
phase(f1>0 & f2<0) = 1;

% Viz
figure(1)
imagesc(xc,yc,phase'), axis image xy
