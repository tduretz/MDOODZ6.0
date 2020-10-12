function main
clear

mode = 2;
xmin = -0.5;
xmax = 0.5;
Lx   = xmax-xmin;
Vx   = 1;
C    = 0.9;
Nx   = 101;
Nt   = 100;
dx   = (xmax-xmin)/(Nx-1);
dt   = C*dx/Vx;
xv   = xmin:dx:xmax;
xc   = xmin+dx/2:dx:xmax-dx/2;
xce  = xmin-dx/2:dx:xmax+dx/2;
% Set up markers
Nxm = 2;                         % 2 markers per cell
dxm = dx/Nxm;                    % Initial marker spacing
xm  = xmin+dxm/2:dxm:xmax-dxm/2;
% Define density on markers
rhom = ones(size(xm));
rhom(xm>-0.1 & xm<0.1) = 2;
% Interpolate density on markers to cell centers
rhoc  = markers2centers(rhom, xc,  xm, dx );
rhom  = centers2markers(rhoc, xce, xm, dx );
drhom = zeros(size(xm));

if mode == 1 % Standard
    for it=1:Nt
        
        xm          = xm + dt*Vx;
        xm(xm<xmin) = xm(xm<xmin) + Lx;
        xm(xm>xmax) = xm(xm>xmax) - Lx;
        rhoc        = markers2centers(rhom, xc,  xm, dx );
        
        figure(1), clf
        hold on
        plot(xc, rhoc, 'o')
        plot(xm, rhom, '*')
        drawnow
    end
    
elseif mode == 2 % New
    for it=1:Nt
        
        drhoc       = markers2centers(drhom, xc,  xm, dx );
        rhoc        = rhoc + drhoc;
        xm          = xm + dt*Vx;
        xm(xm<xmin) = xm(xm<xmin) + Lx;
        xm(xm>xmax) = xm(xm>xmax) - Lx;
        rhom0       = centers2markers(rhoc, xce, xm, dx );
        drhom       = (rhom-rhom0);
        
        figure(1), clf
        hold on
        plot(xc, rhoc, 'o')
        plot(xm, rhom, '*')
        drawnow
        
    end
end

end

function rhoc = markers2centers(rhom, xc, xm, dx )
nmark = length(xm);
rhoc  = zeros(size(xc));
ic    = int16((xm - xc(1))/ dx) + 1;
sumW  = zeros(size(xc));
for k=1:nmark
    icell       = ic(k);
    dst         = abs(xm(k) - xc(icell));
    w           = (1-dst/dx);
    rhoc(icell) = rhoc(icell) + w*rhom(k);
    sumW(icell) = sumW(icell) + w;
end
rhoc = rhoc./sumW;
end

function rhom = centers2markers(rhoc, xce, xm, dx )
Nx    = length(xce)-1;
nmark = length(xm);
rhom  = zeros(size(xm));
iW    = int16((xm - xce(1) - dx/2)/ dx) + 1;
for k=1:nmark
    icellW = iW(k);
    if icellW == 1
        West = Nx-1;
        East = icellW;
    elseif icellW == Nx
        West = Nx-1;
        East = 2;
    else
        West = icellW-1;
        East = West+1;
    end
    dst     = xm(k) - xce(icellW);
    w       = 1 - dst/dx;
    rhom(k) = w*rhoc(West) + (1-w)*rhoc(East);
end
end