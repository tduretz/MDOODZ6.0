function main
clear 
% close all

xmin = -1/2;
xmax = 1/2;
L    = (xmax - xmin);
Nx   = 11;
dx   =  L / (Nx-1);
xg   = xmin:dx:xmax;
xe   = xmin-dx/2:dx:xmax+dx/2;

Nxm_cell  = 20;
dxm       = dx/Nxm_cell;
xm        = xmin+dxm/2:dxm:xmax-dxm/2;

n         = 5;
Tm        = sin(-n*2*pi.*(xm+dx/3)/L);
Np        = length(xm);

iCm       = zeros(size(Tm));
Tg        = zeros(size(xg));
Wg        = zeros(size(xg));
sig       = dx/4;
p         = 10;
for k=1:Np
    
    iC = fix((xm(k)-xe(1))/dx)+1; % corresponding cell
    
    iCm(k) = iC;
    
    iW = iC-1;
    iE = iC+1;
    
    if  iW > 0 % W
        wW = RBF_Weight(xm(k),xg(iW),sig,p);
        Tg(iW) = Tg(iW) + wW*Tm(k);
        Wg(iW) = Wg(iW) + wW;
    end
    wC = RBF_Weight(xm(k),xg(iC),sig,p);
    Tg(iC) = Tg(iC) + wC*Tm(k);
    Wg(iC) = Wg(iC) + wC;
    
    if  iE < Nx+1 % E
        wE = RBF_Weight(xm(k),xg(iE),sig,p);
        Tg(iE) = Tg(iE) + wE*Tm(k);
        Wg(iE) = Wg(iE) + wE;
    end
    
end

iCm

Tg = Tg./Wg;


figure(1), clf
subplot(211)
hold on
plot(xm, xm*0, '*b')
plot(xg, xg*0, '+k')
subplot(212)
hold on
plot(xm, Tm, '-.')
plot(xg, Tg, 'or')


end

function w = RBF_Weight(xm,xg,sig,p)

    w = exp( - abs(xm-xg)/sig^2);
%     w = 1 / abs(xm-xg)^p;

end
