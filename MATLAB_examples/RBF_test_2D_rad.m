function main
clear
% close all

radial_search = 0;


xmin = -1;
xmax = 1;
ymin = -1;
ymax = 1;
L    = (xmax - xmin);
H    = (ymax - ymin);
Nx   = 52;
Ny   = 52;
dx   =  L / (Nx-1);
dy   =  H / (Ny-1);
xg   = xmin:dx:xmax;
xe   = xmin-dx/2:dx:xmax+dx/2;
yg   = ymin:dy:ymax;
ye   = ymin-dy/2:dy:ymax+dy/2;

Nxm_cell  = 4;
Nym_cell  = 4;
dxm       = dx/Nxm_cell;
xm        = xmin+dxm/2:dxm:xmax-dxm/2;
dym       = dy/Nym_cell;
ym        = ymin+dym/2:dym:ymax-dym/2;
[xm2, ym2]=meshgrid(xm,ym);
xm        = xm2(:);
ym        = ym2(:);
Np        = length(xm);

% sinus pertubartion
n         = 5;
Tm        = sin(-n*2*pi.*(xm+dx/3)/L).*sin(-n*2*pi.*(ym+dy/3)/L);

% square
Tm        = ones(size(xm));
Tm( xm> -0.3 & xm<0.3 & ym>-0.3 & ym <0.3 ) = 0;

iCm       = zeros(size(Tm));
Tg        = zeros(Nx,Ny);
Wg        = zeros(Nx,Ny);
sig       = sqrt(dx^2+dy^2)*1;
p         = 0.1;

rc        = 1.0*sig;

for k=1:Np
    
    if radial_search == 1
        
        for i=1:Nx
            for j=1:Ny
                
                r = sqrt((xg(i)-xm(k))^2 +  (yg(j)-ym(k))^2);
                if r < rc
                    w = RBF_Weight_2D(xm(k),xg(i),ym(k),yg(j),sig,p);
                    Tg(i,j) = Tg(i,j) + w*Tm(k);
                    Wg(i,j) = Wg(i,j) + w;
                end
                
            end
        end
        
    else
        
        iC = fix((xm(k)-xe(1))/dx)+1; % corresponding cell
        jC = fix((ym(k)-ye(1))/dy)+1; % corresponding cell
                
        iW = iC-1;
        iE = iC+1;
        jS = jC-1;
        jN = jC+1;
        
        % W
        if  iW > 0
            i = iW; j= jC;
            w = RBF_Weight_2D(xm(k),xg(i),ym(k),yg(j),sig,p);
            Tg(i,j) = Tg(i,j) + w*Tm(k);
            Wg(i,j) = Wg(i,j) + w;
        end
        
        
        % SW
        if  iW > 0 && jS > 0
            i = iW; j= jS;
            w = RBF_Weight_2D(xm(k),xg(i),ym(k),yg(j),sig,p);
            Tg(i,j) = Tg(i,j) + w*Tm(k);
            Wg(i,j) = Wg(i,j) + w;
        end
        
        % S
        if  jS > 0
            i = iC; j= jS;
            w = RBF_Weight_2D(xm(k),xg(i),ym(k),yg(j),sig,p);
            Tg(i,j) = Tg(i,j) + w*Tm(k);
            Wg(i,j) = Wg(i,j) + w;
        end
        
        % SE
        if  iE < Nx+1 && jS > 0
            i = iE; j= jS;
            w = RBF_Weight_2D(xm(k),xg(i),ym(k),yg(j),sig,p);
            Tg(i,j) = Tg(i,j) + w*Tm(k);
            Wg(i,j) = Wg(i,j) + w;
        end
        
        i = iC; j= jC;
        w = RBF_Weight_2D(xm(k),xg(i),ym(k),yg(j),sig,p);
        Tg(i,j) = Tg(i,j) + w*Tm(k);
        Wg(i,j) = Wg(i,j) + w;
        
        % NW
        if  iW > 0 && jN < Ny+1
            i = iW; j= jN;
            w = RBF_Weight_2D(xm(k),xg(i),ym(k),yg(j),sig,p);
            Tg(i,j) = Tg(i,j) + w*Tm(k);
            Wg(i,j) = Wg(i,j) + w;
        end
        
        % N
        if  jN < Ny+1
            i = iC; j= jN;
            w = RBF_Weight_2D(xm(k),xg(i),ym(k),yg(j),sig,p);
            Tg(i,j) = Tg(i,j) + w*Tm(k);
            Wg(i,j) = Wg(i,j) + w;
        end
        
        % NE
        if  iE < Nx+1 && jN < Ny+1
            i = iE; j= jN;
            w = RBF_Weight_2D(xm(k),xg(i),ym(k),yg(j),sig,p);
            Tg(i,j) = Tg(i,j) + w*Tm(k);
            Wg(i,j) = Wg(i,j) + w;
        end
        
        % E
        if  iE < Nx+1 
            i = iE; j= jC;
            w = RBF_Weight_2D(xm(k),xg(i),ym(k),yg(j),sig,p);
            Tg(i,j) = Tg(i,j) + w*Tm(k);
            Wg(i,j) = Wg(i,j) + w;
        end
        
    end
end

Tg = Tg./Wg;

figure(1), clf
subplot(211)
scatter(xm, ym, 30, Tm, 'filled'), axis xy image, colorbar
subplot(212)
imagesc(xg, yg, Tg), axis xy image, colorbar


end

function w = RBF_Weight_2D(xm,xg,ym,yg,sig,p)

d = sqrt((xg-xm)^2 + (yg-ym)^2);
%     w = exp( -d/sig^2);
w = 1 / d^p;

end
