 function [ ev ] = M2Di_EP9_centroids2vertices( ec )
nx    = size(ec,1);
ny    = size(ec,2);
wSW   = ones(nx+1,ny+1); wSW(end,:) = 2; wSW(:,end) = 2; wSW(end,end) = 4; wSW(1  ,:) = 0; wSW(:,  1) = 0;
wSE   = ones(nx+1,ny+1); wSE(  1,:) = 2; wSE(:,end) = 2; wSE(  1,end) = 4; wSE(end,:) = 0; wSE(:,  1) = 0; 
wNW   = ones(nx+1,ny+1); wNW(end,:) = 2; wNW(:,  1) = 2; wNW(end,  1) = 4; wNW(1  ,:) = 0; wNW(:,end) = 0;
wNE   = ones(nx+1,ny+1); wNE(  1,:) = 2; wNE(:,  1) = 2; wNE(  1,  1) = 4; wNE(end,:) = 0; wNE(:,end) = 0;
eSW = zeros(nx+1,ny+1);  eSW(2:end-0,2:end-0) = ec;
eSE = zeros(nx+1,ny+1);  eSE(1:end-1,2:end-0) = ec;
eNW = zeros(nx+1,ny+1);  eNW(2:end-0,1:end-1) = ec;
eNE = zeros(nx+1,ny+1);  eNE(1:end-1,1:end-1) = ec;
ev  = 0.25*(wSW.*eSW + wSE.*eSE + wNW.*eNW + wNE.*eNE );
end