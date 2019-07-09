function [D] = M2Di2_EvaluateConsistenTangent_DP_VM_incomp(plc,dgc,Txxc,Tyyc,Txyc,Jiic,dQdtxxc,dQdtyyc,dQdtxyc,phic,etac,plv,dgv,Txxv,Tyyv,Txyv,Jiiv,dQdtxxv,dQdtyyv,dQdtxyv,phiv,etav,Newton)
dQdtxx = dQdtxxc; dQdtyy = dQdtyyc; dQdtxy = dQdtxyc;
dFdtxx = dQdtxxc; dFdtyy = dQdtyyc; dFdtxy = dQdtxyc;  dFdP = -sin(phic);
D11 = 2*etac; D22 = 2*etac; D33 = 1*etac;

txx = Txxc; tyy = Tyyc; txy = Txyc; J2 = Jiic;
d2Qdsxxdsxx = 0.5000000000e0 .* J2 .^ (-0.1e1 ./ 0.2e1) - 0.2500000000e0 .* txx .^ 2 .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsxxdsyy = -0.2500000000e0 .* txx .* tyy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsxxdsxy = -0.5000000000e0 .* txx .* txy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsyydsxx = -0.2500000000e0 .* txx .* tyy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsyydsyy = 0.5000000000e0 .* J2 .^ (-0.1e1 ./ 0.2e1) - 0.2500000000e0 .* tyy .^ 2 .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsyydsxy = -0.5000000000e0 .* tyy .* txy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsxydsxx = -0.5000000000e0 .* txx .* txy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsxydsyy = -0.5000000000e0 .* tyy .* txy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsxydsxy = 0.1000000000e1 .* J2 .^ (-0.1e1 ./ 0.2e1) - 0.1000000000e1 .* txy .^ 2 .* J2 .^ (-0.3e1 ./ 0.2e1);

dlam = dgc;
if Newton == 1, dlam = 0*dlam; end
M11 = 1 + dlam .* D11 .* d2Qdsxxdsxx;
M12 = dlam .* D11 .* d2Qdsxxdsyy;
M13 = dlam .* D11 .* d2Qdsxxdsxy;
M21 = dlam .* D11 .* d2Qdsyydsxx;
M22 = 1 + dlam .* D11 .* d2Qdsyydsyy;
M23 = dlam .* D11 .* d2Qdsyydsxy;
M31 = dlam .* D33 .* d2Qdsxydsxx;
M32 = dlam .* D33 .* d2Qdsxydsyy;
M33 = 1 + dlam .* D33 .* d2Qdsxydsxy;
M14 = zeros(size(D11));
M24 = zeros(size(D11));
M34 = zeros(size(D11));
M41 = zeros(size(D11));
M42 = zeros(size(D11));
M43 = zeros(size(D11));
M44 =  ones(size(D11));
[ Mi11,Mi12,Mi13,Mi14, Mi21,Mi22,Mi23,Mi24, Mi31,Mi32,Mi33,Mi34, Mi41,Mi42,Mi43,Mi44  ] = M2Di2_inverse_4x4( M11,M12,M13,M14, M21,M22,M23,M24, M31,M32,M33,M34, M41,M42,M43,M44  );

% Consistent tangent - deviator
den = (dFdtxx.*Mi11+dFdtyy.*Mi21+dFdtxy.*Mi31).*D11.*dQdtxx + (dFdtxx.*Mi12+dFdtyy.*Mi22+dFdtxy.*Mi32).*D11.*dQdtyy + (dFdtxx.*Mi13+dFdtyy.*Mi23+dFdtxy.*Mi33).*D33.*dQdtxy;
den(plc==0) = 1;
D.D11c = Mi11 .* D11 .* (1 - 1 ./ den .* (dQdtxx .* dFdtxx .* Mi11 .* D11 + dQdtxx .* dFdtyy .* Mi21 .* D11 + dQdtxx .* dFdtxy .* Mi31 .* D11)) - Mi12 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi11 .* D11 + dQdtyy .* dFdtyy .* Mi21 .* D11 + dQdtyy .* dFdtxy .* Mi31 .* D11) - Mi13 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi11 .* D11 + dQdtxy .* dFdtyy .* Mi21 .* D11 + dQdtxy .* dFdtxy .* Mi31 .* D11);
D.D12c = -Mi11 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi12 .* D11 + dQdtxx .* dFdtyy .* Mi22 .* D11 + dQdtxx .* dFdtxy .* Mi32 .* D11) + Mi12 .* D11 .* (1 - 1 ./ den .* (dQdtyy .* dFdtxx .* Mi12 .* D11 + dQdtyy .* dFdtyy .* Mi22 .* D11 + dQdtyy .* dFdtxy .* Mi32 .* D11)) - Mi13 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi12 .* D11 + dQdtxy .* dFdtyy .* Mi22 .* D11 + dQdtxy .* dFdtxy .* Mi32 .* D11);
D.D13c = -Mi11 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi13 .* D33 + dQdtxx .* dFdtyy .* Mi23 .* D33 + dQdtxx .* dFdtxy .* Mi33 .* D33) - Mi12 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi13 .* D33 + dQdtyy .* dFdtyy .* Mi23 .* D33 + dQdtyy .* dFdtxy .* Mi33 .* D33) + Mi13 .* D33 .* (1 - 1 ./ den .* (dQdtxy .* dFdtxx .* Mi13 .* D33 + dQdtxy .* dFdtyy .* Mi23 .* D33 + dQdtxy .* dFdtxy .* Mi33 .* D33));
D.D14c = -Mi11 .* D11 ./ den .* dQdtxx .* dFdP - Mi12 .* D11 ./ den .* dQdtyy .* dFdP - Mi13 .* D33 ./ den .* dQdtxy .* dFdP;
D.D21c = Mi21 .* D11 .* (1 - 1 ./ den .* (dQdtxx .* dFdtxx .* Mi11 .* D11 + dQdtxx .* dFdtyy .* Mi21 .* D11 + dQdtxx .* dFdtxy .* Mi31 .* D11)) - Mi22 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi11 .* D11 + dQdtyy .* dFdtyy .* Mi21 .* D11 + dQdtyy .* dFdtxy .* Mi31 .* D11) - Mi23 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi11 .* D11 + dQdtxy .* dFdtyy .* Mi21 .* D11 + dQdtxy .* dFdtxy .* Mi31 .* D11);
D.D22c = -Mi21 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi12 .* D11 + dQdtxx .* dFdtyy .* Mi22 .* D11 + dQdtxx .* dFdtxy .* Mi32 .* D11) + Mi22 .* D11 .* (1 - 1 ./ den .* (dQdtyy .* dFdtxx .* Mi12 .* D11 + dQdtyy .* dFdtyy .* Mi22 .* D11 + dQdtyy .* dFdtxy .* Mi32 .* D11)) - Mi23 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi12 .* D11 + dQdtxy .* dFdtyy .* Mi22 .* D11 + dQdtxy .* dFdtxy .* Mi32 .* D11);
D.D23c = -Mi21 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi13 .* D33 + dQdtxx .* dFdtyy .* Mi23 .* D33 + dQdtxx .* dFdtxy .* Mi33 .* D33) - Mi22 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi13 .* D33 + dQdtyy .* dFdtyy .* Mi23 .* D33 + dQdtyy .* dFdtxy .* Mi33 .* D33) + Mi23 .* D33 .* (1 - 1 ./ den .* (dQdtxy .* dFdtxx .* Mi13 .* D33 + dQdtxy .* dFdtyy .* Mi23 .* D33 + dQdtxy .* dFdtxy .* Mi33 .* D33));
D.D24c = -Mi21 .* D11 ./ den .* dQdtxx .* dFdP - Mi22 .* D11 ./ den .* dQdtyy .* dFdP - Mi23 .* D33 ./ den .* dQdtxy .* dFdP;
D.D31c = Mi31 .* D11 .* (1 - 1 ./ den .* (dQdtxx .* dFdtxx .* Mi11 .* D11 + dQdtxx .* dFdtyy .* Mi21 .* D11 + dQdtxx .* dFdtxy .* Mi31 .* D11)) - Mi32 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi11 .* D11 + dQdtyy .* dFdtyy .* Mi21 .* D11 + dQdtyy .* dFdtxy .* Mi31 .* D11) - Mi33 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi11 .* D11 + dQdtxy .* dFdtyy .* Mi21 .* D11 + dQdtxy .* dFdtxy .* Mi31 .* D11);
D.D32c = -Mi31 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi12 .* D11 + dQdtxx .* dFdtyy .* Mi22 .* D11 + dQdtxx .* dFdtxy .* Mi32 .* D11) + Mi32 .* D11 .* (1 - 1 ./ den .* (dQdtyy .* dFdtxx .* Mi12 .* D11 + dQdtyy .* dFdtyy .* Mi22 .* D11 + dQdtyy .* dFdtxy .* Mi32 .* D11)) - Mi33 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi12 .* D11 + dQdtxy .* dFdtyy .* Mi22 .* D11 + dQdtxy .* dFdtxy .* Mi32 .* D11);
D.D33c = -Mi31 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi13 .* D33 + dQdtxx .* dFdtyy .* Mi23 .* D33 + dQdtxx .* dFdtxy .* Mi33 .* D33) - Mi32 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi13 .* D33 + dQdtyy .* dFdtyy .* Mi23 .* D33 + dQdtyy .* dFdtxy .* Mi33 .* D33) + Mi33 .* D33 .* (1 - 1 ./ den .* (dQdtxy .* dFdtxx .* Mi13 .* D33 + dQdtxy .* dFdtyy .* Mi23 .* D33 + dQdtxy .* dFdtxy .* Mi33 .* D33));
D.D34c = -Mi31 .* D11 ./ den .* dQdtxx .* dFdP - Mi32 .* D11 ./ den .* dQdtyy .* dFdP - Mi33 .* D33 ./ den .* dQdtxy .* dFdP;
D.D41c = 0;
D.D42c = 0;
D.D43c = 0;
D.D44c = zeros(size(plc));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dQdtxx = dQdtxxv; dQdtyy = dQdtyyv; dQdtxy = dQdtxyv;
dFdtxx = dQdtxxv; dFdtyy = dQdtyyv; dFdtxy = dQdtxyv; dFdP = -sin(phiv);
D11 = 2*etav; D22 = 2*etav; D33 = 1*etav; D44 = ones(size(etav));

txx = Txxv; tyy = Tyyv; txy = Txyv; J2 = Jiiv;
d2Qdsxxdsxx = 0.5000000000e0 .* J2 .^ (-0.1e1 ./ 0.2e1) - 0.2500000000e0 .* txx .^ 2 .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsxxdsyy = -0.2500000000e0 .* txx .* tyy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsxxdsxy = -0.5000000000e0 .* txx .* txy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsyydsxx = -0.2500000000e0 .* txx .* tyy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsyydsyy = 0.5000000000e0 .* J2 .^ (-0.1e1 ./ 0.2e1) - 0.2500000000e0 .* tyy .^ 2 .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsyydsxy = -0.5000000000e0 .* tyy .* txy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsxydsxx = -0.5000000000e0 .* txx .* txy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsxydsyy = -0.5000000000e0 .* tyy .* txy .* J2 .^ (-0.3e1 ./ 0.2e1);
d2Qdsxydsxy = 0.1000000000e1 .* J2 .^ (-0.1e1 ./ 0.2e1) - 0.1000000000e1 .* txy .^ 2 .* J2 .^ (-0.3e1 ./ 0.2e1);

dlam = dgv;
if Newton==1, dlam = 0*dlam; end
M11 = 1 + dlam .* D11 .* d2Qdsxxdsxx;
M12 = dlam .* D11 .* d2Qdsxxdsyy;
M13 = dlam .* D11 .* d2Qdsxxdsxy;
M21 = dlam .* D11 .* d2Qdsyydsxx;
M22 = 1 + dlam .* D11 .* d2Qdsyydsyy;
M23 = dlam .* D11 .* d2Qdsyydsxy;
M31 = dlam .* D33 .* d2Qdsxydsxx;
M32 = dlam .* D33 .* d2Qdsxydsyy;
M33 = 1 + dlam .* D33 .* d2Qdsxydsxy;
M14 = zeros(size(D11));
M24 = zeros(size(D11));
M34 = zeros(size(D11));
M41 = zeros(size(D11));
M42 = zeros(size(D11));
M43 = zeros(size(D11));
M44 = ones(size(D11));
[ Mi11,Mi12,Mi13,Mi14, Mi21,Mi22,Mi23,Mi24, Mi31,Mi32,Mi33,Mi34, Mi41,Mi42,Mi43,Mi44  ] = M2Di2_inverse_4x4( M11,M12,M13,M14, M21,M22,M23,M24, M31,M32,M33,M34, M41,M42,M43,M44  );

den = (dFdtxx.*Mi11+dFdtyy.*Mi21+dFdtxy.*Mi31).*D11.*dQdtxx+(dFdtxx.*Mi12+dFdtyy.*Mi22+dFdtxy.*Mi32).*D11.*dQdtyy+(dFdtxx.*Mi13+dFdtyy.*Mi23+dFdtxy.*Mi33).*D33.*dQdtxy;
den(plv==0) = 1;
D.D11v = Mi11 .* D11 .* (1 - 1 ./ den .* (dQdtxx .* dFdtxx .* Mi11 .* D11 + dQdtxx .* dFdtyy .* Mi21 .* D11 + dQdtxx .* dFdtxy .* Mi31 .* D11)) - Mi12 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi11 .* D11 + dQdtyy .* dFdtyy .* Mi21 .* D11 + dQdtyy .* dFdtxy .* Mi31 .* D11) - Mi13 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi11 .* D11 + dQdtxy .* dFdtyy .* Mi21 .* D11 + dQdtxy .* dFdtxy .* Mi31 .* D11);
D.D12v = -Mi11 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi12 .* D11 + dQdtxx .* dFdtyy .* Mi22 .* D11 + dQdtxx .* dFdtxy .* Mi32 .* D11) + Mi12 .* D11 .* (1 - 1 ./ den .* (dQdtyy .* dFdtxx .* Mi12 .* D11 + dQdtyy .* dFdtyy .* Mi22 .* D11 + dQdtyy .* dFdtxy .* Mi32 .* D11)) - Mi13 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi12 .* D11 + dQdtxy .* dFdtyy .* Mi22 .* D11 + dQdtxy .* dFdtxy .* Mi32 .* D11);
D.D13v = -Mi11 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi13 .* D33 + dQdtxx .* dFdtyy .* Mi23 .* D33 + dQdtxx .* dFdtxy .* Mi33 .* D33) - Mi12 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi13 .* D33 + dQdtyy .* dFdtyy .* Mi23 .* D33 + dQdtyy .* dFdtxy .* Mi33 .* D33) + Mi13 .* D33 .* (1 - 1 ./ den .* (dQdtxy .* dFdtxx .* Mi13 .* D33 + dQdtxy .* dFdtyy .* Mi23 .* D33 + dQdtxy .* dFdtxy .* Mi33 .* D33));
D.D14v = -Mi11 .* D11 ./ den .* dQdtxx .* dFdP - Mi12 .* D11 ./ den .* dQdtyy .* dFdP - Mi13 .* D33 ./ den .* dQdtxy .* dFdP;
D.D21v = Mi21 .* D11 .* (1 - 1 ./ den .* (dQdtxx .* dFdtxx .* Mi11 .* D11 + dQdtxx .* dFdtyy .* Mi21 .* D11 + dQdtxx .* dFdtxy .* Mi31 .* D11)) - Mi22 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi11 .* D11 + dQdtyy .* dFdtyy .* Mi21 .* D11 + dQdtyy .* dFdtxy .* Mi31 .* D11) - Mi23 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi11 .* D11 + dQdtxy .* dFdtyy .* Mi21 .* D11 + dQdtxy .* dFdtxy .* Mi31 .* D11);
D.D22v = -Mi21 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi12 .* D11 + dQdtxx .* dFdtyy .* Mi22 .* D11 + dQdtxx .* dFdtxy .* Mi32 .* D11) + Mi22 .* D11 .* (1 - 1 ./ den .* (dQdtyy .* dFdtxx .* Mi12 .* D11 + dQdtyy .* dFdtyy .* Mi22 .* D11 + dQdtyy .* dFdtxy .* Mi32 .* D11)) - Mi23 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi12 .* D11 + dQdtxy .* dFdtyy .* Mi22 .* D11 + dQdtxy .* dFdtxy .* Mi32 .* D11);
D.D23v = -Mi21 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi13 .* D33 + dQdtxx .* dFdtyy .* Mi23 .* D33 + dQdtxx .* dFdtxy .* Mi33 .* D33) - Mi22 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi13 .* D33 + dQdtyy .* dFdtyy .* Mi23 .* D33 + dQdtyy .* dFdtxy .* Mi33 .* D33) + Mi23 .* D33 .* (1 - 1 ./ den .* (dQdtxy .* dFdtxx .* Mi13 .* D33 + dQdtxy .* dFdtyy .* Mi23 .* D33 + dQdtxy .* dFdtxy .* Mi33 .* D33));
D.D24v = -Mi21 .* D11 ./ den .* dQdtxx .* dFdP - Mi22 .* D11 ./ den .* dQdtyy .* dFdP - Mi23 .* D33 ./ den .* dQdtxy .* dFdP;
D.D31v = Mi31 .* D11 .* (1 - 1 ./ den .* (dQdtxx .* dFdtxx .* Mi11 .* D11 + dQdtxx .* dFdtyy .* Mi21 .* D11 + dQdtxx .* dFdtxy .* Mi31 .* D11)) - Mi32 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi11 .* D11 + dQdtyy .* dFdtyy .* Mi21 .* D11 + dQdtyy .* dFdtxy .* Mi31 .* D11) - Mi33 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi11 .* D11 + dQdtxy .* dFdtyy .* Mi21 .* D11 + dQdtxy .* dFdtxy .* Mi31 .* D11);
D.D32v = -Mi31 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi12 .* D11 + dQdtxx .* dFdtyy .* Mi22 .* D11 + dQdtxx .* dFdtxy .* Mi32 .* D11) + Mi32 .* D11 .* (1 - 1 ./ den .* (dQdtyy .* dFdtxx .* Mi12 .* D11 + dQdtyy .* dFdtyy .* Mi22 .* D11 + dQdtyy .* dFdtxy .* Mi32 .* D11)) - Mi33 .* D33 ./ den .* (dQdtxy .* dFdtxx .* Mi12 .* D11 + dQdtxy .* dFdtyy .* Mi22 .* D11 + dQdtxy .* dFdtxy .* Mi32 .* D11);
D.D33v = -Mi31 .* D11 ./ den .* (dQdtxx .* dFdtxx .* Mi13 .* D33 + dQdtxx .* dFdtyy .* Mi23 .* D33 + dQdtxx .* dFdtxy .* Mi33 .* D33) - Mi32 .* D11 ./ den .* (dQdtyy .* dFdtxx .* Mi13 .* D33 + dQdtyy .* dFdtyy .* Mi23 .* D33 + dQdtyy .* dFdtxy .* Mi33 .* D33) + Mi33 .* D33 .* (1 - 1 ./ den .* (dQdtxy .* dFdtxx .* Mi13 .* D33 + dQdtxy .* dFdtyy .* Mi23 .* D33 + dQdtxy .* dFdtxy .* Mi33 .* D33));
D.D34v = -Mi31 .* D11 ./ den .* dQdtxx .* dFdP - Mi32 .* D11 ./ den .* dQdtyy .* dFdP - Mi33 .* D33 ./ den .* dQdtxy .* dFdP;
D.D41v = 0;
D.D42v = 0;
D.D43v = 0;
D.D44v = zeros(size(plv));
end