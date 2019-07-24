% check if Differentiation matrices work properly

Globals3D; 
N = 1;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cube.neu');
StartUp3D;
%%
Ez = sin(pi*x).*sin(pi*y);

DxEz_true = pi.*cos(pi*x).*sin(pi*y);
DyEz_true = sin(pi*x).*pi.*cos(pi*y);
DzEz_true = zeros(Np, K);

DxEz = zeros(Np, K);
DyEz = zeros(Np, K);
DzEz = zeros(Np, K);

% DxEz = rx.*(Dr*Ez) + sx.*(Ds*Ez) + tx.*(Dt*Ez);
% DyEz = ry.*(Dr*Ez) + sy.*(Ds*Ez) + ty.*(Dt*Ez);
% DzEz = rz.*(Dr*Ez) + sz.*(Ds*Ez) + tz.*(Dt*Ez);

for i = 1:K
    Dx = rx(1,i) * Dr + sx(1,i) * Ds + tx(1,i) * Dt;
    Dy = ry(1,i) * Dr + sy(1,i) * Ds + ty(1,i) * Dt;
    Dz = rz(1,i) * Dr + sz(1,i) * Ds + tz(1,i) * Dt;
    
    DxEz(:,i) = Dx * Ez(:,i);
    DyEz(:,i) = Dy * Ez(:,i);
    DzEz(:,i) = Dz * Ez(:,i);
end

normDx = norm(DxEz_true - DxEz);
normDy = norm(DyEz_true - DyEz);
normDz = norm(DzEz_true - DzEz);

figure;
PlotPlain3D(0, Ez)
figure;
PlotPlain3D(0, DxEz_true)
figure;
PlotPlain3D(0, DxEz)
