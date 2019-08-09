Globals3D; 
N = 3;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK268.neu');
StartUp3D;

%% test conversion between U and Fields

Ex = sin(1*pi*x).*sin(1*pi*y);
Ey = sin(2*pi*x).*sin(2*pi*y);
Ez = sin(3*pi*x).*sin(3*pi*y);
Hx = sin(4*pi*x).*sin(4*pi*y);
Hy = sin(5*pi*x).*sin(5*pi*y);
Hz = sin(6*pi*x).*sin(6*pi*y);

U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);
[Hx2, Hy2, Hz2, Ex2, Ey2, Ez2] = UToFields(U);

norm(Hx - Hx2)
norm(Hy - Hy2)
norm(Hz - Hz2)
norm(Ex - Ex2)
norm(Ey - Ey2)
norm(Ez - Ez2)

% Conclusion: the conversions work properly.

%% test conversion for idx

Ex = sin(1*pi*x).*sin(1*pi*y);
Ey = sin(2*pi*x).*sin(2*pi*y);
Ez = sin(3*pi*x).*sin(3*pi*y);
Hx = sin(4*pi*x).*sin(4*pi*y);
Hy = sin(5*pi*x).*sin(5*pi*y);
Hz = sin(6*pi*x).*sin(6*pi*y);

U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

% add some value on one element
idx = findNearestNode([0,0,0]);
Ez(:, idx(2)) = Ez(:, idx(2)) + 1000;

idx = [(1:Np)', idx(2)*ones(Np,1)]; 
idx = idxEH_to_idxU(3, idx);
U(idx) = U(idx) + 1000;

[Hx2, Hy2, Hz2, Ex2, Ey2, Ez2] = UToFields(U);

norm(Hx - Hx2)
norm(Hy - Hy2)
norm(Hz - Hz2)
norm(Ex - Ex2)
norm(Ey - Ey2)
norm(Ez - Ez2)

%% test conversion and reordering
fine_idx = [145, 345, 2, 1002, 53];

Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);

Ex = sin(1*pi*x).*sin(1*pi*y);
Ey = sin(2*pi*x).*sin(2*pi*y);
Ez = sin(3*pi*x).*sin(3*pi*y);
Hx = sin(4*pi*x).*sin(4*pi*y);
Hy = sin(5*pi*x).*sin(5*pi*y);
Hz = sin(6*pi*x).*sin(6*pi*y);

U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);
ReorderLawson;
ReorderBackLawson;
[Hx2, Hy2, Hz2, Ex2, Ey2, Ez2] = UToFields(U);

norm(Hx - Hx2)
norm(Hy - Hy2)
norm(Hz - Hz2)
norm(Ex - Ex2)
norm(Ey - Ey2)
norm(Ez - Ez2)