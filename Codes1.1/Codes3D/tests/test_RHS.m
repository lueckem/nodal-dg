% Compare RHS
% Result: The RHS are the same!!

%clear all
Globals3D; 
N = 3;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK86.neu');
StartUp3D;


fine_idx = [];
InitMatLawsonSparse;
%%

Ex = sin(1*pi*x).*sin(1*pi*y);
Ey = sin(2*pi*x).*sin(2*pi*y);
Ez = sin(3*pi*x).*sin(3*pi*y);
Hx = sin(4*pi*x).*sin(4*pi*y);
Hy = sin(5*pi*x).*sin(5*pi*y);
Hz = sin(6*pi*x).*sin(6*pi*y);

% Ex = 0 * ones(Np,K);
% Ey = 0 * ones(Np,K);
% Ez = 0 * ones(Np,K);
% Hx = 0 * ones(Np,K);
% Hy = 0 * ones(Np,K);
% Hz = 0 * ones(Np,K);
% Hx(:,3) = 1 * ones(Np,1);

U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

% Original code
[rhsHx, rhsHy, rhsHz, rhsEx, rhsEy, rhsEz] = MaxwellRHS3D(Hx,Hy,Hz,Ex,Ey,Ez);
%rhsU2 = FieldsToU(rhsHx, rhsHy, rhsHz, rhsEx, rhsEy, rhsEz);

% Lawson
rhsU = Ccoarse * U;
[rhsHx2, rhsHy2, rhsHz2, rhsEx2, rhsEy2, rhsEz2] = UToFields(rhsU);

% Comparison
norm(rhsHx - rhsHx2)
norm(rhsHy - rhsHy2)
norm(rhsHz - rhsHz2)
norm(rhsEx - rhsEx2)
norm(rhsEy - rhsEy2)
norm(rhsEz - rhsEz2)