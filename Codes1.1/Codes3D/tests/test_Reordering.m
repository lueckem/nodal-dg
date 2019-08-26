% Test if the Reordering is done properly

Globals3D; 
N = 1;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK5.neu');
StartUp3D;

% 5 middle elements
nodeIndex = findNearestNode([0,0,0]);
%fine_idx = [nodeIndex(2),EToE(nodeIndex(2),:)];
%fine_idx = unique(fine_idx);
fine_idx = nodeIndex(2);
InitMatLawsonSparse;

% Ex = sin(1*pi*x).*sin(1*pi*y);
% Ey = sin(2*pi*x).*sin(2*pi*y);
% Ez = sin(3*pi*x).*sin(3*pi*y);
% Hx = sin(4*pi*x).*sin(4*pi*y);
% Hy = sin(5*pi*x).*sin(5*pi*y);
% Hz = sin(6*pi*x).*sin(6*pi*y);
% U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

U = rand(6*Np*K,1);

%% RHS without Reordering
rhsU = (Ccoarse+Cfine) * U;

%% RHS with Reordering
ReorderLawson;
U = (Ccoarse+Cfine) * U;
ReorderBackLawson;

norm(rhsU - U)
