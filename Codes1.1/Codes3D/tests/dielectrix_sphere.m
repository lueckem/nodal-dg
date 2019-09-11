% The domain is [0,300]x[0,300]x[0,300].

Globals3D;

% Polynomial order of approximation 
N = 3;

% Read in Mesh
[Nv, VX, VY, VZ, K, EToV, epsilon] = MeshReaderGambit3DMaterial('sph.neu');

% Initialize solver and construct grid and metric
StartUp3D;

% Source function
source = @(t) 100*exp(-0.005*(t-50)^2)*sin(pi*t/5);
source_coordinates = [150,50,150];

%sample node over time
node_idx = findNearestNode([150,250,150]);

FinalTime = 1000;

% No PML
sigmax = 0*ones(1, K);
sigmay = sigmax; sigmaz = sigmax;

%% Maxwell3DMat
fine_idx = [];
InitMatLawsonSparse;
%%
% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DMat(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);
