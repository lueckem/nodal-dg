Globals3D;

% Polynomial order of approximation 
N = 4;

% Read in Mesh
[Nv, VX, VY, VZ, K, EToV, epsilon] = MeshReaderGambit3DMaterial('cubeK5.neu');

% Initialize solver and construct grid and metric
StartUp3D;

% Source function
source = @(t) 0;
source_coordinates = [0,0,0];

% sample node over time
node_idx = findNearestNode([0.5,0.5,0]);

FinalTime = 6;

% no PML
sigmax = 0*ones(1, K);
sigmay = sigmax; sigmaz = sigmax;

% Build C_fine and C_coarse
fine_idx = 5;
InitMatLawsonSparse;

% initial conditions
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K);
xmode = 1; ymode = 1; 
Ez = sin(xmode*pi*x).*sin(ymode*pi*y);

%% Start Simulation
[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DLawsonKrylov(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);