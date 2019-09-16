% The domain is [0,300]x[0,300]x[0,300].

Globals3D;

% Polynomial order of approximation 
N = 2;

% Read in Mesh
[Nv, VX, VY, VZ, K, EToV, epsilon] = MeshReaderGambit3DMaterial('sphBAD.neu');
%epsilon = ones(K,1);

% Initialize solver and construct grid and metric
StartUp3D;

% Source function
source = @(t) 100*exp(-0.005*(t-40)^2)*sin(pi*(t-40)/100); % one pulse at [0,80]
source_coordinates = [150,50,150];

%sample node over time
node_idx = findNearestNode([150,250,150]);

FinalTime = 300;

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


%% Original
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DPointSource(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);
%% Krylov
fine_idx = [86 128 149 221 331 347 461 490 491 492 493 494       495       496       497       498       499       500       501       502       503       504       505       506];

InitMatLawsonSparse;

%%
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = ones(Np, K)* eps;
[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DLawsonKrylov(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);
