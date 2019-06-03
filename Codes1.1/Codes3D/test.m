%Test the Maxwell3D vs Maxwell3DLawson for a point source in the middle of
%the grid

Globals3D;

% Polynomial order of approximation 
N = 1;

% Read in Mesh
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK86.neu');

% Initialize solver and construct grid and metric
StartUp3D;

% Source function
source = @(t) sin(2*pi*t);
%source = @(t) 0.2 * exp(-0.5*(t-5).^2).*sin(2*pi*t);
source_coordinates = [0,0,0];

%sample node over time
node_idx = findNearestNode([0,0,0]);

FinalTime = 2;

%% Maxwell3D
% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DPointSource(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

figure;
plot(Ez_time(1,:), Ez_time(2,:), 'r');
hold on

%% Lawson
%% Initialize Matrices
% fine part of the mesh
fine_idx = [];

InitMatLawson;
%% Time integration

% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DLawson(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

plot(Ez_time(1,:), Ez_time(2,:), 'b');
legend('Maxwell3D', 'Lawson');


