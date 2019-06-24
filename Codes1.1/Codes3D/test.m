%Test the Maxwell3D vs Maxwell3DLawson for a point source in the middle of
%the grid

Globals3D;

% Polynomial order of approximation 
N = 1;

% Read in Mesh
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cube.neu');

% Initialize solver and construct grid and metric
StartUp3D;

% Source function
%source = @(t) sin(2*pi*t);
%source = @(t) 0.2 * exp(-5*(t-1).^2).*sin(4*pi*t);
source = @(t) 0;
source_coordinates = [0,0,0];

%sample node over time
node_idx = findNearestNode([0.5,0.5,0]);

FinalTime = 0.5;

%% Maxwell3D
% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);

node_source = findNearestNode(source_coordinates);
Ez(:,node_source(2)) = 1;

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DPointSource(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

% figure;
% plot(Ez_time(1,:), Ez_time(2,:), 'r');
% hold on

%% Lawson
%% Initialize Matrices
% fine part of the mesh
fine_idx = [];

InitMatLawsonSparse;
%% Time integration

% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DLawson(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

%plot(Ez_time(1,:), Ez_time(2,:), 'b');
%legend('Maxwell3D', 'Lawson');

%% Maxwell3DMat
fine_idx = [];
InitMatLawsonSparse;
%%
% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);
node_source = findNearestNode(source_coordinates);
Ez(:,node_source(2)) = 1;

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DMat(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

