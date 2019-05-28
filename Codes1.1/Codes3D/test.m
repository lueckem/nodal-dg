%Test the Maxwell3D vs Maxwell3DLawson for a point source in the middle of
%the grid

Globals3D;

% Polynomial order of approximation 
N = 2;

% Read in Mesh
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK86.neu');

% Initialize solver and construct grid and metric
StartUp3D;

% Source function
source = @(t) sin(2*pi*t);
source_coordinates = [0,0,0];

%% Maxwell3D
% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);

% Solve Problem
FinalTime = 1;

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DPointSource(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

%% Lawson
% fine part of the mesh
fine_idx = [1,3];

% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);

% Solve Problem
FinalTime = 1;

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DLawson(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,fine_idx,source,source_coordinates);
%% Plot
figure;
plot(Ez_time)
%hold on
%plot(Ez2)
%legend("standard", "lawson")

