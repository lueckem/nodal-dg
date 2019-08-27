%Test the Maxwell3D vs Maxwell3DLawson for a point source in the middle of
%the grid

% Boundary Conditions instable??

Globals3D;

% Polynomial order of approximation 
N = 2;

% Read in Mesh
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK268.neu');

% Initialize solver and construct grid and metric
StartUp3D;

% Source function
source = @(t) sin(4*pi*t);
%source = @(t) 0.2 * exp(-5*(t-1).^2).*sin(4*pi*t);
%source = @(t) 0;
source_coordinates = [0,0,0];

%sample node over time
node_idx = findNearestNode([0.25,0.25,0]);

FinalTime = 0.2;

%% Maxwell3D
% % zero initial condition 
% Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
% Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);
% %xmode = 1; ymode = 1; 
% %Ez = sin(xmode*pi*x).*sin(ymode*pi*y);
% 
% % node_source = findNearestNode(source_coordinates);
% % Ez(:,node_source(2)) = 1;
% 
% [Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DPointSource(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

%% Lawson
%% Initialize Matrices
% fine part of the mesh = 5 middle elements
% nodeIndex = findNearestNode([0,0,0]);
% fine_idx = [nodeIndex(2),EToE(nodeIndex(2),:)];
% fine_idx = unique(fine_idx);

fine_idx = 1:100;

%fine_idx = [190,199,207,232,246,248]; %Bad elements in cubeK268BAD

InitMatLawsonSparse;

% Time integration

% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = ones(Np, K)* eps;

% load initial conditions from file
%load("init_cond.mat")

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DLawsonKrylov(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

%plot(Ez_time(1,:), Ez_time(2,:), 'b');
%legend('Maxwell3D', 'Lawson');

%% Maxwell3DMat
% fine_idx = [];
% 
% InitMatLawsonSparse;

%%
% zero initial condition 
% Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
% Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);
% 1 element
%node_source = findNearestNode(source_coordinates);
%Ez(:,node_source(2)) = 1;
% mode
%xmode = 1; ymode = 1; 
%Ez = sin(xmode*pi*x).*sin(ymode*pi*y);

% load initial conditions from file
% load("init_cond.mat")
% 
% [Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DMat(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

%% Plotting
% for i = size(Ez_time,2):-1:1
%    if abs(Ez_time(2,i)) > 40
%        Ez_time(:,i) = [];
%    end
% end
% for i = size(Ez_time2,2):-1:1
%    if abs(Ez_time2(2,i)) > 40
%        Ez_time2(:,i) = [];
%    end
% end

% figure;
% hold on;
% plot(Ez_time_base(1,:), Ez_time_base(2,:), 'r');
% plot(Ez_time(1,:), Ez_time(2,:), 'b');
% legend("base", "new");

