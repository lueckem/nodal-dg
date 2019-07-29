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
%source = @(t) sin(2*pi*t);
%source = @(t) 0.2 * exp(-5*(t-1).^2).*sin(4*pi*t);
source = @(t) 0;
source_coordinates = [0,0,0];

%sample node over time
node_idx = findNearestNode([0.5,0.5,0]);

FinalTime = 1;

%% Maxwell3D
% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);
xmode = 1; ymode = 1; 
Ez = sin(xmode*pi*x).*sin(ymode*pi*y);

% node_source = findNearestNode(source_coordinates);
% Ez(:,node_source(2)) = 1;

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DPointSource(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

%% Lawson
%% Initialize Matrices
% fine part of the mesh
fine_idx = [100, 52, 679, 1001, 1, 128, 129, 1002];

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
tic;
InitMatLawsonSparse;
toc;
%%
% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);
% 1 element
%node_source = findNearestNode(source_coordinates);
%Ez(:,node_source(2)) = 1;
% mode
xmode = 1; ymode = 1; 
Ez = sin(xmode*pi*x).*sin(ymode*pi*y);

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DMat(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

%% Plotting
% for i = size(Ez_time,2):-1:1
%    if abs(Ez_time(2,i)) > 40
%        Ez_time(:,i) = [];
%    end
% end
for i = size(Ez_time2,2):-1:1
   if abs(Ez_time2(2,i)) > 40
       Ez_time2(:,i) = [];
   end
end

figure;
plot(Ez_time1(1,:), Ez_time1(2,:), 'b');
hold on;
plot(Ez_time2(1,:), Ez_time2(2,:), 'r');
%plot(Ez_time(1,:), Ez_time(2,:), 'k');
legend("upwind", "central", "Lawson");

