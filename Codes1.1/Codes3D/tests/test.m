Globals3D;

% Polynomial order of approximation 
N = 2;

% Read in Mesh
[Nv, VX, VY, VZ, K, EToV, epsilon] = MeshReaderGambit3DMaterial('cube.neu');
epsilon = ones(K,1);
epsilon(1:100) = 2;

% Initialize solver and construct grid and metric
StartUp3D;

% Source function
%source = @(t) sin(pi*t);
source = @(t) exp(-5*(t-1).^2).*sin(4*pi*t);
%source = @(t) 0;
source_coordinates = [0,0,0];

%sample node over time
node_idx = findNearestNode([0.5,0.5,0]);

FinalTime = 6;

%PML every element that is at the boundary
sigmax = 0*ones(1, K);
%sigmax = [1, 1, 1, 1, 0];
% for i = 1:K
%    for j = 1:4
%       if EToE(i,j) == i
%           sigmax(i) = 1;
%           break;
%       end
%    end
% end
sigmay = sigmax; sigmaz = sigmax;
%% Maxwell3D
% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);

% mode
%xmode = 1; ymode = 1; zmode = 1;
%Ez = sin(xmode*pi*x).*sin(ymode*pi*y);%.*sin(zmode*pi*z);

% node_source = findNearestNode(source_coordinates);
% Ez(:,node_source(2)) = 1;

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DPointSource(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

%% Lawson
%% Initialize Matrices
% fine part of the mesh = 5 middle elements
% nodeIndex = findNearestNode([0,0,0]);
% fine_idx = [nodeIndex(2),EToE(nodeIndex(2),:)];
% fine_idx = unique(fine_idx);

fine_idx = 5;
%fine_idx = [190,199,207,232,246,248]; %Bad elements in cubeK268BAD

InitMatLawsonSparse;

% Time integration

% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = ones(Np, K)* eps;

% mode
xmode = 1; ymode = 1; zmode = 1; 
Ez = sin(xmode*pi*x).*sin(ymode*pi*y);%.*sin(zmode*pi*z);

% load initial conditions from file
%load("init_cond.mat")

[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DLawsonKrylov(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

%plot(Ez_time(1,:), Ez_time(2,:), 'b');
%legend('Maxwell3D', 'Lawson');

%% Maxwell3DMat
fine_idx = [];

InitMatLawsonSparse;

%%

% zero initial condition 
Hx = zeros(Np, K); Hy = zeros(Np, K); Hz = zeros(Np, K);
Ex = zeros(Np, K); Ey = zeros(Np, K); Ez = zeros(Np, K);
% Ez = exp(-20*(x.^2 + y.^2));
% 1 element
% node_source = findNearestNode(source_coordinates);
% Ez(:,node_source(2)) = 1;
% Hx(:,node_source(2)) = 1;
% Hy(:,node_source(2)) = -1;
% mode
% xmode = 1; ymode = 1;
% Ez = sin(xmode*pi*x).*sin(ymode*pi*y);

% load initial conditions from file
%load("init_cond.mat")
% 
[Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DMat(Hx,Hy,Hz,Ex,Ey,Ez,FinalTime,source,source_coordinates);

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

figure;
hold on;
plot(Ez_time_mat02(1,:), Ez_time_mat02(2,:), 'k-o');
plot(Ez_time_krylov06(1,:), Ez_time_krylov06(2,:), 'b-x');
plot(Ez_time_krylov10(1,:), Ez_time_krylov10(2,:), 'r-x');
% plot(Ez_time_krylov15(1,:), Ez_time_krylov15(2,:), 'm-x');
% plot(Ez_time_krylov20(1,:), Ez_time_krylov20(2,:), 'g-x');
% plot(Ez_time_krylov30(1,:), Ez_time_krylov30(2,:), 'b-x');
legend("LSRK, dt=0.2","Lawson, dt=0.6","Lawson, dt=1")%,"Lawson, dt=1.5","Lawson, dt=2","Lawson, dt=3.0")

