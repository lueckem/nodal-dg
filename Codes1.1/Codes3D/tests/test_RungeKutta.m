% compare Runge Kutta

Globals3D; 
N = 2;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('sph.neu');
StartUp3D;
epsilon = ones(K,1);

source = @(t) 100*exp(-0.005*(t-40)^2)*sin(pi*(t-40)/100); % one pulse at [0,80]
source_coordinates = [150,50,150];
idx = findNearestNode(source_coordinates);
idxU = idxEH_to_idxU(3, idx);

% compute time step size
dt = dtscale3D;
%% 
fine_idx = [];
InitMatLawsonSparse;
%%

% Ex = sin(1*pi*x).*sin(1*pi*y);
% Ey = sin(2*pi*x).*sin(2*pi*y);
% Ez = sin(3*pi*x).*sin(3*pi*y);
% Hx = sin(4*pi*x).*sin(4*pi*y);
% Hy = sin(5*pi*x).*sin(5*pi*y);
% Hz = sin(6*pi*x).*sin(6*pi*y);

Ex = rand(Np, K);
Ey = rand(Np, K);
Ez = rand(Np, K);
Hx = rand(Np, K);
Hy = rand(Np, K);
Hz = rand(Np, K);
U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

resHx = zeros(Np,K); resHy = zeros(Np,K); resHz = zeros(Np,K); 
resEx = zeros(Np,K); resEy = zeros(Np,K); resEz = zeros(Np,K);
resU = zeros(2*3*Np*K,1);

figure;
hold on;
xlabel('number of timesteps');
ylabel('norm(difference)');

for n = 1:1
    % inject source
    U(idxU) = U(idxU) + source(n*dt);
    Ez(idx(1), idx(2)) = Ez(idx(1),idx(2)) + source(n*dt);
    
for INTRK = 1:5
    % Original
    % compute right hand side of TM-mode Maxwell's equations
    [rhsHx, rhsHy, rhsHz, rhsEx, rhsEy, rhsEz] = MaxwellRHS3D(Hx,Hy,Hz,Ex,Ey,Ez);
    
    % initiate, increment Runge-Kutta residuals and update fields
    resHx = rk4a(INTRK)*resHx + dt*rhsHx;   Hx = Hx+rk4b(INTRK)*resHx;  	
    resHy = rk4a(INTRK)*resHy + dt*rhsHy;   Hy = Hy+rk4b(INTRK)*resHy;  	
    resHz = rk4a(INTRK)*resHz + dt*rhsHz;   Hz = Hz+rk4b(INTRK)*resHz;        
    resEx = rk4a(INTRK)*resEx + dt*rhsEx;   Ex = Ex+rk4b(INTRK)*resEx;  	  
    resEy = rk4a(INTRK)*resEy + dt*rhsEy;   Ey = Ey+rk4b(INTRK)*resEy;  	
    resEz = rk4a(INTRK)*resEz + dt*rhsEz;   Ez = Ez+rk4b(INTRK)*resEz; 
    
    % Matrix
    resU = rk4a(INTRK)*resU + dt*Ccoarse*U;
    U = U + rk4b(INTRK)*resU;
end
%compare
    U2 = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);
    
    plot(n, norm(U - U2), 'x'); drawnow;
end

% comparison
% [Hx2, Hy2, Hz2, Ex2, Ey2, Ez2] = UToFields(U);
U2 = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

% norm(Hx - Hx2)
% norm(Hy - Hy2)
% norm(Hz - Hz2)
% norm(Ex - Ex2)
% norm(Ey - Ey2)
% norm(Ez - Ez2)
norm(U - U2)