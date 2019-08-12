% Compare speed of original RHS with RHS = C * U, ONLY timestepping

% Result: original 23s, Lawson 6s

Globals3D; 
N = 3;
FinalTime = 1;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK86.neu');
StartUp3D;

% compute time step size
dt = dtscale3D;  % TW: buggy

% correct dt for integer # of time steps
Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps;


fine_idx = [];
InitMatLawsonSparse;
%%

Ex = sin(1*pi*x).*sin(1*pi*y);
Ey = sin(2*pi*x).*sin(2*pi*y);
Ez = sin(3*pi*x).*sin(3*pi*y);
Hx = sin(4*pi*x).*sin(4*pi*y);
Hy = sin(5*pi*x).*sin(5*pi*y);
Hz = sin(6*pi*x).*sin(6*pi*y);

U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

% original Runge Kutta 
resHx = zeros(Np,K); resHy = zeros(Np,K); resHz = zeros(Np,K); 
resEx = zeros(Np,K); resEy = zeros(Np,K); resEz = zeros(Np,K);

tic;
for n = 1:100
for INTRK = 1:5
    % compute right hand side of TM-mode Maxwell's equations
    [rhsHx, rhsHy, rhsHz, rhsEx, rhsEy, rhsEz] = MaxwellRHS3D(Hx,Hy,Hz,Ex,Ey,Ez);
    
    % initiate, increment Runge-Kutta residuals and update fields
    resHx = rk4a(INTRK)*resHx + dt*rhsHx;   Hx = Hx+rk4b(INTRK)*resHx;  	
    resHy = rk4a(INTRK)*resHy + dt*rhsHy;   Hy = Hy+rk4b(INTRK)*resHy;  	
    resHz = rk4a(INTRK)*resHz + dt*rhsHz;   Hz = Hz+rk4b(INTRK)*resHz;        
    resEx = rk4a(INTRK)*resEx + dt*rhsEx;   Ex = Ex+rk4b(INTRK)*resEx;  	  
    resEy = rk4a(INTRK)*resEy + dt*rhsEy;   Ey = Ey+rk4b(INTRK)*resEy;  	
    resEz = rk4a(INTRK)*resEz + dt*rhsEz;   Ez = Ez+rk4b(INTRK)*resEz;        
end
end
toc;

% Lawson Runge Kutta
resU = zeros(2*3*Np*K,1);

tic;
for n = 1:100
for INTRK = 1:5
    resU = rk4a(INTRK)*resU + dt*Ccoarse*U;
    U = U + rk4b(INTRK)*resU;
end
end
toc;
