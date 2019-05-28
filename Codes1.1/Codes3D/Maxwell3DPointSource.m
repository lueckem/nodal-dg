function [Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DPointSource(Hx, Hy, Hz, Ex, Ey, Ez, FinalTime, source, source_coordinates)

% function [Hx,Hy,Hz,Ex,Ey,Ez] =Maxwell3D(Hx, Hy, Hz, Ex, Ey, Ez, FinalTime)
% Purpose  : Integrate 3D Maxwell's until FinalTime starting with
%            initial conditions Hx,Hy,Hz, Ex,Ey,Ez
% the point source given by "source" is injected at the node nearest to
% "source_coordinates" into the Hx field

Globals3D;
Ez_time = [];

% find node to inject the source
min_dist = norm([x(1,1), y(1,1), z(1,1)] - source_coordinates);
idx = [1,1];
for k = 1:K
   for i = 1:Np
       dist = norm([x(i,k), y(i,k), z(i,k)] - source_coordinates);
      if dist < min_dist
          min_dist = dist;
          idx = [i,k];
      end
   end
end

% Runge-Kutta residual storage  
resHx = zeros(Np,K); resHy = zeros(Np,K); resHz = zeros(Np,K); 
resEx = zeros(Np,K); resEy = zeros(Np,K); resEz = zeros(Np,K); 

% compute time step size
dt = dtscale3D;  % TW: buggy

% correct dt for integer # of time steps
Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps;

time = 0; tstep = 1;

while (time<FinalTime) % outer time step loop 
    
  % inject the source
  Hx(idx(1), idx(2)) = Hx(idx(1), idx(2)) + source(time);
    
  for INTRK = 1:5   % inner multi-stage Runge-Kutta loop
    
    % compute right hand side of TM-mode Maxwell's equations
    [rhsHx, rhsHy, rhsHz, rhsEx, rhsEy, rhsEz] = MaxwellRHS3D(Hx,Hy,Hz, Ex, Ey, Ez);
    
    % initiate, increment Runge-Kutta residuals and update fields
    resHx = rk4a(INTRK)*resHx + dt*rhsHx;   Hx = Hx+rk4b(INTRK)*resHx;  	
    resHy = rk4a(INTRK)*resHy + dt*rhsHy;   Hy = Hy+rk4b(INTRK)*resHy;  	
    resHz = rk4a(INTRK)*resHz + dt*rhsHz;   Hz = Hz+rk4b(INTRK)*resHz;        
    resEx = rk4a(INTRK)*resEx + dt*rhsEx;   Ex = Ex+rk4b(INTRK)*resEx;  	  
    resEy = rk4a(INTRK)*resEy + dt*rhsEy;   Ey = Ey+rk4b(INTRK)*resEy;  	
    resEz = rk4a(INTRK)*resEz + dt*rhsEz;   Ez = Ez+rk4b(INTRK)*resEz;        
   end;
   
   time = time+dt;    % Increment time
   tstep = tstep+1;
   
   %store field value over time
   Ez_time = [Ez_time, Ez(1,75)];
end
return;


