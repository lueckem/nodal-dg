function [Hx, Hy, Hz, Ex, Ey, Ez] = Maxwell3DMat(Hx, Hy, Hz, Ex, Ey, Ez, FinalTime, source, source_coordinates)
% Purpose  : Integrate 3D Maxwell's until FinalTime starting with
%            initial conditions Hx,Hy,Hz, Ex,Ey,Ez
% the point source given by "source" is injected at the node nearest to
% "source_coordinates" into the Ez field
%
% In this function the Equation d/dt U = C*U from the Lawson Paper is
% solved using a standard LSRK scheme

Globals3D;
Ez_time = [];

U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

% find node to inject the source
idx = findNearestNode(source_coordinates);
idx = idxEH_to_idxU(3, idx);

% compute time step size
dt = dtscale3D;  % TW: buggy

% correct dt for integer # of time steps
Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps

time = 0; tstep = 1;

nextplottime = 0.1;

while (time<FinalTime) % outer time step loop 
    
  % inject the source
  U(idx) = U(idx) + source(time);
  
  % Runge-Kutta loop	
  resU = zeros(2*3*Np*K,1);
  for INTRK = 1:5   
    resU = rk4a(INTRK)*resU + dt*Ccoarse*U;
    U = U + rk4b(INTRK)*resU;
  end
   
   time = time+dt;    % Increment time
   tstep = tstep+1;
   
   %store field value over time
   Ez_time = [Ez_time, [time;Ez(node_idx(1),node_idx(2))]];
   
   %plot
   if time > nextplottime
       [Hx,Hy,Hz,Ex,Ey,Ez] = UToFields(U);
       nextplottime = nextplottime + 0.2;
       f = figure('visible','off');
       PlotPlain3D(0, Ez); drawnow; pause(0.1);
       filename = "field" + num2str(tstep + ".png");
       saveas(f,filename)
   end
end

% convert U back to field components
[Hx,Hy,Hz,Ex,Ey,Ez] = UToFields(U);

return;

