function [Hx, Hy, Hz, Ex, Ey, Ez] = Maxwell3DMat(Hx, Hy, Hz, Ex, Ey, Ez, FinalTime, source, source_coordinates)
% Purpose  : Integrate 3D Maxwell's until FinalTime starting with
%            initial conditions Hx,Hy,Hz, Ex,Ey,Ez
% the point source given by "source" is injected at the element nearest to
% "source_coordinates" into the Ez field
%
% In this function the Equation d/dt U = C*U from the Lawson Paper is
% solved using a standard LSRK scheme

Globals3D;
Ez_time = [];

% Calculate Data for plotting
% [x_grid, y_grid, sampleTets, sampleWeights] = CalcSamplingData(150);
% 
% f = figure('visible','off');
% PlotPlain3DFast(Ez, x_grid, y_grid, sampleTets, sampleWeights); drawnow; pause(0.1);
% filename = "field" + num2str(1+ ".png");
% saveas(f,filename);
% close;

U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

% find element to inject the source
idx = findNearestNode(source_coordinates);
%idx = [(1:Np)', idx(2)*ones(Np,1)]; 
idx = idxEH_to_idxU(3, idx);

% find element so sample field
idx_sample = idxEH_to_idxU(3, node_idx);

% compute time step size
%dt = dtscale3D;  % TW: buggy
dt = 0.4;
% correct dt for integer # of time steps
Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps

time = 0; tstep = 1;

nextplottime = 25;

% PML fields
P = zeros(2*3*Np*K,1);
resP = zeros(2*3*Np*K,1);
[sigma1, sigma2, sigma3] = constructSigma();

resU = zeros(2*3*Np*K,1);

tic;
while (time<FinalTime) % outer time step loop 
    
  % inject the source
  U(idx) = U(idx) + dt * source(time);
  
  % Runge-Kutta loop	
  for INTRK = 1:5   
    resU = rk4a(INTRK)*resU + dt*Ccoarse*U;
    resU = resU + dt * (blkmult(sigma1, U, Np, K) - P); % add PML terms
    U = U + rk4b(INTRK)*resU;
    
    % integrate PML terms
    resP = rk4a(INTRK)*resP + dt * (blkmult(sigma2,P,Np,K) + blkmult(sigma3,U,Np,K));
    P = P + rk4b(INTRK) * resP;
  end
   
   time = time+dt;    % Increment time
   tstep = tstep+1;
   
   %store field value over time
   Ez_time = [Ez_time, [time;U(idx_sample)]];
   U(idx_sample)
   
   % plot
%    if time > nextplottime
%        [~,~,~,~,~,Ez] = UToFields(U);
%        nextplottime = nextplottime + 25;
%        f = figure('visible','off');
%        PlotPlain3DFast(Ez, x_grid, y_grid, sampleTets, sampleWeights); drawnow; pause(0.01);
%        title(num2str(time));
%        filename = "field" + num2str(tstep + ".png");
%        saveas(f,filename);
%        close;
%   end
end
toc;

% convert U back to field components
[Hx,Hy,Hz,Ex,Ey,Ez] = UToFields(U);
return;