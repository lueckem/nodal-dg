function [Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DLawsonKrylov(Hx, Hy, Hz, Ex, Ey, Ez, FinalTime,source,source_coordinates)

% Purpose  : Integrate 3D Maxwell's until FinalTime starting with
%            initial conditions Hx,Hy,Hz, Ex,Ey,Ez ;
% fine_idx = indices of elements belonging to the fine part of the grid;
% the point source given by "source" is injected at the element nearest to
% "source_coordinates" into the Ez field

% In this version matrix exponentials of the form exp(alpha*Cfine)*v are
% calculated during runtime using krylov subspace methods

Globals3D;
Ez_time = [];

% Calculate Data for plotting
% [x_grid, y_grid, sampleTets, sampleWeights] = CalcSamplingData(0);
% 
% f = figure('visible','off');
% PlotPlain3DFast(Ez, x_grid, y_grid, sampleTets, sampleWeights); drawnow; pause(0.1);
% filename = "field" + num2str(1+ ".png");
% saveas(f,filename);
% close;

U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

% Modelreduction: reorders Cfine, Ccoarse and U
ReorderLawson;

% find element to inject the source
idx = findNearestNode(source_coordinates);
%idx = [(1:Np)', idx(2)*ones(Np,1)]; 
idx = idxEH_to_idxU(3, idx);

% find element so sample field
idx_sample = idxEH_to_idxU(3, node_idx);

% Adapt the indices to the new order using C_idx
[~, idx] = ismember(idx, C_idx);
[~, idx_sample] = ismember(idx_sample, C_idx);

% compute time step size
dt = dtscale3D;  % TW: buggy

% correct dt for integer # of time steps
Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps

time = 0; tstep = 1;

rk4cLawson = [0 rk4c];   % c_0 = 0;
size_Cff = 2*3*Np*length(fine_idx);
Cff = Cfine(1:size_Cff,1:size_Cff);
Ccf = Cfine(size_Cff+1:2*3*Np*K,1:size_Cff);

phi2 = zeros(size(U));

% PML fields
P = zeros(2*3*Np*K,1);
resP = zeros(2*3*Np*K,1);
[sigma1, sigma2, sigma3] = constructSigma();


nextplottime = 0.05;

% outer time step loop
tic
while (time<FinalTime)
    
  % inject the source
   U(idx) = U(idx) + dt * source(time);
  
  % Lawson-LSRK scheme
  phi1 = U;
  for k = 1:5
      phi1star = ExpCfinev((rk4cLawson(k+1)-rk4cLawson(k))*dt, phi1);
      phi2 = rk4a(k) * ExpCfinev((rk4cLawson(k+1)-rk4cLawson(k))*dt, phi2) + dt * Ccoarse * phi1star;
      phi2 = phi2 + dt * (blkmult(sigma1, phi1star, Np, K) - P); % add PML terms
      phi1 = phi1star + rk4b(k) * phi2;
      
      % integrate P
      resP = rk4a(k)*resP + dt * (blkmult(sigma2,P,Np,K) + blkmult(sigma3,phi1star,Np,K));
      P = P + rk4b(k) * resP;  
  end
  U = ExpCfinev((1-rk4cLawson(6))*dt, phi1);
  
   % Increment time
   time = time+dt;
   tstep = tstep + 1;
   
   %plot
%    if time > nextplottime
%        [~,~,~,~,~,Ez] = UToFields(ReorderBackULawson(U));
%        nextplottime = nextplottime + 0.05;
%        f = figure('visible','off');
%        PlotPlain3DFast(Ez, x_grid, y_grid, sampleTets, sampleWeights); drawnow; pause(0.01);
%        title(num2str(time));
%        filename = "field" + num2str(tstep + ".png");
%        saveas(f,filename);
%        close;
%   end
   
   %store field value over time
   Ez_time = [Ez_time, [time;U(idx_sample)]];
end
toc
% convert U back to field components
ReorderBackLawson;
[Hx,Hy,Hz,Ex,Ey,Ez] = UToFields(U);
return;