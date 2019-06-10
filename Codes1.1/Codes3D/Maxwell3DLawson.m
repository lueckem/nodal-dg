function [Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DLawson(Hx, Hy, Hz, Ex, Ey, Ez, FinalTime,source,source_coordinates)

% function [Hx,Hy,Hz,Ex,Ey,Ez] =Maxwell3D(Hx, Hy, Hz, Ex, Ey, Ez, FinalTime)
% Purpose  : Integrate 3D Maxwell's until FinalTime starting with
%            initial conditions Hx,Hy,Hz, Ex,Ey,Ez ;
% fine_idx = indices of elements belonging to the fine part of the grid;
% the point source given by "source" is injected at the node nearest to
% "source_coordinates" into the Ez field

Globals3D;
Ez_time = [];

% compute time step size
dt = dtscale3D;  % TW: buggy

% correct dt for integer # of time steps
Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps

time = 0; tstep = 1;

% initialize the vector of unknowns U. The ordering of the elements has to
% be consistent with the basis functions and the Matrices
U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

% find node to inject the source
idx = findNearestNode(source_coordinates);
idx = idxEH_to_idxU(3, idx);


% compute the constants in the Lawson-LRSK scheme
rk4cLawson = [0 rk4c];   % We need c(0). What should the value be?
%matexp1 = expm((1-rk4cLawson(6))*dt*Cfine);
matexp1 = speye(size(Cfine,1)); % for test

matexp2 = cell(1,5);
c1 = cell(1,5);
c2 = cell(1,5);
for k = 1:5
    %matexp2{k} = expm((rk4cLawson(k+1)-rk4cLawson(k))*dt*Cfine);
    matexp2{k} = speye(size(Cfine,1)); % for test
    c1{k} = rk4a(k) * matexp2{k};
    c2{k} = dt * Ccoarse * matexp2{k};
end


%store field value over time
Ez_time = [Ez_time, [time; U(idxEH_to_idxU(3, node_idx))]];
nextplottime = 0.1;

% outer time step loop 
while (time<FinalTime)
    
  % inject the source
  U(idx) = U(idx) + source(time);
  
  % Lawson-LSRK scheme
  phi1 = U;
  phi2 = zeros(size(U));
  for k = 1:5
      phi2 = c1{k} * phi2 + c2{k} * phi1;
      phi1 = matexp2{k} * phi1 + rk4b(k) * phi2;
  end
  U = matexp1 * phi1;
  
   % Increment time
   time = time+dt;
   tstep = tstep + 1;
   
   %store field value over time
   Ez_time = [Ez_time, [time; U(idxEH_to_idxU(3, node_idx))]];
   
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