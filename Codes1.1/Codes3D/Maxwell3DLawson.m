function [Hx,Hy,Hz,Ex,Ey,Ez] = Maxwell3DLawson(Hx, Hy, Hz, Ex, Ey, Ez, FinalTime, fine_idx)

% function [Hx,Hy,Hz,Ex,Ey,Ez] =Maxwell3D(Hx, Hy, Hz, Ex, Ey, Ez, FinalTime)
% Purpose  : Integrate 3D Maxwell's until FinalTime starting with
%            initial conditions Hx,Hy,Hz, Ex,Ey,Ez
% fine_idx = indices of elements belonging to the fine part of the grid

Globals3D;
Ez_time = [];

% compute time step size
dt = dtscale3D;  % TW: buggy

% correct dt for integer # of time steps
Ntsteps = ceil(FinalTime/dt); dt = FinalTime/Ntsteps;

time = 0; tstep = 1;

% initialize the vector of unknowns U. The ordering of the elements has to
% be consistent with the basis functions and the Matrices
U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

% initialize the needed matrices
InitMatLawson;

P = zeros(Np*K,1);
for k = 1:K
   if ismember(k, fine_idx)
      P((k-1)*Np+1 : k * Np) = ones(Np, 1); 
   end
end
P2 = ones(Np*K, 1) - P;
P = [P;P;P;P;P;P];
P2 = [P2;P2;P2;P2;P2;P2];
P = diag(P);
P2 = diag(P2);
Cfine = Cmat * P;
Ccoarse = Cmat * P2;

% compute the constants in the Lawson-LRSK scheme
rk4cLawson = [0 rk4c];   % We need c(0). What should the value be?
matexp1 = expm((1-rk4cLawson(6))*dt*Cfine);

matexp2 = cell(1,5);
c1 = cell(1,5);
c2 = cell(1,5);
for k = 1:5
    matexp2{k} = expm((rk4cLawson(k+1)-rk4cLawson(k))*dt*Cfine);
    c1{k} = rk4a(k) * matexp2{k};
    c2{k} = dt * Ccoarse * matexp2{k};
end


% outer time step loop 
while (time<FinalTime)
  
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
   Ez_time = [Ez_time, U(2*Np*K+1)];
end

% convert U back to field components
[Hx,Hy,Hz,Ex,Ey,Ez] = UToFields(U);
return;