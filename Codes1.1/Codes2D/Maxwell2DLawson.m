function [Hx,Hy,Ez,time] = Maxwell2DLawson(Hx, Hy, Ez, FinalTime)

% function [Hx,Hy,Ez] = Maxwell2D(Hx, Hy, Ez, FinalTime)
% Purpose :Integrate TM-mode Maxwell's until FinalTime starting with initial conditions Hx,Hy,Ez

Globals2D;
time = 0;

% compute time step size
rLGL = JacobiGQ(0,0,N); rmin = abs(rLGL(1)-rLGL(2));
dtscale = dtscale2D; dt = min(dtscale)*rmin*2/3

% initialize the vector of unknowns U. The ordering of the elements has to
% be consistent with the basis functions and the Matrices
U = [];

% initialize the needed matrices
Cfine = [];
Ccoarse = [];

% compute the constants in the Lawson-LRSK scheme
rk4cLawson = [0 rk4c];   % We need c[0]. What should the value be?
matexp1 = expm((1-rk4cLawson(6))*dt*Cfine);

c1 = zeros(1, 5); % have to use cells !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c2 = zeros(1, 5);
for k = 1:5
    matexp2 = expm((rk4cLawson(k+1)-rk4cLawson(k))*dt*Cfine);
    c1(k) = rk4a(k) * matexp2;
    c2(k) = dt * Ccoarse * matexp2;
end


% outer time step loop 
while (time<FinalTime)
  
  if(time+dt>FinalTime), dt = FinalTime-time; end
  
  phi1 = U;
  phi2 = zeros(size(U));
  for k = 1:5
      phi2 = c1(k) * phi2 + c2(k) * phi1;
      phi1 = %%%%%%
       
  end
   
   % Increment time
   time = time+dt;
end
return