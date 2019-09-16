Globals3D;
N = 5;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK5.neu');
StartUp3D;

%sigmax = [1, 1, 1, 1, 1];
sigmax = zeros(1, K);
sigmay = sigmax; sigmaz = sigmax;

epsilon = ones(K,1);
%epsilon(5) = 3;

% initial cond
U_0 = rand(6*Np*K,1);
P_0 = zeros(6*Np*K,1);
dt = dtscale3D;

%% Maxwell3DMat
fine_idx = [];
InitMatLawsonSparse;
[sigma1, sigma2, sigma3] = constructSigma();

U = U_0;
P = P_0;
resU = zeros(6*Np*K,1);
resP = zeros(6*Np*K,1);
for steps = 1:100
for INTRK = 1:5
    resU = rk4a(INTRK)*resU + dt*Ccoarse*U;
    %resU = resU + dt * (blkmult(sigma1, U, Np, K) - P); % add PML terms
    U = U + rk4b(INTRK)*resU;
    
    % integrate PML terms
    %resP = rk4a(INTRK)*resP + dt * (blkmult(sigma2,P,Np,K) + blkmult(sigma3,U,Np,K));
    %P = P + rk4b(INTRK) * resP;
end
end
U_mat = U;
P_mat = P;

%% Krylov
fine_idx = 5;
U = U_0;
P = P_0;
InitMatLawsonSparse;
ReorderLawson;
P = P(C_idx);

rk4cLawson = [0 rk4c];
size_Cff = 2*3*Np*length(fine_idx);
Cff = Cfine(1:size_Cff,1:size_Cff);
Ccf = Cfine(size_Cff+1:2*3*Np*K,1:size_Cff);
[sigma1, sigma2, sigma3] = constructSigma();

phi1 = U;
phi2 = zeros(6*Np*K,1);
resP = zeros(6*Np*K,1);
for steps = 1:100
for k = 1:5
      phi1star = ExpCfinev((rk4cLawson(k+1)-rk4cLawson(k))*dt, phi1);
      phi2 = rk4a(k) * ExpCfinev((rk4cLawson(k+1)-rk4cLawson(k))*dt, phi2) + dt * Ccoarse * phi1star;
      %phi2 = phi2 + dt * (blkmult(sigma1, phi1star, Np, K) - P); % add PML terms
      phi1 = phi1star + rk4b(k) * phi2;
      
      % integrate P
      %resP = rk4a(k)*resP + dt * (blkmult(sigma2,P,Np,K) + blkmult(sigma3,phi1star,Np,K));
      %P = P + rk4b(k) * resP;  
end
    phi1 = ExpCfinev((1-rk4cLawson(6))*dt, phi1);
end
  U_krylov = phi1;
  U_krylov = ReorderBackULawson(U_krylov);
  P_krylov = ReorderBackULawson(P);
  
  %% comparison
  norm(U_mat - U_krylov)
  norm(P_mat - P_krylov)