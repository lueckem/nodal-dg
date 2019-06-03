% Initialize the Matrices needed for the Lawson scheme

% Mass Matrix M
globalMassMatrix = zeros(3 * Np * K);

for k = 1 : K
    Mk = J(1, k) * MassMatrix;
    Mk = blkdiag(Mk, Mk, Mk);
    globalMassMatrix((k-1)*3*Np+1 : k*3*Np, (k-1)*3*Np+1 : k*3*Np) = Mk;
end


%Stiffness Matrix K
globalStiffnessMatrix = zeros(3 * Np * K);
invM = inv(MassMatrix);

for k = 1:K
    Dx = rx(1,k) * Dr + sx(1,k) * Ds + tx(1,k) * Dt;
    Dy = ry(1,k) * Dr + sy(1,k) * Ds + ty(1,k) * Dt;
    Dz = rz(1,k) * Dr + sz(1,k) * Ds + tz(1,k) * Dt;
    
    Sx = invM * Dx;
    Sy = invM * Dy;
    Sz = invM * Dz;
    
    Sk = blkdiag(Sx, Sy, Sz);
    globalStiffnessMatrix((k-1)*3*Np+1 : k*3*Np, (k-1)*3*Np+1 : k*3*Np) = Sk - SurfaceMassLawson(r,s,t,k);
end
for i = 1:K
    for k = 1:K
        if i~=k
            globalStiffnessMatrix((i-1)*3*Np+1 : i*3*Np, (k-1)*3*Np+1 : k*3*Np) = S_ikPlusLawson(i,k,r,s,t);
        end
    end
end

% epsilon = my = 1
A = blkdiag(globalMassMatrix, globalMassMatrix);
clear globalMassMatrix
B = [zeros(3 * Np * K), globalStiffnessMatrix; -globalStiffnessMatrix, zeros(3 * Np * K)];
clear globalStiffnessMatrix
Cmat = A * B;
clear A B

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

clear Cmat P P2