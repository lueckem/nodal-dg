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
for k = 1:K
    Dx = rx(1,k) * Dr + sx(1,k) * Ds;
    Dy = ry(1,k) * Dr + sy(1,k) * Ds;
    %Dz = ...
    
    invMk = inv(J(1, k) * MassMatrix);
    Sx = invMk * Dx;
    Sy = invMk * Dy;
    Sz = zeros(Np);
    %Sz = invMk * Dz;
    
    Sk = blkdiag(Sx, Sy, Sz);
    globalStiffnessMatrix((k-1)*3*Np+1 : k*3*Np, (k-1)*3*Np+1 : k*3*Np) = Sk - SurfaceMassLawson(N,r,s,t,k);
end
for i = 1:K
    for k = 1:K
        if i~=k
            globalStiffnessMatrix((i-1)*3*Np+1 : i*3*Np, (k-1)*3*Np+1 : k*3*Np) = S_ikPlusLawson(i,k,N,r,s,t);
        end
    end
end