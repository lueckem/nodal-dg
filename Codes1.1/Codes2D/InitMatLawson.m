% Initialize the Matrices needed for the Lawson scheme

% Mass Matrix
globalMassMatrix = zeros(3 * Np * K);
globalStiffnessMatrix = zeros(3 * Np * K);

for k = 1 : K
    Mk = J(1, k) * MassMatrix;  % p. 183
    Sr = inv(Mk);
    
    Mk = blkdiag(Mk, Mk, Mk);
    globalMassMatrix((k-1)*3*Np+1 : k*3*Np, (k-1)*3*Np+1 : k*3*Np) = Mk;
    
    Ss = Sr * Ds;   % we need d/dx and not d/dr ?
    Sr = Sr * Dr;
    Sk = blkdiag(Sr, Ss, zeros(Np));
    globalStiffnessMatrix((k-1)*3*Np+1 : k*3*Np, (k-1)*3*Np+1 : k*3*Np) = Sk;
    
    % stiffness not finished!
end