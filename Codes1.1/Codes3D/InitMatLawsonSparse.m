% Initialize the Matrices needed for the Lawson scheme using sparse
% matrices
blksize = 3 * Np;
Ccoarse = spalloc(2 * blksize * K, 2 * blksize * K, blksize^2 * K * 2); % How much memory do I need?
Cfine = spalloc(2 * blksize * K, 2 * blksize * K, blksize^2 * K * 2); % How much memory do I need?

% Upper Right Block
invM = inv(MassMatrix);
for i = 1:K
    % Stiffness
    Dx = rx(1,i) * Dr + sx(1,i) * Ds + tx(1,i) * Dt;
    Dy = ry(1,i) * Dr + sy(1,i) * Ds + ty(1,i) * Dt;
    Dz = rz(1,i) * Dr + sz(1,i) * Ds + tz(1,i) * Dt;
    Sx = invM * Dx;
    Sy = invM * Dy;
    Sz = invM * Dz;
    S = [zeros(Np), -Sz, Sy; Sz, zeros(Np), -Sx; -Sy, Sx, zeros(Np)];
    
    % Mass
    invMi = inv(J(1, i) * MassMatrix);
    invMi = blkdiag(invMi, invMi, invMi);
    
    % Diagonal entry
    S = invMi * (S - SurfaceMassLawson(r,s,t,i));
    
    % fill the row
    for j = 1:K
        % Check whether to write into Cfine or Ccoarse
        if ismember(j, fine_idx)
            if i == j
                Cfine((i-1)*blksize+1:i*blksize, (j-1)*blksize+1+K*blksize:j*blksize+K*blksize) = S;
            elseif ismember(i, EToE(j,:))
                Cfine((i-1)*blksize+1:i*blksize, (j-1)*blksize+1+K*blksize:j*blksize+K*blksize) = -invMi * S_ikPlusLawson(i,j,r,s,t);
            end
        else
            if i == j
                Ccoarse((i-1)*blksize+1:i*blksize, (j-1)*blksize+1+K*blksize:j*blksize+K*blksize) = S;
            elseif ismember(i, EToE(j,:))
                Ccoarse((i-1)*blksize+1:i*blksize, (j-1)*blksize+1+K*blksize:j*blksize+K*blksize) = -invMi * S_ikPlusLawson(i,j,r,s,t);
            end
        end
    end
end

% Lower Left Block (=Upper Right because eps=my=1)
Ccoarse(K*blksize+1:end, 1:K*blksize) = -Ccoarse(1:K*blksize, K*blksize+1:end);
Cfine(K*blksize+1:end, 1:K*blksize) = -Cfine(1:K*blksize, K*blksize+1:end);




