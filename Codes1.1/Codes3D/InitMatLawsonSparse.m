% Initialize the Matrices needed for the Lawson scheme using sparse
% matrices

blksize = 3 * Np;

% allocate memory
if isempty(fine_idx) == 0
    frac_fine = 0; % fraction of fine elements
else
    frac_fine = length(fine_idx)/K;
end
Ccoarse = spalloc(2 * blksize * K, 2 * blksize * K, round((1-frac_fine) * 78*K*Np^2));
Cfine = spalloc(2 * blksize * K, 2 * blksize * K, round(frac_fine * 78*K*Np^2));

% Inverse of MassMatrix
invM = inv(MassMatrix);

% Upper Right Block
for i = 1:K
    % Stiffness
    Dx = rx(1,i) * Dr + sx(1,i) * Ds + tx(1,i) * Dt;
    Dy = ry(1,i) * Dr + sy(1,i) * Ds + ty(1,i) * Dt;
    Dz = rz(1,i) * Dr + sz(1,i) * Ds + tz(1,i) * Dt;
    S = [zeros(Np), -Dz, Dy; Dz, zeros(Np), -Dx; -Dy, Dx, zeros(Np)];
    
    % Mass
    invMi = invM./J(1, i);
    invMi = blkdiag(invMi, invMi, invMi);
    
    % Diagonal entry
    S = S - invMi * SurfaceMassLawson(r,s,t,i);
    
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

% Lower Left Block (=-Upper Right because eps=my=1)
Ccoarse(K*blksize+1:end, 1:K*blksize) = -Ccoarse(1:K*blksize, K*blksize+1:end);
Cfine(K*blksize+1:end, 1:K*blksize) = -Cfine(1:K*blksize, K*blksize+1:end);

%Upper Left Block
for i=1:K
    invMi = invM./J(1, i);
    invMi = blkdiag(invMi, invMi, invMi);
    
    if ismember(i, fine_idx)
        Cfine((i-1)*blksize+1:i*blksize, (i-1)*blksize+1:i*blksize) = invMi * S_iELawson(i,r,s,t);
    else
        Ccoarse((i-1)*blksize+1:i*blksize, (i-1)*blksize+1:i*blksize) = invMi * S_iELawson(i,r,s,t);
    end
end

%Lower Right Block (=Upper Left because eps=my=1)
Ccoarse(K*blksize+1:end, K*blksize+1:end) = Ccoarse(1:K*blksize, 1:K*blksize);
Cfine(K*blksize+1:end, K*blksize+1:end) = Cfine(1:K*blksize, 1:K*blksize);