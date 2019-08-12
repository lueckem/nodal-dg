% Initialize the Matrices needed for the Lawson scheme using sparse
% matrices

blksize = 3 * Np;

% allocate memory
frac_fine = length(fine_idx)/K;
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
    
    % Mass matrix
    invMi = invM./J(1, i);
    invMi = blkdiag(invMi, invMi, invMi);
    
    % Diagonal entry
    S = S - invMi * SurfaceMassInteriorLawson(r,s,t,i);
    if ismember(i, fine_idx)
        Cfine((i-1)*blksize+1:i*blksize, (i-1)*blksize+1+K*blksize:i*blksize+K*blksize) = S;
    else
        Ccoarse((i-1)*blksize+1:i*blksize, (i-1)*blksize+1+K*blksize:i*blksize+K*blksize) = S;
    end
    
    % Add the S_ikPlus
    for j = 1:4
        k = EToE(i,j); % k is the index of the neighbour element
        if k == i % check for boundary
            continue
        end
        
        if ismember(k, fine_idx)
            Cfine((i-1)*blksize+1:i*blksize, (k-1)*blksize+1+K*blksize:k*blksize+K*blksize) = -invMi * S_ikPlusLawson(i,k,r,s,t);
        else
            Ccoarse((i-1)*blksize+1:i*blksize, (k-1)*blksize+1+K*blksize:k*blksize+K*blksize) = -invMi * S_ikPlusLawson(i,k,r,s,t);
        end
    end
end

% Lower Left Block (=-Upper Right because eps=my=1)
Ccoarse(K*blksize+1:end, 1:K*blksize) = -Ccoarse(1:K*blksize, K*blksize+1:end);
Cfine(K*blksize+1:end, 1:K*blksize) = -Cfine(1:K*blksize, K*blksize+1:end);

% Boundary condition not yet working

% %Upper Left Block
% for i=1:K
%     invMi = invM./J(1, i);
%     invMi = blkdiag(invMi, invMi, invMi);
%     
%     if ismember(i, fine_idx)
%         Cfine((i-1)*blksize+1:i*blksize, (i-1)*blksize+1:i*blksize) = invMi * S_iELawson(i,r,s,t);
%     else
%         Ccoarse((i-1)*blksize+1:i*blksize, (i-1)*blksize+1:i*blksize) = invMi * S_iELawson(i,r,s,t);
%     end
% end

% %Lower Right Block (=Upper Left because eps=my=1)
% Ccoarse(K*blksize+1:end, K*blksize+1:end) = Ccoarse(1:K*blksize, 1:K*blksize);
% Cfine(K*blksize+1:end, K*blksize+1:end) = Cfine(1:K*blksize, 1:K*blksize);