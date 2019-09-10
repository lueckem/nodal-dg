% Initialize the Matrices needed for the Lawson scheme using sparse
% matrices.

% PEC boundary conditions are implemented.

blksize = 3 * Np;

% allocate memory
frac_fine = length(fine_idx)/K;
Ccoarse = spalloc(2 * blksize * K, 2 * blksize * K, round((1-frac_fine) * 78*K*Np^2));
Cfine = spalloc(2 * blksize * K, 2 * blksize * K, round(frac_fine * 78*K*Np^2));

% Inverse of MassMatrix
invM = inv(MassMatrix);

% filling in Cfine and Ccoarse
for i = 1:K
    % Stiffness
    Dx = rx(1,i) * Dr + sx(1,i) * Ds + tx(1,i) * Dt;
    Dy = ry(1,i) * Dr + sy(1,i) * Ds + ty(1,i) * Dt;
    Dz = rz(1,i) * Dr + sz(1,i) * Ds + tz(1,i) * Dt;
    S = [zeros(Np), -Dz, Dy; Dz, zeros(Np), -Dx; -Dy, Dx, zeros(Np)];
    
    % Mass matrix
    invMi = invM./J(1, i);
    invMi = blkdiag(invMi, invMi, invMi);
    
    % Diagonal entry for upper right block
    S_E = S - (invMi./epsilon(i)) * SurfaceMassInteriorLawson(r,s,t,i);
    
    % Diagonal entry for lower left block
    S_H = -(S - invMi * SurfaceMassLawson(r,s,t,i));
    
    if ismember(i, fine_idx)
        Cfine((i-1)*blksize+1:i*blksize, (i-1)*blksize+1+K*blksize:i*blksize+K*blksize) = S_E; %upper right
        Cfine((i-1)*blksize+1+K*blksize:i*blksize+K*blksize, (i-1)*blksize+1:i*blksize) = S_H; %lower left
    else
        Ccoarse((i-1)*blksize+1:i*blksize, (i-1)*blksize+1+K*blksize:i*blksize+K*blksize) = S_E; %upper right
        Ccoarse((i-1)*blksize+1+K*blksize:i*blksize+K*blksize, (i-1)*blksize+1:i*blksize) = S_H; %lower left
    end
    
    % Add the S_ikPlus
    for j = 1:4
        k = EToE(i,j); % k is the index of the neighbour element
        if k == i % check for boundary
            continue
        end
        
        if ismember(k, fine_idx)
            Cfine((i-1)*blksize+1:i*blksize, (k-1)*blksize+1+K*blksize:k*blksize+K*blksize) = -(invMi./epsilon(i)) * S_ikPlusLawson(i,k,r,s,t); %upper right
            Cfine((i-1)*blksize+1+K*blksize:i*blksize+K*blksize, (k-1)*blksize+1:k*blksize) = invMi * S_ikPlusLawson(i,k,r,s,t); %lower left
        else
            Ccoarse((i-1)*blksize+1:i*blksize, (k-1)*blksize+1+K*blksize:k*blksize+K*blksize) = -(invMi./epsilon(i)) * S_ikPlusLawson(i,k,r,s,t); %upper right
            Ccoarse((i-1)*blksize+1+K*blksize:i*blksize+K*blksize, (k-1)*blksize+1:k*blksize) = invMi * S_ikPlusLawson(i,k,r,s,t); %lower left
        end
    end
end

% Boundary condition as proposed in the paper not working.
% The Boundary condition is implemented in SurfaceMassLawson.

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
% 
% %Lower Right Block (=Upper Left because eps=my=1)
% Ccoarse(K*blksize+1:end, K*blksize+1:end) = Ccoarse(1:K*blksize, 1:K*blksize);
% Cfine(K*blksize+1:end, K*blksize+1:end) = Cfine(1:K*blksize, 1:K*blksize);