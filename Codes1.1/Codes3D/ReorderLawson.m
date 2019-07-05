% Reorder Cfine and Ccoarse and U for the model reduction

coarse_idx = 1:K;
coarse_idx(fine_idx) = [];

fine_idx_long = zeros(1,length(fine_idx)*3*Np);
for j = 1:length(fine_idx)
   fine_idx_long((j-1)*3*Np+1:j*3*Np) = (fine_idx(j)-1)*3*Np+1:fine_idx(j)*3*Np; 
end

coarse_idx_long = 1:3*Np*K;
coarse_idx_long(fine_idx_long) = [];

C_idx = [fine_idx_long, (fine_idx_long + K*3*Np), coarse_idx_long, (coarse_idx_long + K*3*Np)];

Cfine = Cfine(:,C_idx);
Ccoarse = Ccoarse(:,C_idx);