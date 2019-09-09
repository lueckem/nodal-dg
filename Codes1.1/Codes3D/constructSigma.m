function [sigma1, sigma2, sigma3] = constructSigma()
% construct column vectors sigma1, sigma2, sigma3 for PML terms.
% rhsU ~ sigma1 .* U.
% rhsP ~ sigma2 .* P + sigma3 .* U.
% sigma is also reordered.
% to save memory, the sigmas are stored as small as possible (3*K) and to
% evaluate the product sigma1 .* U the function blkmult is used


Globals3D;

%sigma1
sigma1 = zeros(3*K, 1);
for i = 1:K
   sigma1(3*(i-1)+1:3*i) = -[-sigmax(i)+sigmay(i)+sigmaz(i);
                            sigmax(i)-sigmay(i)+sigmaz(i);
                            sigmax(i)+sigmay(i)-sigmaz(i)];
end

%sigma2
sigma2 = zeros(3*K, 1);
for i = 1:K
   sigma2(3*(i-1)+1:3*i) = [-sigmax(i);-sigmay(i);-sigmaz(i)];
end

%sigma3
sigma3 = zeros(3*K, 1);
for i = 1:K
   sigma3(3*(i-1)+1:3*i) = [sigmax(i)*sigmax(i) + sigmay(i)*sigmaz(i) - sigmay(i)*sigmax(i) - sigmaz(i)*sigmax(i);
                            sigmay(i)*sigmay(i) + sigmax(i)*sigmaz(i) - sigmaz(i)*sigmay(i) - sigmay(i)*sigmax(i);
                            sigmaz(i)*sigmaz(i) + sigmay(i)*sigmax(i) - sigmaz(i)*sigmax(i) - sigmay(i)*sigmaz(i)];
end


% Reordering
fine_idx_long = zeros(1,length(fine_idx)*3);
for j = 1:length(fine_idx)
   fine_idx_long((j-1)*3+1:j*3) = (fine_idx(j)-1)*3+1:fine_idx(j)*3; 
end

coarse_idx_long = 1:3*K;
coarse_idx_long(fine_idx_long) = [];

perm_idx = [fine_idx_long, coarse_idx_long];

sigma1=sigma1(perm_idx);
sigma2=sigma2(perm_idx);
sigma3=sigma3(perm_idx);
end


