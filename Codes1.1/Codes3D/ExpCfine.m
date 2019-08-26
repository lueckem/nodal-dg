function [expcfine] = ExpCfine(alpha)
% After the model reduction, compute exp(alpha * C_fine)

Globals3D;

size_Cff = 2*3*Np*length(fine_idx);
size_Ccc = 2*3*Np*K-size_Cff;
Cff = Cfine(1:size_Cff,1:size_Cff);
expCff = sparse(expm(alpha .* Cff));

%expCcf = (Cfine(size_Cff+1:2*3*Np*K,1:size_Cff)*(expCff-speye(size_Cff)))/Cff; % Cff singular!!!!!!

% % % % % Manual calculation. Very inefficient and probably inaccurate!!!
expCcf = sparse(size_Cff, size_Cff);
temp = speye(size_Cff, size_Cff);
for n = 1:100
   expCcf = expCcf + alpha^n .* temp ./ factorial(n);
   temp = temp * Cff;
end
expCcf = Cfine(size_Cff+1:2*3*Np*K,1:size_Cff) * expCcf;

expcfine = [expCff, sparse(size_Cff,size_Ccc); expCcf, speye(size_Ccc)];
end

