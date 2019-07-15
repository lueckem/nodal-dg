function [expcfine] = ExpCfine(alpha)
% After the model reduction, compute exp(alpha * C_fine)

Globals3D;

size_Cff = 2*3*Np*length(fine_idx);
size_Ccc = 2*3*Np*K-size_Cff;
invCff = Cfine(1:size_Cff,1:size_Cff);
expCff = sparse(expm(alpha .* invCff));
invCff = inv(invCff);

expCcf = Cfine(size_Cff+1:2*3*Np*K,1:size_Cff)*(expCff-speye(size_Cff))*invCff;

expcfine = [expCff, sparse(size_Cff,size_Ccc); expCcf, speye(size_Ccc)];

end

