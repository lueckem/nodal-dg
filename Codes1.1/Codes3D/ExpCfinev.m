function [expcfinev] = ExpCfinev(alpha, v)
% After the model reduction, compute exp(alpha * C_fine)*v using Krylov
% subspace methods

Globals3D;

if v == 0
   expcfinev = sparse(6*Np*K,1);
   return
end

size_Cff = 2*3*Np*length(fine_idx);
Cff = Cfine(1:size_Cff,1:size_Cff);
Ccf = Cfine(size_Cff+1:2*3*Np*K,1:size_Cff);

expcfinev = [expv(alpha, Cff, v(1:size_Cff));
    Ccf * phiv(alpha, Cff, v(1:size_Cff), zeros(size_Cff, 1)) + v(size_Cff+1:end)];
end