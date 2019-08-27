function [expcfinev] = ExpCfinev(alpha, v)
% After the model reduction, compute exp(alpha * C_fine)*v using Krylov
% subspace methods

Globals3D;

expcfinev = [expv(alpha, Cff, v(1:size_Cff));
    Ccf * phiv(alpha, Cff, v(1:size_Cff), zeros(size_Cff, 1)) + v(size_Cff+1:end)];
end