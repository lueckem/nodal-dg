function [expcfinev] = ExpCfinev(alpha, v)
% After the model reduction, compute exp(alpha * C_fine)*v using Krylov
% subspace methods

Globals3D;

tol = 0.5;
m = 4;

expcfinev = [expv(alpha, Cff, v(1:size_Cff), tol, m);
    Ccf * phiv(alpha, Cff, v(1:size_Cff), zeros(size_Cff, 1), tol, m) + v(size_Cff+1:end)];
end