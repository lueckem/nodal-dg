function [U] = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez)

Globals3D;

U = zeros(2*3*K*Np, 1);
blksize = K*3*Np;

% E field
for k = 1:K
    U((k-1)*3*Np+1 : k*3*Np) = [Ex(:,k); Ey(:,k); Ez(:,k)];
end

% H field
for k = 1:K
    U((k-1)*3*Np+1+blksize : k*3*Np+blksize) = [Hx(:,k); Hy(:,k); Hz(:,k)];
end

end