function [Hx, Hy, Hz, Ex, Ey, Ez] = UToFields(U)

Globals3D;
blksize = K*3*Np;

Ex = zeros(Np, K);
Ey = zeros(Np, K);
Ez = zeros(Np, K);
Hx = zeros(Np, K);
Hy = zeros(Np, K);
Hz = zeros(Np, K);

% E field
for k = 1:K
    Ex(:,k) = U((k-1)*3*Np+1 : (k-1)*3*Np+Np);
    Ey(:,k) = U((k-1)*3*Np+Np+1 : (k-1)*3*Np+2*Np);
    Ez(:,k) = U((k-1)*3*Np+2*Np+1 : k*3*Np);
end

% H field
for k = 1:K
    Hx(:,k) = U((k-1)*3*Np+1+blksize : (k-1)*3*Np+Np+blksize);
    Hy(:,k) = U((k-1)*3*Np+Np+1+blksize : (k-1)*3*Np+2*Np+blksize);
    Hz(:,k) = U((k-1)*3*Np+2*Np+1+blksize : k*3*Np+blksize);
end

end

