function [Hx, Hy, Hz, Ex, Ey, Ez] = UToFields(U)

Globals3D;
blksize = K*Np;

%Ex
for k = 1:K
   Ex(:,k) = U((k-1)*Np+1 : k * Np);
end

%Ey
for k = 1:K
    Ey(:,k) = U(blksize + (k-1)*Np+1 : blksize + k * Np);
end

%Ez
for k = 1:K
    Ez(:,k) = U(2*blksize + (k-1)*Np+1 : 2*blksize + k * Np);
end


%Hx
for k = 1:K
   Hx(:,k) = U(3*blksize + (k-1)*Np+1 : 3*blksize + k * Np);
end

%Hy
for k = 1:K
   Hy(:,k) = U(4*blksize + (k-1)*Np+1 : 4*blksize + k * Np);
end

%Hz
for k = 1:K
   Hz(:,k) = U(5*blksize + (k-1)*Np+1 : 5*blksize + k * Np);

end

