function [U] = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez)

Globals3D;

U = zeros(2*3*K*Np, 1);
blksize = K*Np;

%Ex
for k = 1:K
   U((k-1)*Np+1 : k * Np) = Ex(:,k);
end

%Ey
for k = 1:K
   U(blksize + (k-1)*Np+1 : blksize + k * Np) = Ey(:,k);
end

%Ez
for k = 1:K
   U(2*blksize + (k-1)*Np+1 : 2*blksize + k * Np) = Ez(:,k);
end


%Hx
for k = 1:K
   U(3*blksize + (k-1)*Np+1 : 3*blksize + k * Np) = Hx(:,k);
end

%Hy
for k = 1:K
   U(4*blksize + (k-1)*Np+1 : 4*blksize + k * Np) = Hy(:,k);
end

%Hz
for k = 1:K
   U(5*blksize + (k-1)*Np+1 : 5*blksize + k * Np) = Hz(:,k);
end
end