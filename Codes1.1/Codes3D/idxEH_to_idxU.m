function [idxU] = idxEH_to_idxU(field, idxEH)
%convert idxEH [i,k] to idxU
%field is 1,2,3,4,5,6 for Ex,Ey,Ez,Hx,Hy,Hz

Globals3D;
blksize = K*3*Np;

if field < 4
    idxU = (idxEH(:,2)-1)*3*Np + (field-1)*Np + idxEH(:,1);
else
    idxU = (idxEH(:,2)-1)*3*Np+blksize + (field-4)*Np + idxEH(:,1);
end
end

