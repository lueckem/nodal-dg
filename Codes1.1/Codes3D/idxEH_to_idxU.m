function [idxU] = idxEH_to_idxU(field, idxEH)
%convert idxEH [i,k] to idxU
%field is 1,2,3,4,5,6 for Ex,Ey,Ez,Hx,Hy,Hz

Globals3D;

idxU = (field-1)*K*Np + (idxEH(2)-1)*Np + idxEH(1);
end

