function [S_iE] = S_iELawson(i, R, S, T)
%Calculate the Matrix SUM_k S_ik^E

Globals3D;

S_iE = zeros(3*Np);
 
for face=1:Nfaces
  % check if not outer face
  if EToE(i,face) ~= i
      continue
  end
  
  % process face
  if(face==1); faceR = R(Fmask(:,1)); faceS = S(Fmask(:,1)); end
  if(face==2); faceR = R(Fmask(:,2)); faceS = T(Fmask(:,2)); end
  if(face==3); faceR = S(Fmask(:,3)); faceS = T(Fmask(:,3)); end
  if(face==4); faceR = S(Fmask(:,4)); faceS = T(Fmask(:,4)); end
  
  % surface Jacobian
  sJacobian = sJ(Nfp*(face-1)+1,i);
  
  % normal vector
  n_x = nx(Nfp*(face-1)+1,i);
  n_y = ny(Nfp*(face-1)+1,i);
  n_z = nz(Nfp*(face-1)+1,i);
  
  % build F
  F = zeros(Np);
  VFace = Vandermonde2D(N, faceR, faceS);
  massFace = sJacobian .* inv(VFace*VFace');
  
  for j=1:Nfp
     for l=1:Nfp
         F(Fmask(j,face), Fmask(l,face)) = 0.5 * massFace(j, l);
     end
  end
  
  % build S_ikE
  S_iE = S_iE + [(n_y^2+n_z^2).*F, -n_y*n_x.*F, -n_x*n_z.*F;
                 -n_y*n_x.*F, (n_x^2+n_z^2).*F, -n_y*n_z.*F;
                 -n_x*n_z.*F, -n_y*n_z.*F, (n_x^2+n_y^2).*F];
end
end

