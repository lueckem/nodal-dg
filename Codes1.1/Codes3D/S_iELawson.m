function [S_iE] = S_iELawson(i, R, S, T)
%Calculate the Matrix SUM_k S_ik^E

Globals3D;

S_iE = zeros(3*Np);
Fcell = cell(6,1);
 
for face=1:Nfaces
  % check if outer face
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
  
  for l=1:Nfp
     for j=1:Nfp
         F(Fmask(l,face), Fmask(j,face)) = 0.5 * massFace(l, j); %wrong?
     end
  end
  
  % build S_ikE
  factor = dot([n_y^2 + n_z^2; -n_x*n_y; -n_x*n_z], [n_y^2 + n_z^2; -n_x*n_y; -n_x*n_z]);
  Fcell{1} = factor .* F;
  
  factor = dot([-n_y*n_x; n_x^2 + n_z^2; -n_y*n_z], [-n_y*n_x; n_x^2 + n_z^2; -n_y*n_z]);
  Fcell{2} = factor .* F;
  
  factor = dot([-n_z*n_x; -n_z*n_y; n_x^2 + n_y^2], [-n_z*n_x; -n_z*n_y; n_x^2 + n_y^2]);
  Fcell{3} = factor .* F;
  
  Fcell{4} = -n_y*n_x .* F;
  Fcell{5} = -n_x*n_z .* F;
  Fcell{6} = -n_y*n_z .* F;
  
  S_iE = S_iE + [Fcell{1},Fcell{4},Fcell{5}; Fcell{4},Fcell{2},Fcell{6}; Fcell{5},Fcell{6},Fcell{3}];
end
end

