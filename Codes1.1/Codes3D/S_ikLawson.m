function [S_ik] = S_ikLawson(R,S,T,i,k)
% Calculate the matrix S_ik

Globals3D;

S_ik= zeros(3*Np);

% check if element i and k are neighbours
[bool, face_i] = ismember(k, EToE(i,:));
if ~bool
   return; 
end

% surface Jacobian
sJacobian = sJ((face_i-1)*Nfp+1, i);

% normal
normal = [nx(Nfp*(face_i-1)+1,i);
           ny(Nfp*(face_i-1)+1,i);
           nz(Nfp*(face_i-1)+1,i)];

  % process face
  if(face_i==1); faceR = R(Fmask(:,1)); faceS = S(Fmask(:,1)); end
  if(face_i==2); faceR = R(Fmask(:,2)); faceS = T(Fmask(:,2)); end
  if(face_i==3); faceR = S(Fmask(:,3)); faceS = T(Fmask(:,3)); end
  if(face_i==4); faceR = S(Fmask(:,4)); faceS = T(Fmask(:,4)); end
  
  F = zeros(Np);
  VFace = Vandermonde2D(N, faceR, faceS);
  massFace = sJacobian .* inv(VFace*VFace');
  
  for j=1:Nfp
     for l=1:Nfp
         F(Fmask(j,face_i), Fmask(l,face_i)) = 0.5 .* massFace(j, l);
     end
  end
 
  % insert F into SurfaceMassMatrix
  S_ik(1:Np,Np+1:2*Np) = -normal(3).*F;
  S_ik(1:Np,2*Np+1:3*Np) = normal(2).*F;
  S_ik(Np+1:2*Np,1:Np) = normal(3).*F;
  S_ik(Np+1:2*Np,2*Np+1:3*Np) = -normal(1).*F;
  S_ik(2*Np+1:3*Np,1:Np) = -normal(2).*F;
  S_ik(2*Np+1:3*Np,Np+1:2*Np) = normal(1).*F;
end

