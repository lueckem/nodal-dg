function [SurfaceMassMatrix] = SurfaceMassLawson(N,R,S,T,k)
% calculates Sum(S_ik) for the k-th Element

Globals3D;

SurfaceMassMatrix = zeros(3*Np);

% surface Jacobians for the 4 faces
sJacobian = [sJ(1,k), sJ(Nfp+1,k), sJ(2*Nfp+1, k), sJ(3*Nfp+1,k)];

%normals for the 4 faces
normals = [nx(1,k), nx(Nfp+1,k), nx(2*Nfp+1, k), nx(3*Nfp+1,k);
           ny(1,k), ny(Nfp+1,k), ny(2*Nfp+1, k), ny(3*Nfp+1,k);
           nz(1,k), nz(Nfp+1,k), nz(2*Nfp+1, k), nz(3*Nfp+1,k)];

for face=1:Nfaces
  % process face
  if(face==1); faceR = R(Fmask(:,1)); faceS = S(Fmask(:,1)); end
  if(face==2); faceR = R(Fmask(:,2)); faceS = T(Fmask(:,2)); end
  if(face==3); faceR = S(Fmask(:,3)); faceS = T(Fmask(:,3)); end
  if(face==4); faceR = S(Fmask(:,4)); faceS = T(Fmask(:,4)); end
  
  F = zeros(Np);
  VFace = Vandermonde2D(N, faceR, faceS);
  massFace = sJacobian(face) .* inv(VFace*VFace');
  
  for i=1:Nfp
     for j=1:Nfp
         F(Fmask(i,face), Fmask(j,face)) = 0.5 * massFace(i, j);
     end
  end
  
  normal = normals(:,face);
  
  % insert F into SurfaceMassMatrix
  SurfaceMassMatrix(1:Np,Np+1:2*Np) = SurfaceMassMatrix(1:Np,Np+1:2*Np) - normal(3).*F;
  SurfaceMassMatrix(1:Np,2*Np+1:3*Np) = SurfaceMassMatrix(1:Np,2*Np+1:3*Np) + normal(2).*F;
  SurfaceMassMatrix(Np+1:2*Np,1:Np) = SurfaceMassMatrix(Np+1:2*Np,1:Np) + normal(3).*F;
  SurfaceMassMatrix(Np+1:2*Np,2*Np+1:3*Np) = SurfaceMassMatrix(Np+1:2*Np,2*Np+1:3*Np) - normal(1).*F;
  SurfaceMassMatrix(2*Np+1:3*Np,1:Np) = SurfaceMassMatrix(2*Np+1:3*Np,1:Np) - normal(2).*F;
  SurfaceMassMatrix(2*Np+1:3*Np,Np+1:2*Np) = SurfaceMassMatrix(2*Np+1:3*Np,Np+1:2*Np) + normal(1).*F;
  
end
end
