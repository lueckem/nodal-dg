function [SurfaceMassMatrix] = SurfaceMassLawson(N,R,S,T)
% calculates Sum(S_ik) for the reference Element

Globals3D;

SurfaceMassMatrix = zeros(3*Np);

n1 = [0;0;-1];
n2 = [0;-1;0];
n3 = [1;1;1]./sqrt(3);
n4 = [-1;0;0];
normals = [n1, n2, n3, n4];

for face=1:Nfaces
  % process face
  if(face==1); faceR = R(Fmask(:,1)); faceS = S(Fmask(:,1)); end
  if(face==2); faceR = R(Fmask(:,2)); faceS = T(Fmask(:,2)); end
  if(face==3); faceR = S(Fmask(:,3)); faceS = T(Fmask(:,3)); end
  if(face==4); faceR = S(Fmask(:,4)); faceS = T(Fmask(:,4)); end
  
  F = zeros(Np);
  VFace = Vandermonde2D(N, faceR, faceS);
  massFace = inv(VFace*VFace');
  
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

% plot reference element
% figure;
% plot3(r(fmask1),s(fmask1),t(fmask1),'r.')
% hold on;
% plot3(r(fmask2),s(fmask2),t(fmask2),'g.')
% plot3(r(fmask3),s(fmask3),t(fmask3),'b.')
% plot3(r(fmask4),s(fmask4),t(fmask4),'k.')
% xlabel('r')
% ylabel('s')
% zlabel('t')

