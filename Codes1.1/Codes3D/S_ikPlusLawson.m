function [S_ikP] = S_ikPlusLawson(i, k, N, R, S, T)
% calculate the Matrix S_ik^Plus
% todo: boundary conditions?

Globals3D;
S_ikP = zeros(3*Np);

% check if element i and k have a common face
if ~ismember(k, EToE(i,:))
   return; 
end

% face_i is the face number of element i of the common face
[face_i, ~] = ismember(k, EToF(i,:));
[face_k, ~] = ismember(i, EToF(k,:));

% surface Jacobian
sJacobian = sJ(Nfp*(face_i-1)+1,i);

%normal vector
normal = [nx(Nfp*(face_i-1)+1,i);
           ny(Nfp*(face_i-1)+1,i); 
           nz(Nfp*(face_i-1)+1,i)];
       
%nodes on the face
nodes_i = vmapP(4*Nfp*(i-1)+1 + (face_i-1)*Nfp : 4*Nfp*(i-1)+1 + face_i*Nfp);
nodes_k = vmapM(4*Nfp*(i-1)+1 + (face_i-1)*Nfp : 4*Nfp*(i-1)+1 + face_i*Nfp);

% process face
if(face_i==1); faceR = R(Fmask(:,1)); faceS = S(Fmask(:,1)); end
if(face_i==2); faceR = R(Fmask(:,2)); faceS = T(Fmask(:,2)); end
if(face_i==3); faceR = S(Fmask(:,3)); faceS = T(Fmask(:,3)); end
if(face_i==4); faceR = S(Fmask(:,4)); faceS = T(Fmask(:,4)); end
  
F = zeros(Np);
VFace = Vandermonde2D(N, faceR, faceS);
massFace = sJacobian .* inv(VFace*VFace');

for j = 1:Nfp
   for l = 1:Nfp
       %todo
   end
end
  
  
  

end

