function [S_ikP] = S_ikPlusLawson(i, k, R, S, T)
% calculate the Matrix S_ik^Plus

Globals3D;

% check if element i and k have a common face
[bool, face_i] = ismember(k, EToE(i,:));
if ~bool
   return; 
end

S_ikP = zeros(3*Np);

% face_k is the face number of element k of the common face
face_k = EToF(i,face_i);

% surface Jacobian
sJacobian = sJ(Nfp*(face_i-1)+1,i);

%normal vector
normal = [nx(Nfp*(face_i-1)+1,i);
           ny(Nfp*(face_i-1)+1,i);
           nz(Nfp*(face_i-1)+1,i)];
       
%nodes on the face
nodes_i_plus = vmapP(4*Nfp*(i-1)+1 + (face_i-1)*Nfp : 4*Nfp*(i-1) + face_i*Nfp);
nodes_k_minus = vmapM(4*Nfp*(k-1)+1 + (face_k-1)*Nfp : 4*Nfp*(k-1) + face_k*Nfp);

% process face
if(face_i==1); faceR = R(Fmask(:,1)); faceS = S(Fmask(:,1)); end
if(face_i==2); faceR = R(Fmask(:,2)); faceS = T(Fmask(:,2)); end
if(face_i==3); faceR = S(Fmask(:,3)); faceS = T(Fmask(:,3)); end
if(face_i==4); faceR = S(Fmask(:,4)); faceS = T(Fmask(:,4)); end
  
F = zeros(Np);
VFace = Vandermonde2D(N, faceR, faceS);
massFace = sJacobian .* inv(VFace*VFace');

%build the F matrix
for j = 1:Nfp
   for l = 1:Nfp
       %find the index l2 of the point in element i corresponding to the point
       %l in element k
       [~, l2] = ismember(nodes_k_minus(l), nodes_i_plus);
       
       F(Fmask(j,face_i), Fmask(l,face_i)) = -0.5 * massFace(j, l2);
   end
end

%build S_ikPlus
S_ikP(1:Np,Np+1:2*Np) = -normal(3).*F;
S_ikP(1:Np,2*Np+1:3*Np) = normal(2).*F;
S_ikP(Np+1:2*Np,1:Np) = normal(3).*F;
S_ikP(Np+1:2*Np,2*Np+1:3*Np) = -normal(1).*F;
S_ikP(2*Np+1:3*Np,1:Np) = -normal(2).*F;
S_ikP(2*Np+1:3*Np,Np+1:2*Np) = normal(1).*F;

end

