Globals3D; 
N = 1;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK86.neu');
StartUp3D;
%%
%find two elements that share a face
i = 3;
for j = 1:4
   if EToE(i,j) ~= i
      k =  EToE(i,j);
      break;
   end
end

% find the corresponding face numbers
[~, face_i] = ismember(k, EToE(i,:));
face_k = EToF(i,face_i);

% find the global indices of the nodes on the faces
nodes_i_plus = vmapP(4*Nfp*(i-1)+1 + (face_i-1)*Nfp : 4*Nfp*(i-1) + face_i*Nfp);
nodes_k_minus = vmapM(4*Nfp*(k-1)+1 + (face_k-1)*Nfp : 4*Nfp*(k-1) + face_k*Nfp);

l = 1;
[~, l2] = ismember(nodes_k_minus(l), nodes_i_plus);