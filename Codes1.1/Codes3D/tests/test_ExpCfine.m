% Test ExpCfine for a fine part containing the 5 center elements

Globals3D; 
N = 1;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK86.neu');
StartUp3D;

fine_idx = [];
nodeIndex = findNearestNode([0,0,0]);
fine_idx = [nodeIndex(2),EToE(nodeIndex(2),:)];

U = zeros(3*2*Np*K,1);
InitMatLawsonSparse;
ReorderLawson;

%expcfine = ExpCfine(-2);
%figure; spy(expcfine);