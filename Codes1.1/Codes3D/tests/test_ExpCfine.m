% Test ExpCfine for a fine part containing the center element

Globals3D; 
N = 1;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK5.neu');
StartUp3D;

nodeIndex = findNearestNode([0,0,0]);
fine_idx = nodeIndex(2);

U = zeros(3*2*Np*K,1);
InitMatLawsonSparse;
ReorderLawson;

expcfine = expm(Cfine);
expcfine2 = ExpCfine(1);

norm(expcfine - expcfine2)