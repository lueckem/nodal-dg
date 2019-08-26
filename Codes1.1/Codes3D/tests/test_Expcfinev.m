% Test ExpCfinev versus ExpCfine

Globals3D; 
N = 2;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK268.neu');
StartUp3D;


fine_idx = 1:30;

U = rand(3*2*Np*K,1);
InitMatLawsonSparse;
ReorderLawson;

alpha = -0.1;
expcfine = ExpCfine(alpha);
expcfineU = expcfine * U;
expcfineU2 = ExpCfinev(alpha,U);

norm(expcfineU - expcfineU2)