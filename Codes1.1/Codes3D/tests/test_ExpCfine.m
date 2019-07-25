% Test ExpCfine for a fine part containing the 5 center elements

fine_idx = [];
nodeIndex = findNearestNode([0,0,0]);
fine_idx = [nodeIndex(2),EToE(nodeIndex(2),:)];

U = zeros(3*2*Np*K,1);
InitMatLawsonSparse;
ReorderLawson;

expcfine = ExpCfine(-2);
figure; spy(expcfine);