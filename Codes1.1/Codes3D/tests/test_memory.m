% How much memory for C?

Globals3D;
fine_idx = [];

results = [];

grids = {'cubeK5.neu', 'cubeK6.neu', 'cubeK86.neu', 'cubeK268.neu'};

for Nval = 1:4
    for j = 1:length(grids)
        N = Nval;
        [Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D(grids{j});
        StartUp3D;
        InitMatLawsonSparse;
        results = [results; [Np, K, nnz(Ccoarse)]];
    end
end
%%
N = 2;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cube.neu');
StartUp3D;
InitMatLawsonSparse;
results = [results; [Np, K, nnz(Ccoarse)]];
%% fit data
sf = fit([results(:,1), results(:,2)], results(:,3),'poly21');
figure;
plot(sf,[results(:,1), results(:,2)], results(:,3))

maxdif = 0;
index_maxdif = 0;
for j = 1:size(results,1)
    dif = results(j,3) - sf(results(j,1), results(j,2));
    if dif > maxdif
        maxdif = dif;
        index_maxdif = j;
    end
end

%% my own calculated bound 78*Np^2*K
f = @(x,y) 78.*x.^2.*y;
[X,Y] = meshgrid(4:35,5:20:1585);
Z = f(X,Y);
figure;
surf(X,Y,Z)
hold on;
plot3(results(:,1), results(:,2), results(:,3), 'o');

maxdif = 0;
index_maxdif = 0;
for j = 1:size(results,1)
    dif = results(j,3) - f(results(j,1), results(j,2));
    if dif > maxdif
        maxdif = dif;
        index_maxdif = j;
    end
end