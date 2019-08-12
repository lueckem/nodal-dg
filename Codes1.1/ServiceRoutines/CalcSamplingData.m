function [x_grid, y_grid, sampleTets, sampleWeights] = CalcSamplingData(zloc)
% Calculate the Sampling Data for plotting the x-y-crossection at z=zloc

Globals3D;
min_x = min(x,[],'all');
max_x = max(x,[],'all');
step = (max_x - min_x)/50;
[x_grid, y_grid] = meshgrid(min_x:step:max_x);

sampleTets = zeros(size(x_grid));
sampleWeights = zeros(Np, size(x_grid,1), size(x_grid,2));

for i = 1:size(x_grid,1)
    for j = 1:size(y_grid,1)
        [sampleWeights(:,i,j),sampleTets(i,j)] = Sample3D(x_grid(i,j), y_grid(i,j), zloc);
    end
end

end

