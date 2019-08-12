function PlotPlain3DFast(u, x_grid, y_grid, sampleTets, sampleWeights)
% plots a x-y-crossection of the field u according to griddata

% sample u
field = zeros(size(x_grid,1));

for i = 1:size(x_grid,1)
    for j = 1:size(y_grid,1)
        field(i,j) = dot(sampleWeights(:,i,j), u(:,sampleTets(i,j)));
    end
end

%plot
surf(x_grid, y_grid, field);
end


