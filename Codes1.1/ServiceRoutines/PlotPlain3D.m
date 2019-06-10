function PlotPlain3D(zloc, u)
% plots a x-y-crossection at z = zloc of the field u

Globals3D;

[x_grid, y_grid] = meshgrid(-1:0.05:1);

% sample u
field = zeros(size(x_grid,1));
for i = 1:size(x_grid,1)
    for j = 1:size(y_grid,1)
        [sampleweights,sampletet] = Sample3D(x_grid(i,j), y_grid(i,j), zloc);
        field(i,j) = dot(sampleweights, u(:,sampletet));
    end
end

%plot
surf(x_grid, y_grid, field);


end

