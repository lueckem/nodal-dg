function [idx] = findNearestNode(coordinates)
% finds the index (i,k) of the node nearest to coordinates

Globals3D;

min_dist = norm([x(1,1), y(1,1), z(1,1)] - coordinates);
idx = [1,1];
for k = 1:K
   for i = 1:Np
      dist = norm([x(i,k), y(i,k), z(i,k)] - coordinates);
      if dist < min_dist
          min_dist = dist;
          idx = [i,k];
      end
   end
end
end
