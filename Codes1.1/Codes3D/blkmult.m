function [result] = blkmult(sigma, phi, np, k)
% Calculate the multiplikation of a 3*K vector sigma with a 6*Np*K vector phi. The
% multiplication is pointwise for every (Np)-block
result = zeros(6*np*k,1);

for i = 1:3*k
   result(np*(i-1)+1:np*i) = sigma(i).*phi(np*(i-1)+1:np*i);
   result(np*(i-1)+1+3*k:np*i+3*k) = sigma(i).*phi(np*(i-1)+1+3*k:np*i+3*k);
end
end

