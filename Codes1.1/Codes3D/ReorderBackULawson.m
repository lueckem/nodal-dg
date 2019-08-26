function [Utmp] = ReorderBackULawson(U)
% Reorder U back to the original order after the model reduction

Globals3D;

Utmp = zeros(3*2*Np*K,1);
for j = 1:3*2*Np*K
   Utmp(C_idx(j)) = U(j); 
end
end

