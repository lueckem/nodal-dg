% Reorder Cfine and Ccoarse and U for the model reduction

C_idx = CalcOrderingIdx(Np, K, fine_idx);

% columns
% Cfine = Cfine(:,C_idx);
% Ccoarse = Ccoarse(:,C_idx);
% 
% % rows 
% Cfine = Cfine(C_idx,:);
% Ccoarse = Ccoarse(C_idx,:);
U = U(C_idx);