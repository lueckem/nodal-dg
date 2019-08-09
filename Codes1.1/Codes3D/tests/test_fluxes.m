%% The interior fluxes work!

% compare the fluxes from Maxwell3D with the fluxes from Lawson

Globals3D; 
N = 3;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cubeK86.neu');
StartUp3D;

Ex = sin(1*pi*x).*sin(1*pi*y);
Ey = sin(2*pi*x).*sin(2*pi*y);
Ez = sin(3*pi*x).*sin(3*pi*y);
Hx = sin(4*pi*x).*sin(4*pi*y);
Hy = sin(5*pi*x).*sin(5*pi*y);
Hz = sin(6*pi*x).*sin(6*pi*y);

% Ex = 1 * rand(Np,K);
% Ey = 2 * rand(Np,K);
% Ez = 3 * rand(Np,K);
% Hx = 4 * rand(Np,K);
% Hy = 5 * rand(Np,K);
% Hz = 6 * rand(Np,K);

% Ex = 0 * ones(Np,K);
% Ey = 0 * ones(Np,K);
% Ez = 0 * ones(Np,K);
% Hx = 0 * ones(Np,K);
% Hy = 0 * ones(Np,K);
% Hz = 0 * ones(Np,K);
% Hx(:,3) = 1 * ones(Np,1);

U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

%% Flux from Maxwell3D, no boundaries
% storage for field differences at faces
dHx = zeros(Nfp*Nfaces,K); dHy = dHx; dHz = dHx; 
dEx = zeros(Nfp*Nfaces,K); dEy = dEx; dEz = dEx; 

% form field differences at faces  
dHx(:)  = Hx(vmapP)-Hx(vmapM);  dEx(:)  = Ex(vmapP)-Ex(vmapM);	
dHy(:)  = Hy(vmapP)-Hy(vmapM); 	dEy(:)  = Ey(vmapP)-Ey(vmapM);	
dHz(:)  = Hz(vmapP)-Hz(vmapM);  dEz(:)  = Ez(vmapP)-Ez(vmapM);  

% make silver mueller boundary conditions
% dHx(mapB) = -2*Hx(vmapB);  dEx(mapB) = -2*Ex(vmapB); 
% dHy(mapB) = -2*Hy(vmapB);  dEy(mapB) = -2*Ey(vmapB); 
% dHz(mapB) = -2*Hz(vmapB);  dEz(mapB) = -2*Ez(vmapB);

fluxHx = -ny.*dEz + nz.*dEy;
fluxHy = -nz.*dEx + nx.*dEz;
fluxHz = -nx.*dEy + ny.*dEx;

fluxEx =  ny.*dHz - nz.*dHy;
fluxEy =  nz.*dHx - nx.*dHz;
fluxEz =  nx.*dHy - ny.*dHx;

% calculate fluxes
fluxHx = LIFT*(Fscale.*fluxHx/2);
fluxHy = LIFT*(Fscale.*fluxHy/2);
fluxHz = LIFT*(Fscale.*fluxHz/2);

fluxEx =  LIFT*(Fscale.*fluxEx/2);
fluxEy =  LIFT*(Fscale.*fluxEy/2);
fluxEz =  LIFT*(Fscale.*fluxEz/2);

%% Flux from Lawson, no boundaries
% K matrix
Kmat = sparse(3*Np*K,3*Np*K,30*K*Np^2);
blksize = 3 * Np;
invM = inv(MassMatrix);

for i = 1:K %rows
    invMi = invM./J(1, i);
    invMi = blkdiag(invMi, invMi, invMi);
    
    % diagonal element
    Kmat((i-1)*blksize+1:i*blksize, (i-1)*blksize+1:i*blksize) = -invMi * SurfaceMassInteriorLawson(r,s,t,i);
    
    % the SikP
    for j = 1:4
        k = EToE(i,j);
        if k~=i %check for boundary
            Kmat((i-1)*blksize+1:i*blksize, (k-1)*blksize+1:k*blksize) = -invMi * S_ikPlusLawson(i,k,r,s,t);
        end
    end
end


% % SE matrix
% SE = sparse(3*Np*K,3*Np*K,9*K*Np^2);
% for i = 1:K
%     invMi = invM./J(1, i);
%     invMi = blkdiag(invMi, invMi, invMi);
%     SE((i-1)*blksize+1:i*blksize, (i-1)*blksize+1:i*blksize) = invMi * S_iELawson(i,r,s,t);
% end

Kmat = [sparse(3*Np*K, 3*Np*K), Kmat; -Kmat, sparse(3*Np*K, 3*Np*K)];
fluxU = Kmat * U;
[fluxHx2, fluxHy2, fluxHz2, fluxEx2, fluxEy2, fluxEz2] = UToFields(fluxU);

norm(fluxEx - fluxEx2)
norm(fluxEy - fluxEy2)
norm(fluxEz - fluxEz2)
norm(fluxHx - fluxHx2)
norm(fluxHy - fluxHy2)
norm(fluxHz - fluxHz2)


%% check flux for E

% %find two elements that share a face
% i = 3;
% for j = 1:4
%    if EToE(i,j) ~= i
%       k =  EToE(i,j);
%       break;
%    end
% end
% [~, face_i] = ismember(k, EToE(i,:));
% 
% % calculate the flux on the face between element i and k
% 
% dHx = zeros(Nfp*Nfaces,K); dHy = dHx; dHz = dHx;
% 
% dHx(:)  = Hx(vmapP)-Hx(vmapM);
% dHy(:)  = Hy(vmapP)-Hy(vmapM);
% dHz(:)  = Hz(vmapP)-Hz(vmapM);
% 
% fluxEx =  ny.*dHz - nz.*dHy; fluxEx_i = fluxEx((face_i-1)*Nfp+1:face_i*Nfp,i);
% fluxEy =  nz.*dHx - nx.*dHz; fluxEy_i = fluxEy((face_i-1)*Nfp+1:face_i*Nfp,i);
% fluxEz =  nx.*dHy - ny.*dHx; fluxEz_i = fluxEz((face_i-1)*Nfp+1:face_i*Nfp,i);
% 
% Fscale_i = Fscale((face_i-1)*Nfp+1:face_i*Nfp,i);
% LIFT_i = LIFT(1:Np, (face_i-1)*Nfp+1:face_i*Nfp);
% 
% fluxEx_i =  LIFT_i*(Fscale_i.*fluxEx_i/2);
% fluxEy_i =  LIFT_i*(Fscale_i.*fluxEy_i/2);
% fluxEz_i =  LIFT_i*(Fscale_i.*fluxEz_i/2);
% 
% fluxE = [fluxEx_i; fluxEy_i; fluxEz_i];
% 
% 
% %Lawson
% S_ik = S_ikLawson(r,s,t,i,k);
% S_ikPlus = S_ikPlusLawson(i,k,r,s,t);
% invMi = blkdiag(inv(MassMatrix), inv(MassMatrix), inv(MassMatrix)) ./ J(1,i);
% fluxE2 = -invMi * S_ik * [Hx(:,i);Hy(:,i);Hz(:,i)] - invMi * S_ikPlus * [Hx(:,k);Hy(:,k);Hz(:,k)];
% 
% %fluxE = fluxE2 !!!!
