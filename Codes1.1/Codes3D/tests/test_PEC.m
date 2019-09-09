% compare the fluxes from Maxwell3D with the fluxes from Lawson for a PEC
% boundary condition

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
% 
% Ex = 1 * rand(Np,K);
% Ey = 2 * rand(Np,K);
% Ez = 3 * rand(Np,K);
% Hx = 4 * rand(Np,K);
% Hy = 5 * rand(Np,K);
% Hz = 10 * rand(Np,K);

% Ex = 0 * ones(Np,K);
% Ey = 0 * ones(Np,K);
% Ez = 0 * ones(Np,K);
% Hx = 0 * ones(Np,K);
% Hy = 0 * ones(Np,K);
% Hz = 0 * ones(Np,K);
% Hx(:,3) = 1 * ones(Np,1);

U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

%% Flux from Maxwell3D
% storage for field differences at faces
dHx = zeros(Nfp*Nfaces,K); dHy = dHx; dHz = dHx; 
dEx = zeros(Nfp*Nfaces,K); dEy = dEx; dEz = dEx; 

% form field differences at faces  
dHx(:)  = Hx(vmapP)-Hx(vmapM);  dEx(:)  = Ex(vmapP)-Ex(vmapM);	
dHy(:)  = Hy(vmapP)-Hy(vmapM); 	dEy(:)  = Ey(vmapP)-Ey(vmapM);	
dHz(:)  = Hz(vmapP)-Hz(vmapM);  dEz(:)  = Ez(vmapP)-Ez(vmapM);  

% PEC boundary conditions
dEx(mapB) = -2*Ex(vmapB); 
dEy(mapB) = -2*Ey(vmapB); 
dEz(mapB) = -2*Ez(vmapB);

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

KmatE = sparse(3*Np*K,3*Np*K,30*K*Np^2);
KmatH = sparse(3*Np*K,3*Np*K,30*K*Np^2);
blksize = 3 * Np;
invM = inv(MassMatrix);

% KmatE
for i = 1:K %rows
    invMi = invM./J(1, i);
    invMi = blkdiag(invMi, invMi, invMi);
    
    % diagonal element
    KmatE((i-1)*blksize+1:i*blksize, (i-1)*blksize+1:i*blksize) = -invMi * SurfaceMassInteriorLawson(r,s,t,i);
    
    % the SikP
    for j = 1:4
        k = EToE(i,j);
        if k~=i %check for boundary
            KmatE((i-1)*blksize+1:i*blksize, (k-1)*blksize+1:k*blksize) = -invMi * S_ikPlusLawson(i,k,r,s,t);
        end
    end
end

% KmatH
KmatH = KmatE;
for i = 1:K %rows
    invMi = invM./J(1, i);
    invMi = blkdiag(invMi, invMi, invMi);
    
    % diagonal element
    KmatH((i-1)*blksize+1:i*blksize, (i-1)*blksize+1:i*blksize) = -invMi * SurfaceMassLawson(r,s,t,i);
end

Kmat = [sparse(3*Np*K, 3*Np*K), KmatE; -KmatH, sparse(3*Np*K, 3*Np*K)];
fluxU = Kmat * U;
[fluxHx2, fluxHy2, fluxHz2, fluxEx2, fluxEy2, fluxEz2] = UToFields(fluxU);

norm(fluxEx - fluxEx2)
norm(fluxEy - fluxEy2)
norm(fluxEz - fluxEz2)
norm(fluxHx - fluxHx2)
norm(fluxHy - fluxHy2)
norm(fluxHz - fluxHz2)