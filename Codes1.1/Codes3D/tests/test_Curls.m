% compare the curls from Curl3D with the curls from Lawson

Globals3D; 
N = 1;
[Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D('cube.neu');
StartUp3D;

%%
Ex = sin(1*pi*x).*sin(1*pi*y);
Ey = sin(2*pi*x).*sin(2*pi*y);
Ez = sin(3*pi*x).*sin(3*pi*y);
Hx = sin(4*pi*x).*sin(4*pi*y);
Hy = sin(5*pi*x).*sin(5*pi*y);
Hz = sin(6*pi*x).*sin(6*pi*y);
U = FieldsToU(Hx, Hy, Hz, Ex, Ey, Ez);

% curls from Curl3D
[curlHx,curlHy,curlHz] = Curl3D(Hx,Hy,Hz);
[curlEx,curlEy,curlEz] = Curl3D(Ex,Ey,Ez);

% calc Dmat for MaxwellRHS3D
blksize = 3 * Np;
Dmat = sparse(K*3*Np, K*3*Np, 9*Np^2*K);
for i = 1:K
    Dx = rx(1,i) * Dr + sx(1,i) * Ds + tx(1,i) * Dt;
    Dy = ry(1,i) * Dr + sy(1,i) * Ds + ty(1,i) * Dt;
    Dz = rz(1,i) * Dr + sz(1,i) * Ds + tz(1,i) * Dt;
    S = [zeros(Np), -Dz, Dy; Dz, zeros(Np), -Dx; -Dy, Dx, zeros(Np)];
    Dmat((i-1)*blksize+1:i*blksize, (i-1)*blksize+1:i*blksize) = S;
end
Dmat = [sparse(K*3*Np, K*3*Np), Dmat; -Dmat, sparse(K*3*Np, K*3*Np)];

%curls from Lawson
curlU = Dmat*U;
[curlHx2, curlHy2, curlHz2, curlEx2, curlEy2, curlEz2] = UToFields(curlU);

%comparison
norm(-curlEx-curlHx2)
norm(-curlEy-curlHy2)
norm(-curlEz-curlHz2)
norm(curlHx-curlEx2)
norm(curlHy-curlEy2)
norm(curlHz-curlEz2)

% Conclusion: The Curls are the same!

