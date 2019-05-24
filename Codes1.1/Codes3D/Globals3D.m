% Purpose: declare global variables

global Np Nfp N K
global VX VY VZ
global r s t
global Dr Ds Dt LIFT Drw Dsw Dtw
global Fx Fy Fz nx ny nz Fscale sJ
global vmapM vmapP vmapB mapB Fmask
global BCType mapI mapO mapW mapF mapC mapS mapM mapP mapD mapN
global vmapI vmapO vmapW vmapO vmapC vmapS vmapD vmapN
global rx ry rz sx sy sz tx ty tz J
global rk4a rk4b rk4c
global Nfaces EToE EToF EToV
global V invV
global x y z NODETOL

% Variables for Lawson scheme
global MassMatrix globalMassMatrix globalStiffnessMatrix Cmat

% Variable to store field values over time
global Ez_time

In = 1; Out = 2; Wall = 3; Far = 4; Cyl = 5; Dirichlet = 6; Neuman = 7; Slip = 8;

% Low storage Runge-Kutta coefficients
rk4a = [            0.0 ...
        -567301805773.0/1357537059087.0 ...
        -2404267990393.0/2016746695238.0 ...
        -3550918686646.0/2091501179385.0  ...
        -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
         5161836677717.0/13612068292357.0 ...
         1720146321549.0/2090206949498.0  ...
         3134564353537.0/4481467310338.0  ...
         2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
         1432997174477.0/9575080441755.0 ...
         2526269341429.0/6820363962896.0 ...
         2006345519317.0/3224310063776.0 ...
         2802321613138.0/2924317926251.0];
		 