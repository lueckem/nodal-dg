%calculate facemass matrix for the first face of the reference element
%using the code
%Globals3D;
% 
% faceR = r(Fmask(:,1)); faceS = s(Fmask(:,1));
% VFace = Vandermonde2D(N, faceR, faceS);
% massFace = inv(VFace*VFace');

%solving the integral
r1 = [-1;-1;-1]; r2 = [1;-1;-1]; r3 = [-1; 1; -1]; r4 = [-1;-1;1];
r = [r1, r2, r3, r4];

l1([-1;-1;-1], r)

massFace2 = zeros(3);
for i=1:3
    for j=1:3
        fun = @(x,y) l(i,[x;y;-1],r)*l(j,[x;y;-1],r);
        massFace2(i,j) = integral2(fun,-1,1,-1,1);
    end
end

function [out] = l1(in,r)
out = 1;
for m=2:4
    out = out * (in(1)-r(1,m))/(r(1,1)-r(1,m))*(in(2)-r(2,m))/(r(2,1)-r(2,m))*(in(3)-r(3,m))/(r(3,1)-r(3,m));
end
end

function [out] = l2(in,r)
out = 1;
for m=1:4
    if m ~= 2
        out = out * (in(1)-r(1,m))/(r(1,2)-r(1,m))*(in(2)-r(2,m))/(r(2,2)-r(2,m))*(in(3)-r(3,m))/(r(3,2)-r(3,m));
    end
end
end

function [out] = l3(in,r)
out = 1;
for m=1:4
    if m ~= 3
        out = out * (in(1)-r(1,m))/(r(1,3)-r(1,m))*(in(2)-r(2,m))/(r(2,3)-r(2,m))*(in(3)-r(3,m))/(r(3,3)-r(3,m));
    end
end
end

function [out] = l(index, input, r)
if index == 1
    out = l1(input, r);
elseif index == 2
    out = l2(input, r);
else
    out = l3(input, r);
end
end