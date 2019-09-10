function [Nv, VX, VY, VZ, K, EToV, epsilon] = MeshReaderGambit3DMaterial(FileName)

% function [Nv, VX, VY, VZ, K, EToV] = MeshReaderGambit3D(FileName)
% Purpose  : Read in basic grid information to build grid + material info
% NOTE     : gambit *.neu format is assumed

Fid = fopen(FileName, 'rt');

% read intro 
for i=1:6 
  line = fgetl(Fid);
end

% fine number of nodes and number of elements
dims = fscanf(Fid, '%d');

Nv = dims(1); K = dims(2); Ngroups = dims(3);

for i=1:2 
  line = fgetl(Fid);
end

% read node coordinates
xyz = fscanf(Fid, '%lf', [4, Nv]);
xyz = xyz(2:4, :);
VX = xyz(1,:); VY = xyz(2,:); VZ = xyz(3,:);

for i=1:3 
  line = fgetl(Fid);
end

% read element to node connectivity
EToV = zeros(K, 4);
for k = 1:K
  line   = fgetl(Fid);
  tmpcon = sscanf(line, '%lf');
  EToV(k,1:4) = tmpcon(4:7);
end

for i=1:3 
  line = fgetl(Fid);
end

% Read in the material parameters for the different groups
epsilon = ones(K,1);

for j=1:Ngroups
    % get the epsilon for this group
    Str = line;
    Str(strfind(Str, ':')) = [];
    Key   = 'MATERIAL';
    Index = strfind(Str, Key);
    eps_val = sscanf(Str(Index(1) + length(Key):end), '%g', 1);
    
    for i=1:3
        line = fgetl(Fid);
    end
    
    idx = sscanf(line, '%lf'); %indices of the elements in the group
    epsilon(idx) = eps_val;
    
    line = fgetl(Fid);
end

return