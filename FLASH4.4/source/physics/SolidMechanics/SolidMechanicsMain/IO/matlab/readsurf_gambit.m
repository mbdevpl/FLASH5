function [XYZ,ELEM,nnod,nel]=readsurf_gambit(file)
% Routine to read from Gambit generic mesh, neutral format.

[fid] = fopen(file,'r');

%% Read through Header
nlines=6;
for i=1:nlines
   line=fgetl(fid); 
end
vect = str2num(fgetl(fid));
nnod = vect(1)
nel  = vect(2)

%% Read nodes
nlines=2;
for i=1:nlines
   line=fgetl(fid); 
end

XYZ = zeros(nnod+1,3);

for i=2:nnod+1
   vect=str2num(fgetl(fid));
   XYZ(i,1:3) = vect(2:4);
end

%% Read elem:
nlines=2;
for i=1:nlines
   line=fgetl(fid); 
end

ELEM = zeros(nel,3);

for i=1:nel
   vect=str2num(fgetl(fid));
   ELEM(i,1:3) = vect(4:6) + [1 1 1];
end


return

end