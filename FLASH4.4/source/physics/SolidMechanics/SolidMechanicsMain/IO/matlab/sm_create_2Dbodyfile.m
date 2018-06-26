% Matlab program to generate simple 2D firface grids for rigid bodies,
% cylinder, square, ellipse.
clear all
close all
clc

%% Read input file:
file  = 'cylinder_matlab.dat'; 
[fid] = fopen(file,'r');

lg   =10;
[lin]=fgetl(fid);

[lin]=fgetl(fid);
dim = lin(1:2);

[lin]=fgetl(fid);
% Radius
R  = str2num(lin(1:lg));
% DeltaS
[lin]=fgetl(fid);
ds   =str2num(lin(1:lg));

if isequal(dim,'2D')        
    Lz = 0;  
    dz = 1;
elseif isequal(dim,'3D')
    % Length in Z
    [lin]=fgetl(fid);
    Lz = str2num(lin(1:lg));
  
    % dZ
    [lin]=fgetl(fid);
    dz = str2num(lin(1:lg));
else
    disp('Error: dim must be 2D or 3D')
    return
end
fclose(fid);


%% Circumference length:
Lcirc = 2*pi*R

% Compute dtheta
nseg_theta = ceil(Lcirc/ds)
dseg_theta = Lcirc/nseg_theta
dtheta     = dseg_theta/R

% Compute dzeta
nseg_z     = ceil(Lz/dz)
dseg_z     = 1.0;
if (nseg_z > 0.0)
    dseg_z = Lz/nseg_z;
end

% For 3D Compute dR:
if isequal(dim,'3D')
    nseg_r     = floor(R/ds);
    dseg_r     = R/nseg_r;
end

% Create center point:
i = 1;
XYZ(i,1:3) = [0 0 0];

% Compute Surface points
if isequal(dim,'2D')
    npt_theta = nseg_theta + 1;
    npt_z     = 1;
    zmin      = 0;
else
    npt_theta = nseg_theta; % in 3D we dont repeat the last point.
    npt_z     = nseg_z + 1;
    npt_r     = nseg_r - 1; % From row 2 to row nseg_r-1
    zmin      = XYZ(1,3) - Lz/2;
    zmax      = XYZ(1,3) + Lz/2;
end

%% Loops:
% Nodes
np = npt_z*npt_theta + 1;
if isequal(dim,'3D')
    npt = npt + 2*(1+npt_r*npt_theta); % Add side points
end
for ipz = 1:npt_z
    zj = zmin + (ipz-1)*dseg_z;
    for ipt = 1:npt_theta
        i  = i+1;
        ti = (ipt-1)*dtheta;
        XYZ(i,1:3) = [R*cos(ti) R*sin(ti) zj];    
    end
end

% Nodes on sides of 3D Cylinder
if isequal(dim,'3D')
    
    
    % Low z side
    i=i+1;
    XYZ(i,1:3)=[XYZ(1,1) XYZ(1,2) zmin];
    
    for ipr = 1:npt_r
        r=ipr*dseg_r;
        for ipt = 1:npt_theta
            i=i+1;
            ti = (ipt-1)*dtheta;
            XYZ(i,1:3)=[r*cos(ti) r*sin(ti) zmin];
        end
    end
    
    % High z side
    i=i+1;
    XYZ(i,1:3)=[XYZ(1,1) XYZ(1,2) zmax];    

    for ipr = 1:npt_r
        r=ipr*dseg_r;
        for ipt = 1:npt_theta
            i=i+1;
            ti = (ipt-1)*dtheta;
            XYZ(i,1:3)=[r*cos(ti) r*sin(ti) zmax];
        end
    end
        
end


% Elements
iel    =  1;
eltype = 15;  % Nomenclature as in Gmsh
eltype2=  1;  % 2 node segment.
eltype3=  8;  % 3 node triangle.
NOAELEM=  0;  % This element is not part of the Surface grid  
AELEM  =  1;  % This element is part of the Surface grid
ELEM(iel,:) = [iel eltype 3 0 1 NOAELEM 1];
if isequal(dim,'2D')
    for ipt=2:npt_theta
        iel = iel + 1;
        ELEM(iel,1:8) = [iel eltype2 3 0 2 AELEM ipt ipt+1];
    end
else % 3D
    for ipz = 1:npt_z
        for ipt = 1:npt_theta
            iel = iel + 1;
            
            
        end
    end
end

np  = length(XYZ(:,1));
nel = length(ELEM(:,1));

%% Write File: Enhanced Gmsh ASCI II Format.
file_out = ['Cylinder_' dim '_' num2str(ds,4) '.msh'];

[fod]=fopen(file_out,'w');

[h]=fprintf(fod,'$MeshFormat\n');
[h]=fprintf(fod,'2.2  0  8\n');
[h]=fprintf(fod,'$EndMeshFormat\n');
% Nodes
[h]=fprintf(fod,'$Nodes\n');
[h]=fprintf(fod,'%d \n',np);
if isequal(dim,'2D')
    for ip=1:np
        ngcoord = 3;
        tag = [-10 -10 -10]; % This is a constrained rigid surface node in the x,y and theta directions.
        if ip == 1 % Prescribed coordinates.
            tag = [-1 -1 -1];
        end
        [h]=fprintf(fod,'%d %12.8f %12.8f %12.8f\n',[ip XYZ(ip,:)]); %%d %d %d %d  ngcoord tag
    end
else
    
end
[h]=fprintf(fod,'$EndNodes\n',np);
% Elements
[h]=fprintf(fod,'$Elements\n');
[h]=fprintf(fod,'%d \n',nel);
if isequal(dim,'2D')
    for iel=1:nel
        [h]=fprintf(fod,'%d %d %d %d %d %d %d %d\n',ELEM(iel,:));
    end
else
    
end
[h]=fprintf(fod,'$EndElements\n',np);

% % Data:
% [h]=fprintf(fod,'$NodeData\n');
% NumStringTags=0;
% [h]=fprintf(fod,'%d\n',NumStringTags);
% %number-of-string-tags
% %< "string-tag" >
% %...
% NumRealTags=0;
% [h]=fprintf(fod,'%d\n',NumRealTags);
% %number-of-real-tags
% %< real-tag >
% %...
% NumIntegTags=0;
% [h]=fprintf(fod,'%d\n',NumIntegTags);
% %number-of-integer-tags
% %< integer-tag >
% %...
% %node-number value ...
% %...
% [h]=fprintf(fod,'$EndNodeData\n');
% [h]=fprintf(fod,'$ElementData\n');
% NumStringTags=0;
% [h]=fprintf(fod,'%d\n',NumStringTags);
% % number-of-string-tags
% % < "string-tag" >
% % ...
% NumRealTags=0;
% [h]=fprintf(fod,'%d\n',NumRealTags);
% %     number-of-real-tags
% % < real-tag >
% % ...
% NumIntegTags=0;
% [h]=fprintf(fod,'%d\n',NumIntegTags);
% %     number-of-integer-tags
% % < integer-tag >
% % ...
% %     elm-number value ...
% %     ...
% [h]=fprintf(fod,'$EndElementData\n');
% [h]=fprintf(fod,'$ElementNodeData\n')
% NumStringTags=0;
% [h]=fprintf(fod,'%d\n',NumStringTags);
% % number-of-string-tags
% % < "string-tag" >
% % ...
% NumRealTags=0;
% [h]=fprintf(fod,'%d\n',NumRealTags);
% %     number-of-real-tags
% % < real-tag >
% % ...
% NumIntegTags=0;
% [h]=fprintf(fod,'%d\n',NumIntegTags);
% %     number-of-integer-tags
% % < integer-tag >
% % ...
% %     elm-number number-of-nodes-per-element value ...
% %     ...
% [h]=fprintf(fod,'$EndElementNodeData\n')
% [h]=fprintf(fod,'$InterpolationScheme\n')
% % "name"
% % number-of-element-topologies
% % elm-topology
% % number-of-interpolation-matrices
% % num-rows num-columns value ...
% %     ...
% [h]=fprintf(fod,'$EndInterpolationScheme\n')

fclose(fod);


% Figure:
figure
hold on
plot(XYZ(:,1),XYZ(:,2),'ok')
if isequal(dim,'2D')
    for ipt=2:nel
        H=line(XYZ([ELEM(ipt,1+6) ELEM(ipt,1+7)],1), ...
            XYZ([ELEM(ipt,1+6) ELEM(ipt,1+7)],2), ...
            'Color','r','Marker','x','LineWidth',2);
    end
else
    
end
axis equal
grid on
box on

return