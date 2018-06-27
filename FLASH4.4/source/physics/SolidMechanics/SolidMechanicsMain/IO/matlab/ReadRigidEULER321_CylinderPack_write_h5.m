% Script to write Rigid Body h5 file, to be used by Splash.
% Name of file for simulation: sm_body.1.h5
% The Rigid is oriented using a 321 Euler angle sequence.
% -------------------------------------------------------------------------
close all
clear all
clc

%% Define variables:
IAXIS = 1;
JAXIS = 2;
KAXIS = 3;
IPHI  = 4;
NDIM  = 3; % Dimensions

SM_FALSE = 0;
SM_TRUE  = 1;

% Rigid Body constants:
BODYTYPE_RIGID  = 1;

RB_IDENTITY = 1000;
RB_EULER321 = 1321;
RB_QUATERNN = 1111;

RB_ANNSPHERE = 55;
RB_ANNDISC   = 56;
RB_ANNRBC    = 57;


% Kinematics Constants:       
SM_PK_FIXED    =   0;
SM_PK_HARMONIC = 102;
SM_PK_CONSTVEL = 103;

DOFS_per_node   =  9;  % Degrees of freedom per node: 
                       % x,y,z,ang1,ang2,ang3,wx,wy,wz

ix = 1; ex = 3;
ia = 4; ea = 6;
iw = 7; ew = 9;

%% Define variables:
%NumBods = 1; 

np_el   = 3; % Points per wet surface element: 3 - Triangle

% Output File names:
basedir_out = ['/Users/mvanella/Dropbox/Documents/ADAPTIVE/INS_IB/source/physics/SolidMechanics/SolidMechanicsMain/IO/matlab/'];

writefile = SM_TRUE; %SM_FALSE; 

% Prescribed Kinematics Per Body: 
kinemflag=['fixed'];

%% Cylinder Mesh File:

cyldir = KAXIS; % Cylinder direction. Direction 1.

num_seg= 4; % Number of segments per cylinder
nCyl_1 = num_seg;
nCyl_2 = 1;   % Number of Cylinders in the y direction
nCyl_3 = 1;   % Number of Cylinders in the z direction

NumBods = nCyl_1*nCyl_2*nCyl_3;

D  =    1; % Cylinder Diameter
L_t=  4*D;

L = L_t/nCyl_1; % Cylinder length (in idir direction from -L/2 to L/2)

D1 = L;
D2 = 2*D;
D3 = D2;

n_sections_L     = 12; % Number of points in the length
n_sections_theta = 35; % Number of sections in circumf direction.

dL     =    L/n_sections_L;
dtheta = pi*D/n_sections_theta;

% Now make the mesh
[XYZ,ws_IEN,nnodes,nel]=createCylinderMesh(cyldir,D,L,dL,dtheta);

% % Plot figure:
% figure
% hold on
% trimesh(ws_IEN,XYZ(1:nnodes+1,IAXIS),XYZ(1:nnodes+1,JAXIS),XYZ(1:nnodes+1,KAXIS))
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis equal
% 
% % Test normals out:
% for iel=1:nel
%    x12 = XYZ(ws_IEN(iel,2),:) - XYZ(ws_IEN(iel,1),:);
%    x13 = XYZ(ws_IEN(iel,3),:) - XYZ(ws_IEN(iel,1),:);
%    vcr = cross(x12,x13);
%    xcen = 1/3*(XYZ(ws_IEN(iel,1),:)+XYZ(ws_IEN(iel,2),:)+XYZ(ws_IEN(iel,3),:));
%    sign_normal(iel) = sign(dot(xcen,vcr));
%    if sign_normal(iel) < 0
%        aux=ws_IEN(iel,3); ws_IEN(iel,3)=ws_IEN(iel,2); ws_IEN(iel,2)=aux;
%    end
% end
% 
% for iel=1:nel
%    x12 = XYZ(ws_IEN(iel,2),:) - XYZ(ws_IEN(iel,1),:);
%    x13 = XYZ(ws_IEN(iel,3),:) - XYZ(ws_IEN(iel,1),:);
%    vcr = cross(x12,x13);
%    xcen = 1/3*(XYZ(ws_IEN(iel,1),:)+XYZ(ws_IEN(iel,2),:)+XYZ(ws_IEN(iel,3),:));
%    sign_normal2(iel) = sign(dot(xcen,vcr));
% end
% 
% return
% 


figure; hold on
if ( cyldir == IAXIS)

    nCyl(IAXIS) = nCyl_1;
    nCyl(JAXIS) = nCyl_2;
    nCyl(KAXIS) = nCyl_3;
    
    xstart = -(nCyl_1-1)/2*D1;
    ystart = -(nCyl_2-1)/2*D2;
    zstart = -(nCyl_3-1)/2*D3;
    
    delx = D1;
    dely = D2;
    delz = D3;
    
elseif ( cyldir == KAXIS)
    
    nCyl(IAXIS) = nCyl_2;
    nCyl(JAXIS) = nCyl_3;
    nCyl(KAXIS) = nCyl_1;
    
    xstart = -(nCyl_2-1)/2*D2; 
    ystart = -(nCyl_3-1)/2*D3; 
    zstart = -(nCyl_1-1)/2*D1;
    
    delx = D2;
    dely = D3;
    delz = D1;
    
end
    
    
ibd = 0;
for kcyl=1:nCyl(KAXIS)
    
  zpos = zstart + (kcyl-1)*delz;
    
  for jcyl=1:nCyl(JAXIS)
 
   ypos = ystart + (jcyl-1)*dely;   
      
   for iseg = 1:nCyl(IAXIS)

    xpos = xstart + (iseg-1)*delx;
       
    ibd = ibd + 1;

    %% Rigid body properties:
    grav_vec     = [0 0 -1.];
    grav_flag    =  0; % 0 = No gravity, 1 with gravity
    
    
    %% Inertia properties:
    massDensity = 0; % No mass, all prescribed problem.
    volume      = pi*(D/2)^3 * L;
    mass        = (massDensity)*volume;
    IL          = mass*(D/2)^2/2;
    IR          = 1/12*mass*(3*(D/2)^2+L^2);
    if cyldir == IAXIS
        I_body = [IL 0 0; 0 IR 0; 0 0 IR];
    elseif cyldir == JAXIS
        I_body = [IR 0 0; 0 IL 0; 0 0 IR];        
    else
        I_body = [IR 0 0; 0 IR 0; 0 0 IL];
    end
          
    kx          = 0.;
    ky          = kx;
    kz          = kx;
    stiff  = [kx ky kz];
    
    % Analytical body type:
    flag_forceinside = SM_TRUE; % Force inside.
    annbody_type     = RB_ANNDISC;
    annbody_nparam   = 3;
    annbody_param    = [cyldir D L];

    %% Restraints:
    % Restrained Surface:
    nrsurf      =  1;
    fis_nodes_A = [];
    fix_nodes_A = [2:nnodes+1]; % Restrained nodes in global numbering
    nfix = length(fix_nodes_A);
    
    % Restrained DOFS of Node 1:
    restnode      = 1;
    nrestcoords   = 9; % All six + 3 dofs of node 1 are restrained.
    maxrestparams = 1;
    restnodes     = restnode*ones(1,nrestcoords);
    restdofs      = [1 2 3 4 5 6 7 8 9];
    resttype      = SM_PK_FIXED*ones(1,nrestcoords);
    nparam        = 1*ones(1,nrestcoords);
    param         = zeros(maxrestparams,nrestcoords); % All dofs set to zero
        
    % Change restraint on z displacement 
    if (strcmp(kinemflag(1,:),'fixed'))
        
        % Locations in x,y,z
        param(1,IAXIS) = xpos;
        param(1,JAXIS) = ypos-1.09;
        param(1,KAXIS) = zpos;
                
    end

    %% Export Solid Model to file:
    % Values to write:
    TRANSFORMATION  = RB_EULER321;  % Use Euler 321
    nnp             =    nnodes+1;  % Total number of nodes
    if (ibd ==1)
    nel_load        =         nel;  % nel surface elements
    nel             =           1;  % 1 rigid body element
    end
    ned             =     NDIM;
    max_eltype      =       15;  % Max nodes per element type (that is one node elem)
    max_eltype_load =        2;  % Max nodes per surface element (02 is triangle)
    eltype          =       15*ones(1,nel); % Element type : One node Rigid Body
    eltype_load     =        2*ones(1,nel_load); % Surface Element type : Triangle
    kinematics_idx   =       1;  %
    
    x = XYZ(:,IAXIS); % Positions of all nodes : nnodes + 1
    y = XYZ(:,JAXIS);
    z = XYZ(:,KAXIS);
    
    IEN = [1];  % One rigid body with node 1 as connectivity
    
    
    
    if (writefile == SM_TRUE)
    
    % Write hdf5 file:
    % Mesh:
    hfilename = ['sm_body.' num2str(ibd,'%5.5d') '.h5'];
    hfile=[basedir_out hfilename];
    fprintf(1,'     Body %d hdf5 file started...\n',ibd);
    hdf5write(hfile,'BodyType',int32(BODYTYPE_RIGID))
    hdf5write(hfile,'mesh/x',x,'WriteMode','append')
    hdf5write(hfile,'mesh/y',y,'WriteMode','append')
    hdf5write(hfile,'mesh/z',z,'WriteMode','append')
    hdf5write(hfile,'mesh/nnp',int32(nnp),'WriteMode','append')
    hdf5write(hfile,'mesh/ned',int32(ned),'WriteMode','append')
    
    % Body:
    hdf5write(hfile,'body/IEN',int32(IEN'),'WriteMode','append')
    hdf5write(hfile,'body/nel',int32(nel),'WriteMode','append')
    %hdf5write(hfile,'body/nee',int32(nee),'WriteMode','append')
    hdf5write(hfile,'body/max_eltype',int32(max_eltype),'WriteMode','append')
    hdf5write(hfile,'body/DOFS_per_node',int32(DOFS_per_node),'WriteMode','append')
    
    % Body Properties:
    hdf5write(hfile,'body/eltype',int32(eltype),'WriteMode','append')
    hdf5write(hfile,'body/Mass',mass,'WriteMode','append')
    hdf5write(hfile,'body/I_body',I_body,'WriteMode','append')
    hdf5write(hfile,'body/trmatrix',int32(TRANSFORMATION),'WriteMode','append')
    hdf5write(hfile,'body/Stiffness',stiff,'WriteMode','append')
    hdf5write(hfile,'body/gravity',double(grav_vec),'WriteMode','append')
    hdf5write(hfile,'body/gravity_flag',int32(grav_flag),'WriteMode','append')
    hdf5write(hfile,'body/flag_forceinside',int32(flag_forceinside),'WriteMode','append')
    hdf5write(hfile,'body/annbody_type',int32(annbody_type),'WriteMode','append')    
    hdf5write(hfile,'body/annbody_nparam',int32(annbody_nparam),'WriteMode','append')    
    hdf5write(hfile,'body/annbody_param',double(annbody_param),'WriteMode','append')
    
    % WetSurface:
    hdf5write(hfile,'WetSurface/nel',int32(nel_load),'WriteMode','append')
    hdf5write(hfile,'WetSurface/max_eltype',int32(max_eltype_load),'WriteMode','append')
    hdf5write(hfile,'WetSurface/IEN',int32(ws_IEN'),'WriteMode','append')
    hdf5write(hfile,'WetSurface/eltype',int32(eltype_load),'WriteMode','append')
    
    % RestSurface:
    if( isempty(fix_nodes_A) )
        hdf5write(hfile,'RestSurface/fix_list',int32(0),'WriteMode','append')
    else
        hdf5write(hfile,'RestSurface/fix_list',int32(fix_nodes_A),'WriteMode','append')
    end
    hdf5write(hfile,'RestSurface/nrestsurf',int32(nrsurf),'WriteMode','append')
    hdf5write(hfile,'RestSurface/nfix',int32(nfix),'WriteMode','append')
    hdf5write(hfile,'RestSurface/kinematics_idx',int32(kinematics_idx),'WriteMode','append')
    
    
    % Write Restrained nodes, just node 1:
    hdf5write(hfile,'RestNodes/nrestcoords'  ,int32(nrestcoords),'WriteMode','append')
    hdf5write(hfile,'RestNodes/maxrestparams',int32(maxrestparams),'WriteMode','append')
    hdf5write(hfile,'RestNodes/nodes'        ,int32(restnodes),'WriteMode','append')
    hdf5write(hfile,'RestNodes/restdof'      ,int32(restdofs),'WriteMode','append')
    hdf5write(hfile,'RestNodes/restype'      ,int32(resttype),'WriteMode','append')
    hdf5write(hfile,'RestNodes/nparam'       ,int32(nparam),'WriteMode','append')
    hdf5write(hfile,'RestNodes/param'        ,param,'WriteMode','append')
    
    fprintf(1,'   Body %d hdf5 file done.\n',ibd);

    
    end
    
    trimesh(ws_IEN,XYZ(1:nnodes+1,IAXIS)+xpos,XYZ(1:nnodes+1,JAXIS)+ypos,XYZ(1:nnodes+1,KAXIS)+zpos)

    
   end
  end
end

xlabel('x')
ylabel('y')
zlabel('z')
axis equal


return