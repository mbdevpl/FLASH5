function [XYZ,ELEM,nnodes,nel]=createCylinderMesh(dir,D,L,dL,dtheta)

IAXIS = 1;
JAXIS = 2;
KAXIS = 3;


% Position in the length:
xmin=-L/2;
xmax= L/2;
% Number of points in the length:
npL = ceil(L/dL)+1;
% Number of segments in the length:
nsegL = npL-1;
% Segment length:
dl  = L/nsegL; 

% Radius, circumference and npoints in cylinder:
r      =                 D/2;
circmf =              2*pi*r;
npt    = ceil(circmf/dtheta);

% Check if number is odd -> make even:
if mod(npt,2) ~= 0
    npt=npt+1;
end

% Delta theta:
dth=2*pi/npt;

% 1D arrays:
theta = linspace(0,2*pi-dth,npt);
X_circ = r*cos(theta);
Y_circ = r*sin(theta);
Z_L    = linspace(xmin,xmax,npL);

% Nodes and Elements:
nnodes = npt*npL+1;
nel    = 2*npt*nsegL;
XYZ  = zeros(nnodes,3);
ELEM = zeros(nel,3);

% Nodes
inod = 1;
for iL=1:npL
  for it=1:npt
    inod = inod+1;
    XYZ(inod,:) = [X_circ(it) Y_circ(it) Z_L(iL)];  
  end
end

% Elems
iel   = 0;
itvec = [1:npt 1];
for iL=1:npL-1
  for it=1:npt
      
      iel = iel+1;  
      % First element  
      ELEM(iel,:) = [ (iL-1)*npt+itvec(it) iL*npt+itvec(it+1) iL*npt+itvec(it)  ];
    
      iel = iel+1;
      % Second element
      ELEM(iel,:) = [ (iL-1)*npt+itvec(it) (iL-1)*npt+itvec(it+1) iL*npt+itvec(it+1)  ];
    
  end
end
ELEM = ELEM + 1;

nnodes = nnodes - 1; % Dont count body origin node.

% trimesh(ELEM,XYZ(:,1),XYZ(:,2),XYZ(:,3))
% grid on
% axis equal

disp(['Nodes =' num2str(nnodes)])
disp(['Elems =' num2str(nel)])

if dir == IAXIS
    vec=XYZ(:,JAXIS);
    XYZ(:,JAXIS) = XYZ(:,IAXIS);
    XYZ(:,IAXIS) = XYZ(:,KAXIS);
    XYZ(:,KAXIS) = vec;
elseif dir == JAXIS
    vec=XYZ(:,IAXIS);
    XYZ(:,IAXIS) = XYZ(:,JAXIS);
    XYZ(:,JAXIS) = XYZ(:,KAXIS);
    XYZ(:,KAXIS) = vec;
end

return

end