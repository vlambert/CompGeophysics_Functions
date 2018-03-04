
%% % % % % % % % % % % % % % % % % % % % % % % % % %
% Generate triangular element mesh for rectangle   %
%          with an elliptical crack                %
%       Valere Lambert + Ludo Gil, 2017            %
% % % % % % % % % % % % % % % % % % % % % % % % % %%
clear all;close all

% Make rectangle
R = [3,4,0,1,1,0,0,0,1,1]';

% Make ellipse
ratio = 0.01; % ratio of crack area to bulk
ecc = 20;     % eccentricity
% axis lengths
A2 = sqrt(ratio / (pi *ecc ) ); 
A1 = ecc*A2;
angle = pi/4 ;
E = [4,0.5,0.5,A1,A2,angle]';

% pad to make same size as rectangle
E = [E;zeros(length(R)-length(E),1)];
ns = (char('R1','C1'))';
sf = 'R1-C1';

% Create geometry and generate mesh
geom = [R1,C1];
gd = decsg(geom,sf,ns);
pdem = createpde();
gm = geometryFromEdges(pdem, gd);
generateMesh(pdem,'Hmax',0.2);
pdeplot(pdem)
axis equal

% Save nodal connectivity of mesh
nodeLoc = pdem.Mesh.Nodes;
N_nodes=pdem.Mesh.Elements - 1; % 1 -> 0 for C++
file=fopen('connectivity_0_2_alpha45.dat','w');
fprintf(file,'%6.4f,%6.4f\n',nodeLoc);
fprintf(file,'connectivity\n');
fprintf(file,'%1.i,%1.i,%1.i\n',N_nodes);
fclose(file);

% Calculate area of each element
E_nodes = pdem.Mesh.Elements;
Area = zeros(size(E_nodes,2),1);
for i = 1:size(E_nodes,2)
    p1 = nodeLoc(:,E_nodes(1,i));
    p2 = nodeLoc(:,E_nodes(2,i));
    p3 = nodeLoc(:,E_nodes(3,i));
    v1 = p2 - p1;
    v2 = p3 - p1;
    vd = dot(v1,v2)*v1/dot(v1,v1);
    vd2 = v2 -vd;
    Area(i) = norm(vd2)*norm(v1)/2;
end

% Save area of each element
file=fopen('Areas_0_2_alpha45.dat','w');
fprintf(file,'%6.4f\n',Area);
fclose(file);