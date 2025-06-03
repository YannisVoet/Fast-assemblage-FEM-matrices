function [Mesh] = readGMSH(file_name)

% readGMSH: Reads a GMSH mesh file (.msh) and constructs the coordinate,
% connectivity, boundary nodes, material tag and boundary tag matrices
% INPUT:
% file_name: GMSH mesh file (.msh) 
% OUTPUT:
% Mesh: Structure containing all fields related to the geometry, boundary
% conditions and material tags
%   coord: Matrix of node coordinates
%   connectivites: Connectivity matrix of the system
%   boundary_nodes: Matrix containing the boundary nodes (both displacements
%   and external forces)
%   material_tag: Vector containing the material tag of the elements
%   boundary_tag: Vector containing the boundary tags of the elements 
%   forming boundaries (enables to make the distinction between Dirichlet
%   and Neumann boundary conditions)

fileID = fopen(file_name);
text = textscan(fileID, '%s', 'delimiter', '\n');
fclose(fileID);

text=text{1};
% Construction of the coordinate matrix
nodes_start = findPos(text, '$Nodes') + 1;

fileID = fopen(file_name);
% Node coordinates
coord = textscan(fileID, '%f %f %f %f', 'HeaderLines', nodes_start);
% Elements
connect = textscan(fileID, [repmat('%f', 1, 20) '%*[^\n]'], 'HeaderLines', 3, 'EndOfLine', '\n', 'CollectOutput', 1);
fclose(fileID);

% Coordinate matrix
coord=cell2mat(coord);
% Ignore z component
coord=coord(:,2:end-1);

% Connectivity matrix
connect=cell2mat(connect);

I=sum(~isnan(connect), 2);
v=unique(I);

% Boundary elements
B=connect(I==v(1), 1:v(1));
% Elements
E=connect(I==v(2), 1:v(2));

% Boundary tags
B_tags=B(:,4);
% Element tags
E_tags=E(:,4);

% Boundary elements connectivity matrix
B_connect=B(:,6:end);
% Elements connectivity matrix
E_connect=E(:,6:end);

% Initialize mesh structure
Mesh.coord_mat=coord;
Mesh.BC_nodes=B_connect;
Mesh.connect_mat=E_connect;
Mesh.BC_tag=B_tags;
Mesh.material_tag=E_tags;

% Identification of the type of element in GMSH
% 2 -> T3
% 9 -> T6
% 21 -> T10
% 1 -> line
% 8 -> quadratic curve
% 26 -> cubic curve
end


function[start]=findPos(text, string)
lines = size(text, 1);
start = 1;
for i = 1:lines
    if strcmp(text{i}, string)
        start = i;
        break
    end 
end
end