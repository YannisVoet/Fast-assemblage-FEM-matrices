function[V]=computePMatrix(Mesh, Data, Temp)

% computePMatrix: Function which computes the entries of the global
% projection matrix.
% Restricted to 2D problems and when all elements are of the same type
% INPUT:
% Mesh:     Structure containing all the mesh parameters
% Data:     Structure containing all the data parameters
% Temp:     Temperature vector
% OUTPUT:
% V:        Values of non-zero entries

% Mesh parameters
% Connectivity matrix
connect=Mesh.connect_mat;
% Number of nodes per element
nb_nodes_elem=Mesh.nb_nodes_elem;
% Total number of elements
nb_elem=Mesh.nb_elem;

% Operator properties
c=      Data.Coefficient.c.function;
Phi=    Data.Operator.PMatrix.Phi;
Q=      Data.Operator.PMatrix.Q;
weights=Data.Operator.PMatrix.weights;

% Mesh properties
absDetJK=Data.MeshProperties.absDetJK;
R=Data.MeshProperties.R;

T=connect';
G=T(:);

% Interpolated temperatures at Gauss points
U=reshape(Temp(G), [nb_nodes_elem, nb_elem]);
temp_int=Phi'*U;

lambda=weights'.*c(temp_int).*absDetJK;
S = Q*kr(lambda, R);
V=S(:);
end