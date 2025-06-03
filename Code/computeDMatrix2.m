function[I, J, M]=computeDMatrix2(Mesh, Data, U)

% % computeDMatrix2: Algorithm 2 for the computation of the mass matrix

% Mesh parameters
% Coordinate matrix
coord=Mesh.coord_mat;
% Connectivity matrix
connect=Mesh.connect_mat;
% Finite element order
order=Mesh.order;
% Total number of elements
nb_elem=Mesh.nb_elem;
% Number of nodes per element
nb_nodes_elem=Mesh.nb_nodes_elem;

% Data parameters
% Coefficient
d_coeff=Data.Coefficient.d.function;

% Retrieve functions
[Functions]=getFunctions(order);
% Retrieve basis functions
phi=Functions.phi;

% Retrieve quadrature nodes and weights
switch order
    case 1
        [points, weights] = getQuadrature(2, 'bulk');
    case 2
        [points, weights] = getQuadrature(4, 'bulk');
    case 3
        [points, weights] = getQuadrature(6, 'bulk');
    otherwise
        error('Not yet implemented');
end

% Computation of the indices of nonzero entries
T=connect';
G=T(:);

I=repmat(connect, [1 nb_nodes_elem])';
J=repmat(G', [nb_nodes_elem 1]);

% Evaluation of the basis functions at the Gauss points
Phi=phi(points);

% Coordinates of the nodes of the elements
coord_nodes=coord(G,:);
% Temperatures at nodes
temp=reshape(U(G), [nb_nodes_elem, nb_elem]);
temp_int=Phi'*temp;
% Vertices of the elements
a=coord_nodes(1:nb_nodes_elem:end,:);
b=coord_nodes(2:nb_nodes_elem:end,:);
c=coord_nodes(3:nb_nodes_elem:end,:);
% Determinants of Jacobian matrices
detJK=(b(:,1)-a(:,1)).*(c(:,2)-a(:,2))-(b(:,2)-a(:,2)).*(c(:,1)-a(:,1));
lambda=weights'.*d_coeff(temp_int).*abs(detJK)';

M=kr(Phi, Phi)*lambda;
end