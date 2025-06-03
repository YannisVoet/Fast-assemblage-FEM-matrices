function[I, J, C]=computeCMatrix2(Mesh, Data, U)

% computeCMatrix2: Algorithm 2 for the computation of the conductivity
% matrix

% Mesh parameters
% Coordinate matrix
coord=Mesh.coord_mat;
% Connectivity matrix
connect=Mesh.connect_mat;
% Finite element order
order=Mesh.order;
% Number of nodes per element
nb_nodes_elem=Mesh.nb_nodes_elem;
% Total number of elements
nb_elem=Mesh.nb_elem;

% Data parameters
% c coefficient
c_coeff=Data.Coefficient.c.function;

% Retrieve functions
[Functions]=getFunctions(order);
% Retrieve basis functions
phi=Functions.phi;
% Jacobian matrix
J_phi=Functions.J_phi;

% Retrieve quadrature nodes and weights
switch order
    case 1
        [points, weights] = getQuadrature(1, 'bulk');
    case 2
        [points, weights] = getQuadrature(2, 'bulk');
    case 3
        [points, weights] = getQuadrature(4, 'bulk');
    otherwise
        error('Not yet implemented');
end

% Number of quadrature points
nb_points=length(weights);

% Computation of the indices of nonzero entries
T=connect';
G=T(:);

I=repmat(connect, [1 nb_nodes_elem])';
J=repmat(G', [nb_nodes_elem 1]);

% Evaluation of the basis functions and Jacobian matrices at the Gauss points
Phi=phi(points);
JPhi=J_phi(points);
JPhi=mat2cell(JPhi, nb_nodes_elem, 2*ones(1,nb_points));

Q=cellfun(@(x,y) kron(x,y), JPhi, JPhi, 'UniformOutput', false);
Q=cell2mat(Q);

% Coordinates of the nodes of the elements
coord_nodes=coord(G,:);
% Interpolated temperatures at Gauss points
temp_int=Phi'*reshape(U(G), [nb_nodes_elem, nb_elem]);
% Vertices of the elements
a=coord_nodes(1:nb_nodes_elem:end,:);
b=coord_nodes(2:nb_nodes_elem:end,:);
c=coord_nodes(3:nb_nodes_elem:end,:);
% Determinants of Jacobian matrices
detJK=(b(:,1)-a(:,1)).*(c(:,2)-a(:,2))-(b(:,2)-a(:,2)).*(c(:,1)-a(:,1));

lambda=weights'.*c_coeff(temp_int).*abs(detJK');
V=[(vecnorm(c-a,2,2).^2) -sum((c-a).*(b-a),2) -sum((c-a).*(b-a),2) (vecnorm(b-a,2,2).^2)];
W=1./(detJK').^2.*(V');

C=Q*kr(lambda, W);
end