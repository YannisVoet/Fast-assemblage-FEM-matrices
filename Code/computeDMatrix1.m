function[I, J, M]=computeDMatrix1(Mesh, Data, U)

% computeDMatrix1: Algorithm 1 for the computation of the mass matrix

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

% Initialization
M=zeros(nb_nodes_elem, nb_nodes_elem*nb_elem);

T=connect';
G=T(:);

I=repmat(connect, [1 nb_nodes_elem])';
J=repmat(G', [nb_nodes_elem 1]);

% Evaluation of the basis functions at the Gauss points
Phi=phi(points);


for e = 1:nb_elem
    % Recovering the nodes of element e
    nodes_elem = connect(e, :);
    % Recovering the node coordinates of element e
    coord_nodes=coord(nodes_elem, :);
    % Temperatures at nodes
    temp=U(nodes_elem);
    % Computation of the local stiffness matrix
    M(:,nb_nodes_elem*(e-1)+1:nb_nodes_elem*e) = computeMassLoc(temp, coord_nodes');
end


    function[M_local] = computeMassLoc(temp, coord_nodes)
        
        % Initialization
        M_local=zeros(nb_nodes_elem);
        
        % Retrieve coordinates of vertices
        a=coord_nodes(:,1);
        b=coord_nodes(:,2);
        c=coord_nodes(:,3);
        
        Bk=[b-a c-a];
        detJK=det(Bk);
        
        % Integration using Gaussian quadrature
        for i = 1:length(weights)
            Tp=temp'*Phi(:,i);
            M_local = M_local+weights(i)*d_coeff(Tp)*Phi(:,i)*Phi(:,i)'*abs(detJK);
        end
    end
end