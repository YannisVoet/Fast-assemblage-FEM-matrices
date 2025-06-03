function[I, J, K]=computeCMatrix1(Mesh, Data, Temp)

% computeCMatrix1: Algorithm 1 for the computation of the conductivity
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
% Number of local degrees of freedom
nb_dof_local=Mesh.nb_dof_local;

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

% Initialization
K=zeros(nb_nodes_elem, nb_nodes_elem*nb_elem);

% Computation of the indices of nonzero entries
T=connect';
G=T(:);

I=repmat(connect, [1 nb_nodes_elem])';
J=repmat(G', [nb_nodes_elem 1]);

% Evaluation of the basis functions and Jacobian matrices at the Gauss points
Phi=phi(points);
JPhi=J_phi(points);

for e = 1:nb_elem
    % Recovering the nodes of element e
    nodes_elem = connect(e, :);
    % Recovering the node coordinates of element e
    coord_nodes=coord(nodes_elem, :);
    % Temperatures at nodes
    temp=Temp(nodes_elem);
    % Computation of the local stiffness matrix
    K(:,nb_nodes_elem*(e-1)+1:nb_nodes_elem*e) = computeStiffnessLoc(temp, coord_nodes');
end

    function [K_local] = computeStiffnessLoc(temp, coord_nodes)
        
        % Retrieve the coordinates of the vertices
        a=coord_nodes(:,1);
        b=coord_nodes(:,2);
        c=coord_nodes(:,3);
        
        Bk=[b-a c-a];
        detJK=det(Bk);
        
        % Initialization
        K_local=zeros(nb_nodes_elem);
        
        % Integration using Gaussian quadrature
        for i = 1:length(weights)
            JN= JPhi(:, (i-1)*nb_dof_local+1:i*nb_dof_local)/Bk;
            Tp=temp'*Phi(:,i);
            K_local=K_local+weights(i)*JN*c_coeff(Tp)*JN'*abs(detJK);
        end
    end
end