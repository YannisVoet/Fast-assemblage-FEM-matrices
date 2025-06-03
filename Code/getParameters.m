function[Mesh]=getParameters(Mesh)

% getParameters: Initializes mesh parameters
% INPUT:
% Mesh: Structure containing the mesh geometry
% OUTPUT:
% Mesh: Updated structure containing all mesh parameters needed to run the code

% Retrieve mesh parameters
% Coordinate matrix
coord=Mesh.coord_mat;
% Connectivity matrix
connect=Mesh.connect_mat;
% Boundary nodes
BC_nodes=Mesh.BC_nodes;

% Geometry [m] 
x_min=min(coord(:,1));
x_max=max(coord(:,1));
y_min=min(coord(:,2));
y_max=max(coord(:,2));
% Maximum dimension [m]
L_max=max([abs(x_max-x_min) abs(y_max-y_min)]);
% Number of elements
nb_elem=size(connect,1);
% Number of nodes per element
nb_nodes_elem=size(connect,2);
% Number of nodes
nb_nodes=size(coord,1);
% Local number of degrees of freedom
nb_dof_local=size(coord,2);
% Total number of degrees of freedom
nb_dof=nb_dof_local*nb_nodes;
% Number of boundary nodes
nb_BC_nodes=size(BC_nodes,1);

% Order of polynomials
% It is assumed all elements of the mesh are of same type
order_map=containers.Map([3 6 10], [1 2 3]);
order=order_map(nb_nodes_elem);

% Initialization
% Matrix of equation numbers   
Eq=zeros(nb_elem,nb_dof_local*nb_nodes_elem);
% Matrix of equation numbers for boundary nodes
Eq_BC=zeros(nb_BC_nodes,nb_dof_local*order);
% Matrix containing the sets of degrees of freedom for each node
Dof_set=zeros(nb_dof_local, nb_nodes);

% Iteration over the dimension
for k=1:nb_dof_local
    d1=k:nb_dof_local:nb_dof_local*nb_nodes_elem;
    d2=k:nb_dof_local:nb_dof_local*(order+1);
    Eq(:,d1)=nb_dof_local*connect-nb_dof_local+k;
    Eq_BC(:,d2)=nb_dof_local*BC_nodes-nb_dof_local+k;
    Dof_set(k,:)=nb_dof_local*(1:nb_nodes)-nb_dof_local+k;
end

% Update mesh structure
Mesh.frame=[x_min y_min; x_max y_max];
Mesh.L_max=L_max;
Mesh.nb_elem=nb_elem;
Mesh.nb_nodes_elem=nb_nodes_elem;
Mesh.nb_nodes=nb_nodes;
Mesh.nb_dof_local=nb_dof_local;
Mesh.nb_dof=nb_dof;
Mesh.order=order;
Mesh.Eq=Eq;
Mesh.Eq_BC=Eq_BC;
Mesh.nb_BC_nodes=nb_BC_nodes;
Mesh.Dof_set=Dof_set;
end