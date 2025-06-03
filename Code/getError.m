function[e] = getError(Mesh, Data, Solution, Uex, gradUex1, gradUex2)

% Retrieve mesh parameters
% Coordinate matrix
coord=Mesh.coord_mat;
% Connectivity matrix
connect=Mesh.connect_mat;
% Number of elements
nb_elem=Mesh.nb_elem;
% Number of nodes per element
nb_nodes_elem=Mesh.nb_nodes_elem;
% Number of local degrees of freedom
nb_dof_local=Mesh.nb_dof_local;

T=connect';
G=T(:);

% Retrieve solution
U=Solution.U;

% Retrieve discretization parameters
N=Data.Discretization.N;
delta_t=Data.Discretization.delta_t;
time_vec=Data.Discretization.time_vec;

% Object properties
Phi=    Data.Operator.EMatrix.Phi;
Q=      Data.Operator.EMatrix.Q;
points= Data.Operator.EMatrix.points;
weights=Data.Operator.EMatrix.weights;

% Mesh properties
absDetJK=Data.MeshProperties.absDetJK;
B1=Data.MeshProperties.B1;
B2=Data.MeshProperties.B2;

% Number of quadrature points
nb_points=length(weights);
% Initialization
C=zeros(nb_nodes_elem*nb_dof_local,nb_dof_local*nb_elem);
X=zeros(nb_elem, nb_points, 2);

% Coordinates of the nodes of the elements
coord_nodes=coord(G,:);
% Vertices of the elements
a=coord_nodes(1:nb_nodes_elem:end,:);
b=coord_nodes(2:nb_nodes_elem:end,:);
c=coord_nodes(3:nb_nodes_elem:end,:);
% Interpolated coordinates
P=a(:)+[b(:)-a(:) c(:)-a(:)]*points';
X(:,:,1)=P(1:nb_elem,:);
X(:,:,2)=P(nb_elem+1:end,:);


e=0;
for k= 1:N+1
    
    U_reshape=reshape(U(G, k), [nb_nodes_elem, nb_elem]);
    
    C(:,1:nb_dof_local:end)=kr(U_reshape,B1);
    C(:,2:nb_dof_local:end)=kr(U_reshape,B2);
    
    gradUh=Q*C;
    
    % Gradient of numerical solution
    gradUh1=gradUh(:,1:nb_dof_local:end);
    gradUh2=gradUh(:,2:nb_dof_local:end);
    
    % Gradient of exact solution
    gradU1=gradUex1(X, time_vec(k))';
    gradU2=gradUex2(X, time_vec(k))';
    
    grad=(gradU1-gradUh1).^2+(gradU2-gradUh2).^2;
    
    % Exact solution
    Uexact=Uex(X, time_vec(k))';
    % Numerical solution
    Uint=Phi'*U_reshape;
    
    L2_err=weights*(((Uexact-Uint).^2).*absDetJK);
    H1_err=weights*(grad.*absDetJK);
    
    e_local=sum(L2_err+H1_err);
    
    e = e + delta_t*e_local;
end
e = sqrt(e);
end