function[Solution]=getFlux(Mesh, Data, Solution)

% getFlux: Function used to compute heat fluxes once the temperatures
% are known.
% INPUT:
% Mesh:     Structure containing all mesh parameters
% Data:     Structure containing all data parameters
% Solution: Structure containing the computed solutions
% OUTPUT:
% Solution: Updated solution structure including heat fluxes

% Mesh parameters
% Number of local degrees of freedom
nb_dof_local=Mesh.nb_dof_local;
% Total number of degrees of freedom
nb_dof=Mesh.nb_dof;
% Number of nodes
nb_nodes=Mesh.nb_nodes;
% Data parameters
% Number of sub-intervals in time
N=Data.Discretization.N;
% Retrieve computed solution
U=Solution.U;

% Compute mass matrix
[V]=computeDMatrix(Data);
% Indices of the mass matrix
I=Data.Operator.DMatrix.I;
% Assemblage of the mass matrix
M=sparse(I(:,1), I(:,2), V, nb_nodes, nb_nodes);
% Expansion to full size
M=kron(M, speye(nb_dof_local));

% Sparse Cholesky factorization of the mass matrix
p=dissect(M);
L=chol(M(p,p), 'lower');
% Indices of the projection matrix
I=Data.Operator.PMatrix.I;

% Initialization:
Q=zeros(nb_dof, N+1);

for k=1:N+1
    [V]=computePMatrix(Mesh, Data, U(:,k));
    % Matrix-vector multiplication
    V_new=U(I(:,2),k).*V;
    V=-accumarray(I(:,1), V_new);
    Q(p,k)=L'\(L\V(p,:));
end

% Sort the array of fluxes
Q_sort=zeros(nb_nodes,nb_dof_local*(N+1));
for k=1:nb_dof_local
    Q_sort(:,k:nb_dof_local:(nb_dof_local*(N+1)))=Q(k:nb_dof_local:nb_dof,:);
end

% Update solution structure
Solution.Q_sort=Q_sort;
end