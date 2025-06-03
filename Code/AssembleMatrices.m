function[Vmd, Vtd, Vmc, Vtc, Vma, Vta, Vmbc, Vtbc]=AssembleMatrices(Mesh, Data, U)

% AssembleMatrices: Recomputes the non-zero entries of the matrices which
% have changed
% INPUT:
% Mesh:     Structure containing all the mesh parameters
% Data:     Structure containing all the data parameters
% U:        Vector of temperatures where the matrices and tensors must be evaluated
% OUTPUT:
% Vmd:      Vector containing the non-zeros entries of the D-Matrix
% Vtd:      Vector containing the non-zeros entries of the D-Tensor
% Vmc:      Vector containing the non-zeros entries of the C-Matrix
% Vtc:      Vector containing the non-zeros entries of the C-Tensor
% Vma:      Vector containing the non-zeros entries of the A-Matrix
% Vta:      Vector containing the non-zeros entries of the A-Tensor
% Vmbc:     Vector containing the non-zeros entries of the BC-Matrix
% Vtbc:     Vector containing the non-zeros entries of the BC-Tensor

% Mesh parameters
% Coordinate matrix
coord=Mesh.coord_mat;
% Connectivity matrix
connect=Mesh.connect_mat;
% Total number of elements
nb_elem=Mesh.nb_elem;
% Number of nodes per element
nb_nodes_elem=Mesh.nb_nodes_elem;
% Order
order=Mesh.order;
% Size of the local matrices
s=order+1;

T=connect';
G=T(:);
Coefficient=Data.Coefficient;
U_reshape=reshape(U(G), [nb_nodes_elem, nb_elem]);

% Mesh properties
A=Data.MeshProperties.A;
absDetJK=Data.MeshProperties.absDetJK;


% Recompute the matrices which are temperature dependent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass matrix
if strcmp(Coefficient.d.label, 'temp')
    
    [Q, lambda]=computeLambda(Data, 'd', 'DMatrix', U_reshape, absDetJK);
    [Data.Operator.DMatrix.V]=implicitAssemblage(Q, lambda);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stiffness matrix bulk
if strcmp(Coefficient.c.label, 'temp')
    
    [Q, lambda]=computeLambda(Data, 'c', 'CMatrix', U_reshape, absDetJK);
    lambda=kr(lambda, A);
    [Data.Operator.CMatrix.V]=implicitAssemblage(Q, lambda);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass tensor
if strcmp(Coefficient.der_d.label, 'temp')
    
    [Q, lambda]=computeLambda(Data, 'der_d', 'DTensor', U_reshape, absDetJK);
    [Data.Operator.DTensor.V]=implicitAssemblage(Q, lambda);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stiffness tensor bulk
if strcmp(Coefficient.der_c.label, 'temp')
    
    [Q, lambda]=computeLambda(Data, 'der_c', 'CTensor', U_reshape, absDetJK);
    lambda=kr(lambda, A);
    [Data.Operator.CTensor.V]=implicitAssemblage(Q, lambda);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add contribution from "a" coefficient
if strcmp(Coefficient.a.label, 'temp')
    
    [Q, lambda]=computeLambda(Data, 'a', 'AMatrix', U_reshape, absDetJK);
    [Data.Operator.AMatrix.V]=implicitAssemblage(Q, lambda);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Coefficient.der_a.label, 'temp')
    
    [Q, lambda]=computeLambda(Data, 'der_a', 'ATensor', U_reshape, absDetJK);
    [Data.Operator.ATensor.V]=implicitAssemblage(Q, lambda);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contribution from Robin BC
if Data.Model.radiation
    
    % Initialization
    Vmbc=[];
    Vtbc=[];
    
    % Operator properties
    Phi_mbc=Data.Operator.BCMatrix.Phi;
    Q_mbc=Data.Operator.BCMatrix.Q;
    weights_mbc=Data.Operator.BCMatrix.weights;
    
    Phi_tbc=Data.Operator.BCTensor.Phi;
    Q_tbc=Data.Operator.BCTensor.Q;
    weights_tbc=Data.Operator.BCTensor.weights;
    
    % Retrieve Data
    N_BC=Data.N_BC;
    % Number of declared Neumann boundaries
    nb_N_BC=length(N_BC);
    
    for m=1:nb_N_BC
        % Retrieve data
        data=N_BC{m};
        type=data.type;
        edges=data.edges;
        
        if strcmp(type,'convection_radiation')
            for flag=edges
                % Convection coefficient
                hc=data.convection_coeff;
                % Radiation coefficient
                hr=data.emissivity*Data.Constants.sigma;
                % Retriveve edges making up the boundary
                Neum_edges=Mesh.BC_nodes(Mesh.BC_tag==flag,:);
                % Number of Neumann boundary edges
                nb_Neum_edges=size(Neum_edges, 1);
                
                T=Neum_edges';
                Gbc=T(:);
                
                % Coordinates of the nodes of the elements
                coord_nodes=coord(Gbc,:);
                % Temperatures at nodes
                temp=reshape(U(Gbc), [s, nb_Neum_edges]);
                temp_int_mbc=temp'*Phi_mbc;
                temp_int_tbc=temp'*Phi_tbc;
                % Endpoints of the edges
                a=coord_nodes(1:(order+1):end,:);
                b=coord_nodes(2:(order+1):end,:);
                bk=b-a;
                norm_bk=vecnorm(bk,2,2);
                
                lambdac=hc*norm_bk.*weights_mbc';
                lambdar=hr*norm_bk.*(temp_int_mbc).^3.*weights_mbc';
                
                lambda=3*hr*norm_bk.*(temp_int_tbc).^2.*weights_tbc';
                
                Kc=Q_mbc*lambdac';
                Kr=Q_mbc*lambdar';
                
                Sr=Q_tbc*lambda';
                
                % Assemblage
                Vmbc=[Vmbc; Kc(:)+Kr(:)];
                Vtbc=[Vtbc; Sr(:)];
            end
        end
    end
    
    Data.Operator.BCMatrix.V=Vmbc;
    Data.Operator.BCTensor.V=Vtbc;
end

Vmd=Data.Operator.DMatrix.V;
Vmc=Data.Operator.CMatrix.V;
Vma=Data.Operator.AMatrix.V;
Vmbc=Data.Operator.BCMatrix.V;

Vtd=Data.Operator.DTensor.V;
Vtc=Data.Operator.CTensor.V;
Vta=Data.Operator.ATensor.V;
Vtbc=Data.Operator.BCTensor.V;

    function[Q, lambda]=computeLambda(struct, coeff_name, mat_name, Ur, det)
        
        % computeLambda: Computes the matrix of coefficients
        % INPUT:
        % coeff_name:   Name of the coefficient for which the matrix or
        %               tensor representation must be computed
        % mat_name:     Name of the corresponding matrix
        % Ur:           Reordered temperatures
        % det:          Determinants of Jacobian matrices
        
        % Operator properties
        coeff=  struct.Coefficient.(coeff_name).function;
        Phi=    struct.Operator.(mat_name).Phi;
        Q=      struct.Operator.(mat_name).Q;
        weights=struct.Operator.(mat_name).weights;
        
        % Interpolated temperatures at Gauss points
        temp_int=Phi'*Ur;
        lambda=weights'.*coeff(temp_int).*det;
        
    end

    function[V]=implicitAssemblage(Q, lambda)
        
        % implicitAssemblage: Computes the non-zero entries of tensors and
        % matrices
        % INPUT:
        % Q:        Matrix containing the dependency on the basis functions
        % lambda:   Matrix containing the coefficients
        % OUTPUT:
        % V:        Vector of all non-zero entries
        
        M=Q*lambda;
        % Values of non-zero entries
        V=M(:);
        
    end

end