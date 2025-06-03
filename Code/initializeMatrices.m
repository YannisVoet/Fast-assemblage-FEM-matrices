function[Data, varargout]=initializeMatrices(Mesh, Data, U)

% initializeMatrices: Initializes the matrices and corrects the set of
% indices to begin the computations
% INPUT:
% Mesh:     Structure containing all the mesh parameters
% Data:     Structure containing all the data parameters
% U:        Vector of temperatures where the matrices and tensors must be evaluated
% OUTPUT:
% Data:     Updated data structure containing all the parameters
%           implicitly defining the operators 
% OPTIONAL OUTPUT:
% mu:       Multiplication mode for mutiplication between tensors and
%           vectors

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

% Initialize indices
Im=Data.Operator.DMatrix.I;
Imbc=Data.Operator.BCMatrix.I;

if strcmp(Data.Model.type, 'nonlinear')
    It=Data.Operator.DTensor.I;
end

% Initialize the matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass matrix
[Q, lambda]=computeLambda(Data, 'd', 'DMatrix', U_reshape, absDetJK);
[Data.Operator.DMatrix.V]=implicitAssemblage(Q, lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stiffness matrix bulk
[Q, lambda]=computeLambda(Data, 'c', 'CMatrix', U_reshape, absDetJK);
lambda=kr(lambda, A);
[Data.Operator.CMatrix.V]=implicitAssemblage(Q, lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass tensor
if ~strcmp(Coefficient.der_d.label, 'zero')
    
    [Q, lambda]=computeLambda(Data, 'der_d', 'DTensor', U_reshape, absDetJK);
    [Data.Operator.DTensor.V]=implicitAssemblage(Q, lambda);
    
else
    Data.Operator.DTensor.V=[];
    Data.Operator.DTensor.I=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stiffness tensor bulk
if ~strcmp(Coefficient.der_c.label, 'zero')
    
    [Q, lambda]=computeLambda(Data, 'der_c', 'CTensor', U_reshape, absDetJK);
    lambda=kr(lambda, A);
    [Data.Operator.CTensor.V]=implicitAssemblage(Q, lambda);
    
else
    Data.Operator.CTensor.V=[];
    Data.Operator.CTensor.I=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add contribution from "a" coefficient
if ~strcmp(Coefficient.a.label, 'zero')
    
    [Q, lambda]=computeLambda(Data, 'a', 'AMatrix', U_reshape, absDetJK);
    [Data.Operator.AMatrix.V]=implicitAssemblage(Q, lambda);
    
else
    Data.Operator.AMatrix.V=[];
    Data.Operator.AMatrix.I=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(Coefficient.der_a.label, 'zero')
    
    [Q, lambda]=computeLambda(Data, 'der_a', 'ATensor', U_reshape, absDetJK);
    [Data.Operator.ATensor.V]=implicitAssemblage(Q, lambda);
    
else
    Data.Operator.ATensor.V=[];
    Data.Operator.ATensor.I=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contribution from Robin BC
if Data.Model.convection || Data.Model.radiation
    
    % Initialization
    Vmbc=[];
    
    % Operator properties
    Phi_mbc=Data.Operator.BCMatrix.Phi;
    Q_mbc=Data.Operator.BCMatrix.Q;
    weights_mbc=Data.Operator.BCMatrix.weights;
    
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
                % Endpoints of the edges
                a=coord_nodes(1:(order+1):end,:);
                b=coord_nodes(2:(order+1):end,:);
                bk=b-a;
                norm_bk=vecnorm(bk,2,2);
                
                lambdac=hc*norm_bk.*weights_mbc';
                lambdar=hr*norm_bk.*(temp_int_mbc).^3.*weights_mbc';
                
                Kc=Q_mbc*lambdac';
                Kr=Q_mbc*lambdar';
                
                % Assemblage
                Vmbc=[Vmbc; Kc(:)+Kr(:)];
            end
        end
    end
    
    Data.Operator.BCMatrix.V=Vmbc;
    
else
    Data.Operator.BCMatrix.V=[];
    Data.Operator.BCMatrix.I=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contribution from Robin BC
if Data.Model.radiation
    
    % Initialization
    Vtbc=[];
    
    % Operator properties
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
                temp_int_tbc=temp'*Phi_tbc;
                % Endpoints of the edges
                a=coord_nodes(1:(order+1):end,:);
                b=coord_nodes((order+1):(order+1):end,:);
                bk=b-a;
                norm_bk=vecnorm(bk,2,2);
                
                lambda=3*hr*norm_bk.*(temp_int_tbc).^2.*weights_tbc';
                
                Sr=Q_tbc*lambda';
                
                % Assemblage
                Vtbc=[Vtbc; Sr(:)];
            end
        end
    end
    
    Data.Operator.BCTensor.V=Vtbc;
    
else
    Data.Operator.BCTensor.V=[];
    Data.Operator.BCTensor.I=[];
end

if strcmp(Data.Model.type, 'nonlinear')
    % Rework the indices for computations
    % Compute indices corresponding to multiplication between
    % tensors and vectors
    % mu:       Multiplication mode
    % res:      Vector defining how the other modes of the tensor should be
    %           arranged after the product is carried out
    mu=[2];
    res=[1 3];
    
    varargout{1}=mu;
    
    Data.G=G;

    
    % Concatenation
    % Vector
    Ivec = [Data.G;
            Data.Gbc];
    
    Ivec = [Ivec;
        (1:Mesh.nb_nodes)'];
    
    % Matrix
    Imat=[Im;
         Imbc];
    
    [Data.Operation.Iv, ~, Data.Operation.Rv]=unique(Ivec, 'rows');
    [Data.Operation.Im, ~, Data.Operation.Rm]=unique(Imat, 'rows');
end




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