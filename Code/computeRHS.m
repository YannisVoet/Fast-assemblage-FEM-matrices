function[F]=computeRHS(Mesh, Data, time)

% computeRHS: Function which assembles the right-hand side vector
% INPUT:
% Mesh:             Structure containing the mesh parameters
% Data:             Structure containing the data parameters
% time:             Time at which the RHS vector must be evaluated
% OUTPUT:
% F:                Right-hand side vector


% Assemble the global vector of body forces
[Bs]=computeSource(Mesh, Data, time);
% Assemble the global vector of surface tractions
[Bn]=computeNeumann(Mesh, Data, time);

% Sum the contribution from internal heat generation and Neumann boundary
% conditions
F=Bs+Bn;

    function[Bs]=computeSource(Mesh, Data, time)
        
        % computeSource: Assembles the component corresponding to the
        % function f. Restricted to 2D problems and when all elements are
        % of the same type
        % INPUT:
        % Mesh:     Structure containing all the mesh parameters
        % Data:     Structure containing all the data parameters
        % time:     Time at which the component must be assembled
        % OUTPUT:
        % Bs:       Global vector resulting from the function f
        
        % Retrieve mesh data
        % Total number of nodes
        nb_nodes=Mesh.nb_nodes;
        % Coordinate matrix
        coord=Mesh.coord_mat;
        % Connectivity matrix
        connect=Mesh.connect_mat;
        % Number of elements
        nb_elem=Mesh.nb_elem;
        % Number of nodes per element
        nb_nodes_elem=Mesh.nb_nodes_elem;
        
        % Object properties
        f=      Data.Coefficient.f.function;
        Phi=    Data.Operator.FVector.Phi;
        I=      Data.Operator.FVector.I;
        points= Data.Operator.FVector.points;
        weights=Data.Operator.FVector.weights;
        
        % Mesh properties
        absDetJK=Data.MeshProperties.absDetJK;
        
        % Number of quadrature points
        nb_points=length(weights);
        % Initialization
        X=zeros(nb_elem, nb_points, 2);
        
        T=connect';
        G=T(:);
        
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
        
        % Determinants of Jacobian matrices
        lambda=absDetJK'.*f(X, time).*weights;
        
        Bs=Phi*lambda';
        
        % Assemblage
        Bs=accumarray(I(:), Bs(:), [nb_nodes 1]);
        
    end

    function[Bn]=computeNeumann(Mesh, Data, time)
        
        % computeNeumann: Assembles the component corresponding to the
        % Neumann boundary conditions. Restricted to 2D problems and when
        % all elements are of the same type
        % INPUT:
        % Mesh:     Structure containing all the mesh parameters
        % OUTPUT:
        % Bn:       Global vector resulting from the Neumann boundary
        %           conditions
        
        
        % Retrieve mesh parameters
        % Coordinate matrix
        coord=Mesh.coord_mat;
        % Order
        order=Mesh.order;
        % Total number of nodes
        nb_nodes=Mesh.nb_nodes;
        % Size of the local matrices
        s=order+1;
        
        % Initialization
        Bn=zeros(nb_nodes,1);
        
        % Operator properties
        Phi=    Data.Operator.BCMatrix.Phi;
        points= Data.Operator.BCMatrix.points;
        weights=Data.Operator.BCMatrix.weights;
        
        % Number of quadrature points
        nb_points=length(weights);
        
        % Retrieve Data
        N_BC=Data.N_BC;
        % Number of declared Neumann boundaries
        nb_N_BC=length(N_BC);
        
        for m=1:nb_N_BC
            % Retrieve data
            data=N_BC{m};
            type=data.type;
            edges=data.edges;
            
            
            switch type
                case 'flux'
                    for flag=edges
                        % Flux function
                        f=data.function;
                        % Retriveve edges making up the boundary
                        Neum_edges=Mesh.BC_nodes(Mesh.BC_tag==flag,:);
                        % Number of Neumann boundary edges
                        nb_Neum_edges=size(Neum_edges, 1);
                        
                        % Initialization
                        X=zeros(nb_Neum_edges, nb_points, 2);
                        
                        C=Neum_edges';
                        G=C(:);
                        I=C;
                        
                        % Coordinates of the nodes of the elements
                        coord_nodes=coord(G,:);
                        % Endpoints of the edges
                        a=coord_nodes(1:s:end,:);
                        b=coord_nodes(2:s:end,:);
                        bk=b-a;
                        norm_bk=vecnorm(bk,2,2);
                        
                        % Interpolated coordinates
                        P=a(:)+(b(:)-a(:))*points';
                        X(:,:,1)=P(1:nb_Neum_edges,:);
                        X(:,:,2)=P(nb_Neum_edges+1:end,:);
                        
                        lambda=norm_bk.*f(X, time).*weights';
                        Tl=Phi*lambda';
                        
                        % Assemblage
                        Bn=Bn+accumarray(I(:), Tl(:), [nb_nodes 1]);
                    end
                    
                case 'convection_radiation'
                    for flag=edges
                        % Ambient temperature
                        Ta=data.ambient_temp;
                        % Convection coefficient
                        hc=data.convection_coeff;
                        % Radiation coefficient
                        hr=data.emissivity*Data.Constants.sigma;
                        % Retriveve edges making up the boundary
                        Neum_edges=Mesh.BC_nodes(Mesh.BC_tag==flag,:);
                        % Number of Neumann boundary edges
                        nb_Neum_edges=size(Neum_edges, 1);
                        
                        % Initialization
                        X=zeros(nb_Neum_edges, nb_points, 2);
                        
                        C=Neum_edges';
                        G=C(:);
                        I=C;
                        
                        % Coordinates of the nodes of the elements
                        coord_nodes=coord(G,:);
                        % Endpoints of the edges
                        a=coord_nodes(1:s:end,:);
                        b=coord_nodes(2:s:end,:);
                        bk=b-a;
                        norm_bk=vecnorm(bk,2,2);
                        
                        % Interpolated coordinates
                        P=a(:)+(b(:)-a(:))*points';
                        X(:,:,1)=P(1:nb_Neum_edges,:);
                        X(:,:,2)=P(nb_Neum_edges+1:end,:);
                        
                        lambdac=hc*norm_bk.*Ta(X, time).*weights';
                        lambdar=hr*norm_bk.*(Ta(X, time).^4).*weights';
                        
                        Fc=Phi*lambdac';
                        Fr=Phi*lambdar';
                        
                        % Assemblage
                        Tl=accumarray(I(:), Fc(:)+Fr(:), [nb_nodes 1]);
                        Bn=Bn+Tl;
                        
                    end
            end
        end
    end
end