function[Data]=preProcess(Mesh, Data)

% Define quadrature table
Data.Operator.DMatrix.quad=[2 4 6];
Data.Operator.CMatrix.quad=[1 2 4];
Data.Operator.AMatrix.quad=[2 4 6];
Data.Operator.BCMatrix.quad=[3 6 8];
Data.Operator.FVector.quad=[5 5 5];
% Pre-compute indices
[Data]=setParameters(Mesh, Data, 'DMatrix');
[Data]=setParameters(Mesh, Data, 'CMatrix');
[Data]=setParameters(Mesh, Data, 'AMatrix');
[Data]=setParameters(Mesh, Data, 'BCMatrix');
[Data]=setParameters(Mesh, Data, 'FVector');

if strcmp(Data.Model.type, 'nonlinear')
    % Define quadrature table
    Data.Operator.DTensor.quad=[3 6 9];
    Data.Operator.CTensor.quad=[1 4 7];
    Data.Operator.ATensor.quad=[3 6 9];
    Data.Operator.BCTensor.quad=[3 6 8];
    % Pre-compute indices
    [Data]=setParameters(Mesh, Data, 'DTensor');
    [Data]=setParameters(Mesh, Data, 'CTensor');
    [Data]=setParameters(Mesh, Data, 'ATensor');
    [Data]=setParameters(Mesh, Data, 'BCTensor');
end

if any(contains(Data.Model.postprocess, 'flux'))
    % Define quadrature table
    Data.Operator.PMatrix.quad=[1 3 5];
    % Pre-compute indices
    [Data]=setParameters(Mesh, Data, 'PMatrix');
end

if any(contains(Data.Model.postprocess, 'error'))
    % Define quadrature table
    Data.Operator.EMatrix.quad=[2 4 6];
    % Pre-compute indices
    [Data]=setParameters(Mesh, Data, 'EMatrix');
end

% Pre-compute mesh constants
[Data]=setMeshProperties(Mesh, Data);

    function[Data]=setParameters(Mesh, Data, name)
        
        % Mesh parameters
        % Connectivity matrix
        connect=Mesh.connect_mat;
        % Finite element order
        order=Mesh.order;
        % Size of the local matrices for integration on the boundary
        s=order+1;
        % Total number of elements
        nb_elem=Mesh.nb_elem;
        % Number of nodes per element
        nb_nodes_elem=Mesh.nb_nodes_elem;
        % Number of local degrees of freedom
        nb_dof_local=Mesh.nb_dof_local;
        % Matrix linking an element to the degrees of freedom associated to it.
        Eq=Mesh.Eq;
        
        % Data parameters
        % Retrieve functions
        [Functions]=getFunctions(order);
        % Retrieve boundary basis functions
        phi_boundary=Functions.phi_boundary;
        % Retrieve basis functions
        phi=Functions.phi;
        % Jacobian matrix
        J_phi=Functions.J_phi;
        
        
        if strcmp(name, 'DMatrix') || strcmp(name, 'CMatrix') || strcmp(name, 'AMatrix')
            
            % Retrieve quadrature nodes and weights
            [points, weights] = getQuadrature(Data.Operator.(name).quad(order), 'bulk');
            
            % Number of quadrature points
            nb_points=length(weights);
            
            % Compute the indices
            T=connect';
            G=T(:);
            
            I=repmat(connect, [1 nb_nodes_elem])';
            J=repmat(G', [nb_nodes_elem 1]);
            I=[I(:) J(:)];
            
            % Initialization
            Phi=zeros(nb_nodes_elem, nb_points);
            Jr=zeros(nb_nodes_elem^2, nb_points*nb_dof_local^2);
            
            % Loop over the Gauss points
            for i = 1:nb_points
                Phi(:,i)=phi(points(i,:));
                Jphi=J_phi(points(i,:));
                Jr(:,(i-1)*nb_dof_local^2+1:i*nb_dof_local^2)=kron(Jphi, Jphi);
            end
            
            switch name
                
                case 'DMatrix'
                    Jr=kr(Phi, Phi);
                    
                case 'CMatrix'
                    Jr=Jr;
                    
                case 'AMatrix'
                    Jr=kr(Phi, Phi);
                    
            end
            
            % Store constant quantities
            Data.Operator.(name).Phi=Phi;
            Data.Operator.(name).Q=Jr;
            Data.Operator.(name).I=I;
            Data.Operator.(name).points=points;
            Data.Operator.(name).weights=weights;
            
            
        elseif strcmp(name, 'DTensor') || strcmp(name, 'CTensor') || strcmp(name, 'ATensor')
            
            % Retrieve quadrature nodes and weights
            [points, weights] = getQuadrature(Data.Operator.(name).quad(order), 'bulk');
            
            % Number of quadrature points
            nb_points=length(weights);
            
            T=connect';
            G=T(:);
            
            I=repmat(connect, [1 nb_nodes_elem^2])';
            J=repmat(G, [1 nb_nodes_elem])';
            J=reshape(J, [nb_nodes_elem^2 nb_elem]);
            J=repmat(J, [nb_nodes_elem, 1]);
            K=repmat(G', [nb_nodes_elem^2 1]);
            I=[I(:) J(:) K(:)];
            
            % Initialization
            Phi=zeros(nb_nodes_elem, nb_points);
            Jr=zeros(nb_nodes_elem^3, nb_points*nb_dof_local^2);
            
            % Loop over the Gauss points
            for i = 1:nb_points
                Phi(:,i)=phi(points(i,:));
                Jphi=J_phi(points(i,:));
                Jr(:,(i-1)*nb_dof_local^2+1:i*nb_dof_local^2)=kron(Phi(:,i), kron(Jphi, Jphi));
            end
            
            switch name
                case 'DTensor'
                    Jr=kr(Phi, kr(Phi, Phi));
                    
                case 'CTensor'
                    Jr=Jr;
                    
                case 'ATensor'
                    Jr=kr(Phi, kr(Phi, Phi));
            end
            
            % Store constant quantities
            Data.Operator.(name).Phi=Phi;
            Data.Operator.(name).Q=Jr;
            Data.Operator.(name).I=I;
            Data.Operator.(name).points=points;
            Data.Operator.(name).weights=weights;
            

        elseif contains(name, 'BC')
            
            % Retrieve quadrature nodes and weights
            [points, weights] = getQuadrature(Data.Operator.(name).quad(order), 'boundary');
            
            % Number of quadrature points
            nb_points=length(weights);
            
            % Retrieve Data
            N_BC=Data.N_BC;
            % Number of declared Neumann boundaries
            nb_N_BC=length(N_BC);
            
            % Initialization
            I=[];
            Gbc=[];
            nb_edges=0;
            Phi=zeros(s, nb_points);
            
            % Loop over the Gauss points
            for i = 1:nb_points
                Phi(:,i)=phi_boundary(points(i));
            end
            
            
            for m=1:nb_N_BC
                % Retrieve data
                data=N_BC{m};
                type=data.type;
                edges=data.edges;
                
                if strcmp(type,'convection_radiation')
                    for flag=edges
                        % Retriveve edges making up the boundary
                        Neum_edges=Mesh.BC_nodes(Mesh.BC_tag==flag,:);
                        % Number of Neumann boundary edges
                        nb_Neum_edges=size(Neum_edges, 1);
                        
                        T=Neum_edges';
                        G=T(:);
                        Gbc=[Gbc; G];
                        nb_edges=nb_edges+nb_Neum_edges;
                        
                        switch name
                            case 'BCMatrix'
                                Il=repmat(Neum_edges, [1 s])';
                                Jl=repmat(G', [s 1]);
                                I=[I; [Il(:) Jl(:)]];
                                
                            case 'BCTensor'
                                
                                Il=repmat(Neum_edges, [1 s^2])';
                                Jl=repmat(G, [1 s])';
                                Jl=reshape(Jl, [s^2 nb_Neum_edges]);
                                Jl=repmat(Jl, [s, 1]);
                                Kl=repmat(G', [s^2 1]);
                                I=[I; [Il(:) Jl(:) Kl(:)]];
                        end
                    end
                    
                end
            end
        
        switch name
            case 'BCMatrix'
                Jr=kr(Phi, Phi);
                
            case 'BCTensor'
                Jr=kr(Phi, kr(Phi, Phi));
                
        end
        
        % Store constant quantities
        Data.Operator.(name).Phi=Phi;
        Data.Operator.(name).Q=Jr;
        Data.Operator.(name).I=I;
        Data.Operator.(name).points=points;
        Data.Operator.(name).weights=weights;
        Data.Gbc=Gbc;
        Data.nb_edges=nb_edges;
        
        elseif strcmp(name, 'FVector')
            
            % Retrieve quadrature nodes and weights
            [points, weights] = getQuadrature(Data.Operator.(name).quad(order), 'bulk');
            
            % Number of quadrature points
            nb_points=length(weights);
            
            % Compute the indices
            T=connect';
            G=T(:);
            
            % Initialization
            Phi=zeros(nb_nodes_elem, nb_points);
            
            % Loop over the Gauss points
            for i = 1:nb_points
                Phi(:,i)=phi(points(i,:));
            end
            
            % Store constant quantities
            Data.Operator.(name).Phi=Phi;
            Data.Operator.(name).I=G;
            Data.Operator.(name).points=points;
            Data.Operator.(name).weights=weights;
            
        elseif strcmp(name, 'PMatrix')
            
            % Retrieve quadrature nodes and weights
            [points, weights] = getQuadrature(Data.Operator.(name).quad(order), 'bulk');
            
            % Number of quadrature points
            nb_points=length(weights);
            
            % Compute the indices
            T=connect';
            G=T(:);
            
            I=repmat(Eq, [1 nb_nodes_elem])';
            J=repmat(G', [nb_nodes_elem*nb_dof_local 1]);
            I=[I(:) J(:)];
            
            % Initialization
            n1=nb_nodes_elem^2*nb_dof_local;
            n2=nb_dof_local^2*nb_nodes_elem^2;
            Jr=zeros(n1, nb_points*n2);
            Phi=zeros(nb_nodes_elem, nb_points);
            
            % Loop over the Gauss points
            for i = 1:nb_points
                Phi(:,i)=phi(points(i,:));
                Jphi=J_phi(points(i,:));
                Jr(:,(i-1)*n2+1:i*n2)=kron(Phi(:,i)', kron(Jphi, speye(nb_nodes_elem*nb_dof_local)));
            end
            
            Jr=sparse(Jr);
            
            % Store constant quantities
            Data.Operator.(name).Phi=Phi;
            Data.Operator.(name).Q=Jr;
            Data.Operator.(name).I=I;
            Data.Operator.(name).points=points;
            Data.Operator.(name).weights=weights;
            
        elseif strcmp(name, 'EMatrix')
            
            % Retrieve quadrature nodes and weights
            [points, weights] = getQuadrature(Data.Operator.(name).quad(order), 'bulk');
            
            % Number of quadrature points
            nb_points=length(weights);
            
            % Initialization
            Phi=zeros(nb_nodes_elem, nb_points);
            Jr=zeros(nb_points, nb_dof_local*nb_nodes_elem);
            
            for i=1:nb_points
                Phi(:,i)=phi(points(i,:));
                Jphi=J_phi(points(i,:))';
                Jr(i,:)=Jphi(:);
            end
            
            % Store constant quantities
            Data.Operator.(name).Phi=Phi;
            Data.Operator.(name).Q=Jr;
            Data.Operator.(name).points=points;
            Data.Operator.(name).weights=weights;
        end
    end

    function[Data]=setMeshProperties(Mesh, Data)
        
        % Mesh parameters
        % Coordinate matrix
        coord=Mesh.coord_mat;
        % Connectivity matrix
        connect=Mesh.connect_mat;
        % Number of nodes per element
        nb_nodes_elem=Mesh.nb_nodes_elem;
        % Total number of elements
        nb_elem=Mesh.nb_elem;
        % Number of local degrees of freedom
        nb_dof_local=Mesh.nb_dof_local;
        
        T=connect';
        G=T(:);
        
        % Coordinates of the nodes of the elements
        coord_nodes=coord(G,:);
        % Vertices of the elements
        a=coord_nodes(1:nb_nodes_elem:end,:);
        b=coord_nodes(2:nb_nodes_elem:end,:);
        c=coord_nodes(3:nb_nodes_elem:end,:);
        % Determinants of Jacobian matrices
        detJK=(b(:,1)-a(:,1)).*(c(:,2)-a(:,2))-(b(:,2)-a(:,2)).*(c(:,1)-a(:,1));
        
        V=[(vecnorm(c-a,2,2).^2) -sum((c-a).*(b-a),2) -sum((c-a).*(b-a),2) (vecnorm(b-a,2,2).^2)];
        A=1./(detJK').^2.*(V');
        
        % Store constant mesh properties
        Data.MeshProperties.absDetJK=abs(detJK');
        Data.MeshProperties.A=A;
        
        % R-matrix for the projection matrix
        if any(contains(Data.Model.postprocess, 'flux'))
            JK_invT=zeros(nb_dof_local, nb_dof_local*nb_elem);
            JK_invT(:,1:nb_dof_local:nb_dof_local*nb_elem)=1./detJK'.*[c(:,2)-a(:,2) -(c(:,1)-a(:,1))]';
            JK_invT(:,2:nb_dof_local:nb_dof_local*nb_elem)=1./detJK'.*[-(b(:,2)-a(:,2)) b(:,1)-a(:,1)]';
            [R]=product(eye(nb_nodes_elem), JK_invT, nb_dof_local);
            Data.MeshProperties.R=R;
        end
        
        % B-matrix for computation of H1 error
        if any(contains(Data.Model.postprocess, 'error'))
            
            a=a'; b=b'; c=c';
            
            B1=[c(2,:)-a(2,:);
                -(b(2,:)-a(2,:))].*(1./detJK');
            B2=[-(c(1,:)-a(1,:));
                b(1,:)-a(1,:)].*(1./detJK');
            
            Data.MeshProperties.B1=B1;
            Data.MeshProperties.B2=B2;
        end
        
        
        
        function[R]=product(A, B, nc)
            
            % nc: Number of rows and columns of the matrices Bk
            [m, n]=size(A);
            [p, q]=size(B);
            
            R=zeros(m*p*n*nc, q/nc);
            
            for k=1:n
                y=repmat(A(:,k),[1 q]);
                R((k-1)*m*p*nc+1:k*m*p*nc,:)=reshape(kr(y, B), [m*p*nc, q/nc]);
            end
        end
        
        
    end
end