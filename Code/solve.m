function [Solution] = solve(Mesh, Data)

% solve: Numerical approximation of the solution to the system of 
% differential equations M*d1_x+K*x=f
% INPUT:
% Mesh:     Structure containing all the mesh parameters
% Data:     Structure containing all the data parameters
% OUTPUT:
% Solution: Structure containing the solution and its time derivative

% Initialize solution and initial conditions
[Solution]=evaluateInitialConditions(Mesh, Data);
% Compute boundary conditions
[Solution]=computeBC(Mesh, Data, Solution);
% Solves the problem using the Newmark scheme
[Solution, Data]=scheme(Mesh, Data, Solution);
% Include post-processing if needed
if any(contains(Data.Model.postprocess, 'flux'))
    [Solution]=getFlux(Mesh, Data, Solution);
end

    function[Solution]=evaluateInitialConditions(Mesh, Data)
        
        % evaluateInitialConditions: Function which initializes the solution
        % matrices with the initial conditions
        % INPUT:
        % Mesh:     Structure containing all the mesh parameters
        % Data:     Structure containing all the data parameters
        % OUTPUT:
        % Solution: Structure containing the solution and its time derivative
        
        % Retrieve mesh parameters
        % Coordinate matrix
        coord=Mesh.coord_mat;
        % Total number of nodes
        nb_nodes=Mesh.nb_nodes;
        
        % Retrieve data parameters
        % Number of sub-intervals in time
        N=Data.Discretization.N;
        % Initialization
        U=zeros(nb_nodes,N+1);
        d1_U=zeros(nb_nodes,N+1);
        
        % Initial conditions
        u0=Data.u0;
        
        % Evaluation
        for k=1:nb_nodes
            U(k,1)=u0(coord(k,:));
        end
        
        % Create solution structure
        Solution.U=U;
        Solution.d1_U=d1_U;
    end

    function[Solution]=computeBC(Mesh, Data, Solution)
        
        % computeBC: Function which computes the solution at Dirichlet BC nodes
        % INPUT:
        % Mesh:     Structure containing all the mesh parameters
        % Data:     Structure containing all the data parameters
        % Solution: Structure containing the solution and its time derivative 
        %           initialized with the initial conditions
        % OUTPUT:
        % Solution: Structure containing the solution and its time derivative 
        %           completed with the evaluation at Dirichlet nodes
        
        % Retrieve mesh parameters
        % Coordinate matrix
        coord=Mesh.coord_mat;
        
        % Retrieve data
        % Boundary condition structure
        D_BC=Data.D_BC;
        % Set of nodes associated to Dirichlet BC
        set_D_nodes=Data.set_D_nodes;
        % Discretization parameters
        Discretization=Data.Discretization;
        % Number of sub-intervals in time
        N=Discretization.N;
        % Vector containing the discrete times
        time_vec=Discretization.time_vec;
        % Number of boundaries with Dirichlet BC
        nb_D_BC=length(D_BC);
        
        % Retrieve solutions
        U=Solution.U;
        d1_U=Solution.d1_U;
        
        for i=1:nb_D_BC
            data=D_BC{i};
            % Retrieve function, first time derivative
            g=data.function;
            d1_g=data.derivative;
            nodes=set_D_nodes{i};
            
            % Evaluation
            for j=2:N+1
                for k=1:length(nodes)
                    U(nodes(k),j)=g(coord(nodes(k),:), time_vec(j));
                    d1_U(nodes(k),j)=d1_g(coord(nodes(k),:), time_vec(j));
                end
            end
        end
        
        % Update solution
        Solution.U=U;
        Solution.d1_U=d1_U;
    end

    function[Solution, Data]=scheme(Mesh, Data, Solution)
        
        % scheme: Numerical approximation of the solution to the system of
        % differential equations M*d1_x+K*x=f using the theta method
        % INPUT:
        % Mesh:     Structure containing all the mesh parameters
        % Data:     Structure containing all the data parameters
        % Solution: Structure containing the solution and its time derivative
        % OUTPUT:
        % Solution: Structure containing the solution and its time derivative
        
        % Mesh parameters
        nb_nodes=Mesh.nb_nodes;
        
        % Discretization parameters
        Discretization=Data.Discretization;
        % Discretized time vector
        T=Discretization.time_vec;
        % Time-step
        delta_t=Discretization.delta_t;
        % Number of sub-intervals in time
        N=Discretization.N;     
        
        % Retrieve Interior degrees of freedom
        Nodes_F=Data.Nodes_F;
        Nodes_D=Data.Nodes_D;
        
        % Retrieve solutions
        U=Solution.U;
        d1_U=Solution.d1_U;
        
        
        if strcmp(Data.Model.type, 'linear')
            
            % Initialize indices
            [Data]=preProcess(Mesh, Data);
            
            % Initialize the matrices
            [Data]=initializeMatrices(Mesh, Data, U(:,1));

            Imd=Data.Operator.DMatrix.I;
            Imbc=Data.Operator.BCMatrix.I;
            
            Vmd=Data.Operator.DMatrix.V;
            Vmc=Data.Operator.CMatrix.V;
            Vma=Data.Operator.AMatrix.V;
            Vmbc=Data.Operator.BCMatrix.V;
            
            Vmk=[sumObject(Vmc, Vma); Vmbc];
            Imk=[Imd; Imbc];
            
            % Assemblage of the mass matrix
            M=sparse(Imd(:,1), Imd(:,2), Vmd, nb_nodes, nb_nodes);
            % Assemblage of the stiffness matrix
            K=sparse(Imk(:,1), Imk(:,2), Vmk, nb_nodes, nb_nodes);
            
            % Partitioning of K and M in block matrices
            K11=K(Nodes_F,Nodes_F);
            K12=K(Nodes_F,Nodes_D);
            
            M11=M(Nodes_F,Nodes_F);
            M12=M(Nodes_F,Nodes_D);
            
            
            switch Data.Solver.ODE.name
                case 'explicit_euler'
                    B=M11-delta_t*K11;
                    % Nested dissection algorithm
                    p=dissect(M11);
                    R=chol(M11(p,p));
                    
                    switch Data.Model.vec
                        case {'space_time', 'time'}
                            for k=1:N
                                f=computeRHS(Mesh, Data, T(k));
                                RHS=delta_t*(f(Nodes_F)-M12*d1_U(Nodes_D,k)-K12*U(Nodes_D,k))+B*U(Nodes_F,k);
                                U(Nodes_F(p),k+1)=R\(R'\RHS(p));
                            end
                        otherwise % Right-hand side vector is constant
                            f=computeRHS(Mesh, Data, 1);
                            for k=1:N
                                RHS=delta_t*(f(Nodes_F)-M12*d1_U(Nodes_D,k)-K12*U(Nodes_D,k))+B*U(Nodes_F,k);
                                U(Nodes_F(p),k+1)=R\(R'\RHS(p));
                            end
                    end
                    
                case 'implicit_euler'
                    B=M11+delta_t*K11;
                    % Nested dissection algorithm
                    p=dissect(B);
                    R=chol(B(p,p));
                    
                    switch Data.Model.vec
                        case {'space_time', 'time'}
                            for k=1:N
                                f=computeRHS(Mesh, Data, T(k));
                                RHS=delta_t*(f(Nodes_F)-M12*d1_U(Nodes_D,k+1)-K12*U(Nodes_D,k+1))+M11*U(Nodes_F,k);
                                U(Nodes_F(p),k+1)=R\(R'\RHS(p));
                            end
                        otherwise % Right-hand side vector is constant
                            f=computeRHS(Mesh, Data, 1);
                            for k=1:N
                                RHS=delta_t*(f(Nodes_F)-M12*d1_U(Nodes_D,k+1)-K12*U(Nodes_D,k+1))+M11*U(Nodes_F,k);
                                U(Nodes_F(p),k+1)=R\(R'\RHS(p));
                            end
                    end
                    
                case 'crank_nicolson'
                    B=M11+1/2*delta_t*K11;
                    C=M11-1/2*delta_t*K11;
                    % Nested dissection algorithm
                    p=dissect(B);
                    R=chol(B(p,p));
                    
                    switch Data.Model.vec
                        case {'space_time', 'time'}
                            for k=1:N
                                f1=computeRHS(Mesh, Data, T(k));
                                f2=computeRHS(Mesh, Data, T(k+1));
                                RHS=1/2*delta_t*(f2(Nodes_F)+f1(Nodes_F)-M12*(d1_U(Nodes_D,k+1)+d1_U(Nodes_D,k))-K12*(U(Nodes_D,k+1)+U(Nodes_D,k)))+C*U(Nodes_F,k);
                                U(Nodes_F(p),k+1)=R\(R'\RHS(p));
                            end
                        otherwise
                            f=computeRHS(Mesh, Data, 1);
                            for k=1:N
                                RHS=1/2*delta_t*(2*f(Nodes_F)-M12*(d1_U(Nodes_D,k+1)+d1_U(Nodes_D,k))-K12*(U(Nodes_D,k+1)+U(Nodes_D,k)))+C*U(Nodes_F,k);
                                U(Nodes_F(p),k+1)=R\(R'\RHS(p));
                            end
                    end
                    
            end
        else % Nonlinear solver
            % Time integration scheme: Implicit Euler method
            % Retrieve parameters for the nonlinear solver
            % Maximum number of iterations for the Newton method
            max_iter=Data.Solver.System.max_iter;
            % Tolerance for the Newton method
            tol=Data.Solver.System.tolerance;
            
            % Initialize indices
            [Data]=preProcess(Mesh, Data);
            
            % Initialize the matrices
            [Data]=initializeMatrices(Mesh, Data, U(:,1));
            
            % Right-hand side vector
            f=computeRHS(Mesh, Data, 1);
            
            Operation=Data.Operation;
            
            for k=1:N
                % Initilization
                U0=U(:,k);
                Ui=zeros(nb_nodes,1);
                d1_Ui=d1_U(:,k+1);
                Ui(Nodes_D)=U(Nodes_D,k+1);
                Ui(Nodes_F)=U0(Nodes_F);
                beta=zeros(max_iter,1);
                
                if strcmp(Data.Model.vec, 'space_time') || strcmp(Data.Model.vec, 'time')
                    f=computeRHS(Mesh, Data, T(k+1));
                end
                
                for m=1:max_iter
                    
                    U1=zeros(nb_nodes,1);
                    U2=zeros(nb_nodes,1);
                    
                    U1(Nodes_F)=Ui(Nodes_F)-U0(Nodes_F);
                    U1(Nodes_D)=delta_t*d1_Ui(Nodes_D);
                    
                    U2(Nodes_F)=delta_t*Ui(Nodes_F);
                    U2(Nodes_D)=delta_t*Ui(Nodes_D);
                    
                    % Reassemble the matrices which need to be updated
                    [Vmd, Vtd, Vmc, Vtc, Vma, Vta, Vmbc, Vtbc]=AssembleMatrices(Mesh, Data, Ui);
                    
                    % Matrix-vector operations
                    [vec1]=multiplyMatVec(Mesh, Data, Vmd, U1, 'bulk');
                    [vec2]=multiplyMatVec(Mesh, Data, sumObject(Vmc, Vma), U2, 'bulk');
                    [vec3]=multiplyMatVec(Mesh, Data, Vmbc, U2, 'bc');
                    
                    % Sum of contributions
                    G=[sumObject(vec1, vec2); vec3; -delta_t*f];
                    G=accumarray(Operation.Rv, G);
                    G=accumarray(Operation.Iv, G, [nb_nodes, 1]);
                    G=G(Nodes_F);
                    
                    % Tensor-vector operations
                    [mat1]=multiplyTenVec(Mesh, Data, Vtd, U1, 'bulk');
                    [mat2]=multiplyTenVec(Mesh, Data, sumObject(Vtc, Vta), U2, 'bulk');
                    [mat3]=multiplyTenVec(Mesh, Data, Vtbc, U2, 'bc');
                    
                    % Sum different contributions
                    m1=sumObject(mat1, mat2);
                    m2=delta_t*sumObject(Vmc, Vma);
                    m2=sumObject(Vmd, m2);
                    m1=sumObject(m1, m2);
                    m2=sumObject(mat3, delta_t*Vmbc);
                    
                    % Jacobian matrix
                    dG=[m1; m2];
                    dG=accumarray(Operation.Rm, dG);
                    dG=accumarray(Operation.Im, dG, [nb_nodes, nb_nodes], [], [], true);
                    dG=dG(Nodes_F, Nodes_F);
%                     delta_U=-(dG\G);
                    
                    [delta_U, ~, ~]=gmres_restart(@(x) dG*x, -G, zeros(length(Data.Nodes_F), 1));
%                     delta_U=-gmres(dG, G, 100);
%                     delta_U=-bicgstab(dG, G, 1e-6, 100);
                    
                    Ui(Nodes_F)=Ui(Nodes_F)+delta_U;
                    beta(m)=norm(delta_U);
                    
                    if beta(m)<tol*norm(Ui)
                        break
                    end
                end
                
                if beta(m)>tol*norm(Ui)
                    disp(['The Newton method did not converge within the specified number of iterations. Norm of increment was ' num2str(beta(m))]);
                end
                U(:,k+1)=Ui;
            end
        end
        
        % Update solution structure
        Solution.U=U;
        
        function[vec]=linearOperator(Mesh, Data, m1, m2, x)
            
            u=zeros(Mesh.nb_nodes,1);
            u(Data.Nodes_F)=x;
            [v1]=multiplyMatVec(Mesh, Data, m1, u, 'bulk');
            [v2]=multiplyMatVec(Mesh, Data, m2, u, 'bc');
            vec=accumarray([Data.G; Data.Gbc], [v1; v2], [Mesh.nb_nodes, 1]);
            vec=vec(Data.Nodes_F);
            
        end
        
        function[vec]=sumObject(V1, V2)
            
            if isempty(V1)
                vec=V2;
            elseif isempty(V2)
                vec=V1;
            else
                vec=V1+V2;
            end
            
        end
        
        function[vec]=multiplyMatVec(Mesh, Data, V, U, ind)
            
            if isempty(V)
                vec=[];
            else
                if strcmp(ind, 'bulk')
                    nne=Mesh.nb_nodes_elem;
                    ne=Mesh.nb_elem;
                    I=Data.G;
                else
                    nne=Mesh.order+1;
                    ne=Data.nb_edges;
                    I=Data.Gbc;
                end
                
                U=U(I);
                V=reshape(V, [nne, nne, ne]);
                U=reshape(U, [nne, 1, ne]);
                vec=pagemtimes(V, 'transpose', U, 'none');
                vec=reshape(vec, [], 1);
                
            end
            
        end
        
        function[mat]=multiplyTenVec(Mesh, Data, V, U, ind)
            
            % multiplyTenVec: Computes the 1st mode multiplication between an implicitly
            % defined third order tensor and a vector
            
            if isempty(V)
                mat=[];
            else
                if strcmp(ind, 'bulk')
                    nne=Mesh.nb_nodes_elem;
                    ne=Mesh.nb_elem;
                    I=Data.G;
                else
                    nne=Mesh.order+1;
                    ne=Data.nb_edges;
                    I=Data.Gbc;
                end
                
                U=U(I);
                V=reshape(V, [nne, nne^2, ne]);
                U=reshape(U, [nne, 1, ne]);
                mat=pagemtimes(V, 'transpose', U, 'none');
                mat=reshape(mat, [], 1);

            end
        end
        
        
        
    end

end