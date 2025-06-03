% Experiment 1: Numerical experiments designed to compare the performance 
% of a standard algorithm with the new algorithm for the fast assemblage 
% of finite element matrices. Warning: this experiment might take a while.

% Limitation: 2D
clc
clear variables
close all

names={'p1_0_05', 'p1_0_025', 'p1_0_0125', 'p1_0_0063', 'p1_0_0031', 'p1_0_0016'};
% names={'p2_0_05', 'p2_0_025', 'p2_0_0125', 'p2_0_0063', 'p2_0_0031', 'p2_0_0016'};

% Number of launches
N=5;
n=length(names);
sizes=zeros(n,1);

timeDMatrix1=zeros(n,1);
timeDMatrix2=zeros(n,1);

timeDMatrix1_full=zeros(n,1);
timeDMatrix2_full=zeros(n,1);

timeCMatrix1=zeros(n,1);
timeCMatrix2=zeros(n,1);

timeCMatrix1_full=zeros(n,1);
timeCMatrix2_full=zeros(n,1);

for k=1:N
    for i=1:n
        % Read GMSH mesh file (.msh)
        [Mesh] = readGMSH(['../GMSH/Meshes/experiment1_' names{i} '.msh']);
        
        %% Parameter prescription
        % Mesh parameters
        [Mesh]=getParameters(Mesh);
        % Initialize model
        [Data]=setModel('model', 'transient_heat', 'submodel', 'standard', 'type', 'linear', 'postprocess', 'none');
        
        %% Set material parameters
        % Thermal conductivity [W/(m*K)]
        kc=1;
        % Specific mass [kg/m^3]
        rho=1000;
        % Specific heat [J/(kg*K)]
        c=1;
        % Stefan Boltzman constant W/(m^2*K^4)
        sigma=5.670373e-8;
        % Convection coefficient [W/(m^2*K)]
        hc=1;
        % Emissivity [-]
        e=0.5;
        % Ambient temperature [ºC] or [K]
        Ta=300;
        % Thickness [m]
        t=1;
        
        [Data]=setCoefficient(Data, 'd', @(T) t*c*rho);
        [Data]=setCoefficient(Data, 'c', @(T) t*kc);
        [Data]=setCoefficient(Data, 'a', @(T) 0);
        
        [Data]=setCoefficient(Data, 'der_d', @(T) 0);
        [Data]=setCoefficient(Data, 'der_c', @(T) 0);
        [Data]=setCoefficient(Data, 'der_a', @(T) 0);
        
        [Data]=setCoefficient(Data, 'f', @(x,t) 0);
        
        %% Set boundary conditions
        
        BC={{'Type', 'Convection_Radiation', 'Edges', 1, 'Convection_Coeff', 1, 'Emissivity', 0, 'Ambient_Temp', 0}};
        
        p2={'BC',      BC};
        [Data]=setBoundaryConditions(Mesh, Data, p2);
        
        %% Set initial conditions
        % Value of the solution at time 0
        [Data]=setInitialConditions(Data, @(x) 1000);
        % Set numerical parameters
        % Time discretization (number of subintervals)
        % Time interval [s]
        [Data]=setNumericalParameters(Data, 'sub_Intervals', 1800, 'interval', [0, 1800]);
        
        %% Assemblage of the matrices
        sizes(i)=Mesh.nb_nodes;
        
        tic
        [I, J, V]=computeDMatrix1(Mesh, Data, ones(Mesh.nb_nodes,1));
        timeDMatrix1(i)=timeDMatrix1(i)+toc;
        M=sparse(I(:), J(:), V(:), Mesh.nb_nodes, Mesh.nb_nodes);
        timeDMatrix1_full(i)=timeDMatrix1_full(i)+toc;
        tic
        [I, J, V]=computeDMatrix2(Mesh, Data, ones(Mesh.nb_nodes,1));
        timeDMatrix2(i)=timeDMatrix2(i)+toc;
        M=sparse(I(:), J(:), V(:), Mesh.nb_nodes, Mesh.nb_nodes);
        timeDMatrix2_full(i)=timeDMatrix2_full(i)+toc;
        tic
        [I, J, V]=computeCMatrix1(Mesh, Data, ones(Mesh.nb_nodes,1));
        timeCMatrix1(i)=timeCMatrix1(i)+toc;
        K=sparse(I(:), J(:), V(:), Mesh.nb_nodes, Mesh.nb_nodes);
        timeCMatrix1_full(i)=timeCMatrix1_full(i)+toc;
        tic
        [I, J, V]=computeCMatrix2(Mesh, Data, ones(Mesh.nb_nodes,1));
        timeCMatrix2(i)=timeCMatrix2(i)+toc;
        K=sparse(I(:), J(:), V(:), Mesh.nb_nodes, Mesh.nb_nodes);
        timeCMatrix2_full(i)=timeCMatrix2_full(i)+toc;
        
        
    end
end
% Compute average time
timeDMatrix1=timeDMatrix1/N;
timeDMatrix2=timeDMatrix2/N;
timeCMatrix1=timeCMatrix1/N;
timeCMatrix2=timeCMatrix2/N;

timeDMatrix1_full=timeDMatrix1_full/N;
timeDMatrix2_full=timeDMatrix2_full/N;
timeCMatrix1_full=timeCMatrix1_full/N;
timeCMatrix2_full=timeCMatrix2_full/N;
%% Visualization without the sparse function call

% Importing data from FreeFem++
Times_FF=importdata('../Data/Times_p1_FF++.dat', ' ', 0);

timeDMatrix_FF=Times_FF(:,1);
timeCMatrix_FF=Times_FF(:,2);

xf=[sizes', fliplr(sizes')];

DMatrix1yf=[timeDMatrix1', fliplr(timeDMatrix1_full')];
DMatrix2yf=[timeDMatrix2', fliplr(timeDMatrix2_full')];

CMatrix1yf=[timeCMatrix1', fliplr(timeCMatrix1_full')];
CMatrix2yf=[timeCMatrix2', fliplr(timeCMatrix2_full')];

figure
loglog(sizes, timeDMatrix1, '--ob', 'LineWidth', 2)
hold on; grid on;
loglog(sizes, timeDMatrix1_full, '-ob', 'LineWidth', 2)
loglog(sizes, timeDMatrix2, '--sr', 'LineWidth', 2)
loglog(sizes, timeDMatrix2_full, '-sr', 'LineWidth', 2)
loglog(sizes, timeDMatrix_FF, '-xk', 'LineWidth', 2)
loglog(sizes, sizes*(3e-7), '--k', 'LineWidth', 2)
fill(xf, DMatrix1yf, 'b', 'FaceAlpha', 0.3)
fill(xf, DMatrix2yf, 'r', 'FaceAlpha', 0.3)
title('Global mass matrix assembly')
legend('Algorithm 1, Formation', 'Algorithm 1, Assemblage', 'Algorithm 2, Formation', 'Algorithm 2, Assemblage', 'FreeFEM++', 'Order 1')
xlabel('Matrix size')
ylabel('Computation time [s]');
set(gca,'FontSize',15)

figure
loglog(sizes, timeCMatrix1, '--ob', 'LineWidth', 2)
hold on; grid on;
loglog(sizes, timeCMatrix1_full, '-ob', 'LineWidth', 2)
loglog(sizes, timeCMatrix2, '--sr', 'LineWidth', 2)
loglog(sizes, timeCMatrix2_full, '-sr', 'LineWidth', 2)
loglog(sizes, timeCMatrix_FF, '-xk', 'LineWidth', 2)
loglog(sizes, sizes*(3e-7), '--k', 'LineWidth', 2)
fill(xf, CMatrix1yf, 'b', 'FaceAlpha', 0.3)
fill(xf, CMatrix2yf, 'r', 'FaceAlpha', 0.3)
title('Global conductivity matrix assembly')
legend('Algorithm 1, Formation', 'Algorithm 1, Assemblage', 'Algorithm 2, Formation', 'Algorithm 2, Assemblage', 'FreeFEM++', 'Order 1')
xlabel('Matrix size')
ylabel('Computation time [s]');
set(gca,'FontSize',15)

fprintf('Speedup factor for the formation of element matrices: \n\n')
DMatrixAlgorithm1=timeDMatrix1./timeDMatrix2;
CMatrixAlgorithm1=timeCMatrix1./timeCMatrix2;
DMatrixFreeFEM=timeDMatrix_FF./timeDMatrix2;
CMatrixFreeFEM=timeCMatrix_FF./timeCMatrix2;

formation=table(DMatrixAlgorithm1, DMatrixFreeFEM, CMatrixAlgorithm1, CMatrixFreeFEM);
disp(formation)

fprintf('Speedup factor for the assemblage of finite element matrices: \n\n')
DMatrixAlgorithm1=timeDMatrix1_full./timeDMatrix2_full;
CMatrixAlgorithm1=timeCMatrix1_full./timeCMatrix2_full;
DMatrixFreeFEM=timeDMatrix_FF./timeDMatrix2_full;
CMatrixFreeFEM=timeCMatrix_FF./timeCMatrix2_full;

assemblage=table(DMatrixAlgorithm1, DMatrixFreeFEM, CMatrixAlgorithm1, CMatrixFreeFEM);
disp(assemblage)

fprintf('Time fraction for sparse function: \n\n')
DMatrixAlgorithm1=(1-timeDMatrix1./timeDMatrix1_full)*100;
DMatrixAlgorithm2=(1-timeDMatrix2./timeDMatrix2_full)*100;
CMatrixAlgorithm1=(1-timeCMatrix1./timeCMatrix1_full)*100;
CMatrixAlgorithm2=(1-timeCMatrix2./timeCMatrix2_full)*100;

fraction=table(DMatrixAlgorithm1, DMatrixAlgorithm2, CMatrixAlgorithm1, CMatrixAlgorithm2);
disp(fraction)