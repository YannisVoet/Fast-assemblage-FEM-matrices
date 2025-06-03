function[Data]=setBoundaryConditions(Mesh, Data, varargin)

% setBoundaryConditions: Function which initializes the boundary data
% INPUT:
% Mesh:     Structure containing all the mesh parameters
% Data:     Structure containing all the data parameters
% OPTIONAL INPUT:
% If no boundary conditions are specified, zero flux conditions (Natural
% Neumann boundary condition) are assumed on all the boundary
% Types of boundary conditions:
%   Temperature: Assign a temperature to an edge or several edges
%       Edges:      Edge label(s) where the temperature should be prescribed
%       Function:   Function of the prescribed temperature
%       Derivative: Time derivative of the prescribed temperature
%   Flux: Assign a flux to an edge or several edges
%       Edges:      Edge label(s) where the flux should be prescribed
%       Function:   Function of the prescribed flux
%   Convection_Radiation: Assign a convecive and/or radiative flux to an
%   edge or several edges
%       Edges:              Edge label(s) where the boundary condition should be prescribed
%       Convection_Coeff:   Convection coefficient
%       Emissivity:         Emissivity value
%       Ambient_Temp:       Function for the ambient temperature
% OUTPUT:
% Data: Updated data structure with boundary conditions


%% Set boundary conditions
if mod(length(varargin{1}), 2)~=0
    error('Missing a name or value')
end

% Retrieve boundary nodes
BC_tags=unique(Mesh.BC_tag);

% Define available boundary conditions type
available_types={'temperature', 'flux', 'convection_radiation'};
% Dirichlet type
D_types={'temperature'};
% Neumann types
N_types=setdiff(available_types, D_types);
% Count to detect duplicate boundaries
count=[];
% Count number of Dirichlet and Neumann boundaries specified
nb_D_BC=0;
nb_N_BC=0;

% Retrieve labels and values entered by the user
labels_in=varargin{1}(1:2:end);
values_in=varargin{1}(2:2:end);

% Set default Parameters
Data.BC=[];
% Boolean indicators specifying if the model considers convection/radiation
Data.Model.convection=0;
Data.Model.radiation=0;
% Boolean indicator specifying if the right-hand side vector is time
% dependent
Data.Model.vec='const';

% Checking for boundary conditions
c1=strcmp(labels_in, 'BC');

if any(c1)
    BC=values_in{c1};
    nb_BC=length(BC);
    Data.BC=cell(nb_BC,1);
    
    for i=1:nb_BC
        bc=BC{i};
        l=length(bc);
        
        if mod(l,2)~=0
            error('Missing a name or value')
        end
        
        for j=1:2:l
            arg=lower(bc{j});
            val=bc{j+1};
            
            switch arg
                
                case 'type'
                    if any(ismember(available_types, lower(val)))
                        Data.BC{i}.(arg)=lower(val);
                        
                        if any(ismember(D_types, lower(val)))
                            nb_D_BC=nb_D_BC+1;
                        else
                            nb_N_BC=nb_N_BC+1;
                        end
                    else
                        error('Unrecognized boundary condition')
                    end
                    
                case 'edges'
                    
                    edges=intersect(BC_tags, val);
                    
                    if isempty(edges)
                        error('One of the edges specified does not exist')
                    elseif ~isempty(intersect(count, val)) || length(unique(val))~=length(val)
                        error('Duplicate edges detected')   
                    end
                    count=[count val];
                    Data.BC{i}.(arg)=val;
                    
                case {'function', 'ambient_temp', 'derivative'}
                    
                    if isa(val,'double')
                        Data.BC{i}.label='const';
                        Data.BC{i}.(arg)= @(x,t) val;
                        
                    elseif isa(val,'function_handle')
                        Data.BC{i}.label=detectDependency(val);
                        Data.BC{i}.(arg)=val;
                    else
                        error('Input parameters must be either constants or function handles')
                    end
                    
                    
                case 'convection_coeff'
                    
                    if val<0
                        error('Unphysical value of convection coefficient')
                    elseif val>0
                        Data.Model.convection=1;
                    end
                    Data.BC{i}.(arg)=val;
                    
                case 'emissivity'
                    
                    if val<0 || val >1
                        error('Unphysical value of emissivity')
                    elseif val>0
                        Data.Model.radiation=1;
                        Data.Model.type='nonlinear';
                        Data.Model.bc='nonlinear';
                    end
                    Data.BC{i}.(arg)=val;
                    
                otherwise
                    error('Unrecognized argument')
            end
            
        end
        
        
    end
    
end

%% Post-Processing: Partitioning the structure BC into Dirichlet BC and Neumann BC and checking for missing data
BC=Data.BC;
l=length(BC);

% Initialize cells
D_BC=cell(nb_D_BC,1);
N_BC=cell(nb_N_BC,1);
% Set of Dirichlet and Neumann tags
D_tags=[];
N_tags=[];
% Initialize counters
d=1;
n=1;

for i=1:l
    bc=BC{i};
    type=bc.type;
    
    if ~isfield(bc, 'edges')
        error('The edges must be specified')
    end
    
    if any(ismember(D_types, type))
        D_tags=[D_tags bc.edges];
        D_BC{d}.type=type; 
        D_BC{d}.edges=bc.edges;
        
        if isfield(bc, 'function')
            D_BC{d}.function=bc.function;
            D_BC{d}.label=bc.label;
        else
            warning('A temperature function has not been prescribed and will be assumed zero')
            D_BC{d}.function=@(x,t) 0;
            D_BC{d}.label='const';
        end
        
        if isfield(bc, 'derivative')
            D_BC{d}.derivative=bc.derivative;
        else
            error('The time derivative has not been prescribed')
        end
        
        d=d+1;
        
    else
        N_tags=[N_tags bc.edges];
        N_BC{n}.type=type;
        N_BC{n}.edges=bc.edges;
        
        if strcmp(type, 'flux')
            
            if isfield(bc, 'function')
                N_BC{n}.function=bc.function;
                N_BC{n}.label=bc.label;
            else
                warning('A flux function has not been prescribed and will be assumed zero')
                N_BC{n}.function=@(x,t) 0;
                N_BC{n}.label='const';
            end
            
            
        elseif strcmp(type, 'convection_radiation')
            
            if isfield(bc, 'convection_coeff')
                N_BC{n}.convection_coeff=bc.convection_coeff;
            else
                warning('The convection coefficient has not been prescribed and will be assumed to be 10 [W/(m^2*K)]')
                N_BC{n}.convection_coeff=10;
            end
            
            if isfield(bc, 'emissivity')
                N_BC{n}.emissivity=bc.emissivity;
            else
                warning('The emissivity has not been prescribed and will be assumed zero')
                N_BC{n}.emissivity=0;
            end
            
            if isfield(bc, 'ambient_temp')
                N_BC{n}.ambient_temp=bc.ambient_temp;
                N_BC{n}.label=bc.label;
            else
                warning('The ambient temperature has not been prescribed and will be assumed zero')
                N_BC{n}.ambient_temp=@(x,t) 0;
                N_BC{n}.label='const';
            end
            
            
        end
        
        % RHS vector becomes time dependent if the flux function is or the 
        % ambient temperature is
        [Data.Model.vec]=updateLabel(Data.Model.vec, N_BC{n}.label);
        n=n+1;
    end
end

% RHS vector becomes time dependent if the user defined function f is
[Data.Model.vec]=updateLabel(Data.Model.vec, Data.Coefficient.f.label);

Data.D_BC=D_BC;
Data.N_BC=N_BC;

Data.D_tags=D_tags;
Data.N_tags=N_tags;

%% Process coefficients to detect nonlinearities

coefficients={'d', 'c', 'a'};

for k=1:length(coefficients)
    if strcmp(Data.Coefficient.(coefficients{k}).label, 'temp')
        Data.Model.type='nonlinear';
    end
    
end

%% Sorting nodes
% If there is a node belonging to conflicting boundary conditions, it is
% set using the following priority list: Dirichlet, Neumann, Interior
BC_nodes=Mesh.BC_nodes;
BC_tag=Mesh.BC_tag;
nb_nodes=Mesh.nb_nodes;
Dof_set=Mesh.Dof_set;

set_D_nodes=cell(nb_D_BC,1);
set_N_nodes=cell(nb_N_BC,1);

for k=1:nb_D_BC
    nodes=BC_nodes(any(BC_tag==D_BC{k}.edges,2),:);
    set_D_nodes{k}=unique(nodes);
end
for k=1:nb_N_BC
    nodes=BC_nodes(any(BC_tag==N_BC{k}.edges,2),:);
    set_N_nodes{k}=unique(nodes);
end
Nodes_D=unique(cell2mat(set_D_nodes));
Nodes_N=setdiff(unique(BC_nodes), Nodes_D);
Nodes_I=setdiff(1:nb_nodes, union(Nodes_D, Nodes_N));

% Correct the set for Neumann nodes
for k=1:nb_N_BC
    set_N_nodes{k}=setdiff(set_N_nodes{k}, Nodes_D);
end

% Define the sets of nodes
% Interior nodes
Ni=Dof_set(:,Nodes_I);
% Dirichlet nodes
Nd=Dof_set(:,Nodes_D);
% Neumann nodes
Nn=Dof_set(:,Nodes_N);

Ni=Ni(:);
Nd=Nd(:);
Nn=Nn(:);

% Update data
Data.Nodes_I=Nodes_I;
Data.Nodes_N=Nodes_N;
Data.Nodes_D=Nodes_D;
Data.Nodes_F=union(Nodes_I, Nodes_N);

Data.Ni=Ni;
Data.Nd=Nd;
Data.Nn=Nn;

Data.set_D_nodes=set_D_nodes;
Data.set_N_nodes=set_N_nodes;

%% Auxiliary function

    function[label]=detectDependency(func)
        
        string=func2str(func);
        
        if ~strcmp(string(1:6), '@(x,t)')
            error('Function handle must contain space and time variables')
        end
        string=string(7:end);
        
        if contains(string, 'x') && contains(string, 't')
            label='space_time';
        elseif contains(string, 'x')
            label='space';
        elseif contains(string, 't')
            label='time';
        elseif contains(string(1), '0')
            label='zero';
        else
            label='const';
        end
        
    end

    function[label_out]=updateLabel(label_in, label)
        
        dependency={'zero', 'const', 'space', 'time', 'space_time'};
        [B_in, I_in]=ismember(dependency, label_in);
        [B, I]=ismember(dependency, label);
        
        indx_in=find(I_in);
        indx=find(I);
        
        if indx > indx_in
            label_out=dependency{B};
        else
            label_out=dependency{B_in};
        end
    end

end
