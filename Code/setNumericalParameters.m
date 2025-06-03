function[Data]=setNumericalParameters(Data, varargin)

% setNumericalParameters: Initializes the numerical parameters of the
% solver
% INPUT:
% Data: Structure containing all the data parameters
% OPTIONAL INPUT:
% Solver:           Name of the ODE solver to be used
%                   Default for linear models: Crank-Nicolson
%                   Default for nonlinear models: Implicit Euler
% System:           Name of the solver for nonlinear systems of equations
%                   Default: Newton
% Interval:         Time interval [s]
%                   Default: [0 5] [s]
% Sub_Intervals:    Number of sub-intervals in time
%                   Default: 100
% Tolerance:        Tolerance for the solution of the nonlinear systems
%                   Default: 1e-7
% Max_Iter:         Maximum number of iterations for the solution of the
%                   nonlinear systems
%                   Default: 30
% OUTPUT:
% Data:             Updated data structure with solver parameters

if mod(length(varargin), 2)~=0
    error('Missing a name or value')
end
labels_in=varargin(1:2:end);
values_in=varargin(2:2:end);

% Set default parameters
% Number of sub-intervals in time
Discretization.N=100;
% Time interval [s]
Discretization.I=[0 5];

% Solver parameters
ODE_solvers_linear={ 'crank_nicolson', 'explicit_euler', 'implicit_euler'};
ODE_solvers_nonlinear={'implicit_euler'};
System_solvers={'newton'};

if strcmp(Data.Model.type, 'linear')
    Solver.ODE.name=ODE_solvers_linear{1};
else
    Solver.ODE.name=ODE_solvers_nonlinear{1};
end
Solver.System.name=System_solvers{1};
Solver.System.tolerance=1e-7;
Solver.System.max_iter=30;

for i=1:length(labels_in)
    arg=lower(labels_in{i});
    val=values_in{i};
    
    switch arg
        case 'solver'
            if strcmp(Data.Model.type, 'linear')
                if any(ismember(ODE_solvers_linear, lower(val)))
                    Solver.ODE.name=lower(val);
                else
                    error('Unrecognized ODE solver or unsupported for linear systems of ODEs')
                end
            else
                if any(ismember(ODE_solvers_nonlinear, lower(val)))
                    Solver.ODE.name=lower(val);
                else
                    error('Unrecognized ODE solver or unsupported for nonlinear systems of ODEs')
                end
            end
            
        case 'system'
            if any(ismember(System_solvers, lower(val)))
                Solver.System.name=lower(val);
            else
                error('Unrecognized nonlinear system solver')
            end
            
        case 'sub_intervals'
            if val<2 || ~isa(val, 'double')
                warning('Erroneous number of sub-intervals, changing to default value')
            else
                Discretization.N=val;
            end
            
            
        case 'interval'
            if length(val)~=2
                error('Time interval must have a starting time and an ending time')
            elseif val(2)<=val(1)
                error('Starting time must be smaller than ending time')
            else
                Discretization.I=val;
            end
            
        case 'tolerance'
            if val<0 || val>1
                warning('Erroneous tolerance value, changing to default value')
            else
                Solver.System.tolerance=val;
            end
            
        case 'max_iter'
            if val<2 || ~isa(val, 'double')
                warning('Erroneous maximum number of iterations, changing to default value')
            else
                Solver.System.max_iter=val;
            end
            
        otherwise
            error('Unrecognized argument')
    end
    
end

N=Discretization.N;
I=Discretization.I;
% Time step [s]
Discretization.delta_t=(I(2)-I(1))/N;
% Time vector [s]
Discretization.time_vec=linspace(I(1), I(2), N+1);

% Update data structure
Data.Solver=Solver;
Data.Discretization=Discretization;
end