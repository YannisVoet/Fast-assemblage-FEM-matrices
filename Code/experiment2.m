% Experiment 2: Example 2 of SAFIR validation examples. 
% Finite element simulation of a nonlinear transient heat transfer problem. 
% Warning: this experiment might take a few minutes.

% Limitations: 2D
clc
clear variables
close all

% Read GMSH mesh file (.msh)
[Mesh] = readGMSH('../GMSH/Meshes/experiment2_p2_0_018.msh');

%% Parameter prescription
% Mesh parameters
[Mesh]=getParameters(Mesh);
% Initialize model
[Data]=setModel('model', 'transient_heat', 'submodel', 'standard', 'postprocess', 'flux');

%% Set material parameters
T1=0;       k1=1.5;
T2=200;     k2=0.7;
T3=1000;    k3=0.5;
% Thermal conductivity [W/(m*K)]
k=@(T) (k1+(k2-k1)/(T2-T1)*(T-T1)).*(T>=T1 & T<=T2)+(k2+(k3-k2)/(T3-T2)*(T-T2)).*(T>T2 & T<=T3)+k3*(T>T3);
% Derivative of the thermal conductivity
der_k = @(T) (k2-k1)/(T2-T1)*(T>=T1 & T<=T2)+(k3-k2)/(T3-T2)*(T>T2 & T<=T3)+0*(T>T3);
% Specific mass [kg/m^3]
rho=2400;
% Specific heat [J/(kg*K)]
c=1000;
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
[Data]=setCoefficient(Data, 'c', @(T) t*k(T));
[Data]=setCoefficient(Data, 'a', @(T) 0);

[Data]=setCoefficient(Data, 'der_d', @(T) 0);
[Data]=setCoefficient(Data, 'der_c', @(T) t*der_k(T));
[Data]=setCoefficient(Data, 'der_a', @(T) 0);

[Data]=setCoefficient(Data, 'f', @(x,t) 0);

%% Set boundary conditions

BC={{'Type', 'Convection_Radiation', 'Edges', [1 2 3 4], 'Convection_Coeff', 10, 'Emissivity', 0.8, 'Ambient_Temp', 1000}};

p2={'BC',      BC};
[Data]=setBoundaryConditions(Mesh, Data, p2);

%% Set initial conditions
% Value of the solution at time 0
[Data]=setInitialConditions(Data, @(x) 0);
% Set numerical parameters
% Time discretization (number of subintervals)
% Time interval [s]
[Data]=setNumericalParameters(Data, 'sub_Intervals', 1080, 'interval', [0, 10800]);

%% Computations
[Solution] = solve(Mesh, Data);

%% Visualization of the solution
% snapshot: selected time-step at which the solution should be viewed
%   0 -> none (default)
%   any time step in the range defined
% view_video: view the video of the solution
%   1 -> yes (default)
%   0 -> no
% save_video: choose whether to save the video or not
%   1 -> yes
%   0 -> no (default)
% plot_type: choose which solution to visualize
%   temp (default)
%   flux
% metric: which type of temperature should be visualized
%   standard (default)
% metric: which type of flux should be visualized
%   standard (default)
plotSolution(Mesh, Data, Solution, {200, 0, 0}, {'temp', 'flux'}, {'standard', 'standard'});
plotTemperature(Mesh, Data, Solution, [0.1; 0.1])