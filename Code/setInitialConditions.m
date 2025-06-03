function[Data]=setInitialConditions(Data, u0)

% setInitialConditions: Initializes the initial conditions of the problem
% INPUT:
% u0:   Value of the solution for time t=0 defined in the open set
% OUTPUT:
% Data: Updated data structure

Data.u0=u0;
end