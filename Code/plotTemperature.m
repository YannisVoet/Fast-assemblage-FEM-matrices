function[]=plotTemperature(Mesh, Data, Solution, coord_want)

% plotTemperature: Plots the temperature curves at a specific point
% Mesh:             Structure containing the mesh parameters
% Data:             Structure containing the data parameters
% Solution:         Structure containing the solution
% coord_want:       Coordinates of a point where the temperature curve
%                   should be plotted


% Retrieve mesh parameters
% Coordinate matrix
coord=Mesh.coord_mat';

% Retrieve data
% Time vector
time_vec=Data.Discretization.time_vec;

% Retrieve solution
% Temperatures
U=Solution.U;

[~, indx]=min(vecnorm(coord-coord_want, 2, 1));
coord_select=coord(:, indx);
U_select=U(indx,:);

figure
plot(time_vec, U_select, '-b', 'displayname', 'MATLAB')
hold on

% Numerical results from Experiment 2
T=load('../Data/experiment2_mathworks.mat');
T_mathworks=T.T_mathworks;
plot(time_vec, T_mathworks, '-m', 'displayname', 'MATLAB PDE toolbox')


title(['Temperature at point (' num2str(coord_select(1)) ', ' num2str(coord_select(2)) ')'])
xlabel('Time [s]')
ylabel('Temperature [ºC]')
legend show
grid on
end