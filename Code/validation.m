% validation: Equivalent of experiment 2 using MATLAB’s PDE toolbox 
% (the toolbox must be installed beforehand)


clc
close all
clear variables

model = createpde('thermal','transient');

width = 0.2; 
height = 0.2;
gdm = [3 4 0 width width 0 0 0 height height]';
g = decsg(gdm, 'S1', ('S1')');


geometryFromEdges(model,g);
figure; 
pdegplot(model,'EdgeLabels','on'); 
axis([-.1 1.1 -.1 1.1]);
title 'Geometry With Edge Labels Displayed';

msh=generateMesh(model,'Hmax',0.005);
figure 
pdeplot(model); 
axis equal
title 'Block With Finite Element Mesh Displayed'

T1=0;       k1=1.5;
T2=200;     k2=0.7;
T3=1000;    k3=0.5;

k = @(~,state) (k1+(k2-k1)/(T2-T1)*(state.u-T1)).*(state.u>=T1 & state.u<=T2)+(k2+(k3-k2)/(T3-T2)*(state.u-T2)).*(state.u>T2 & state.u<=T3)+k3*(state.u>T3);
thermalProperties(model,'ThermalConductivity',k,...
                        'MassDensity',2400,...
                        'SpecificHeat',1000);

model.StefanBoltzmannConstant = 5.670373E-8; 
thermalBC(model,'Edge',[1 2 3 4],'ConvectionCoefficient', 10, 'Emissivity',0.8,'AmbientTemperature',1000)

thermalIC(model, 0)

endTime = 10800;
tlist = 0:10:endTime;

R = solve(model,tlist);
T = R.Temperature;

getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
[~,nid] = getClosestNode( msh.Nodes, 0.1, 0.1 );

figure
plot(tlist, T(nid,:), '-b', 'displayname', 'MATLAB PDE toolbox'); 
grid on; hold on;

title('Temperature at the center')
xlabel('Time [s]')
ylabel('Temperature [ºC]')
legend('location', 'northwest')