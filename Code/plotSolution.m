function[]=plotSolution(Mesh, Data, Solution, varargin)

% plotSolution: Function which plots the displacements and/or stresses 
% computed
% INPUT:
% Mesh:     Structure containing the mesh parameters
% Data:     Structure containing the data parameters
% Solution: Structure containing the solution computed
% OPTIONAL INPUT
% snapshot: selected time-step at which the solution should be viewed
%   0 -> none (default)
%   any time step in the range defined
% view_video: view the video of the solution
%   1 -> yes (default)
%   0 -> no
% save_video: choose whether to save the video or not
%   1 -> yes
%   0 -> no (default)
% plot_arg: choose which solution to visualize
%   dynamic_solid module
%       disp (default)
%       stress
%   transient_heat module
%       temp (default)
%       flux
% plot_type: decide which type of metric should be used for visualization
%   displacement metric
%       standard -> structure with amplified displacements (default)
%       norm
%       absx
%       absy
%   temperature metric
%       standard (default)
%   stress metric
%       von_mises (default)
%       sigma_xx
%       sigma_yy
%       sigma_xy
%       frobenius
%   flux metric
%       standard (default)


%% Initialize Parameters

switch Data.Model.name
    case 'dynamic_solid'
        plot_argv={'disp', 'stress'};
        n=length(plot_argv);
        plot_typev=cell(n,1);
        plot_typev{1}={'standard', 'norm', 'absx', 'absy'};
        plot_typev{2}={'von_mises', 'sigma_xx', 'sigma_yy', 'sigma_xy', 'frobenius'};
        
    case 'transient_heat'
        plot_argv={'temp', 'flux'};
        n=length(plot_argv);
        plot_typev=cell(n,1);
        plot_typev{1}={'standard'};
        plot_typev{2}={'standard'};
        
end

switch nargin
    case 3 % Set all default parameters
        % Default visualization parameters
        Param.visua_param.snapshot=0;
        Param.visua_param.view_video=1;
        Param.visua_param.write_video=0;
        % Default plot argument
        Param.plot_param.arg=plot_argv(1);
        % Default plot type
        Param.plot_param.type=plot_typev{1}(1);
        
    case 4 % Set default plot argument and plot type      
        [Param]=setParameters(Data, varargin, 1);
        % Default plot argument
        Param.plot_param.arg=plot_argv(1);
        % Default plot type
        Param.plot_param.type=plot_typev{1}(1);
        
    case 5 % Set default plot type
        [Param]=setParameters(Data, varargin, 2);
        % Default plot type
        
        plot_arg=Param.plot_param.arg;
        for t=1:length(plot_arg)
            index=strcmp(plot_arg{t}, plot_argv);
            plot_type{t}=plot_typev{index}{1};
        end
        Param.plot_param.type=plot_type;
        
    case 6 % Check all parameters entered
        [Param]=setParameters(Data, varargin, 3);
end


%% Plotting interface
% Import data
N=Data.Discretization.N;


plot_arg=Param.plot_param.arg;
plot_type=Param.plot_param.type;
n=length(plot_arg);
Plot=cell(n,1);
Axis=cell(n,1);
Colorbar=cell(n,1);
Dynamic=cell(n,1);

fig=figure;
% fig.Units='normalized';
% fig.Position=[0 0 1 1];

% fig.Units='normalized';
% fig.Position=[0 1 1/2 1/2];

for k=1:n
    switch plot_arg{k}
        case {'disp', 'temp'}
            [PlotParam, AxisParam, ColorbarParam, DynamicParam]=plotPatch(Mesh, Data, Solution, plot_type{k});
            subplot(1,n,k)
            axis=gca;
            p=patch(axis);
            c=colorbar(axis);
            [p]=    setObject(p, PlotParam);
            [axis]= setObject(axis, AxisParam);
            [c]=    setObject(c, ColorbarParam);
            
            
            
        case 'stress'
            [PlotParam, AxisParam, ColorbarParam, DynamicParam]=plotSurf(Param, Mesh, Data, Solution, plot_type{k});
            subplot(1,n,k)
            axis=gca;
            p=surf(axis);
            c=colorbar(axis);
            [p]=    setObject(p, PlotParam);
            [axis]= setObject(axis, AxisParam);
            [c]=    setObject(c, ColorbarParam);
            
            
            
        case 'flux'
            [PlotParam, AxisParam, ColorbarParam, DynamicParam]=plotQuiver(Mesh, Data, Solution, plot_type{k});
            subplot(1,n,k)
            axis=gca;
            p=quiver(axis, PlotParam{1,2}, PlotParam{2,2});
            c=colorbar(axis);
            [p]=    setObject(p, PlotParam);
            [axis]= setObject(axis, AxisParam);
            [c]=    setObject(c, ColorbarParam);
            
    end
    
    Plot{k}=p;
    Axis{k}=axis;
    Colorbar{k}=c;
    Dynamic{k}=DynamicParam;
end



if Param.visua_param.view_video % View the solution over the entire simulation
    
    for k=1:n
        [Plot{k}, Axis{k}]=updateObject(Plot{k}, Axis{k}, Dynamic{k}, 1, N);
    end
    
    frames(1)=getframe(fig);
    for step = 2:N+1
        pause(0.05)
        for k=1:n
            [Plot{k}, Axis{k}]=updateObject(Plot{k}, Axis{k}, Dynamic{k}, step, N);
            refreshdata(Axis{k});
        end
        frames(step)=getframe(fig);
    end
   
else % View solution at a particular time snapshot
    step=Param.visua_param.snapshot;
        for k=1:n
            [Plot{k}, Axis{k}]=updateObject(Plot{k}, Axis{k}, Dynamic{k}, step, N);
        end
end


if Param.visua_param.write_video
%     vid_duration=25; % seconds
    vid_duration=10; % seconds
    video=VideoWriter('Output.mp4', 'MPEG-4');
    video.FrameRate=(N+1)/vid_duration;
    open(video);
    writeVideo(video,frames);
    close(video);
end


%% Setting all parameters
    function[Param]=setParameters(Data, arguments, label)
        
        
        if label > 0 % Check if numerical parameters are valid
            visua_param=arguments{1};
            
            if length(visua_param)~=3
                error('All visualization parameters must be specified')
            end
            
            if visua_param{1} < 0 || visua_param{1} > Data.Discretization.N+1
                error(['Time snapshot invalid: it must be an integer between ' num2str(0) ' and ' num2str(Data.Discretization.N+1)])
            else
                Param.visua_param.snapshot=visua_param{1};
            end
            
            if any(visua_param{2}<0) || any(visua_param{2}>1)
                error('View video parameter must be a boolean')
            else
                Param.visua_param.view_video=visua_param{2};
            end
            
            if any(visua_param{3}<0) || any(visua_param{3}>1)
                error('Write video parameter must be a boolean')
            else
                Param.visua_param.write_video=visua_param{3};
            end
            
            
            if label > 1 % Check if plot arguments are valid
                plot_arg=arguments{2};
                
                if any(~ismember(plot_arg, plot_argv))
                    error('Plot argument invalid: see available plot arguments')
                else
                    Param.plot_param.arg=plot_arg;
                end
                
                if label > 2 % Check if plot types are valid
                    plot_type=arguments{3};
                    
                    if length(plot_arg)~=length(plot_type)
                        error('All plot argument and plot type pairs must be specified')
                    end
                    
                    for j=1:length(plot_type)
                        index=strcmp(plot_arg{j}, plot_argv);
                        if any(~ismember(plot_type{j}, plot_typev{index}))
                            error('Plot type invalid: see available plot types')
                        end
                    end
                    Param.plot_param.type=plot_type;
                    
                end
                
            end
        end
    end

%% Parameters for patch plots
    function[PlotParam, AxisParam, ColorbarParam, DynamicParam]=plotPatch(Mesh, Data, Solution, type)
        
        
        % Retrieve mesh parameters
        % Coordinate matrix
        coord=Mesh.coord_mat;
        % Connectivity matrix
        connect=Mesh.connect_mat;
        % Total number of nodes
        nb_nodes=Mesh.nb_nodes;
        % Frame surrounding the structure
        frame=Mesh.frame;
        % Largest dimension
        L_max=Mesh.L_max;
        % Number of local degrees of freedom
        nb_dof_local=Mesh.nb_dof_local;
        % Total number of degrees of freedom
        nb_dof=Mesh.nb_dof;
        
        % Retrieve data parameters
        % Number of sub-intervals in time
        N=Data.Discretization.N;
        
        % Retrieve solution
        U=Solution.U;
        
        
        switch Data.Model.name
            
            case 'dynamic_solid'
                
                % Initilization
                Vertices=zeros(nb_nodes,2,N+1);
                CData=cell(N+1,1);
                
                
                switch type
                    case 'standard'
                        
                        % Maximum displacement
                        disp_max=max(max(abs(U)));
                        % Amplification factor
                        fact_ampl=L_max/(2*disp_max);
                        % Retrieve degrees of freedom
                        d1=1:nb_dof_local:nb_dof;
                        d2=2:nb_dof_local:nb_dof;
                        % Aplified displacements
                        ux_ampl=fact_ampl*U(d1,:);
                        uy_ampl=fact_ampl*U(d2,:);
                        % Deformed state
                        def_x=coord(:,1)+ux_ampl;
                        def_y=coord(:,2)+uy_ampl;
                        % Vertices in deformed state
                        Vertices(:,1,:)=def_x;
                        Vertices(:,2,:)=def_y;
                        XLim=[min(min(def_x));
                            max(max(def_x))];
                        YLim=[min(min(def_y));
                            max(max(def_y))];
                        CLim=[0,1];
                        Colormap=[];
                        ColorbarStatus='off';
                        ColorbarLabel='';
                        
                        Title='Amplified displacements';
                        EdgeColor='black';
                        FaceColor='none';
                        
                    case 'norm'
                        
                        % Retrieve degrees of freedom
                        d1=1:nb_dof_local:nb_dof;
                        d2=2:nb_dof_local:nb_dof;
                        Ux=U(d1,:);
                        Uy=U(d2,:);
                        U_norm=sqrt(Ux.^2+Uy.^2);
                        
                        XLim=frame(:,1);
                        YLim=frame(:,2);
                        CLim=[0 max(U_norm, [], 'all')];
                        Colormap=jet;
                        ColorbarStatus='on';
                        ColorbarLabel='Displacement [m]';
                        Title='Norm of displacements';
                        EdgeColor='none';
                        FaceColor='interp';
                        
                        for m=1:N+1
                            Vertices(:,1,m)=coord(:,1);
                            Vertices(:,2,m)=coord(:,2);
                            CData{m}=U_norm(:,m);
                        end
                        
                    case 'absx'
                        
                        % Retrieve degrees of freedom
                        d1=1:nb_dof_local:nb_dof;
                        Ux_abs=abs(U(d1,:));
                        
                        XLim=frame(:,1);
                        YLim=frame(:,2);
                        CLim=[0 max(Ux_abs, [], 'all')];
                        Colormap=jet;
                        ColorbarStatus='on';
                        ColorbarLabel='Displacement [m]';
                        Title='Displacements in absolute value along x';
                        EdgeColor='none';
                        FaceColor='interp';
                        
                        for m=1:N+1
                            Vertices(:,1,m)=coord(:,1);
                            Vertices(:,2,m)=coord(:,2);
                            CData{m}=Ux_abs(:,m);
                        end
                        
                    case 'absy'
                        
                        % Retrieve degrees of freedom
                        d2=2:nb_dof_local:nb_dof;
                        Uy_abs=abs(U(d2,:));
                        
                        XLim=frame(:,1);
                        YLim=frame(:,2);
                        CLim=[0 max(Uy_abs, [], 'all')];
                        Colormap=jet;
                        ColorbarStatus='on';
                        ColorbarLabel='Displacement [m]';
                        Title='Displacements in absolute value along y';
                        EdgeColor='none';
                        FaceColor='interp';
                        
                        for m=1:N+1
                            Vertices(:,1,m)=coord(:,1);
                            Vertices(:,2,m)=coord(:,2);
                            CData{m}=Uy_abs(:,m);
                        end
                        
                        
                end
                
                
                % Set static parameters
                PlotParam = {'Vertices',        coord;...
                             'Faces',           connect(:,1:3);...
                             'FaceVertexCData', CData{1};...
                             'FaceColor',       FaceColor;...
                             'EdgeColor',       EdgeColor;...
                             'LineWidth',       1};
                
                AxisParam = {'XLim',            XLim;...
                             'YLim',            YLim;...
                             'CLim',            CLim;...
                             'Colormap',        Colormap;...
                             'DataAspectRatio', [1 1 1];...
                             'XLabel.String',          'x coordinate [m]';...
                             'YLabel.String',          'y coordinate [m]'};
                
                ColorbarParam = {'Visible',         ColorbarStatus;...
                                 'Label.String',    ColorbarLabel;...
                                 'Label.FontSize',  11};
                
                % Set dynamic parameters
                DynamicParam={'Vertices',   num2cell(Vertices, [1,2]);...
                              'CData',      CData;...
                              'Title',      Title};
                          
            case 'transient_heat'
                
                % Initilization
                CData=cell(N+1,1);
                
                switch type
                    case 'standard'
                        
                        XLim=frame(:,1);
                        YLim=frame(:,2);
                        CLim=[min(U, [], 'all') max(U, [], 'all')];
                        Colormap=jet;
                        ColorbarStatus='on';
                        ColorbarLabel='Temperature [ºC]';
                        % ColorbarLabel='Temperature [K]';
                        Title='Temperature';
                        EdgeColor='none';
                        FaceColor='interp';
                        
                        for m=1:N+1
                            CData{m}=U(:,m);
                        end
                        
                end
                
                
                % Set static parameters
                PlotParam = {'Vertices',        coord;...
                             'Faces',           connect(:,1:3);...
                             'FaceVertexCData', CData{1};...
                             'FaceColor',       FaceColor;...
                             'EdgeColor',       EdgeColor;...
                             'LineWidth',       1};
                
                AxisParam = {'XLim',            XLim;...
                             'YLim',            YLim;...
                             'CLim',            CLim;...
                             'Colormap',        Colormap;...
                             'DataAspectRatio', [1 1 1];...
                             'XLabel.String',          'x coordinate [m]';...
                             'YLabel.String',          'y coordinate [m]'};
                
                ColorbarParam = {'Visible',         ColorbarStatus;...
                                 'Label.String',    ColorbarLabel;...
                                 'Label.FontSize',  11};
                
                % Set dynamic parameters
                DynamicParam={'CData',      CData;...
                              'Title',      Title};
                   
        end
        
        
    end

%% Parameters for surf plots
    function[PlotParam, AxisParam, ColorbarParam, DynamicParam]=plotSurf(Param, Mesh, Data, Solution, type)
  
        % Retrieve mesh parameters
        % Coordinate matrix of interpolation nodes
        coord_inter=Mesh.coord_inter;
        % Frame surrounding the structure
        frame=Mesh.frame;
        
        % Retrieve data parameters
        % Number of sub-intervals in time
        N=Data.Discretization.N;
        % Number of grid points along each direction for visualization
        S=Data.Discretization.S;
        
        % Retrieve solution
        sigma=Solution.sigma;
        % Total number of interpolation nodes
        nb_nodes_inter=size(coord_inter,1);
        
        % Initialization
        XLim_stress = frame(:,1);
        YLim_stress = frame(:,2);
        [Xm,Ym]=meshgrid(linspace(XLim_stress(1), XLim_stress(2), S(1)), linspace(YLim_stress(1), YLim_stress(2), S(2)));
        
        switch type
            case 'von_mises'
                
                % Von Mises stress
                sigma_vonMises=sqrt(sigma(1,:,:).^2-sigma(1,:,:).*sigma(2,:,:)+sigma(2,:,:).^2+3*sigma(3,:,:).^2);
                sigma_select=reshape(sigma_vonMises, nb_nodes_inter,N+1);
                Title_stress='Von Mises stress field';
                
            case 'sigma_xx'
                
                % Abolute value of sigma_xx
                sigma_xx=abs(sigma(1,:,:));
                sigma_select=reshape(sigma_xx, nb_nodes_inter,N+1);
                Title_stress='|\sigma_{xx}|';
                
            case 'sigma_yy'
                
                % Absolute value of sigma_yy
                sigma_yy=abs(sigma(3,:,:));
                sigma_select=reshape(sigma_yy, nb_nodes_inter,N+1);
                Title_stress='|\sigma_{yy}|';
                
            case 'sigma_xy'
                
                % Absolute value of sigma_xy
                sigma_xy=abs(sigma(2,:,:));
                sigma_select=reshape(sigma_xy, nb_nodes_inter,N+1);
                Title_stress='|\sigma_{xy}|';
                
            case 'frobenius'
                
                % Frobenius norm of sigma
                sigma_fro=sqrt(sigma(1,:,:).^2+2*sigma(2,:,:).^2+sigma(3,:,:).^2);
                sigma_select=reshape(sigma_fro, nb_nodes_inter,N+1);
                Title_stress='Frobenius norm of \sigma';
                       
        end
        
        
        % If the video is viewed, compute all stresses for all time
        % steps and set the maximum to the maximum among interpolated
        % values
        Zm=zeros(S(2), S(1), N+1);
        if Param.visua_param.view_video
            
            
            for m=1:N+1
                Vm=scatteredInterpolant(coord_inter(:,1), coord_inter(:,2), sigma_select(:,m));
                Zm(:,:,m)=Vm(Xm,Ym);
            end
            
            CLim_stress=[min(min(min(Zm))), max(max(max(Zm)))];

        else % Only a time snapshot. Set the maximum stress to the "raw" maximum stress over all time steps
            
            Zm=zeros(S(2), S(1), 1);
            time_step=Param.visua_param.snapshot;
            Vm=scatteredInterpolant(coord_inter(:,1), coord_inter(:,2), sigma_select(:,time_step));
            Zm(:,:, time_step)=Vm(Xm,Ym);
            CLim_stress = [min(min(sigma_select)) max(max(sigma_select))];
        end
    
    
      
        % Set static parameters
        PlotParam = {'XData',     Xm;...
                     'YData',     Ym;...
                     'ZData',     Zm(:,:,1);...
                     'EdgeColor',    'interp'};
        
        AxisParam = {'XLim',            XLim_stress;...
                     'YLim',            YLim_stress;...
                     'CLim',            CLim_stress;...
                     'View',            [0 90];...
                     'Colormap',        jet;...
                     'DataAspectRatio', [1 1 1];...
                     'XLabel.String',   'x coordinate [m]';...
                     'YLabel.String',   'y coordinate [m]'};
        
        ColorbarParam = {'Visible',         'on';...
                         'Label.String',    'Stress [Pa]';...
                         'Label.FontSize',  11};
        
        % Set dynamic parameters
        DynamicParam={'ZData',   num2cell(Zm, [1,2]);...
                      'Title',      Title_stress};
  
    end

%% Parameters for quiver plots
    function[PlotParam, AxisParam, ColorbarParam, DynamicParam]=plotQuiver(Mesh, Data, Solution, type)
        
        % Retrieve mesh parameters
        % Coordinate matrix
        coord=Mesh.coord_mat;
        % Frame surrounding the structure
        frame=Mesh.frame;
        
        % Retrieve data parameters
        % Number of sub-intervals in time
        N=Data.Discretization.N;
        
        % Retrieve solution
        Q_sort=Solution.Q_sort;
        
        % Initialization
        U=cell(1,N+1);
        V=cell(1,N+1);
        
        switch type
            case 'standard'
                % Standard visualization of fluxes
                Title='Flux field';
                
                for m=1:N+1
                    U{m}=Q_sort(:,2*m-1);
                    V{m}=Q_sort(:,2*m);
                end
                
        end
        
        
        % Set static parameters
        PlotParam = {'XData',        coord(:,1);...
                     'YData',        coord(:,2);...
                     'LineWidth',       1};
        
        AxisParam = {'XLim',            frame(:,1);...
                     'YLim',            frame(:,2);...
                     'DataAspectRatio', [1 1 1];...
                     'XLabel.String',          'x coordinate [m]';...
                     'YLabel.String',          'y coordinate [m]'};
        
        ColorbarParam = {'Visible',         'off';...
                         'Label.String',    '';...
                         'Label.FontSize',  11};
        
        % Set dynamic parameters
        DynamicParam={'UData',      U;...
                      'VData',      V;...
                      'Title',      Title};
        
        
    end






%% Setting and updating graphic properties
    function[object]=setObject(object, param)
        
        argument=param(:,1);
        value=param(:,2);
        
        for i=1:length(argument)         
            if contains(argument{i}, '.')
                field=extractBetween(argument{i}, 1, '.');
                subfield=extractAfter(argument{i}, '.');
                object.(field{1}).(subfield)=value{i};
            else 
                object.(argument{i})=value{i};
            end
        end
        
    end

    function[plot, axis]=updateObject(plot, axis, param, step, N)
        
        argument=param(:,1);
        value=param(:,2);
        
        for i=1:length(argument)-1
            plot.(argument{i})=value{i}{step};
        end
        
        axis.Title.String=[value{end} ' at step ' num2str(step) ' out of ' num2str(N+1)];
    end
end