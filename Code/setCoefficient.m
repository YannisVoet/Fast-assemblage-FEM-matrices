function[Data]=setCoefficient(Data, varargin)

% setCoefficient: Initializes the model coefficients
% OPTIONAL INPUT:
% d:        "d" coefficient of the PDE
% c:        "c" coefficient of the PDE
% a:        "a" coefficient of the PDE
% der_d:    Derivative of the "d" coefficient with respect to temperature
% der_c:    Derivative of the "c" coefficient with respect to temperature
% der_a:    Derivative of the "a" coefficient with respect to temperature
% f:        Function "f" of the PDE
% OUTPUT:
% Data:     Data structure with initialized coefficients

if mod(length(varargin), 2)~=0
    error('Missing a name or value')
end
labels_in=varargin(1:2:end);
values_in=varargin(2:2:end);

available_coeff={'d', 'c', 'a', 'der_d', 'der_c', 'der_a', 'f'};

for i=1:length(labels_in)
    arg=lower(labels_in{i});
    val=values_in{i};
    
    if any(ismember(available_coeff, arg))
        
        
        switch arg
            case {'d', 'c', 'a', 'der_d', 'der_c', 'der_a'} % Only temperature dependent
                if isa(val,'double')
                    
                    if val ~= 0
                        Data.Coefficient.(arg).label='const';
                    else
                        Data.Coefficient.(arg).label='zero';
                    end
                    Data.Coefficient.(arg).function=@(T) val;
                    
                elseif isa(val,'function_handle')
                    Data.Coefficient.(arg).label=detectDependency1(val);
                    Data.Coefficient.(arg).function=val;
                else
                    error('Input parameters must be either constants or function handles')
                end
                
            case 'f'
                if isa(val,'double') % Only space-time dependent
                    Data.Coefficient.f.label='const';
                    Data.Coefficient.f.function=@(x,t) val;
                elseif isa(val,'function_handle')
                    Data.Coefficient.f.label=detectDependency2(val);
                    Data.Coefficient.f.function=@(x,t) val(x,t);
                else
                    error('Input parameters must be either constants or function handles')
                end
                
        end
        
    else
        error('Unrecognized coefficient')
    end
end


    function[label]=detectDependency1(func)
        
        string=func2str(func);
        
        if ~strcmp(string(1:4), '@(T)')
            error('Coefficient can only be constant or temperature dependent')
        end
        
        if contains(string(5:end), 'T')
            label='temp';
        elseif contains(string(5), '0')
            label='zero';
        else
            label='const';
        end
    end

    function[label]=detectDependency2(func)
        
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
end