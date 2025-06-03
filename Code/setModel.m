function[Data]=setModel(varargin)

% setModel: Initializes the numerical model
% OPTIONAL INPUT:
% model:        Specifies which model should be used
%               Default: transient_heat
% submodel:     Specifies which submodel from the specified model should be
%               used
%               Default: standard
% type:         Specifies whether the model is linear or nonlinear
%               Default: linear
% postprocess:  Quantity which should be computed in addition to the
%               temperatures
%               Default: none
% OUTPUT:
% Data:         Data structure with initialized model

if mod(length(varargin), 2)~=0
    error('Missing a name or value')
end
labels_in=varargin(1:2:end);
values_in=varargin(2:2:end);

models={'transient_heat'};
submodels={'standard'};
types={'linear', 'nonlinear'};
postprocess={'none', 'flux', 'error'};

% Set default Parameters
Data.Model.name=models{1};
Data.Model.submodel=submodels{1};
Data.Model.type=types{1};
Data.Model.postprocess=postprocess{1};

% Define constants
% Stefan Boltzman constant W/(m^2*K^4)
Data.Constants.sigma=5.670373e-8;

for i=1:length(labels_in)
    arg=lower(labels_in{i});
    val=lower(values_in{i});
    
    switch arg
        case 'model'
            if any(ismember(models, val))
                Data.Model.name=val;
            else
                error('Unrecognized model')
            end
            
        case 'submodel'
            if any(ismember(submodels, val))
                Data.Model.submodel=val;
            else
                error('Unrecognized submodel')
            end
            
        case 'type'
            if any(ismember(types, val))
                Data.Model.type=val;
            else
                error('Unrecognized type')
            end
            
        case 'postprocess'
            if any(ismember(postprocess, val))
                Data.Model.postprocess=val;
            else
                error('Unrecognized quantity')
            end
            
            
        otherwise
            error('Unrecognized argument')
    end
    
end
end