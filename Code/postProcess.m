function[Solution, varargout]=postProcess(Mesh, Data, Solution, varargin)

% postProcess: Computes quantities in addition to the temperatures. These
% might include fluxes or the H1 error in case the exact solution is known.
% This last case is used to verify experimentially the convergence order
% INPUT:
% Mesh:     Structure containing all the mesh parameters
% Data:     Structure containing all the data parameters
% Solution: Structure containing the solution and its time derivative
% OPTIONAL INPUT:
% types:    Additional quantities which should be computed
% Uex:      Exact solution
% gradUex1: First component of the gradient of the exact solution
% gradUex2: Second component of the gradient of the exact solution
% OUTPUT:
% Solution: Updated solution structure containing post-processed quantities
% OPTIONAL OUTPUT:
% error:    Computed H1 error

if nargin >3
    types=varargin{1};
else
    types={'flux'};
end

% Pre-process the data to compute invariant matrices
Data.Model.postprocess=types;
[Data]=preProcess(Mesh, Data);

for k=1:length(types)
    type=types{k};
    
    switch type
        case 'flux'
            [Solution]=getFlux(Mesh, Data, Solution);
            
        case 'error'
            if length(varargin)<4
                error('Exact solution, first and second component of the gradient must be provided')
            end
            [varargout{1}]=getError(Mesh, Data, Solution, varargin{2}, varargin{3}, varargin{4});
            
        otherwise
            error('Unrecognized quantity')
            
    end
end
end