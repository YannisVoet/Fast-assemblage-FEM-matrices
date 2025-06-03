function[V]=computeDMatrix(Data)

% computeDMatrix: Function which computes the entries of the
% global consistent mass matrix (without "d" coefficient). Restricted to
% 2D problems and when all elements are of the same type
% INPUT:
% Mesh:     Structure containing all the data parameters
% OUTPUT:
% V:        Values of non-zero entries

% Operator properties
Q=      Data.Operator.DMatrix.Q;
weights=Data.Operator.DMatrix.weights;

% Mesh properties
absDetJK=Data.MeshProperties.absDetJK;

lambda=absDetJK'.*weights;
M=Q*lambda';
V=M(:);
end