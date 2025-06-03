function[Functions]=getFunctions(order)

% getFunctions: Returns the basis functions, their jacobian matrices and
% the basis functions needed to integrate over the boundary
% INPUT:
% order:        Finite element order
% OUTPUT:
% Functions:    Structure containing all the functions


switch order
    case 1
Functions.phi=@(x) [1-x(:,1)-x(:,2) x(:,1) x(:,2)]';
 
Functions.J_phi=@(x) reshape([-ones(1,size(x,1)) zeros(1,size(x,1)) ones(1,size(x,1)) ...
                              -ones(1,size(x,1)) ones(1,size(x,1)) zeros(1,size(x,1))]', 3, 2*size(x,1));
                                                 
Functions.phi_boundary=@(x) [1-x x]';

    case 2
Functions.phi=@(x) [(1-x(:,1)-x(:,2)).*(1-2*x(:,1)-2*x(:,2)) x(:,1).*(2*x(:,1)-1) x(:,2).*(2*x(:,2)-1) 4*x(:,1).*(1-x(:,1)-x(:,2)) 4*x(:,1).*x(:,2) 4*x(:,2).*(1-x(:,1)-x(:,2))]';
                 
Functions.J_phi=@(x) reshape([-3+4*x(:,2)+4*x(:,1) -1+4*x(:,1) 0*x(:,1) 4-4*x(:,2)-8*x(:,1) 4*x(:,2) -4*x(:,2) ...
                              -3+4*x(:,1)+4*x(:,2) 0*x(:,1) -1+4*x(:,2) -4*x(:,1) 4*x(:,1) 4-4*x(:,1)-8*x(:,2)]', 6, 2*size(x,1));

Functions.phi_boundary=@(x) [(2*x-1).*(x-1) x.*(2*x-1) 4*x.*(1-x)]';
               
    otherwise
        error('Unimplemented!');
end
end

