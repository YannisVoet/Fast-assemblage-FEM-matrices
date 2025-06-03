function[x, beta, count]=gmres_restart(linearOperator, b, x0)

% gmres_restart: Basic implementation of the Generalized minimal residual method
% INPUT:
% linearOperator:   Function handle for a matrix-vector multiplication
% b:                Right-hand side
% x0:               Starting vector
% OUTPUT:
% x:                Approximate solution
% beta:             Vector containing the errors
% count:            Number of iterations

x=x0;
beta_error=1;
count=0;

maxiter=40;
maxiter_tot=200;
reorth_tol=0.7;
tol_error=1e-7;


while beta_error>tol_error && count <= maxiter_tot
    [beta, x, H, k] = GMRES(linearOperator, b, x, maxiter, reorth_tol, tol_error);
    beta_error=beta(end);
    count=count+k;
end

if beta_error>tol_error && count > maxiter_tot
    disp(['The restarted GMRES method did not converge after ' num2str(count-k) ' iterations'])
end
end