function[beta, x, H, k] = GMRES(linearOperator, b, x0, maxiter, reorth_tol, tol_error)

% GMRES: Basic implementation of the Generalized minimal residual method
% INPUT:
% linearOperator:   Function handle for a matrix-vector multiplication
% b:                Right-hand side
% x0:               Starting vector
% maxiter:          Maximum number of iterations
% reorth_tol:       Reorthogonalization tolerance
% error_tol:        Error tolerance
% OUTPUT:
% beta:             Vector containing the errors
% x:                Approximate solution
% H:                Upper Hessenberg matrix
% k:                Number of iterations

r0=b-linearOperator(x0);
n=length(r0);
U=zeros(n, maxiter);
H=zeros(maxiter+1, maxiter);
beta=zeros(maxiter+1,1);
beta(1)=norm(r0);
U(:,1)=r0/beta(1);


for k=1:maxiter
    w=linearOperator(U(:,k));
    h=U(:,1:k)'*w;
    u_tilde=w-U(:,1:k)*h;
    
    if norm(u_tilde-w)<reorth_tol
        h_hat=U(:,1:k)'*u_tilde;
        h=h+h_hat;
        u_tilde=u_tilde-U(:,1:k)*h_hat;
    end
    
    h_tilde=norm(u_tilde);
    U(:,k+1)=u_tilde/h_tilde;
    H(1:(k+1) ,k)=[h; h_tilde];
    
    % Be carefull: parentheses are crucial
    y=beta(1)*(H(1:k+1,1:k)\eye(k+1,1));
    beta(k+1)=norm(beta(1)*eye(k+1,1)-H(1:k+1,1:k)*y);
    
    if beta(k+1)<tol_error*beta(1)
        break
    end
    
end
x=x0+U(:,1:k)*y;
H=H(1:k,1:k);
end