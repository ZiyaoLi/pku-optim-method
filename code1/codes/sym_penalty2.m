function [ f,g,G ] = sym_penalty2( n,Hesse )
% SYM_PENALTY2 calculates the function and gradient symbolic expression
% f(x;n) and g(x;n), with x = [x1 x2 ... xn] the n-dimensional symbol.
% 
% Input:
% m:        The function argument and scale n, scalar in N with n >= 2.
% Hesse:    Whether the symbolic expression of the function's Hesse matrix 
%           should be calculated, logical.
% 
% Output:
% f:        The calculated symbolic expression of f.
% g:        The calculated symbolic expression of g.
% G:        The calculated symbolic expression of G.
% 
% Call
% [ f,g ] = sym_penalty2( n,0 )
% [ f,g,G ] = sym_penalty2( n,1 )

% Date:     Apr 11th, 2018
% Creator:  BroC


a = 1e-5;
X = sym('x',[1,n]);
r_bf = @(x,i) sqrt(a)*(exp(x(i)/10)+exp(x(i-1)/10)-exp(i/10)-exp((i-1)/10));
r_aft = @(x,i) sqrt(a)*(exp(x(i-n+1)/10)+exp(-0.1));
r1 = X(1) - 0.2;
r2n = X.^2*(n:-1:1)' - 1;
f = r1^2 + r2n^2;
for i=2:n
    r = feval(r_bf,X,i);
    f = f + r^2;
end
for i=n+1:2*n-1
    r = feval(r_aft,X,i);
    f = f + r^2;
end
for i=1:n
    g(i,1) = diff(f, X(i));
end
if Hesse
    for i=1:n
        for j=1:n
            G(i,j) = diff(g(i,1),X(j));
        end
    end
end
end

