function [ f,g,G ] = sym_box3d( m,Hesse )
% SYM_BOX3D calculates the function and gradient symbolic expression f(x;m) 
% and g(x;m), with x = [x1 x2 x3] the three-dimensional symbol.
% 
% Input:
% m:        The function argument m, scalar in N with m >= 3.
% Hesse:    Whether the symbolic expression of the function's Hesse matrix 
%           should be calculated, logical.
% 
% Output:
% f:        The calculated symbolic expression of f.
% g:        The calculated symbolic expression of g.
% G:        The calculated symbolic expression of G.
% 
% Call
% [ f,g ] = sym_box3d( m,0 )
% [ f,g,G ] = sym_box3d( m,1 )

% Date:     Apr 11th, 2018
% Creator:  BroC


X = sym('x', [1,3]);
r = @(x,t) exp(-t*x(1))-exp(-t*x(2))-(exp(-t)-exp(-10*t))*x(3);
f = 0;
for i=1:m
    ri = feval(r, X, 0.1*i);
    f = f + ri^2;
end
for i=1:3
    g(i,1) = diff(f, X(i));
end
if Hesse
    for i=1:3
        for j=1:3
            G(i,j) = diff(g(i,1),X(j));
        end
    end
end
end


