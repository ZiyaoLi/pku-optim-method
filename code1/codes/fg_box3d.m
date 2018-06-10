function [ f,g ] = fg_box3d( x,m )
% FG_BOX3D.M function calculates the Box Three-Dimensional Function value
% and gradient on given point x and argument m.
%
% Input:
% x:        x in f, n-vector.
% m:        Argument, scalar \in N with m >= 3.
% 
% Output:
% f:        The calculated function value f(x), scalar.
% g:        The calculated gradient vector g(x), n-vector.
% 
% Call:
% [f, g] = fg_box3d(x,m)

% Date:     Apr 13th, 2018
% Creator:  BroC


% function handle of r_t(x)
rf = @(x,t) exp(-t*x(1))-exp(-t*x(2))-(exp(-t)-exp(-10*t))*x(3);

% initialization
f = 0;
g = [0 0 0]';

for i=1:m
    r = feval(rf, x, 0.1*i);
    f = f + r^2;
    g(1) = g(1) - 0.2 * i * r * exp(-0.1 * i * x(1));
    g(2) = g(2) + 0.2 * i * r * exp(-0.1 * i * x(2));
    g(3) = g(3) - 2 * r * (exp(-0.1 * i) - exp(-i));
end
end


