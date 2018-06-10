function [ f,g ] = fg_gencube( x, calc_g )
% FG_BOX3D.M function calculates the Generalized CUBE Function value
% and gradient on given point x and argument m.
%
% Input:
% x:        x in f, n-vector.
% calc_g:   Whether should g be calculated, logical.
% 
% Output:
% f:        The calculated function value f(x), scalar.
% g:        The calculated gradient vector g(x), n-vector.
% 
% Call:
% [f, g] = fg_gencube(x)

% Date:     Apr 13th, 2018
% Creator:  BroC


n = length(x);

% initialization
r1 = x(1) - 1;
f = r1 * r1 / 100;
g = false;
if calc_g
    g = zeros(n, 1);
    g(1) = 2 * r1 / 100;
    for i=2:n
        % pre-load elements to save calculating time of accessing memory.
        xi = x(i);
        xi_1 = x(i - 1);
        r = xi - xi_1 * xi_1 * xi_1;
        f = f + r * r;
        g(i) = g(i) + 2 * r;
        g(i - 1) = g(i - 1) - 6 * r * xi_1 * xi_1;
    end
    g = g .* 100;
else
    for i=2:n
        % pre-load elements to save calculating time of accessing memory.
        xi = x(i);
        xi_1 = x(i - 1);
        r = xi - xi_1 * xi_1 * xi_1;
        f = f + r * r;
    end
end

f = f * 100;

end


