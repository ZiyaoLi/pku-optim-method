function [ f,g ] = fg_genbard( x, calc_g, y )
% FG_BOX3D.M function calculates the Generalized BARD Function value
% and gradient on given point x and argument m.
%
% Input:
% x:        x in f, n-vector.
% calc_g:   Whether should g be calculated, logical.
% y:        Argument, in the reference CUTEr the default:
%             y = [14 18 22 25 29 32 35 39 37 58 73 96 134 210 439] ./ 100;
% 
% Output:
% f:        The calculated function value f(x), scalar.
% g:        The calculated gradient vector g(x), n-vector.
% 
% Call:
% [f, g] = fg_genbard(x,y)

% Date:     Jun 5th, 2018
% Creator:  BroC


n = length(x);
t = n - 2;
m = length(y);

% y = [14 18 22 25 29 32 35 39 37 58 73 96 134 210 439] ./ 100;
% initialization
f = 0;
g = false;
w = min(1:m, m:-1:1);
v = m:-1:1;

if calc_g
    g = zeros(n, 1);
    for j=1:t
        % pre-load elements to save calculating time of accessing memory.
        xj = x(j);
        xj1 = x(j + 1);
        xj2 = x(j + 2);
        gj = g(j);
        gj1 = g(j + 1);
        gj2 = g(j + 2);
        for i=1:m
            vi = v(i);
            wi = w(i);
            s = 1 / (vi * xj1 + wi * xj2);
            r = y(i) - xj - i * s;
            f = f + r * r;
            gj = gj - 2 * r;
            gj1 = gj1 + 2 * r * i * vi * s * s;
            gj2 = gj2 + 2 * r * i * wi * s * s;
        end
        g(j) = gj;
        g(j + 1) = gj1;
        g(j + 2) = gj2;
    end
else
    for j=1:t
        % pre-load elements to save calculating time of accessing memory.
        xj = x(j);
        xj1 = x(j + 1);
        xj2 = x(j + 2);
        for i=1:m
            vi = v(i);
            wi = w(i);
            s = 1 / (vi * xj1 + wi * xj2);
            r = y(i) - xj - i * s;
            f = f + r * r;
        end
    end
end


