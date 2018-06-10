function [ f,g ] = fg_s303( x, calc_g )
% FG_BOX3D.M function calculates the s303-305 Function value
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
% [f, g] = fg_s303(x,y)

% Date:     Jun 5th, 2018
% Creator:  BroC


n = length(x);

% initialization
s = sum((1 : n)' .* x) / 2;
s2 = s * s;
f = x' * x + s2 + s2 * s2;
if calc_g
    g = 2 .* x + (s + 2 * s * s2) .* (1 : n)';
end
end


