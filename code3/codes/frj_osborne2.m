function [ f,r,J ] = frj_osborne2( x, calc_j )
% FRJ_OSBORNE1.M function calculates the Osborne 2 Function value,
% residuals and Jacobian matrix on given point x.
%
% Input:
% x:        x in f, 11-vector.
% calc_j:   Whether should J be calculated, logical.
% 
% Output:
% f:        The calculated function value f(x), scalar.
% r:        The calculated residuals r(x), m-vector
% J:        The calculated Jacobian matrix J(x), (m*n)-matrix.
% 
% Call:
% [f r] = frj_osborne2(x, false)
% [f r J] = frj_osborne2(x, true)

% Date:     Jun 8th, 2018
% Creator:  BroC


n = length(x);
if n ~= 11
    error('wrong variable length: x must be an 11-vector.')
end

% initialization
m = 65;
r = zeros(m, 1);
if calc_j
    J = zeros(m, n);
end

y = [1.366, 1.191, 1.112, 1.013, 0.991, ...
     0.885, 0.831, 0.847, 0.786, 0.725, ...
     0.746, 0.679, 0.608, 0.655, 0.616, ...
     0.606, 0.602, 0.626, 0.651, 0.724, ...
     0.649, 0.649, 0.694, 0.644, 0.624, ...
     0.661, 0.612, 0.558, 0.533, 0.495, ...
     0.500, 0.423, 0.395, 0.375, 0.372, ...
     0.391, 0.396, 0.405, 0.428, 0.429, ...
     0.523, 0.562, 0.607, 0.653, 0.672, ...
     0.708, 0.633, 0.668, 0.645, 0.632, ...
     0.591, 0.559, 0.597, 0.625, 0.739, ...
     0.710, 0.729, 0.720, 0.636, 0.581, ...
     0.428, 0.292, 0.162, 0.098, 0.054];

for i=1:m
    t = (i - 1) / 10;
    t9 = t - x(9);
    t10 = t - x(10);
    t11 = t - x(11);
    t92 = t9 * t9;
    t102 = t10 * t10;
    t112 = t11 * t11;
    e5 = exp(-t * x(5));
    e96 = exp(-t92 * x(6));
    e107 = exp(-t102 * x(7));
    e118 = exp(-t112 * x(8));
    x1e = x(1) * e5;
    x2e = x(2) * e96;
    x3e = x(3) * e107;
    x4e = x(4) * e118;
    r(i) = y(i) - x1e - x2e - x3e - x4e;
    if calc_j
        J(i, 1) = -e5;
        J(i, 2) = -e96;
        J(i, 3) = -e107;
        J(i, 4) = -e118;
        J(i, 5) = x1e * t;
        J(i, 6) = x2e * t92;
        J(i, 7) = x3e * t102;
        J(i, 8) = x4e * t112;
        J(i, 9) = -2 * x2e * x(6) * t9;
        J(i, 10) = -2 * x3e * x(7) * t10;
        J(i, 11) = -2 * x4e * x(8) * t11;
    end
end

f = r' * r;

end

