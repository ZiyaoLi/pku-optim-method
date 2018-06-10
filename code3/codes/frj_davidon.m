function [ f,r,J ] = frj_davidon( x, calc_j, m )
% FRJ_DAVIDON.M function calculates the Davidon Function value,
% residuals and Jacobian matrix on given point x.
%
% Input:
% x:        x in f, 3-vector.
% calc_j:   Whether should J be calculated, logical.
% m:        the number of residuals generated, scalar, default=10
% 
% Output:
% f:        The calculated function value f(x), scalar.
% r:        The calculated residuals r(x), m-vector
% J:        The calculated Jacobian matrix J(x), (m*n)-matrix.
% 
% Call:
% [f r] = frj_davidon(x, false)
% [f r] = frj_davidon(x, false, m)
% [f r J] = frj_davidon(x, true)
% [f r J] = frj_davidon(x, true, m)

% Date:     Jun 8th, 2018
% Creator:  BroC


if nargin == 2 || isempty(m)
    m = 20;
end

n = length(x);
if n ~= 4
    error('wrong variable length: x must be a 4-vector.')
end

% initialization
r = zeros(m, 1);
if calc_j
    J = zeros(m, n);
end

eta = (1 : m) ./ 5;
eeta = exp(eta);
sineta = sin(eta);
coseta = cos(eta);

for i=1:m
    r1 = x(1) + eta(i) * x(2) - eeta(i);
    r2 = x(3) + sineta(i) * x(4) -coseta(i);
    r(i) = r1 * r1 + r2 * r2;
    if calc_j
        J(i, 1) = 2 * r1;
        J(i, 2) = 2 * r1 * eta(i);
        J(i, 3) = 2 * r2;
        J(i, 4) = 2 * r2 * sineta(i);
    end
end

f = r' * r;

end

