function [ f,r,J ] = frj_js( x, calc_j, m )
% FRJ_JS.M function calculates the Jennrich and Sampson Function value, 
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
% [f r] = frj_js(x, false)
% [f r] = frj_js(x, false, m)
% [f r J] = frj_js(x, true)
% [f r J] = frj_js(x, true, m)

% Date:     Jun 8th, 2018
% Creator:  BroC


if nargin == 2 || isempty(m)
    m = 10;
end

n = length(x);
if n ~= 2
    error('wrong variable length: x must be a 2-vector.')
end

% initialization
r = zeros(m, 1);
if calc_j
    J = zeros(m, n);
end

for i=1:m
    e1 = -exp(i * x(1));
	e2 = -exp(i * x(2));
    r(i) = 2 + 2 * i + e1 + e2;
    if calc_j
        J(i, 1) = i * e1;
        J(i, 2) = i * e2;
    end
end

f = r' * r;

end

