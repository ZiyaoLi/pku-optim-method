function [ f,r,J ] = frj_fr( x, calc_j )
% FRJ_FR.M function calculates the Freudenstein and Roth Function value,
% residuals and Jacobian matrix on given point x.
%
% Input:
% x:        x in f, 3-vector.
% calc_j:   Whether should J be calculated, logical.
% 
% Output:
% f:        The calculated function value f(x), scalar.
% r:        The calculated residuals r(x), m-vector
% J:        The calculated Jacobian matrix J(x), (m*n)-matrix.
% 
% Call:
% [f r] = frj_fr(x, false)
% [f r J] = frj_fr(x, true)

% Date:     Jun 8th, 2018
% Creator:  BroC


n = length(x);
if n ~= 2
    error('wrong variable length: x must be a 2-vector.')
end

% initialization
m = 2;
r = zeros(m, 1);
if calc_j
    J = zeros(m, n);
end

x1 = x(1);
x2 = x(2);
t11 = 5 - x2;
t12 = t11 * x2 - 2;
r(1) = t12 * x2 + x1 - 13;
t21 = x2 + 1;
t22 = t21 * x2 - 14;
r(2) = t22 * x2 + x1 - 29;
if calc_j
    J(1, 1) = 1;
    J(1, 2) = t12 + x2 * (t11 - x2);
    J(2, 1) = 1;
    J(2, 2) = t22 + x2 * (t21 + x2);
end

f = r' * r;

end

