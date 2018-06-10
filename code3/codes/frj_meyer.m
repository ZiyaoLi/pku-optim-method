function [ f,r,J ] = frj_meyer( x, calc_j )
% FRJ_MEYER.M function calculates the Meyer Function value,
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
% [f r] = frj_meyer(x, false)
% [f r J] = frj_meyer(x, true)

% Date:     Jun 8th, 2018
% Creator:  BroC


n = length(x);
if n ~= 3
    error('wrong variable length: x must be a 3-vector.')
end

% initialization
m = 16;
r = zeros(m, 1);
if calc_j
    J = zeros(m, n);
end

y = [34780, 28610, 23650, 19630, 16370, ...
     13720, 11540,  9744,  8261,  7030, ...
      6005,  5147,  4427,  3820,  3307, ...
      2872];

for i=1:m
    t = 45 + 5 * i;
    t3 = t + x(3);
    e23 = exp(x(2) / t3);
    x1e = x(1) * e23;
    r(i) = y(i) - x1e;
    if calc_j
        J(i, 1) = -e23;
        J(i, 2) = -x1e / t3;
        J(i, 3) = x1e * x(2) / t3 / t3;
    end
end

f = r' * r;

end

