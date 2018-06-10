function [ f,r,J ] = frj_osborne1( x, calc_j )
% FRJ_OSBORNE1.M function calculates the Osborne 1 Function value,
% residuals and Jacobian matrix on given point x.
%
% Input:
% x:        x in f, 5-vector.
% calc_j:   Whether should J be calculated, logical.
% 
% Output:
% f:        The calculated function value f(x), scalar.
% r:        The calculated residuals r(x), m-vector
% J:        The calculated Jacobian matrix J(x), (m*n)-matrix.
% 
% Call:
% [f r] = frj_osborne1(x, false)
% [f r J] = frj_osborne1(x, true)

% Date:     Jun 7th, 2018
% Creator:  BroC


n = length(x);
if n ~= 5
    error('wrong variable length: x must be a 5-vector.')
end

% initialization
m = 33;
r = zeros(m, 1);
if calc_j
    J = zeros(m, n);
end

y = [0.844, 0.908, 0.932, 0.936, 0.925, ...
     0.908, 0.881, 0.850, 0.818, 0.784, ...
     0.751, 0.718, 0.685, 0.658, 0.628, ...
     0.603, 0.580, 0.558, 0.538, 0.522, ...
     0.506, 0.490, 0.478, 0.467, 0.457, ...
     0.448, 0.438, 0.431, 0.424, 0.42, ...
     0.414, 0.411, 0.406];

for i=1:m
    t_ = 10 * (1 - i);
    etx4 = exp(t_ * x(4));
    etx5 = exp(t_ * x(5));
    r(i) = y(i) - x(1) - x(2) * etx4 - x(3) * etx5;
    if calc_j
        J(i, 1) = -1;
        J(i, 2) = -etx4;
        J(i, 3) = -etx5;
        J(i, 4) = -x(2) * t_ * etx4;
        J(i, 5) = -x(3) * t_ * etx5;
    end
end

f = r' * r;

end

