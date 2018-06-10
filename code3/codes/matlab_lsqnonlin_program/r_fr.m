function [ r ] = r_fr( x )

% Date:     Jun 10th, 2018
% Creator:  BroC

m = 2;
r = zeros(m, 1);

r(1) = 13 - x(1) - ((5 - x(2)) * x(2) - 2) * x(2);
r(2) = 29 - x(1) - ((x(2) + 1) * x(2) - 14) * x(2);

end

