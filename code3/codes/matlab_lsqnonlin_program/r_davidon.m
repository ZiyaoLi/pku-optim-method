function [ r ] = r_davidon( x  )

% Date:     Jun 10th, 2018
% Creator:  BroC

m = 20;
r = zeros(m, 1);

eta = (1 : m) ./ 5;
eeta = exp(eta);
sineta = sin(eta);
coseta = cos(eta);

for i=1:m
    r1 = x(1) + eta(i) * x(2) - eeta(i);
    r2 = x(3) + sineta(i) * x(4) -coseta(i);
    r(i) = r1 * r1 + r2 * r2;
end

end

