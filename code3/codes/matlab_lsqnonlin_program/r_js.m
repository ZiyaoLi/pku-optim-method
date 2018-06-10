function [ r ] = r_js( x )

% Date:     Jun 10th, 2018
% Creator:  BroC

m = 20;
r = zeros(m, 1);

for i=1:m
    r(i) = 2 + 2 * i - exp(i * x(1)) - exp(i * x(2));
end

end

