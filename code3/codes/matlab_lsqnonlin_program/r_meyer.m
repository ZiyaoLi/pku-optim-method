function [ r ] = r_meyer( x )

% Date:     Jun 10th, 2018
% Creator:  BroC

n = length(x);
if n ~= 3
    error('wrong variable length: x must be a 3-vector.')
end

m = 16;
r = zeros(m, 1);

y = [34780, 28610, 23650, 19630, 16370, ...
     13720, 11540,  9744,  8261,  7030, ...
      6005,  5147,  4427,  3820,  3307, ...
      2872];

for i=1:m
    r(i) = y(i) - x(1) * exp(x(2) / (45 + 5 * i + x(3)));
end

end

