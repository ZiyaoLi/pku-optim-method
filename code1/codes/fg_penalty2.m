function [ f,g ] = fg_penalty2( x,n )
% FG_PENALTY2.M function calculates the Penalty II Function value and 
% gradient on given point x and argument n.
% 
% Input:
% x:        x in f, n-vector.
% n:        Argument and problem scale, scalar in N with n >= 2.
% 
% Output:
% f:        The calculated function value f(x), scalar.
% g:        The calculated gradient vector g(x), n-vector.
% 
% Call:
% [f, g] = fg_penalty2(x,n)

% Date:     Apr 13th, 2018
% Creator:  BroC



a_sqrt = sqrt(1e-5);
r_bf = @(x,i) a_sqrt*(exp(x(i)/10)+exp(x(i-1)/10)-exp(i/10)-exp((i-1)/10));
r_aft = @(x,i) a_sqrt*(exp(x(i-n+1)/10)+exp(-0.1));
r1 = x(1) - 0.2;
r2n = (x' .^ 2) * (n:-1:1)' - 1;
f = r1^2 + r2n^2;
g = (4 * r2n) .* x .* (n:-1:1)';
g(1) = g(1) + (2 * r1);
for i=2:n
    r = feval(r_bf,x,i);
    f = f + r^2;
    g(i) = g(i) + r * a_sqrt / 5 * exp(x(i) / 10);
    g(i-1) = g(i-1) + r * a_sqrt / 5 * exp(x(i-1) / 10);
end
for i=n+1:2*n-1
    r = feval(r_aft,x,i);
    f = f + r^2;
    g(i-n+1) = g(i-n+1) + r * a_sqrt / 5 * exp(x(i-n+1) / 10);
end

end

