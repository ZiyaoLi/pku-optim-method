function [ u ] = u_bound( m,n,type )
% U_BOUND.M formulates pre-defined kinds of boundary constraints for a 
% 2D grid-differential problem. The number of rows and columns should be
% provided.
% 
% Input:
% m,n:      The number of rows and columns of the grid. Note that the
%           increasing rows indicates the increasing x, and columns, the y.
% type:     The pre-defined types of boundary constraints.
%               0:  A linear-increasing boundary. The calculated
%                   surface shall be a plane.
%               1:  Upper-bound z = 1 + sin [0, pi]
%                   Lower-bound z = 1 + sin [0, pi]
%                   Left-bound  z = 1 - sin [0, pi]
%                   Right-bound z = 1 - sin [0, pi]
%               2:  Upper-bound z = 1 + sin [0, pi]
%                   Lower-bound z = 1 + sin [0, pi]
%                   Left-bound  z = 1 + sin [0, pi]
%                   Right-bound z = 1 + sin [0, pi]
%               3:  Upper-bound z = 1 + sin [0, pi]
%                   Lower-bound z = 1
%                   Left-bound  z = 1 + sin [0, pi]
%                   Right-bound z = 1
%               4:  Upper-bound z = 1 + sin [0, 2pi]
%                   Lower-bound z = 1 + sin [0, 2pi]
%                   Left-bound  z = 1 + sin [0, 2pi]
%                   Right-bound z = 1 + sin [0, 2pi]
%               5:  Upper-bound z = 1 + sin [0, 3pi]
%                   Lower-bound z = 1 + sin [0, 3pi]
%                   Left-bound  z = 1 - sin [0, 3pi]
%                   Right-bound z = 1 - sin [0, 3pi]
%               6:  Upper-bound z = - sin [0, pi]
%                   Lower-bound z = 2 + sin [0, pi]
%                   Left-bound  z = 1 + sin [0.5pi, 1.5pi]
%                   Right-bound z = 1 - sin [0.5pi, 1.5pi]
% 
% Output:
% u:        The established boundary points, struct.
%               u.up, right, down, left are the regarding four boundaries
%               formulated. Note that the last point of every boundary is
%               saved in the regarding next boundary, i.e. up(n+1) is in
%               right(1), etc.
% 
% Call:
% [u] = u_bound(m,n,type)

% Date:     Apr 14th, 2018
% Creator:  BroC


switch type
    case 0
        u.up = (n+1:2*n)./n;
        u.right = (2*n+1:-(n/m):n+2)./n;
        u.down = (2:n+1)./n;
        u.left = (n:-(n/m):1)./n;
    case 1
        u.up = 1 + sin((0:n-1)./n*pi);
        u.right = 1 - sin((0:m-1)./m*pi);
        u.down = 1 + fliplr(sin((0:n-1)./n*pi));
        u.left = 1 - fliplr(sin((0:m-1)./m*pi));
    case 2
        u.up = 1 + sin((0:n-1)./n*pi);
        u.right = 1 + sin((0:m-1)./m*pi);
        u.down = 1 + fliplr(sin((0:n-1)./n*pi));
        u.left = 1 + fliplr(sin((0:m-1)./m*pi));
    case 3
        u.up = 1 + sin((0:n-1)./n*pi);
        u.right = ones(1,m);
        u.down = 1 + fliplr(sin((0:n-1)./n*pi));
        u.left = ones(1,m);
    case 4
        u.up = 1 + sin((0:n-1)./n*2*pi);
        u.right = 1 + sin((0:m-1)./m*2*pi);
        u.down = 1 + fliplr(sin((0:n-1)./n*2*pi));
        u.left = 1 + fliplr(sin((0:m-1)./m*2*pi));
    case 5
        u.up = 1 + sin((0:n-1)./n*3*pi);
        u.right = 1 - sin((0:m-1)./m*3*pi);
        u.down = 1 + fliplr(sin((0:n-1)./n*3*pi));
        u.left = 1 - fliplr(sin((0:m-1)./m*3*pi));
    case 6
        u.up = - sin((0:n-1)./n*pi);
        u.right = 1 - sin((0:m-1)./m*pi + pi/2);
        u.down = 2 + fliplr(sin((0:n-1)./n*pi));
        u.left = 1 + fliplr(sin((0:m-1)./m*pi + pi/2));
end
end


