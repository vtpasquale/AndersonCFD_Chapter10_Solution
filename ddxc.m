function dfdx = ddxc(f,dx)
% Calculate central-difference derivative wrt x on even-spaced 2D grid 
% 
% INPUTS
% f  = [ny,nx double] array with function values
% dx = [double] grid spacing
%
% OUTPUTS
% dfdx = [ny,nx double] array with function derivatives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

[ny,nx] = size(f);
inX = 2:nx-1; % interior x

% initialize
dfdx = zeros(ny,nx);

% central difference interior points  [O(dx^2)]
dfdx(:,inX) = (f(:,inX+1) - f(:,inX-1))./(2*dx);

% % one-sided difference boundary points  [O(dx)]
dfdx(:,1) = (f(:,2) - f(:,1))./(dx);
dfdx(:,nx) = (f(:,nx) - f(:,nx-1))./(dx);

% one-sided difference boundary points [O(dx^2)]
% dfdx(:,1) = (-f(:,3) + 4*f(:,2) - 3*f(:,1))./(2*dx);
% dfdx(:,nx) = (3*f(:,nx) - 4*f(:,nx-1) + f(:,nx-2))./(2*dx);