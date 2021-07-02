function dfdx = ddxf(f,dx)
% Calculates forward-difference derivative wrt x on even-spaced 2D grid
% 
% INPUTS
% f  = [ny,nx,: double] array with function values
% dx = [double] grid spacing
%
% OUTPUTS
% dfdx = [ny,nx,: double] array with function derivatives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

[ny,nx,n3] = size(f);
inX = 1:nx-1;
dfdx = zeros(ny,nx,n3);

% forward difference interior points
dfdx(:,inX,:) = (f(:,inX+1,:) - f(:,inX,:))./dx;

% Set boundary point derivatives equal to adjacent points. Boundary points 
% are essentially unused in the solution, but they are provided for array 
% size consistency.
dfdx(:,nx,:) = dfdx(:,nx-1,:);