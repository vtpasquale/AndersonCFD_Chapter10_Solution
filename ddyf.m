function dfdy = ddyf(f,dy)
% Calculates forward-difference derivative wrt y on even-spaced 2D grid
% 
% INPUTS
% f  = [ny,nx,: double] array with function values
% dy = [double] grid spacing
%
% OUTPUTS
% dfdy = [ny,nx,: double] array with function derivatives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

[ny,nx,n3] = size(f);
inY = 1:ny-1;
dfdy = zeros(ny,nx,n3);

% forward difference interior points
dfdy(inY,:,:) = (f(inY+1,:,:) - f(inY,:,:))./dy;

% Set boundary point derivatives equal to adjacent points. Boundary points 
% are essentially unused in the solution, but they are provided for array 
% size consistency.
dfdy(ny,:,:) = dfdy(ny-1,:,:);