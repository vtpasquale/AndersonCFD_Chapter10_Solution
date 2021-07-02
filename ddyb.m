function dfdy = ddyb(f,dy)
% Calculates backward-difference derivative wrt y on even-spaced 2D grid
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
inY = 2:ny; 
dfdy = zeros(ny,nx,n3);

% backward difference interior points
dfdy(inY,:,:) = (f(inY,:,:) - f(inY-1,:,:))./dy;

% Set boundary point derivatives equal to adjacent points. Boundary points 
% are essentially unused in the solution, but they are provided for array 
% size consistency.
dfdy(1,:,:) = dfdy(2,:,:);