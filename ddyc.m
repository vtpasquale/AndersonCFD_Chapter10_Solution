function dfdy = ddyc(f,dy)
% Calculates central-difference derivative wrt y on even-spaced 2D grid 
% 
% INPUTS
% f  = [ny,nx double] array with function values
% dy = [double] grid spacing
%
% OUTPUTS
% dfdy = [ny,nx double] array with function derivatives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

[ny,nx] = size(f);
inY = 2:ny-1; % interior y

% initialize
dfdy = zeros(ny,nx);

% central difference interior points  [O(dy^2)]
dfdy(inY,:) = (f(inY+1,:) - f(inY-1,:))./(2*dy);

% % one-sided difference boundary points [O(dy)]
dfdy(1,:) = (f(2,:) - f(1,:))./(dy);
dfdy(ny,:) = (f(ny,:) - f(ny-1,:))./(dy);

% one-sided difference boundary points [O(dy^2)]
% dfdy(1,:) = (-f(3,:) + 4*f(2,:) - 3*f(1,:))./(2*dy);
% dfdy(ny,:) = (3*f(ny,:) - 4*f(ny-1,:) + f(ny-2,:))./(2*dy);
