function U=calculateU(primitives)
% Calculate solution array from primitive variables
%
% INPUTS
% primitives = [Primitives] Domain primitives
%
% OUTPUTS
% U = [ny,ny,4 double] Solution array (See Anderson Eq. 10.4b)
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

U = zeros(size(primitives.u,1),size(primitives.u,2),4);
U(:,:,1) = primitives.r;
U(:,:,2) = U(:,:,1).*primitives.u; % r*u
U(:,:,3) = U(:,:,1).*primitives.v; % r*v
U(:,:,4) = primitives.Et;
