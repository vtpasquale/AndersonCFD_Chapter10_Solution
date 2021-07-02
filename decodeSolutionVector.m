function primitives=decodeSolutionVector(U)
% Calculate primitive variables from solution vector
%
% INPUTS
% U = [ny,ny,4 double] Solution array (See Anderson Eq. 10.4b)
%
% OUTPUTS
% primitives = [Primitives] Domain primitives
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.

% Anthony Ricciardi
% July 2021

r = U(:,:,1);
u = U(:,:,2)./r;
v = U(:,:,3)./r;
Et = U(:,:,4);
e = Et./r - .5*(u.^2 + v.^2);
cv =  Primitives.R/(Primitives.gm-1);
T = e./cv;
p = r.*Primitives.R.*T;
primitives = Primitives(u,v,p,T);