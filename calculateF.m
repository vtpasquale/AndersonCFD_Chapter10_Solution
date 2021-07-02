function F = calculateF(primitives,dx,dy,direction)
% Calculate F flux array from primitive variables
%
% INPUTS
% primitives = [Primitives] Domain primitives
% dx         = [double] x grid spacing
% dy         = [double] y grid spacing
% direction  = [char] Finite difference direction = 'forward' or 'backward'
%
% OUTPUTS
% F = [ny,ny,4 double] Solution array (See Anderson Eq. 10.4d [typo in term 3])
%
% ny and nx are the number of x and y grid points, respectively. It is
% assumed that the 2D grid is created using the meshgrid() function.
%
% From Anderson page 454:
% To maintain second-order accuracy, the x-derivative terms appearing in E
% are differenced in the opposite direction to that used for dE/dx, while
% the y-derivative terms are approximated using central differences.
% Likewise, the y-derivative terms appearing in F are differenced in the
% opposite direction to that used for dF/dy, while the x-derivative terms
% in F are central differenced. 

% Anthony Ricciardi
% July 2021

% primitives
[u,v,p,T] = primitives.deal();
[mu,lambda,k] = primitives.getMuLambdaK();

r = primitives.r;
Et = primitives.Et;

% gradients
switch direction
    case 'forward'
        dudy = ddyf(u,dy);
        dvdy = ddyf(v,dy);
        dTdy = ddyf(T,dy);
    case 'backward'
        dudy = ddyf(u,dy);
        dvdy = ddyf(v,dy);
        dTdy = ddyf(T,dy);
    otherwise
        error('direction not allowed')
end
dudx = ddxc(u,dx);
dvdx = ddxc(v,dx);
% dTdx = ddxc(T,dx);

% stresses and heat fluxes
% txx = lambda.*(dudx + dvdy)+ 2*mu.*(dudx);
tyy = lambda.*(dudx + dvdy)+ 2*mu.*(dvdy);
txy = mu.*(dudy + dvdx);
% qx = -k.*dTdx;
qy = -k.*dTdy;

F = zeros(size(r,1),size(r,2),4);
F(:,:,1) = r.*v;
F(:,:,2) = r.*u.*v - txy;
F(:,:,3) = r.*v.^2 + p - tyy;
F(:,:,4) = (Et+p).*v- u.*txy - v.*tyy + qy;