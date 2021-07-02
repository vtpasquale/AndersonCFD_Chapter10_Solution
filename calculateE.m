function E = calculateE(primitives,dx,dy,direction)
% Calculate E flux array from primitive variables
%
% INPUTS
% primitives = [Primitives] Domain primitives
% dx         = [double] x grid spacing
% dy         = [double] y grid spacing
% direction  = [char] Finite difference direction = 'forward' or 'backward'
%
% OUTPUTS
% E = [ny,ny,4 double] Solution array (See Anderson Eq. 10.4c)
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
        dudx = ddxf(u,dx);
        dvdx = ddxf(v,dx);
        dTdx = ddxf(T,dx);
    case 'backward'
        dudx = ddxb(u,dx);
        dvdx = ddxb(v,dx);
        dTdx = ddxb(T,dx);
    otherwise
        error('direction not allowed')
end
dudy = ddyc(u,dy);
dvdy = ddyc(v,dy);
% dTdy = ddyc(T,dy);

% stresses and heat fluxes
txx = lambda.*(dudx + dvdy)+ 2*mu.*(dudx);
% tyy = lambda.*(dudx + dvdy)+ 2*mu.*(dvdy);
txy = mu.*(dudy + dvdx);
qx = -k.*dTdx;
% qy = -k.*dTdy;

E = zeros(size(r,1),size(r,2),4);
E(:,:,1) = r.*u;
E(:,:,2) = r.*u.^2 + p - txx;
E(:,:,3) = r.*u.*v - txy;
E(:,:,4) = (Et+p).*u- u.*txx - v.*txy + qx;
