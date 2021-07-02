function [txx,tyy,txy,qx,qy]=postStressAndHeatFlux(primitives,x,y)
% For postprocessing stresses and heat fluxes. Central differences are used
% to calculate all internal gradients - this is different than the gradient
% calculations used for the MacCormack processing.
%
% INPUTS
% primitives = [Primitives] Domain primitives
% x          = [ny,nx double] Grid point x locations (must be uniform spacing)
% x          = [ny,nx double] Grid point y locations (must be uniform spacing)
%
% OUPUTS - [ny,nx double] stresses and heat fluxes

% Anthony Ricciardi
% July 2021

dx = x(1,2)-x(1,1);
dy = y(2,1)-y(1,1);

% Extract primitives
[u,v,~,T] = primitives.deal();
[mu,lambda,k] = primitives.getMuLambdaK();


% velocity gradients & stresses
dudx = ddxc(u,dx);
dudy = ddyc(u,dy);
dvdx = ddxc(v,dx);
dvdy = ddyc(v,dy);
txx = lambda.*(dudx + dvdy)+ 2*mu.*(dudx);
tyy = lambda.*(dudx + dvdy)+ 2*mu.*(dvdy);
txy = mu.*(dudy + dvdx);

% thermal gradients and heat fluxes
dTdx = ddxc(T,dx);
dTdy = ddyc(T,dy);
qx = -k.*dTdx;
qy = -k.*dTdy;