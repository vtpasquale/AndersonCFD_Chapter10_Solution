% This script sets up the solution to supersonic flow over a flat plate as 
% described in Ref. [1] Chapter 10. Figures 10.10-10.15 are replicated.
% Minor deviations are expected because the author does not get into 
% details such as order of accuracy used for boundary derivative 
% computations and selected Courant number used for the reference solutions.
%
% Anthony Ricciardi
% July 2021
%
% [1] John D. Anderson, Jr., Computational Fluid Dynamics, McGraw-Hill 
% Education; 1st edition. 1995.
%
clear all; close all; clc

%% Set up 
% Plate length
lhori = 0.00001; % m

% Courant number
K = 0.8;

% grid size and max iterations
nx = 70;
ny = 70; 
maxiter = 10000;

% Inflow conditions
Minf = 4;
pinf = 101325; % Pa
Tinf = 288.16; % Kelvin
inflow = Primitives(0,0,pinf,Tinf);
inflow.u =Minf*inflow.a; % m/s
Reinf = inflow.calculateReynoldsNumber(lhori);

% boundary layer size & vertical domain size
delta = 5*lhori/sqrt(Reinf);
lvert = 5*delta;

% grid
[x,y] = meshgrid(linspace(0,lhori,nx),linspace(0,lvert,ny));

%% Set initial conditions
% Set all intial values to inflow values. Conditions on boundaries updated
% by solveMacCormack().
primitives = Primitives(inflow.u*ones(ny,nx),...
                        inflow.v*ones(ny,nx),...
                        inflow.p*ones(ny,nx),...
                        inflow.T*ones(ny,nx));

%% Solve two wall temperature conditions

Tw_Tinf =  1.0; % constant wall temperature
constantTw = solveMacCormack(primitives,inflow,Tw_Tinf,K,x,y,maxiter);

Tw_Tinf = -1.0; % adiabatic
adiabaticTw = solveMacCormack(primitives,inflow,Tw_Tinf,K,x,y,maxiter);

%% Figure 10.10a - pressure along plate surface
figure(1)
plot(1:nx,constantTw.p(1,:)./pinf,'o-',...
     1:nx,adiabaticTw.p(1,:)./pinf,'v-',...
     'MarkerSize',3)
xlabel('Horizontal Grid Number')
ylabel('p/p_{\infty}')
title('Fig. 10.10(a) - Plate Surface Pressure')
legend('Constant T_w','Adiabatic')
grid on

%% Figure 10.10b - pressure at outflow
yBar = y(:,end)./x(1,end) * sqrt(Reinf);

figure(2)
plot(constantTw.p(:,end)./pinf,yBar,'o-',...
     adiabaticTw.p(:,end)./pinf,yBar,'v-',...
     'MarkerSize',3)
xlabel('p/p_{\infty}')
ylabel('$\bar{y}$','interpreter','latex')
title('Fig. 10.10(b) - Outflow Pressure')
legend('Constant T_w','Adiabatic')
grid on

%% Figure 10.12(a) - temperature at outflow
figure(4)
plot(constantTw.T(:,end)./Tinf,yBar,'o-',...
     adiabaticTw.T(:,end)./Tinf,yBar,'v-',...
     'MarkerSize',3)
xlabel('T/T_{\infty}')
ylabel('$\bar{y}$','interpreter','latex')
title('Fig. 10.12(a) - Outflow Temperature')
legend('Constant T_w','Adiabatic')
grid on

%% Figure 10.14(a)- u at outflow
figure(5)
plot(constantTw.u(:,end)./inflow.u,yBar,'o-',...
     adiabaticTw.u(:,end)./inflow.u,yBar,'v-',...
     'MarkerSize',3)
xlabel('u/u_{\infty}')
ylabel('$\bar{y}$','interpreter','latex')
title('Fig. 10.14(a) - Outflow X Velocity')
legend('Constant T_w','Adiabatic','location','NorthWest')
grid on

%% Figure 10.15 - M at outflow
figure(6)
plot(constantTw.u(:,end)./constantTw.a(:,end),yBar,'o-',...
     adiabaticTw.u(:,end)./adiabaticTw.a(:,end),yBar,'v-',...
     'MarkerSize',3)
xlabel('M_{local}')
ylabel('$\bar{y}$','interpreter','latex')
title('Fig. 10.15 - Outflow Local Mach Number')
legend('Constant T_w','Adiabatic','location','NorthWest')
grid on

%% Bonus Plot
figure(7)
contourf(x,y./delta,sqrt(adiabaticTw.u.^2+adiabaticTw.v.^2)./adiabaticTw.a,6)
% contourf(x,y./delta,constantTw.Et)
% contourf(x,y./delta,constantTw.e)
% contourf(x,y./delta,constantTw.r.*constantTw.u)
% contourf(x,y./delta,.5*constantTw.r.*(constantTw.u.^2+constantTw.v.^2))
% contourf(x,y./delta,constantTw.u.^2+constantTw.v.^2)
xlabel('x (m)')
ylabel('y/\delta')
title('Mach Number Contours - Adiabatic Wall')
view(0,90)
colorbar

%%
% [txx,tyy,txy,qx,qy]=postStressAndHeatFlux(adiabaticTw,x,y);
% contourf(x,y./delta,qy)
% xlabel('x (m)')
% ylabel('y/\delta')
% view(0,90)
% colorbar

