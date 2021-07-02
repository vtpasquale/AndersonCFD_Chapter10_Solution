% Script to graphically check finite difference functions
%
% Anthony Ricciardi
% July 2021	
clear all; close all; clc

% % preview grid
[x,y] = meshgrid(linspace(0,2*pi,70),linspace(0,2*pi,70));
[nx,ny] = size(x);
dx = x(1,2)-x(1,1);
dy = y(2,1)-y(1,1);

% analytic function and derivatives
z = pi*sin(4*x-2) - 3/2*pi*cos(3*y-1);
dzdx = 4*pi*cos(4*x-2);
dzdy = 3*3/2*pi*sin(3*y-1);

%% central difference
dzdx_fd = ddxc(z,dx);
dzdy_fd = ddyc(z,dy);

% plot results
figure(1)
surf(x,y,dzdx)
hold on
surf(x,y,dzdx_fd)
xlabel('x')
ylabel('y')
hold off

figure(2)
surf(x,y,dzdy)
hold on
surf(x,y,dzdy_fd)
xlabel('x')
ylabel('y')
hold off

%% forward difference
dzdx_ffd = ddxf(z,dx);
dzdy_ffd = ddyf(z,dy);

figure(3)
surf(x,y,dzdx)
hold on
surf(x,y,dzdx_ffd)
xlabel('x')
ylabel('y')
hold off

figure(4)
surf(x,y,dzdy)
hold on
surf(x,y,dzdy_ffd)
xlabel('x')
ylabel('y')
hold off

%% backward difference
dzdx_bfd = ddxb(z,dx);
dzdy_bfd = ddyb(z,dy);

figure(5)
surf(x,y,dzdx)
hold on
surf(x,y,dzdx_bfd)
xlabel('x')
ylabel('y')
hold off

figure(6)
surf(x,y,dzdy)
hold on
surf(x,y,dzdy_bfd)
xlabel('x')
ylabel('y')
hold off