function [x,y,xc,yc,nx,ny,eps] = fiber(n,r,side,dx,dy);
n = 1;
r = 0.1;
dx = 0.02;
dy = 0.02;

nx = round((sum(r))/dx);
ny = round((sum(r))/dy);

xc = (1:nx)'*dx - dx/2;
yc = (1:ny)*dy - dy/2;
x = (0:nx)'*dx;
y = (0:ny)*dy;
eps = ones(nx,ny);

