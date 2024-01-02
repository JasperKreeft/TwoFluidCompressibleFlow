function [NX,NY,dx,dy,x,y]=grids

X = 4.1;
Y = 1.0;

NX = 161; %81;
NY = 40;% 20;

dx = X/NX;
dy = Y/NY;

x = ((1:NX)-.5)*dx;
y = Y-((1:NY)-.5)*dy;