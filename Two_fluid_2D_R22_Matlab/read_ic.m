%-------------------------------------------%
%                                           %
% Initial conditions and settings           %
%                                           %
%-------------------------------------------%

function [dx,dy,dt,T_final,W,Q,Nout]=read_ic()

global NX NY NT
global gamma1 gamma2

%-------------------------------------------%
%       Generate grid                       %
%-------------------------------------------%

X = 16.0e-2;
Y = 8.90e-2;
T_final = 1e0;

dx = X/NX;
dy = Y/NY;
dt = T_final/NT;

gamma1 = 1.4e0;
gamma2 = 1.249e0;

Nout = 1;

load 0.mat

Q = Prim2Cons(W);