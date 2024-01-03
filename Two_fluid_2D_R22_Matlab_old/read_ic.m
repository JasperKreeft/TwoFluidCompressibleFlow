%-------------------------------------------%
%                                           %
% Initial conditions and settings           %
%                                           %
%-------------------------------------------%

function [dx,dy,dt,T_final,W,Q,gamma1,gamma2,Nout]=read_ic() %#ok<STOUT>

global NX NY NT

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

%Convert from primary (W) to conservation (Q)
Q = zeros(NY,NX,8);
for i=1:NY
    for j=1:NX
        Q(i,j,3) = W(i,j,3);
        Q(i,j,4) = W(i,j,3)*W(i,j,4);
        Q(i,j,5) = W(i,j,3)*W(i,j,5);
        Q(i,j,6) = W(i,j,6)*(W(i,j,8)/(gamma1-1e0)+(1e0-W(i,j,8))/(gamma2-1e0))+0.5e0*W(i,j,3)*(W(i,j,4)^(2e0)+W(i,j,5)^(2e0));
        Q(i,j,7) = W(i,j,3)*W(i,j,7);
        Q(i,j,8) = W(i,j,6)*W(i,j,8)/(gamma1-1e0)+1e0/2e0*W(i,j,3)*(W(i,j,4)^(2e0)+W(i,j,5)^(2e0))*W(i,j,7);
    end
end