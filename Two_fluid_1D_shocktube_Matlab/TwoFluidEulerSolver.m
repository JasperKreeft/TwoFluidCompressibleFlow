%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Compressible Two-Fluid Flow Shocktube Solver.                           %
%                                                                         %
% The code is described in:                                               %
%                                                                         %
% Jasper Kreeft & Barry Koren,                                            %
% Journal of Computational Physics volume 229 (2010) 6220-6242            %
%                                                                         %
% Written by Jasper Kreeft (2007)                                         %
% Updated by Jasper Kreeft (2010)                                         %
% Updated by Jasper Kreeft (2015)                                         %
%                                                                         %
% Contact: j.j.kreeft@gmail.com                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

WindowSize

global N L

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options
% nr = 3; N = 200; T = 5; time = 3; flux = 3; BC_L = 1; BC_R = 1; limiter = 3; ani = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh

L  = 0.5;
dx = L/N;

if limiter==0
    CFL = 0.95;
else
    CFL = 0.48;
end

x = ((1:N)-1/2)*dx-L/2;

y = [0 0.2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial conditions

InitialCondition

Q = Prim2Cons(W);
q = Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t(n)   = 0;
Lambda = 20;

Qf = zeros(5,2*N);

if ani==1
    han = waitbar(0,'Please wait...','Position',[500 500 288 60]);
elseif ani==3 || ani==4
%     figure('units','normalized','outerposition',[0.2 0.2 0.6 0.6])
    figure(winFull{:})
end

while t(n)<T

n      = n+1;
dt     = CFL*dx/Lambda;
t(n)   = t(n-1)+dt;
if ani==1
    waitbar(t(n)/T,han)
end
Lambda = [];

disp(['n=' num2str(n) ', t=' num2str(t(n))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Third-Order Runge-Kutta scheme                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dQ = zeros(5,N,3);
for stage = 1:time

w  = Cons2Prim(q);
Wf = limitersPrim(w,limiter);
if flux==1
%     [F,sL,sR,Lambda] = ExactFlux(N,Wf,gamma1,gamma2,pi1,pi2,Lambda);
    [F,sL,sR,Lambda] = ExactFlux(Wf,Lambda);
elseif flux==2
    [F,sL,sR,Lambda] = OsherFlux(Wf,Lambda);
elseif flux==3
    [F,sL,sR,Lambda] = LinearOsherFlux(Wf,Lambda);
end
S = sourceintegral(w,Wf);
F = BoundaryConditions(F,Wf,BC_L,BC_R);

dQ(1:4,:,stage) = -dt/dx*(F(1:4,2:N+1)-F(1:4,1:N));
dQ(5,:,stage)   = -dt/dx*(F(5,2:N+1)-F(5,1:N))+dt/dx*(sR+S+sL);

q = Q + RK_coeff(time,(stage-1)*3+1)*dQ(:,:,1) ...
      + RK_coeff(time,(stage-1)*3+2)*dQ(:,:,2) ...
      + RK_coeff(time,(stage-1)*3+3)*dQ(:,:,3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = q;

plotten

end %time step

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

postprocessen(W,t,x)