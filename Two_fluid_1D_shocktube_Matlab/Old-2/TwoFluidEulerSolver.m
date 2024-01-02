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
% Updated by Jasper Kreeft (2017)                                         %
%                                                                         %
% Contact: j.j.kreeft@gmail.com                                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

WindowSize

global N

%% Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% options
nr = 3; N = 200; T = 10.12; time = 2; flux = 1; BC_L = 1; BC_R = 1; limiter = 2; ani = 3; movie = 0;

%% Domain and mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L  = 0.5;
dx = L/N;

x = ((1:N)-1/2)*dx-L/2;

y = [0 0.2];

%% Initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

InitialCondition

Q = Prim2Cons(W);

q = Q;


%% Time loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while t(n)<T

n    = n+1;
dt   = CFL*dx/Lambda;
t(n) = t(n-1)+dt;

if ani==1; waitbar(t(n)/T,han); end

Lambda = [];

if ~rem(n,10)
disp(['n=' num2str(n,'%03d') ', t=' num2str(t(n),'%4.3f') ', dt=' num2str(dt,'%3.1e')]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Third-Order Runge-Kutta scheme                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dQ = zeros(5,N,time);
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

    q = Q;
    for i=1:time
        q = q + RK_coeff{time}(stage,i)*dQ(:,:,i);
    end

end

%% New Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = q;

%% Visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotten

end %time step

%% Post-processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if movie
    close(writerObj);
    set(gcf,'visible','on')
end

W = Cons2Prim(Q);

postprocessen(W,t,x)

%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%