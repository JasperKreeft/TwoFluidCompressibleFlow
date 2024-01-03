clear all
close all
clc

%Two-Fluid flow shocktube solver

%-------------------------------------------%
%-------------------------------------------%
%Program: Shocktube solver for two-fluid    %
%         flow problems                     %
%                                           %
%        author:  Jasper Kreeft             %
%          j.j.kreeft@student.tudelft.nl    %
%                                           %
% Written 23-04-2007, changed 23-04-2007    %
%-------------------------------------------%
%-------------------------------------------%

global NX NY NT

NX = 200;
NY = 100;
NT = 8000;

dQ1 = zeros(NY,NX,8);
dQ2 = zeros(NY,NX,8);
dQ3 = zeros(NY,NX,8);

% ic
% IC_Bubble
[dx,dy,dt,T_final,W,Q,gamma1,gamma2,Nout]=read_ic();

disp('initial OK')

I = 5;
n = I*(Nout-1);
t = n*dt;

%-----------------------------------------------------------------%
%               Third-Order Runge-Kutta scheme                    %
%-----------------------------------------------------------------%

%dt     = T_final/NT
while t<T_final
t      = t+dt;
n      = n+1;
disp([n,t])


% Stage 1
QRK3 = Q;

% Stage 2
[WRK3,Wfx,Wfy]=limiterPrim(dx,dy,QRK3,gamma1,gamma2);
disp('limiter OK')
% [FBC,GBC]=boundaryconditions(NX,NY,W,Wfx,Wfy,gamma1,gamma2);
[FBC,GBC]=linboundaryconditions(NX,NY,W,Wfx,Wfy,gamma1,gamma2);
disp('BoundaryConditions OK')
% CALL osherflux_x(NX,NY,FBC,Wfx,gamma1,gamma2,F,sL,sR)
[F,sL,sR]=linosherflux_x(NX,NY,FBC,Wfx,gamma1,gamma2);
disp('flux_x OK')
% CALL osherflux_y(NX,NY,GBC,Wfy,gamma1,gamma2,G,sA,sB)
[G,sA,sB]=linosherflux_y(NX,NY,GBC,Wfy,gamma1,gamma2);
disp('flux_y OK')
[S]=sourceintegral(NX,NY,WRK3,Wfx,Wfy,gamma1,gamma2);
disp('source OK')

for k=1:5
    for i=1:NY
        for j=1:NX
            dQ1(i,j,k+2)  = -dt/dx*(F(i,j+1,k)-F(i,j,k))-dt/dy*(G(i+1,j,k)-G(i,j,k));
            QRK3(i,j,k+2) = Q(i,j,k+2)+dQ1(i,j,k+2);
        end
    end
end

for i=1:NY
    for j=1:NX
        dQ1(i,j,8)=-dt/dx*(F(i,j+1,6)-F(i,j,6))-dt/dy*(G(i+1,j,6)-G(i,j,6))+dt/dx*(sR(i,j)+sL(i,j+1))+dt/dy*(sB(i+1,j)+sA(i,j))+dt/(dx*dy)*S(i,j);
        QRK3(i,j,8)= Q(i,j,8)+dQ1(i,j,8);
    end
end
disp('time stage 2 OK')

% Stage 3

[WRK3,Wfx,Wfy]=limiterPrim(dx,dy,QRK3,gamma1,gamma2);
disp('limiter OK')
% [FBC,GBC]=boundaryconditions(NX,NY,W,Wfx,Wfy,gamma1,gamma2);
[FBC,GBC]=linboundaryconditions(NX,NY,W,Wfx,Wfy,gamma1,gamma2);
disp('BoundaryConditions OK')
% CALL osherflux_x(NX,NY,FBC,Wfx,gamma1,gamma2,F,sL,sR)
[F,sL,sR]=linosherflux_x(NX,NY,FBC,Wfx,gamma1,gamma2);
disp('flux_x OK')
% CALL osherflux_y(NX,NY,GBC,Wfy,gamma1,gamma2,G,sA,sB)
[G,sA,sB]=linosherflux_y(NX,NY,GBC,Wfy,gamma1,gamma2);
disp('flux_y OK')
[S]=sourceintegral(NX,NY,WRK3,Wfx,Wfy,gamma1,gamma2);
disp('source OK')

for k=1:5
    for i=1:NY
        for j=1:NX
            dQ2(i,j,k+2)  = -dt/dx*(F(i,j+1,k)-F(i,j,k))-dt/dy*(G(i+1,j,k)-G(i,j,k));
            QRK3(i,j,k+2) = Q(i,j,k+2)+1D0/4D0*(dQ1(i,j,k+2)+dQ2(i,j,k+2));
        end
    end
end

for i=1:NY
    for j=1:NX
        dQ2(i,j,8)= -dt/dx*(F(i,j+1,6)-F(i,j,6))-dt/dy*(G(i+1,j,6)-G(i,j,6))+dt/dx*(sR(i,j)+sL(i,j+1))+dt/dy*(sB(i+1,j)+sA(i,j))+dt/(dx*dy)*S(i,j);
        QRK3(i,j,8)= Q(i,j,8)+1D0/4D0*(dQ1(i,j,8)+dQ2(i,j,8));
    end
end
disp('time stage 3 OK')

% Stage 4

[WRK3,Wfx,Wfy]=limiterPrim(dx,dy,QRK3,gamma1,gamma2);
disp('limiter OK')
% [FBC,GBC]=boundaryconditions(NX,NY,W,Wfx,Wfy,gamma1,gamma2);
[FBC,GBC]=linboundaryconditions(NX,NY,W,Wfx,Wfy,gamma1,gamma2);
disp('BoundaryConditions OK')
% CALL osherflux_x(NX,NY,FBC,Wfx,gamma1,gamma2,F,sL,sR)
[F,sL,sR]=linosherflux_x(NX,NY,FBC,Wfx,gamma1,gamma2);
disp('flux_x OK')
% CALL osherflux_y(NX,NY,GBC,Wfy,gamma1,gamma2,G,sA,sB)
[G,sA,sB]=linosherflux_y(NX,NY,GBC,Wfy,gamma1,gamma2);
disp('flux_y OK')
[S]=sourceintegral(NX,NY,WRK3,Wfx,Wfy,gamma1,gamma2);
disp('source OK')


for k=1:5
    for i=1:NY
        for j=1:NX
            dQ3(i,j,k+2)  = -dt/dx*(F(i,j+1,k)-F(i,j,k))-dt/dy*(G(i+1,j,k)-G(i,j,k));
            Q(i,j,k+2) = Q(i,j,k+2)+1D0/6D0*(dQ1(i,j,k+2)+dQ2(i,j,k+2)+4D0*dQ3(i,j,k+2));
        end
    end
end

for i=1:NY
    for j=1:NX
        dQ3(i,j,8)= -dt/dx*(F(i,j+1,6)-F(i,j,6))-dt/dy*(G(i+1,j,6)-G(i,j,6))+dt/dx*(sR(i,j)+sL(i,j+1))+dt/dy*(sB(i+1,j)+sA(i,j))+dt/(dx*dy)*S(i,j);
        Q(i,j,8)   = Q(i,j,8)+1D0/6D0*(dQ1(i,j,8)+dQ2(i,j,8)+4D0*dQ3(i,j,8));
    end
end
disp('time stage 4 OK')

%Q = QRK3

%Convert from conservation (Q) to primary (W) variables
for i=1:NY
    for j=1:NX
        W(i,j,3)=Q(i,j,3);
        W(i,j,4)=Q(i,j,4)/Q(i,j,3);
        W(i,j,5)=Q(i,j,5)/Q(i,j,3);
        W(i,j,6)=(gamma1-gamma2)*(Q(i,j,8)-((Q(i,j,7)*(Q(i,j,4)^(2D0)+Q(i,j,5)^(2D0)))/(2D0*Q(i,j,3)^(2D0))))+(gamma2-1D0)*(Q(i,j,6)-(Q(i,j,4)^(2D0)+Q(i,j,5)^(2D0))/(2D0*Q(i,j,3)));
        W(i,j,7)=Q(i,j,7)/Q(i,j,3);
        W(i,j,8)=((gamma1-1D0)*(Q(i,j,8)-Q(i,j,7)*(Q(i,j,4)^(2D0)+Q(i,j,5)^(2D0))/(2D0*Q(i,j,3)^(2D0))))/((gamma1-gamma2)*(Q(i,j,8)-Q(i,j,7)*(Q(i,j,4)^(2D0)+Q(i,j,5)^(2D0))/(2D0*Q(i,j,3)^(2D0)))+(gamma2-1D0)*(Q(i,j,6)-(Q(i,j,4)^(2D0)+Q(i,j,5)^(2D0))/(2D0*Q(i,j,3))));
    end
end
disp('Cons-to-Prim OK')

if n==I*Nout
    disp(['Saving Matlab-readable output file:',num2str(Nout),'.mat'])
    save(num2str(Nout),'W');
    Nout = Nout+1;
end

end %time step