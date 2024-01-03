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
[dx,dy,dt,T_final,W,Q,Nout]=read_ic();



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
[WRK3,Wfx,Wfy]=limiterPrim(dx,dy,QRK3);         disp('limiter OK')
[FBC,GBC]=boundaryconditions(WRK3,Wfx,Wfy);     disp('BoundaryConditions OK')
% [FBC,GBC]=linboundaryconditions(WRK3,Wfx,Wfy);  disp('BoundaryConditions OK')
% [F,sL,sR]=osherflux_x(FBC,Wfx);
[F,sL,sR]=linosherflux_x(FBC,Wfx);              disp('flux_x OK')
% [G,sA,sB]=osherflux_y(GBC,Wfy);
[G,sA,sB]=linosherflux_y(GBC,Wfy);              disp('flux_y OK')
S=sourceintegral(WRK3,Wfx,Wfy);                 disp('source OK')

for k=1:5
    dQ1(:,:,k)  = -dt/dx*(F(:,2:NX+1,k)-F(:,1:NX,k))-dt/dy*(G(2:NY+1,:,k)-G(1:NY,:,k));
    QRK3(:,:,k) = Q(:,:,k+2)+dQ1(:,:,k);
end
dQ1(:,:,6)  = -dt/dx*(F(:,2:NX+1,6)-F(:,1:NX,6))-dt/dy*(G(2:NY+1,:,6)-G(1:NY,:,6))+dt/dx*(sR(:,1:NX)+sL(:,2:NX+1))+dt/dy*(sB(2:NY+1,:)+sA(1:NY,:))+dt/(dx*dy)*S;
QRK3(:,:,6) = Q(:,:,8)+dQ1(:,:,6);
% !!!!!!!!!!!!
for k=1:6; Q(:,:,k+2) = QRK3(:,:,k); end

% % Stage 3
% 
% [WRK3,Wfx,Wfy]=limiterPrim(dx,dy,QRK3);         disp('limiter OK')
% % [FBC,GBC]=boundaryconditions(W,Wfx,Wfy);
% [FBC,GBC]=linboundaryconditions(W,Wfx,Wfy);     disp('BoundaryConditions OK')
% % [F,sL,sR]=osherflux_x(FBC,Wfx);
% [F,sL,sR]=linosherflux_x(FBC,Wfx);              disp('flux_x OK')
% % [G,sA,sB]=osherflux_y(GBC,Wfy);
% [G,sA,sB]=linosherflux_y(GBC,Wfy);              disp('flux_y OK')
% S=sourceintegral(WRK3,Wfx,Wfy);                 disp('source OK')
% 
% for k=1:5
%     dQ2(:,:,k)  = -dt/dx*(F(:,2:NX+1,k)-F(:,1:NX,k))-dt/dy*(G(2:NY+1,:,k)-G(1:NY,:,k));
%     QRK3(:,:,k) = Q(:,:,k+2)+1/4*(dQ1(:,:,k)+dQ2(:,:,k));
% end
% dQ2(:,:,6)  = -dt/dx*(F(:,2:NX+1,6)-F(:,1:NX,6))-dt/dy*(G(2:NY+1,:,6)-G(1:NY,:,6))+dt/dx*(sR(:,1:NX)+sL(:,2:NX+1))+dt/dy*(sB(2:NY+1,:)+sA(1:NY,:))+dt/(dx*dy)*S;
% QRK3(:,:,6) = Q(:,:,8)+1/4*(dQ1(:,:,6)+dQ2(:,:,6));
% 
% 
% % Stage 4
% 
% [WRK3,Wfx,Wfy]=limiterPrim(dx,dy,QRK3);         disp('limiter OK')
% % [FBC,GBC]=boundaryconditions(W,Wfx,Wfy);
% [FBC,GBC]=linboundaryconditions(W,Wfx,Wfy);     disp('BoundaryConditions OK')
% % [F,sL,sR]=osherflux_x(FBC,Wfx);
% [F,sL,sR]=linosherflux_x(FBC,Wfx);              disp('flux_x OK')
% % [G,sA,sB]=osherflux_y(GBC,Wfy);
% [G,sA,sB]=linosherflux_y(GBC,Wfy);              disp('flux_y OK')
% S=sourceintegral(WRK3,Wfx,Wfy);                 disp('source OK')
% 
% 
% for k=1:5
%     dQ3(:,:,k) = -dt/dx*(F(:,2:NX+1,k)-F(:,1:NX,k))-dt/dy*(G(2:NY+1,:,k)-G(1:NY,:,k));
%     Q(:,:,k+2) = Q(:,:,k+2)+1/6*(dQ1(:,:,k)+dQ2(:,:,k)+4*dQ3(:,:,k));
% end
% dQ3(:,:,6) = -dt/dx*(F(:,2:NX+1,6)-F(:,1:NX,6))-dt/dy*(G(2:NY+1,:,6)-G(1:NY,:,6))+dt/dx*(sR(:,1:NX)+sL(:,2:NX+1))+dt/dy*(sB(2:NY+1,:)+sA(1:NY,:))+dt/(dx*dy)*S;
% Q(:,:,8)   = Q(:,:,8)+1/6*(dQ1(:,:,6)+dQ2(:,:,6)+dQ3(:,:,k));


%Convert from conservation (Q) to primary (W) variables
W = Cons2Prim(Q);

if n==I*Nout
    disp(['Saving Matlab-readable output file:',num2str(Nout),'.mat'])
    save(num2str(Nout),'W');
    Nout = Nout+1;
end

figure(1)
plot(W(50,:,1),W(50,:,3),'.-')

end %time step