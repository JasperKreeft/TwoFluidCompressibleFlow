clear % all
close all
clc

addpath(genpath("solver"))

BC = [2 1 2 5];

limiter = 2;

flux = 2;
time = 2;

CFL = 0.45;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     grid

[NX,NY,dx,dy,x,y]=grids;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial conditions
n = 1;

[W,gamma] = initialconditions(NX,NY);

Q = Prim2Cons(W,gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t(n)   = 0;
Lambda = 20;

dQ1 = zeros(NY,NX,4);
dQ2 = zeros(NY,NX,4);
dQ3 = zeros(NY,NX,4);

error = 1;

while error>1e-10 && n<1000

n      = n+1;
dt     = CFL*dx/Lambda;
t(n)   = t(n-1)+dt;
Lambda = [];

q = Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Third-Order Runge-Kutta scheme                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Stage 1

[w]        = Cons2Prim(q,gamma);
[Wfx,Wfy]  = limiters(NX,NY,w,limiter,eps);
if flux==1
    [F,Lambda] = ExactFlux_x(NX,NY,Wfx,gamma,Lambda);
    [G,Lambda] = ExactFlux_y(NX,NY,Wfy,gamma,Lambda);
elseif flux==2
    [F,Lambda] = OsherFlux_x(NX,NY,Wfx,gamma,Lambda);
    [G,Lambda] = OsherFlux_y(NX,NY,Wfy,gamma,Lambda);
end
[F,G]      = BoundaryConditions(F,G,Wfx,Wfy,gamma,NX,NY,BC);

for k=1:4
    dQ1(:,:,k) = -dt/dx*(F(:,2:(NX+1),k)-F(:,1:NX,k))-dt/dy*(G(2:(NY+1),:,k)-G(1:NY,:,k));
    q(:,:,k)   = Q(:,:,k)+dQ1(:,:,k);
end

if time==2
% Stage 2

[w]        = Cons2Prim(q,gamma);
[Wfx,Wfy]  = limiters(NX,NY,w,limiter,eps);
if flux==1
        [F,Lambda] = ExactFlux_x(NX,NY,Wfx,gamma,Lambda);
        [G,Lambda] = ExactFlux_y(NX,NY,Wfy,gamma,Lambda);
elseif flux==2
        [F,Lambda] = OsherFlux_x(NX,NY,Wfx,gamma,Lambda);
        [G,Lambda] = OsherFlux_y(NX,NY,Wfy,gamma,Lambda);
end
[F,G]      = BoundaryConditions(F,G,Wfx,Wfy,gamma,NX,NY,BC);

for k=1:4
    dQ2(:,:,k) = -dt/dx*(F(:,2:(NX+1),k)-F(:,1:NX,k))-dt/dy*(G(2:(NY+1),:,k)-G(1:NY,:,k));
    q(:,:,k)   = Q(:,:,k)+1/4*dQ1(:,:,k)+1/4*dQ2(:,:,k);
end

% Stage 3

[w]        = Cons2Prim(q,gamma);
[Wfx,Wfy]  = limiters(NX,NY,w,limiter,eps);
if flux==1
    [F,Lambda] = ExactFlux_x(NX,NY,Wfx,gamma,Lambda);
    [G,Lambda] = ExactFlux_y(NX,NY,Wfy,gamma,Lambda);
elseif flux==2
    [F,Lambda] = OsherFlux_x(NX,NY,Wfx,gamma,Lambda);
    [G,Lambda] = OsherFlux_y(NX,NY,Wfy,gamma,Lambda);
end
[F,G]      = BoundaryConditions(F,G,Wfx,Wfy,gamma,NX,NY,BC);

for k=1:4
    dQ3(:,:,k) = -dt/dx*(F(:,2:(NX+1),k)-F(:,1:NX,k))-dt/dy*(G(2:(NY+1),:,k)-G(1:NY,:,k));
    q(:,:,k)   = Q(:,:,k)+1/6*dQ1(:,:,k)+1/6*dQ2(:,:,k)+2/3*dQ3(:,:,k);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qold = Q;
Q = q;

error = max(max(abs((Qold(:,:,1)-Q(:,:,1))./Qold(:,:,1))));

% figure(1)
% semilogy(n,error,'.r')
% hold on
% % ylim([1e-17,1])
% pause(0.01)

W = Cons2Prim(Q,gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot(W(5,:,1),'-')
% hold on
% plot(W(1,:,1),'-g')
% plot(W(10,:,1),'-r')
% hold off
% imagesc(W(:,:,4))
% surf(W(:,:,1))
% pcolor(W(:,:,3))
% shading interp
% colorbar

if rem(n,50)==0
M=sqrt(W(:,:,2).^2+W(:,:,3).^2)./sqrt(gamma*W(:,:,4)./W(:,:,1));
pcolor(x,y,M)
shading interp
axis equal
axis([0 4.1 0 1])
% colorbar('NorthOutside')
colorbar('horiz')
title('Mach nr')
colormap jet
set(gca,'clim',[1.9 2.9])
pause(0.01)
end

end %time step

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
M=sqrt(W(:,:,2).^2+W(:,:,3).^2)./sqrt(gamma*W(:,:,4)./W(:,:,1));
contourf(x,y,M,20)
axis equal
axis([0 4.1 0 1])
colorbar('NorthOutside')
% colorbar('horiz')
title('Mach nr')