clear all
close all
clc
% warning('off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
N = 200;

% T = 3e-5*sqrt(1e5); % full Epoxy-Spinel
T = 0.2; % Sod
% T = 2e-4*sqrt(1e5);

time = 2;
BC_L = 1;   % alleen open bc's mogelijk voor mixture two fluid met Stiffened gas
BC_R = 1;
limiter = 2;
ani = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L  = 1;
dx = L/N;

CFL = 0.45;

x = ((1:N)-1/2)*dx;

y = [0 0.2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial conditions
n = 1;

% gamma1 = 1.4;
% gamma2 = 4.4;
% pi1    = 0.0;
% pi2 = 6e3;
% 
% W(1,1:(N/2)) = 525;
% W(2,1:(N/2)) = 0.0;
% W(3,1:(N/2)) = 1e4;
% W(4,1:(N/2)) = 25/525;
% W(5,1:(N/2)) = 0.5;
% 
% 
% W(1,(N/2+1):N) = 525;
% W(2,(N/2+1):N) = 0.0;
% W(3,(N/2+1):N) = 1;
% W(4,(N/2+1):N) = 25/525;
% W(5,(N/2+1):N) = 0.5;

% % Epoxy-Spinel
% gamma1 = 2.43;
% gamma2 = 1.62;
% pi1    = 5.3e4;
% pi2    = 141.45e4;
% 
% W(1,1:(N*6/10)) = 2171;
% W(2,1:(N*6/10)) = 0.0;
% W(3,1:(N*6/10)) = 2e6;
% W(4,1:(N*6/10)) = 0.3250;
% W(5,1:(N*6/10)) = 0.5954;
% 
% W(1,(N*6/10+1):N) = 2171;
% W(2,(N*6/10+1):N) = 0.0;
% W(3,(N*6/10+1):N) = 1;
% W(4,(N*6/10+1):N) = 0.3250;
% W(5,(N*6/10+1):N) = 0.5954;


% Mixture Sod problem
gamma1 = 1.6;
gamma2 = 1.4;
pi1    = 0;
pi2    = 0;

W(1,1:(N*5/10)) = 1.0;
W(2,1:(N*5/10)) = 0.0;
W(3,1:(N*5/10)) = 1.0;
W(4,1:(N*5/10)) = 0.5;
W(5,1:(N*5/10)) = 0.5;

W(1,(N*5/10+1):N) = 0.125;
W(2,(N*5/10+1):N) = 0.0;
W(3,(N*5/10+1):N) = 0.1;
W(4,(N*5/10+1):N) = 0.5;
W(5,(N*5/10+1):N) = 0.5;

% % Strong Interface
% gamma1 = 1.6;
% gamma2 = 1.4;
% pi1    = 0;
% pi2    = 0;
% 
% W(1,1:(N*5/10)) = 1000.0;
% W(2,1:(N*5/10)) = 1.0;
% W(3,1:(N*5/10)) = 1.0;
% W(4,1:(N*5/10)) = 1.0;
% W(5,1:(N*5/10)) = 1.0;
% 
% W(1,(N*5/10+1):N) = 1.0;
% W(2,(N*5/10+1):N) = 1.0;
% W(3,(N*5/10+1):N) = 1.0;
% W(4,(N*5/10+1):N) = 0.0;
% W(5,(N*5/10+1):N) = 0.0;

% % Water-air problem (Murrone 5.2.1)
% gamma1 = 1.4;
% gamma2 = 4.4;
% pi1    = 0;
% pi2    = 6e3;
% 
% W(1,1:(N*5/10)) = 525;
% W(2,1:(N*5/10)) = 0.0;
% W(3,1:(N*5/10)) = 1e4;
% W(4,1:(N*5/10)) = 0.5*50/525;
% W(5,1:(N*5/10)) = 0.5;
% 
% W(1,(N*5/10+1):N) = 525;
% W(2,(N*5/10+1):N) = 0.0;
% W(3,(N*5/10+1):N) = 1;
% W(4,(N*5/10+1):N) = 0.5*1000/525;
% W(5,(N*5/10+1):N) = 0.5;

% % no-reflection problem
% gamma1 = 1.667;
% gamma2 = 1.2;
% pi1    = 0;
% pi2    = 0;
% 
% W(1,1:(N*5/10)) = 3.1748;
% W(2,1:(N*5/10)) = 9.4350;
% W(3,1:(N*5/10)) = 100;
% W(4,1:(N*5/10)) = 1.0;
% W(5,1:(N*5/10)) = 1.0;
% 
% W(1,(N*5/10+1):N) = 1.0;
% W(2,(N*5/10+1):N) = 0.0;
% W(3,(N*5/10+1):N) = 1.0;
% W(4,(N*5/10+1):N) = 0.0;
% W(5,(N*5/10+1):N) = 0.0;

% % Strong Interface
% gamma1 = 1.4;
% gamma2 = 1.4;
% pi1    = 0;
% pi2    = 6e-3;
% 
% 
% rho1 = 1000; rho2 = 50; alphal = 1-1e-1; alphar = 1e-1;
% W(1,1:(N*5/10)) = alphal*rho1+(1-alphal)*rho2;
% W(2,1:(N*5/10)) = 1000.0;
% W(3,1:(N*5/10)) = 1.0;
% W(4,1:(N*5/10)) = alphal*rho1/(alphal*rho1+(1-alphal)*rho2);
% W(5,1:(N*5/10)) = alphal;
% 
% W(1,(N*5/10+1):N) = alphar*rho1+(1-alphar)*rho2;
% W(2,(N*5/10+1):N) = 1000.0;
% W(3,(N*5/10+1):N) = 1.0;
% W(4,(N*5/10+1):N) = alphar*rho1/(alphar*rho1+(1-alphar)*rho2);
% W(5,(N*5/10+1):N) = alphar;

Q = Prim2Cons(W,gamma1,gamma2,pi1,pi2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t(n)   = 0; %global t N
Lambda = 20;

Qf = zeros(5,2*N);

if ani==1
    han = waitbar(0,'Please wait...','Position',[500 500 288 60]);
elseif ani==3 || ani==4
%     figure('Position',[1,1024,1280,1024])
end

while t(n)<T

n      = n+1;
dt     = CFL*dx/Lambda;
t(n)   = t(n-1)+dt;
if ani==1
    waitbar(t(n)/T,han)
end
Lambda = [];

% disp([n t(n)]);

q = Q;

l = limiter;
if n<5
    limiter = min(1,l);
else
    limiter = l;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Third-Order Runge-Kutta scheme                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stage 1

[w]   = Cons2Prim(q,gamma1,gamma2,pi1,pi2);
[Wf]  = limiters(N,w,limiter); 
% [F,sL,sR,Lambda] = OsherFlux(N,Wf,gamma1,gamma2,pi1,pi2,Lambda); 
% [F,sL,sR,Lambda] = LinearOsherFlux(N,Wf,gamma1,gamma2,pi1,pi2,Lambda); 
[F,sL,sR,Lambda] = exactriemann(N,Wf,gamma1,gamma2,pi1,pi2,Lambda); 
[S]   = sourceintegral(N,w,Wf,gamma1,gamma2,pi1,pi2);
% [S]   = sourceintegral2(N,w,Wf,gamma1,gamma2,pi1,pi2);  % Averaging at cell-face
[F]   = BoundaryConditions(F,Wf,gamma1,gamma2,pi1,pi2,N,BC_L,BC_R); 

dQ1(1:4,:) = -dt/dx*(F(1:4,2:N+1)-F(1:4,1:N));
dQ1(5,:)   = -dt/dx*(F(5,2:N+1)-F(5,1:N))+dt/dx*(sR+S+sL);
q          = Q+dQ1;

if time==2
% Stage 2

[w]   = Cons2Prim(q,gamma1,gamma2,pi1,pi2); 
[Wf]  = limiters(N,w,limiter); 
% [F,sL,sR,Lambda] = OsherFlux(N,Wf,gamma1,gamma2,pi1,pi2,Lambda); 
% [F,sL,sR,Lambda] = LinearOsherFlux(N,Wf,gamma1,gamma2,pi1,pi2,Lambda); 
[F,sL,sR,Lambda] = exactriemann(N,Wf,gamma1,gamma2,pi1,pi2,Lambda); 
[S]   = sourceintegral(N,w,Wf,gamma1,gamma2,pi1,pi2);
% [S]   = sourceintegral2(N,w,Wf,gamma1,gamma2,pi1,pi2);  % Averaging at cell-face
[F]   = BoundaryConditions(F,Wf,gamma1,gamma2,pi1,pi2,N,BC_L,BC_R);

dQ2(1:4,:) = -dt/dx*(F(1:4,2:N+1)-F(1:4,1:N));
dQ2(5,:)   = -dt/dx*(F(5,2:N+1)-F(5,1:N))+dt/dx*(sR+S+sL);
q          = Q+1/4*(dQ1+dQ2);


% Stage 3

[w]   = Cons2Prim(q,gamma1,gamma2,pi1,pi2); 
[Wf]  = limiters(N,w,limiter); 
% [F,sL,sR,Lambda] = OsherFlux(N,Wf,gamma1,gamma2,pi1,pi2,Lambda);
% [F,sL,sR,Lambda] = LinearOsherFlux(N,Wf,gamma1,gamma2,pi1,pi2,Lambda);
[F,sL,sR,Lambda] = exactriemann(N,Wf,gamma1,gamma2,pi1,pi2,Lambda);
[S]   = sourceintegral(N,w,Wf,gamma1,gamma2,pi1,pi2);
% [S2]   = sourceintegral2(N,w,Wf,gamma1,gamma2,pi1,pi2);  % Averaging at cell-face
[F]   = BoundaryConditions(F,Wf,gamma1,gamma2,pi1,pi2,N,BC_L,BC_R);

dQ3(1:4,:) = -dt/dx*(F(1:4,2:N+1)-F(1:4,1:N));
dQ3(5,:)   = -dt/dx*(F(5,2:N+1)-F(5,1:N))+dt/dx*(sR+S+sL);
q          = Q+1/6*(dQ1+dQ2+4*dQ3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% speeds(:,:,n) = speedy;
Q = q;

W = Cons2Prim(Q,gamma1,gamma2,pi1,pi2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ani==2
    figure(1)
    plot(x,W(1,:),'+-b','linewidth',2)
    grid on
    title('density')
    xlim([0 1])
    pause(0.05)
elseif ani==3
    figure(1)
    subplot(2,2,1)
    plot(x,W(1,:),'.-','linewidth',1)
    title('mixture density')
    grid on
    xlim([0 1])
%     axis([0 1 0 1200])
%     set(gca,'ytick',0:120:1200)
    subplot(2,2,2)
    plot(x,W(2,:),'.-r','linewidth',1)%*sqrt(10^5)
    grid on
    title('velocity')
%     axis([0 1 -200 800])
%     set(gca,'ytick',-200:200:800)
%     axis([0 1 -500 4500])
%     set(gca,'ytick',-500:500:4500)
    subplot(2,2,3)
    plot(x,W(3,:),'.-g','linewidth',1)%*10^5
    grid on
    title('pressure')
%     axis([0 1 0 10^9])
%     set(gca,'ytick',(0:10)*10^8)
    xlim([0 1])
    subplot(2,2,4)
    plot(x,W(5,:),'.-m','linewidth',1)
    title('epoxy volume fraction')
    grid on
%     axis([0 1 0.46 0.62])
%     set(gca,'ytick',0.46:0.02:0.62)
    pause(0.01)
elseif ani==4
    z = zeros(2,N);
    z(1,:) = W(1,:);
    z(2,:) = W(1,:);
    figure(1)
    pcolor(x,y,z)
    shading interp
    axis equal
    axis off
    colorbar('horiz')
    colormap(hot)
    pause(0.1)
end

end %time step
plotjes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%