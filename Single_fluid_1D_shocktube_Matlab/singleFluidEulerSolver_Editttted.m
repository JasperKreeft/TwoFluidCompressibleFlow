clear all
close all
clc

WindowSize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
options;
% nr = 4; N = 200; T = 5.2; time = 2; flux = 3; BC_L = 1; BC_R = 1;
% limitedvariables = 1; limiter = 2; ani = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generating grid

L   = 1.0; % Size domain

eta_nodes = linspace(0,L,N+1);
x_nodes = eta_nodes+0.00*sin(2*pi*eta_nodes);
dx = diff(x_nodes);
x_cellcenter = (x_nodes(1:end-1)+x_nodes(2:end))/2;

y = [0 0.2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions

initialCondition

Q = Prim2Cons(W,gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if time==1
    if nr==4
        CFL = 0.49;
    else
        CFL = 0.9;
    end
elseif time==2
    if nr==4
        CFL = 0.49;
    else
        CFL = .9;
    end
end

t(n)   = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ani==1
    han = waitbar(0,'Please wait...');
elseif ani==3 || ani==4
%     figure('Position',[1,1,1280,1024])
    figure(winFull{:})
end


Lambda = 20;
Qf=zeros(3,2*N);

while t(n)<T

n      = n+1;
dt     = CFL*min(dx)/Lambda; % !!!!!!!!!!!!! CFL*min(dx./Lambda);
t(n)   = t(n-1)+dt;
if ani==1
    waitbar(t(n)/T,han)
end
Lambda = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Third-Order Runge-Kutta scheme                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dX = ones(3,1)*dx;

% Stage 0
q = Q; %#ok<NASGU>

% Stage 1
FluxFunction

dQ1 = -dt./dX.*(F(:,2:N+1)-F(:,1:N));
q   = Q+dQ1;

if time==2

% Stage 2
FluxFunction

dQ2 = -dt./dX.*(F(:,2:N+1)-F(:,1:N));
q   = Q+1/4*(dQ1+dQ2); %#ok<NASGU>

% Stage 3
FluxFunction

dQ3 = -dt./dX.*(F(:,2:N+1)-F(:,1:N));
q   = Q+1/6*(dQ1+dQ2+4*dQ3);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Q = q;

W = Cons2Prim(Q,gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end %time step
if ani==2
    figure(1)
    plot(x_cellcenter,W(1,:),'+-b','linewidth',2)
    grid on
    title('density')
    pause(0.1)
elseif ani==3
    figure(1)
    subplot(2,2,1)
    plot(x_cellcenter,W(1,:),'.-','linewidth',1)
    title('density')
    grid on
    xlim([0 1])
    if     nr==1; ylim([0 1.1]);
    elseif nr==2; ylim([0 1100]);
    elseif nr==3; ylim([0 1.1]);
    elseif nr==4; ylim([0 11]);
    elseif nr==5; ylim([0 1.1]);
    end
    subplot(2,2,2)
    plot(x_cellcenter,W(2,:),'.-','linewidth',1)
    grid on
    title('velocity')
    xlim([0 1])
    if     nr==1; ylim([0.9 1.1]);
    elseif nr==2; ylim([0.9 1.1]);
    elseif nr==3; ylim([0 1.1]);
    elseif nr==4; ylim([0 2]);
    elseif nr==5; ylim([-2.5 2.5]);
    end
    subplot(2,2,3)
    plot(x_cellcenter,W(3,:),'.-','linewidth',1)
    grid on
    title('pressure')
    xlim([0 1])
    if     nr==1; ylim([0.9 1.1]);
    elseif nr==2; ylim([0.9 1.1]);
    elseif nr==3; ylim([0 1.1]);
    elseif nr==4; ylim([0 11]);
    elseif nr==5; ylim([0 0.5]);
    end
    subplot(2,2,4)
    e = W(3,:)./((gamma-1)*W(1,:));
    plot(x_cellcenter,e,'.-','linewidth',1)
    title('internal energy')
    grid on
    xlim([0 1])
    if     nr==1; ylim([0 30]);
    elseif nr==2; ylim([0 3]);
    elseif nr==3; ylim([1.5 3]);
    elseif nr==4; ylim([0 5]);
    elseif nr==5; ylim([0 1.1]);
    end
    pause(0.05)
elseif ani==4
    xp = [ x_cellcenter ; x_cellcenter ];
    yp = y'*ones(1,N);
    zp = zeros(2,N);
    zp(1,:) = W(1,:);
    zp(2,:) = W(1,:);
    figure(1)
    pcolor(xp,yp,zp)
    shading interp
    axis equal
    axis off
    cbar = colorbar('horiz');
    cbar.Label.String = 'Density';
    set(gca,'clim',[0 1]);
    colormap jet % cool
    title(['Time = ' num2str(t(n),'%4.3f') ' sec'])
    pause(0.01)
end

end %time step
toc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Jplus = zeros(length(t),N);
% Jmin  = zeros(length(t),N);
% 
% for k=1:length(t)
%     Jplus(k,:) = W(2,:,k)+2/(gamma-1)*sqrt(gamma*W(3,:,k)./W(1,:,k));
%     Jmin(k,:)  = W(2,:,k)-2/(gamma-1)*sqrt(gamma*W(3,:,k)./W(1,:,k));
% end
% 
% figure
% contour(x_cellcenter,t,Jplus,30)
% hold on
% contour(x_cellcenter,t,Jmin,30)
% colormap([0 0 0])
% xlabel('location x')
% ylabel('time t')
% 
% 
% z = zeros(length(t),N);
% 
% for j=1:length(t)
%     for i=1:N
%         z(j,i) = W(1,i,j);
%     end
% end
% 
% figure
% contourf(x_cellcenter,t,z,50)
% xlabel('location x')
% ylabel('time t')
% title('x-t diagram with contourlines of the density \rho')
% 
% figure
% pcolor(x_cellcenter,t,z)
% shading interp
% % set(gca,'clim',[0 1])
% xlabel('location x')
% ylabel('time t')
% title('x-t diagram with contourlines of the density \rho')
% 
% if abs(W(1,1,1)-W(1,1,end))<1e-8 && abs(W(1,end,1)-W(1,end,end))<1e-8
%     postprocessen(W,t,x_cellcenter)
% end