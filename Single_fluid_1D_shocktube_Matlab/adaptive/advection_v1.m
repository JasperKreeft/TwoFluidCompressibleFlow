% ### STABILITY #####

% m-file with the solution for Problem 1 of the practical assignment
% in Week 5 of the course AE3-030.

clear all
clf%ose all
clc
% figure('Position',[100,100,800,600])
% set(gca,'FontSize',14)
% xlabel('x','FontSize',16)
% ylabel('u','FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings
a     = 0.25;    % velocity
N     = 30;      % Number of grid cells
alpha = 1;    % CFL number
Tend  = 2;       % end time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_nodes = linspace(0,1,N+1);
dx(1) = 1/N;
dx(2) = dx(1)/2;
dx(3) = dx(2)/2;
x = (x_nodes(1:end-1)+x_nodes(2:end))/2;% x_cellcenter

cellnr = 1:N;
position = x;
level = ones(1,N);
leftnb = [N 1:N-1];
rightnb = [2:N 1];

nr = N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initial solution (u0)
for i=1:N
    if x(i)<=0.5
        u0(i,1) = 1/4-1/4*cos(4*pi*x(i));
    else
        u0(i,1) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condition
u = u0;

tlevel = 1;
t = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Marching in time
while t<Tend
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dt = alpha*dx(max(level))/abs(a); % time step
    t = t+dt;

% --  Advance in time with Forward Euler --
    % inner region
    for i=1:length(cellnr)
       u_new(i)= u(i) - a*dt/dx(level(i))*( u(i) - u(cellnr==leftnb(i)) );
    end
    u = u_new;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:length(cellnr)
        Cdp(i) = max(abs(u(i)-u(cellnr==leftnb(i)))/(.5+2^(level(i)-level(cellnr==leftnb(i))-1)),...
                     abs(u(cellnr==rightnb(i))-u(i))/(.5+2^(level(i)-level(cellnr==leftnb(i))-1)));
    end

    cellstosplit = sort((level<3).*(Cdp>0.02).*cellnr);
    
    for k=cellstosplit(end:-1:1)
        if k~=0
            nrofcells = length(cellnr);
            rowtosplit = (1:nrofcells)*(cellnr==k)';

            cellnr   = [cellnr nr+1 nr+2];
            position = [position position(rowtosplit)-dx(level(rowtosplit)+1)/2 position(rowtosplit)+dx(level(rowtosplit)+1)/2];
            level    = [level level(rowtosplit)+1 level(rowtosplit)+1];
            leftnb   = [leftnb leftnb(rowtosplit) nr+1];
            rightnb  = [rightnb nr+2 rightnb(rowtosplit)];

            rightnb(cellnr==leftnb(rowtosplit)) = nr+1;
            leftnb(cellnr==rightnb(rowtosplit)) = nr+2;

            u(end+1) = ((1/2*dx(level(cellnr==leftnb(rowtosplit)))+1/4*dx(level(rowtosplit)))*u(rowtosplit)+1/4*dx(level(rowtosplit))*u(cellnr==leftnb(rowtosplit)))/(1/2*dx(level(cellnr==leftnb(rowtosplit)))+1/2*dx(level(rowtosplit)));
            u(end+1) = ((1/2*dx(level(cellnr==rightnb(rowtosplit)))+1/4*dx(level(rowtosplit)))*u(rowtosplit)+1/4*dx(level(rowtosplit))*u(cellnr==rightnb(rowtosplit)))/(1/2*dx(level(cellnr==rightnb(rowtosplit)))+1/2*dx(level(rowtosplit)));

            cellnr(rowtosplit) = []; position(rowtosplit) = []; level(rowtosplit) = []; leftnb(rowtosplit) = []; rightnb(rowtosplit) = [];
            u(rowtosplit) = [];

            nr = nr+2;
        end
    end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot solutions in figure
%     legend_strings = [];
plot(x,u0,'k--'); %legend_strings = [legend_strings;'u0   '];
    hold on
    grid on
%     plot(x,u_exact,'k.--','Linewidth',2); %legend_strings = [legend_strings;'exact'];
    plot(position,u,'o','Linewidth',2);     %legend_strings = [legend_strings;'FE   '];
    plot(position,zeros(1,length(position)),'r+','linewidth',2)
%     legend(legend_strings)
    ylim([-0.2 0.8])
    xlabel('x','FontSize',16)
    ylabel('u','FontSize',16)
    pause%(0.05)
%     F(tlevel) = getframe;
    hold off
%     tlevel= tlevel+1;

end % end timestep