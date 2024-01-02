function postprocessen(W,t,x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate exact solution                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rho   = W(1,:,end);
u     = W(2,:,end);
p     = W(3,:,end);

t     = t(end);
gamma = 1.4;

%Initial conditions
rho1   = rho(1);
u1     = u(1);
p1     = p(1);

rho4   = rho(end);
u4     = u(end);
p4     = p(end);

%initial speed of sound
a1   = sqrt(gamma*p1/rho1);
a4   = sqrt(gamma*p4/rho4);

eps = 1e-12;

%initial guess final state pressure
pf01 = (((gamma-1)/2*(u1-u4)+a1+a4)/(a1/p1^((gamma-1)/(2*gamma))+a4/p4^((gamma-1)/(2*gamma))))^((2*gamma)/(gamma-1));
pf02 = (rho1*a1*p4+rho4*a4*p1-rho1*a1*rho4*a4*(u4-u1))/(rho1*a1+rho4*a4);

i = 0;
pfnew = 10;
pfold = 0;

while abs(pfnew-pfold)>1e-8
    i = i + 1;

    if i==1
        pf = max(eps,pf02);
    elseif i==2
        pf = max(eps,pf01);
    else
        pf = max(eps,pfnew);
    end

    if pf>=p1%left == 1 %shock
        m1 = rho1*a1*sqrt(1+(gamma+1)/(2*gamma)*(pf-p1)/p1);
    elseif pf<p1 %left == 2 % expansion
        m1 = rho1*a1*(gamma-1)/(2*gamma)*(1-pf/p1)/(1-(pf/p1)^((gamma-1)/(gamma*2)));
    end

    if pf>=p4 %right == 1 %shock
        m4 = rho4*a4*sqrt(1+(gamma+1)/(2*gamma)*(pf-p4)/p4);
    elseif pf<p4 %right == 2 %expansion
        m4 = rho4*a4*(gamma-1)/(2*gamma)*(1-pf/p4)/(1-(pf/p4)^((gamma-1)/(gamma*2)));
    end

    %pressure of final state
    G = (m1*p4+m4*p1-m1*m4*(u4-u1))/(m1+m4) - pf;

    if i>1
        pfnew = pf - G*(pf-pfold+eps)/(G-Gold+eps);
    end

    %save data
    Gold  = G;
    pfold = pf;
end
pf = pfnew;

%velocity of final state
if pf>=p1 %left == 1 %shock
    m1 = rho1*a1*sqrt(1+(gamma+1)/(2*gamma)*(pf-p1)/p1);
elseif pf<p1 %left == 2 %expansion
    m1 = rho1*a1*(gamma-1)/(2*gamma)*(1-pf/p1)/(1-(pf/p1)^((gamma-1)/(gamma*2)));
end
    
if pf>=p4 %right == 1 % shock 
    m4 = rho4*a4*sqrt(1+(gamma+1)/(2*gamma)*(pf-p4)/p4);
elseif pf<p4 %right == 2 %expansion
    m4 = rho4*a4*(gamma-1)/(2*gamma)*(1-pf/p4)/(1-(pf/p4)^((gamma-1)/(gamma*2)));    
end

uf = (m1*u1+m4*u4-(p4-p1))/(m1+m4);


%density
if pf>=p1 %left == 1 %shock
    rho2 = rho1*(1+(gamma+1)/(gamma-1)*pf/p1)/((gamma+1)/(gamma-1)+pf/p1);
elseif pf<p1 %left == 2 %expansion
    rho2 = rho1*(pf/p1)^(1/gamma);
end

if pf>=p4 %right ==1 %shock
    rho3 = rho4*(1+(gamma+1)/(gamma-1)*pf/p4)/((gamma+1)/(gamma-1)+pf/p4);
elseif pf<p4 %right == 2 %expansion
    rho3 = rho4*(pf/p4)^(1/gamma);
end

%speed of sound
a2 = sqrt(gamma*pf/rho2);
a3 = sqrt(gamma*pf/rho3);


%direction characteristics
if pf>=p1
    dxdts1 = u1 - a1*sqrt(1+(gamma+1)/(2*gamma)*(pf-p1)/p1);
elseif pf<p1
    dxdtfirst1 = u1 - a1;
    dxdtlast1  = uf - a2;
end

if pf>=p4
    dxdts4 = u4 + a4*sqrt(1+(gamma+1)/(2*gamma)*(pf-p4)/p4);
elseif pf<p4
    % dxdtfirst4 = u4 + a4;
    % dxdtlast4  = uf + a3;
    dxdtlast4  = u4 + a4;
    dxdtfirst4 = uf + a3;
end

%contact discontinuity
dxdt0 = uf;

N = length(x);
dx = 1/N;

Nexp = 1e3;

x0 = x(1)-dx/2;
x6 = x(end)+dx/2;

if pf>=p1
    x1 = dxdts1 * t+(x6-x0)/2;
    x2 = dxdts1 * t+(x6-x0)/2;
    xleft = x2;
elseif pf<p1
    x1 = dxdtfirst1 * t+(x6-x0)/2;
    x2  = dxdtlast1 * t+(x6-x0)/2;
    xleft(:,1) = x1+(1:Nexp)*(x2-x1)/Nexp;
end

x3 = dxdt0 * t+(x6-x0)/2;

if pf>=p4
    x4 = dxdts4 * t+(x6-x0)/2;
    x5 = dxdts4 * t+(x6-x0)/2;
    xright = x4;
elseif pf<p4
    x4 = dxdtfirst4 * t+(x6-x0)/2;
    x5  = dxdtlast4 * t+(x6-x0)/2;
    xright(:,1) = x4+(1:Nexp)*(x5-x4)/Nexp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact solution per cell                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



rhoexact   = zeros(N,1);
uexact     = zeros(N,1);
pexact     = zeros(N,1);

% left initial state
for j=1:floor(x1/dx+1/2)
    rhoexact(j)   = rho1;
    uexact(j)     = u1;
    pexact(j)     = p1;
end

rhosol   = [rho1;rho1];
usol     = [u1;u1];
psol     = [p1;p1];

% left wave state
if pf<p1
    for j=ceil(x1/dx+1/2):ceil(x2/dx+1/2)
        xexp          = j*dx-dx/2;
        uexact(j)     = 2/(gamma+1)*((xexp-(x6-x0)/2)/t+a1)+(gamma-1)/(gamma+1)*u1;
        aexpleft      = -(gamma-1)/(gamma+1)*((xexp-(x6-x0)/2)/t-u1)+2/(gamma+1)*a1;
        rhoexact(j)   = rho1*(aexpleft/a1)^(2/(gamma-1));
        pexact(j)     = p1*(aexpleft/a1)^((2*gamma)/(gamma-1));
    end
        usol(3:(Nexp+2),1)     = 2/(gamma+1)*((xleft-(x6-x0)/2)/t+a1)+(gamma-1)/(gamma+1)*u1;
        asol(1:Nexp,1)         = -(gamma-1)/(gamma+1)*((xleft-(x6-x0)/2)/t-u1)+2/(gamma+1)*a1;
        psol(3:(Nexp+2),1)     = p1*(asol(1:Nexp)/a1).^((2*gamma)/(gamma-1));
        rhosol(3:(Nexp+2),1)   = rho1*(asol(1:Nexp)/a1).^(2/(gamma-1));
else
    rhosol(3,1)   = rho2;
    usol(3,1)     = uf;
    psol(3,1)     = pf;
end

%  final state 2
for j=ceil(x2/dx+1/2):ceil(x3/dx+1/2)
    rhoexact(j)   = rho2;
    uexact(j)     = uf;
    pexact(j)     = pf;
end

rhosol((length(rhosol)+1):(length(rhosol)+2),1)       = [rho2;rho2];
usol((length(usol)+1):(length(usol)+2),1)             = [uf;uf];
psol((length(psol)+1):(length(psol)+2),1)             = [pf;pf];

% final state 3
for j=ceil(x3/dx+1/2):floor(x4/dx+1/2)
    rhoexact(j)   = rho3;
    uexact(j)     = uf;
    pexact(j)     = pf;
end

rhosol((length(rhosol)+1):(length(rhosol)+2),1)       = [rho3;rho3];
usol((length(usol)+1):(length(usol)+2),1)             = [uf;uf];
psol((length(psol)+1):(length(psol)+2),1)             = [pf;pf];

if pf<p4
    for j=ceil(x4/dx+1/2):floor(x5/dx+1/2)
        xexp          = j*dx-dx/2;
        uexact(j)     = 2/(gamma+1)*((xexp-(x6-x0)/2)/t-a4)+(gamma-1)/(gamma+1)*u4;
        aexpright     = (gamma-1)/(gamma+1)*((xexp-(x6-x0)/2)/t-u4)+2/(gamma+1)*a4;
        rhoexact(j)   = rho4*(aexpright/a4)^(2/(gamma-1));
        pexact(j)     = p4*(aexpright/a4)^((2*gamma)/(gamma-1));
    end

    usol((length(usol)+1):(length(usol)+Nexp),1)       = 2/(gamma+1)*((xright-(x6-x0)/2)/t-a4)+(gamma-1)/(gamma+1)*u4;
    asol(1:Nexp,1)                                     = (gamma-1)/(gamma+1)*((xright-(x6-x0)/2)/t-u4)+2/(gamma+1)*a4;
    psol((length(psol)+1):(length(psol)+Nexp),1)       = p4*(asol(1:Nexp)/a4).^((2*gamma)/(gamma-1));
    rhosol((length(rhosol)+1):(length(rhosol)+Nexp),1) = rho4*(asol(1:Nexp)/a4).^(2/(gamma-1));
else
    rhosol((length(rhosol)+1),1) = rho3;
    usol((length(usol)+1),1)     = uf;
    psol((length(psol)+1),1)     = pf;
end

for j=ceil(x5/dx+1/2):N
    rhoexact(j) = rho4;
    uexact(j)   = u4;
    pexact(j)   = p4;
end

rhosol((length(rhosol)+1):(length(rhosol)+2),1) = [rho4;rho4];
usol((length(usol)+1):(length(usol)+2),1)       = [u4;u4];
psol((length(psol)+1):(length(psol)+2),1)       = [p4;p4];

xplot = [x0;x1;x1;xleft;x3;x3;xright;x5;x5;x6];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Errors                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L1(1) = 1/length(rho)*sum(abs(rhoexact-rho'));
L1(2) = 1/length(u)*sum(abs(uexact-u'));
L1(3) = 1/length(p)*sum(abs(pexact-p'));
L1(4) = 1/length(p)*sum(abs(pexact./(rhoexact*(gamma-1))-(p./(rho*(gamma-1)))'));

L2(1) = sqrt(1/length(rho)*sum((rhoexact-rho').^2));
L2(2) = sqrt(1/length(u)*sum((uexact-u').^2));
L2(3) = sqrt(1/length(p)*sum((pexact-p').^2));
L2(4) = sqrt(1/length(p)*sum((pexact./(rhoexact*(gamma-1))-(p./(rho*(gamma-1)))').^2));

Linf(1) = max(abs(rhoexact-rho'));
Linf(2) = max(abs(uexact-u'));
Linf(3) = max(abs(pexact-p'));
Linf(4) = max(abs(pexact./(rhoexact*(gamma-1))-(p./(rho*(gamma-1)))'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figures                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
% subplot(2,2,1)
hold on
plot(xplot,rhosol,'r','linewidth',1)
plot(x,rho,'o','MarkerSize',7,'linewidth',1)
set(gca,'fontsize',14)
grid on
xlim([0 1])
% ylim([0 1.2])
Y = ylim;
title('density')
xlabel('x')
ylabel('\rho')
legend('exact','approx')
text(0.5,Y(2)*19/20,strcat('L_1 = ',num2str(L1(1),'%e')),'fontsize',14)
text(0.5,Y(2)*18/20,strcat('L_2 = ',num2str(L2(1),'%e')),'fontsize',14)
text(0.5,Y(2)*17/20,strcat('L_\infty = ',num2str(Linf(1),'%e')),'fontsize',14)


figure
% subplot(2,2,2)
hold on
plot(xplot,usol,'r','linewidth',1)
plot(x,u,'o','MarkerSize',7,'linewidth',1)
set(gca,'fontsize',14)
grid on 
xlim([0 1])
% ylim([-0.2 2])
Y = ylim;
title('velocity')
xlabel('x')
ylabel('u')
legend('exact','approx')
text(0.5,Y(2)*19/20,strcat('L_1 = ',num2str(L1(2),'%e')),'fontsize',14)
text(0.5,Y(2)*18/20,strcat('L_2 = ',num2str(L2(2),'%e')),'fontsize',14)
text(0.5,Y(2)*17/20,strcat('L_\infty = ',num2str(Linf(2),'%e')),'fontsize',14)

figure
% subplot(2,2,3)
hold on
plot(xplot,psol,'r','linewidth',1)
plot(x,p,'o','MarkerSize',7,'linewidth',1)
set(gca,'fontsize',14)
grid on
xlim([0 1])
% ylim([0 1.2])
Y = ylim;
title('pressure')
xlabel('x')
ylabel('p')
legend('exact','approx')
text(0.5,Y(2)*19/20,strcat('L_1 = ',num2str(L1(3),'%e')),'fontsize',14)
text(0.5,Y(2)*18/20,strcat('L_2 = ',num2str(L2(3),'%e')),'fontsize',14)
text(0.5,Y(2)*17/20,strcat('L_\infty = ',num2str(Linf(3),'%e')),'fontsize',14)

figure
% subplot(2,2,4)
hold on
plot(xplot,psol./(rhosol*(gamma-1)),'r','linewidth',1)
plot(x,p./(rho*(gamma-1)),'o','MarkerSize',7,'linewidth',1)
set(gca,'fontsize',14)
grid on
xlim([0 1])
% ylim([0 1.2])
Y = ylim;
title('internal energy')
xlabel('x')
ylabel('e')
legend('exact','approx')
text(0.5,Y(2)*19/20,strcat('L_1 = ',num2str(L1(4),'%e')),'fontsize',14)
text(0.5,Y(2)*18/20,strcat('L_2 = ',num2str(L2(4),'%e')),'fontsize',14)
text(0.5,Y(2)*17/20,strcat('L_\infty = ',num2str(Linf(4),'%e')),'fontsize',14)