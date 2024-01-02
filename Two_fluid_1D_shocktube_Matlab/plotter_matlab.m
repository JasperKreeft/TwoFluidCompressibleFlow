% clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Fortran file                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output = load('output.0443');

rho   = W(1,:);
u     = W(2,:);
p     = W(3,:);
beta  = W(4,:);
alpha = W(5,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate exact solution                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


t      = t(end);
gamma1 = 1.667;
gamma4 = 1.2;

%Initial conditions


rho1 = W(1,1);
u1   = W(2,1);
p1   = W(3,1);

rho4 = W(1,end);
u4   = W(2,end);
p4   = W(3,end);

%initial speed of sound
a1   = sqrt(gamma1*p1/rho1);
a4   = sqrt(gamma4*p4/rho4);

ptest = 1;
ptestold = 0;

TOL = 1e-10;

%initial guess final state pressure
gamma = min(gamma1,gamma4);
pf01 = (((gamma-1)/2*(u1-u4)+a1+a4)/(a1/p1^((gamma-1)/(2*gamma))+a4/p4^((gamma-1)/(2*gamma))))^((2*gamma)/(gamma-1));
pf02 = (rho1*a1*p4+rho4*a4*p1-rho1*a1*rho4*a4*(u4-u1))/(rho1*a1+rho4*a4);

i = 0;
pfnew = 10;
pfold = 0;

while abs(pfnew-pfold)>1e-8
i = i + 1;

if i==1
pf = max(TOL,pf02);
elseif i==2
pf = max(TOL,pf01);
else
pf = max(TOL,pfnew);
end

if pf>=p1%left == 1 %shock
m1 = rho1*a1*sqrt(1+(gamma1+1)/(2*gamma1)*(pf-p1)/p1);
elseif pf<p1 %left == 2 % expansion
m1 = rho1*a1*(gamma1-1)/(2*gamma1)*(1-pf/p1)/(1-(pf/p1)^((gamma1-1)/(gamma1*2)));
else
m1 = 0;
end

if pf>=p4 %right == 1 %shock
m4 = rho4*a4*sqrt(1+(gamma4+1)/(2*gamma4)*(pf-p4)/p4);
elseif pf<p4 %right == 2 %expansion
m4 = rho4*a4*(gamma4-1)/(2*gamma4)*(1-pf/p4)/(1-(pf/p4)^((gamma4-1)/(gamma4*2)));
else
m4 = 0;
end

%pressure of final state
G = (m1*p4+m4*p1-m1*m4*(u4-u1))/(m1+m4) - pf;

if i>1
pfnew = pf - G*(pf-pfold+TOL)/(G-Gold+TOL);
end

%save data
Gold  = G;
pfold = pf;
end
pf = pfnew;

%velocity of final state
if pf>=p1 %left == 1 %shock
m1 = rho1*a1*sqrt(1+(gamma1+1)/(2*gamma1)*(pf-p1)/p1);
elseif pf<p1 %left == 2 %expansion
m1 = rho1*a1*(gamma1-1)/(2*gamma1)*(1-pf/p1)/(1-(pf/p1)^((gamma1-1)/(gamma1*2)));
else
m1 = 0;
end
    
if pf>=p4 %right == 1 % shock 
m4 = rho4*a4*sqrt(1+(gamma4+1)/(2*gamma4)*(pf-p4)/p4);
elseif pf<p4 %right == 2 %expansion
m4 = rho4*a4*(gamma4-1)/(2*gamma4)*(1-pf/p4)/(1-(pf/p4)^((gamma4-1)/(gamma4*2)));    
else
m4 = 0;
end

uf = (m1*u1+m4*u4-(p4-p1))/(m1+m4);


%density
if pf>=p1 %left == 1 %shock
rho2 = rho1*(1+(gamma1+1)/(gamma1-1)*pf/p1)/((gamma1+1)/(gamma1-1)+pf/p1);
elseif pf<p1 %left == 2 %expansion
rho2 = rho1*(pf/p1)^(1/gamma1);
else
rho2 = rho1;
end

if pf>=p4 %right ==1 %shock
rho3 = rho4*(1+(gamma4+1)/(gamma4-1)*pf/p4)/((gamma4+1)/(gamma4-1)+pf/p4);
elseif pf<p4 %right == 2 %expansion
rho3 = rho4*(pf/p4)^(1/gamma4);
else
rho3 = rho4;
end

%speed of sound
a2 = sqrt(gamma1*pf/rho2);
a3 = sqrt(gamma4*pf/rho3);


%direction characteristics
if pf>=p1
dxdts1 = u1 - a1*sqrt(1+(gamma1+1)/(2*gamma1)*(pf-p1)/p1);
elseif pf<p1
dxdtfirst1 = u1 - a1;
dxdtlast1  = uf - a2;
else
end

if pf>=p4
dxdts4 = u4 + a4*sqrt(1+(gamma4+1)/(2*gamma4)*(pf-p4)/p4);
elseif pf<p4
% dxdtfirst4 = u4 + a4;
% dxdtlast4  = uf + a3;
dxdtlast4  = u4 + a4;
dxdtfirst4 = uf + a3;
else
end

%contact discontinuity
dxdt0 = uf;

N = length(x);
dx = 0.5/N;

Nexp = 1e3;

x0 = x(1)-dx/2;
    
if pf>=p1
x1 = dxdts1 * t;
x2 = dxdts1 * t;
xleft = x2;
elseif pf<p1
x1 = dxdtfirst1 * t;
x2  = dxdtlast1 * t;
xleft(:,1) = x1+(1:Nexp)*(x2-x1)/Nexp;
else
end

x3    = dxdt0 * t;

if pf>=p4
x4 = dxdts4 * t;
x5 = dxdts4 * t;
xright = x4;
elseif pf<p4
x4 = dxdtfirst4 * t;
x5  = dxdtlast4 * t;
xright(:,1) = x4+(1:Nexp)*(x5-x4)/Nexp;
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact solution per cell                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x6 = x(end)+dx/2;

rhoexact   = zeros(N,1);
uexact     = zeros(N,1);
pexact     = zeros(N,1);
betaexact  = zeros(N,1);
alphaexact = zeros(N,1);

for j=1:floor((x1+0.25)/dx)
    rhoexact(j)   = rho1;
    uexact(j)     = u1;
    pexact(j)     = p1;
    betaexact(j)  = 1;
    alphaexact(j) = 1;
end

rhosol   = [rho1;rho1];
usol     = [u1;u1];
psol     = [p1;p1];
betasol  = [1;1];
alphasol = [1;1];

if pf<p1
for j=ceil((x1+0.25)/dx):floor((x2+0.25)/dx)
    xexp          = j*dx-0.25;
    uexact(j)     = 2/(gamma1+1)*(xexp/t+a1)+(gamma1-1)/(gamma1+1)*u1;
    aexpleft      = -(gamma1-1)/(gamma1+1)*(xexp/t-u1)+2/(gamma1+1)*a1;
    rhoexact(j)   = rho1*(aexpleft/a1)^(2/(gamma1-1));
    pexact(j)     = p1*(aexpleft/a1)^((2*gamma1)/(gamma1-1));
    betaexact(j)  = 1;
    alphaexact(j) = 1;
end
    usol(3:(Nexp+2),1)     = 2/(gamma1+1)*((xleft)/t+a1)+(gamma1-1)/(gamma1+1)*u1;
    asol(1:Nexp,1)         = -(gamma1-1)/(gamma1+1)*((xleft)/t-u1)+2/(gamma1+1)*a1;
    psol(3:(Nexp+2),1)     = p1*(asol(1:Nexp)/a1).^((2*gamma1)/(gamma1-1));
    rhosol(3:(Nexp+2),1)   = rho1*(asol(1:Nexp)/a1).^(2/(gamma1-1));
    betasol(3:(Nexp+2),1)  = ones(Nexp,1);
    alphasol(3:(Nexp+2),1) = ones(Nexp,1);
else
    rhosol(3,1)   = rho2;
    usol(3,1)     = uf;
    psol(3,1)     = pf;
    betasol(3,1)  = 1;
    alphasol(3,1) = 1;
end

for j=ceil((x2+0.25)/dx):floor((x3+0.25)/dx)
    rhoexact(j)   = rho2;
    uexact(j)     = uf;
    pexact(j)     = pf;
    betaexact(j)  = 1;
    alphaexact(j) = 1;
end

rhosol((length(rhosol)+1):(length(rhosol)+2),1)       = [rho2;rho2];
usol((length(usol)+1):(length(usol)+2),1)             = [uf;uf];
psol((length(psol)+1):(length(psol)+2),1)             = [pf;pf];
betasol((length(betasol)+1):(length(betasol)+2),1)    = [1;1];
alphasol((length(alphasol)+1):(length(alphasol)+2),1) = [1;1];

for j=ceil((x3+0.25)/dx):floor((x4+0.25)/dx)
    rhoexact(j)   = rho3;
    uexact(j)     = uf;
    pexact(j)     = pf;
    betaexact(j)  = 0;
    alphaexact(j) = 0;
end

rhosol((length(rhosol)+1):(length(rhosol)+2),1)       = [rho3;rho3];
usol((length(usol)+1):(length(usol)+2),1)             = [uf;uf];
psol((length(psol)+1):(length(psol)+2),1)             = [pf;pf];
betasol((length(betasol)+1):(length(betasol)+2),1)    = [0;0];
alphasol((length(alphasol)+1):(length(alphasol)+2),1) = [0;0];

if pf<p4
for j=ceil((x4+0.25)/dx):floor((x5+0.25)/dx)
    xexp          = j*dx-0.25;
    uexact(j)     = 2/(gamma4+1)*(xexp/t-a4)+(gamma4-1)/(gamma4+1)*u4;
    aexpright     = (gamma4-1)/(gamma4+1)*(xexp/t-u4)+2/(gamma4+1)*a4;
    rhoexact(j)   = rho4*(aexpright/a4)^(2/(gamma4-1));
    pexact(j)     = p4*(aexpright/a4)^((2*gamma4)/(gamma4-1));
    betaexact(j)  = 0;
    alphaexact(j) = 0;
end



    usol((length(usol)+1):(length(usol)+Nexp),1)             = 2/(gamma4+1)*(xright/t-a4)+(gamma4-1)/(gamma4+1)*u4;
    asol(1:Nexp,1)                                           = (gamma4-1)/(gamma4+1)*(xright/t-u4)+2/(gamma4+1)*a4;
    psol((length(psol)+1):(length(psol)+Nexp),1)             = p4*(asol(1:Nexp)/a4).^((2*gamma4)/(gamma4-1));
    rhosol((length(rhosol)+1):(length(rhosol)+Nexp),1)       = rho4*(asol(1:Nexp)/a4).^(2/(gamma4-1));
    betasol((length(betasol)+1):(length(betasol)+Nexp),1)    = zeros(Nexp,1);
    alphasol((length(alphasol)+1):(length(alphasol)+Nexp),1) = zeros(Nexp,1);
else
    rhosol((length(rhosol)+1),1)     = rho3;
    usol((length(usol)+1),1)         = uf;
    psol((length(psol)+1),1)         = pf;
    betasol((length(betasol)+1),1)   = 0;
    alphasol((length(alphasol)+1),1) = 0;
end

for j=ceil((x5+0.25)/dx):N
    rhoexact(j)   = rho4;
    uexact(j)     = u4;
    pexact(j)     = p4;
    betaexact(j)  = 0;
    alphaexact(j) = 0;
end
rhosol((length(rhosol)+1):(length(rhosol)+2),1)       = [rho4;rho4];
usol((length(usol)+1):(length(usol)+2),1)             = [u4;u4];
psol((length(psol)+1):(length(psol)+2),1)             = [p4;p4];
betasol((length(betasol)+1):(length(betasol)+2),1)    = [0;0];
alphasol((length(alphasol)+1):(length(alphasol)+2),1) = [0;0];


xplot = [x0;x1;x1;xleft;x3;x3;xright;x5;x5;x6];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figures                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
% subplot(2,2,1)
hold on
plot(xplot,rhosol,'r')%,'linewidth',1)
plot(x,rho,'o','MarkerSize',5,'linewidth',1)
% set(gca,'fontsize',14)
% grid on
xlim([-0.25 0.25])
% ylim([0.2 1.2])
xlabel('x')
ylabel('\rho')

figure(2)
% subplot(2,2,2)
hold on
plot(xplot,usol,'r')%,'linewidth',1)
plot(x,u,'o','MarkerSize',5,'linewidth',1)
% set(gca,'fontsize',14)
% grid on 
xlim([-0.25 0.25])
% ylim([-0.2 2])
xlabel('x')
ylabel('u')

figure(3)
% subplot(2,2,3)
hold on
plot(xplot,psol,'r')%,'linewidth',1)
plot(x,p,'o','MarkerSize',5,'linewidth',1)
% set(gca,'fontsize',14)
% grid on
xlim([-0.25 0.25])
% ylim([0 11])
xlabel('x')
ylabel('p')

figure(4)
% subplot(2,2,4)
subplot(1,2,1)
hold on
plot(xplot,alphasol,'r')%,'linewidth',1)
plot(x,alpha,'o','MarkerSize',5,'linewidth',1)
% set(gca,'fontsize',14)
% legend('\beta','\alpha','exact',3)
% grid on
xlim([0 0.25])
ylim([-0.2 1.2])
xlabel('x')
ylabel('\alpha')

subplot(1,2,2)
hold on
plot(xplot,betasol,'r')%,'linewidth',1)
plot(x,beta,'o','MarkerSize',5,'linewidth',1)
% set(gca,'fontsize',14)
xlim([0 0.25])
ylim([-0.2 1.2])
xlabel('x')
ylabel('\beta')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Errors                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L1(1) = 1/length(rho)*sum(abs(rhoexact-rho));
L1(2) = 1/length(u)*sum(abs(uexact-u));
L1(3) = 1/length(p)*sum(abs(pexact-p));
L1(4) = 1/length(beta)*sum(abs(betaexact-beta));
L1(5) = 1/length(alpha)*sum(abs(alphaexact-alpha));

L2(1) = sqrt(1/length(rho)*sum((rhoexact-rho).^2));
L2(2) = sqrt(1/length(u)*sum((uexact-u).^2));
L2(3) = sqrt(1/length(p)*sum((pexact-p).^2));
L2(4) = sqrt(1/length(beta)*sum((betaexact-beta).^2));
L2(5) = sqrt(1/length(alpha)*sum((alphaexact-alpha).^2));

Linf(1) = max(abs(rhoexact-rho));
Linf(2) = max(abs(uexact-u));
Linf(3) = max(abs(pexact-p));
Linf(4) = max(abs(betaexact-beta));
Linf(5) = max(abs(alphaexact-alpha));