function postprocessen(W,t,x)

global gamma1 gamma2

L = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate exact solution                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rho   = W(1,:);
u     = W(2,:);
p     = W(3,:);
beta  = W(4,:);
alpha = W(5,:);
gamma = 1./(alpha/gamma1+(1-alpha)/gamma2);
e     = p./((gamma-1).*rho);

t     = t(end);

%Initial conditions
rho1   = rho(1);
u1     = u(1);
p1     = p(1);
beta1  = beta(1);
alpha1 = alpha(1);
e1     = e(1);
gammaL = 1/(alpha1/gamma1+(1-alpha1)/gamma2);

rho4   = rho(end);
u4     = u(end);
p4     = p(end);
beta4  = beta(end);
alpha4 = alpha(end);
e4     = e(end);
gammaR = 1/(alpha4/gamma1+(1-alpha4)/gamma2);


%initial speed of sound
a1   = sqrt(gammaL*p1/rho1);
a4   = sqrt(gammaR*p4/rho4);

eps = 1e-12;

%initial guess final state pressure
gamma0 = min(gammaL,gammaR);
pf01 = (((gamma0-1)/2*(u1-u4)+a1+a4)/(a1/p1^((gamma0-1)/(2*gamma0))+a4/p4^((gamma0-1)/(2*gamma0))))^((2*gamma0)/(gamma0-1));
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
        m1 = rho1*a1*sqrt(1+(gammaL+1)/(2*gammaL)*(pf-p1)/p1);
    elseif pf<p1 %left == 2 % expansion
        m1 = rho1*a1*(gammaL-1)/(2*gammaL)*(1-pf/p1)/(1-(pf/p1)^((gammaL-1)/(gammaL*2)));
    end

    if pf>=p4 %right == 1 %shock
        m4 = rho4*a4*sqrt(1+(gammaR+1)/(2*gammaR)*(pf-p4)/p4);
    elseif pf<p4 %right == 2 %expansion
        m4 = rho4*a4*(gammaR-1)/(2*gammaR)*(1-pf/p4)/(1-(pf/p4)^((gammaR-1)/(gammaR*2)));
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
    m1 = rho1*a1*sqrt(1+(gammaL+1)/(2*gammaL)*(pf-p1)/p1);
elseif pf<p1 %left == 2 %expansion
    m1 = rho1*a1*(gammaL-1)/(2*gammaL)*(1-pf/p1)/(1-(pf/p1)^((gammaL-1)/(gammaL*2)));
end
    
if pf>=p4 %right == 1 % shock 
    m4 = rho4*a4*sqrt(1+(gammaR+1)/(2*gammaR)*(pf-p4)/p4);
elseif pf<p4 %right == 2 %expansion
    m4 = rho4*a4*(gammaR-1)/(2*gammaR)*(1-pf/p4)/(1-(pf/p4)^((gammaR-1)/(gammaR*2)));
end

uf = (m1*u1+m4*u4-(p4-p1))/(m1+m4);


%density
if pf>=p1 %left == 1 %shock
    rho2 = rho1*(1+(gammaL+1)/(gammaL-1)*pf/p1)/((gammaL+1)/(gammaL-1)+pf/p1);
elseif pf<p1 %left == 2 %expansion
    rho2 = rho1*(pf/p1)^(1/gammaL);
end

if pf>=p4 %right ==1 %shock
    rho3 = rho4*(1+(gammaR+1)/(gammaR-1)*pf/p4)/((gammaR+1)/(gammaR-1)+pf/p4);
elseif pf<p4 %right == 2 %expansion
    rho3 = rho4*(pf/p4)^(1/gammaR);
end

%speed of sound
a2 = sqrt(gammaL*pf/rho2);
a3 = sqrt(gammaR*pf/rho3);

e2 = pf/(((1/(alpha1/gamma1+(1-alpha1)/gamma2))-1)*rho2);
e3 = pf/(((1/(alpha4/gamma1+(1-alpha4)/gamma2))-1)*rho3);

%direction characteristics
if pf>=p1
    dxdts1 = u1 - a1*sqrt(1+(gammaL+1)/(2*gammaL)*(pf-p1)/p1);
elseif pf<p1
    dxdtfirst1 = u1 - a1;
    dxdtlast1  = uf - a2;
end

if pf>=p4
    dxdts4 = u4 + a4*sqrt(1+(gammaR+1)/(2*gammaR)*(pf-p4)/p4);
elseif pf<p4
    % dxdtfirst4 = u4 + a4;
    % dxdtlast4  = uf + a3;
    dxdtlast4  = u4 + a4;
    dxdtfirst4 = uf + a3;
end

%contact discontinuity
dxdt0 = uf;

N = length(x);
dx = L/N;

Nexp = 1e3;

x0 = x(1)-dx/2;
x6 = x(end)+dx/2;

if pf>=p1
    x1 = dxdts1 * t+(x6+x0)/2;
    x2 = dxdts1 * t+(x6+x0)/2;
    xleft = x2;
elseif pf<p1
    x1 = dxdtfirst1 * t+(x6+x0)/2;
    x2  = dxdtlast1 * t+(x6+x0)/2;
    xleft(:,1) = x1+(1:Nexp)*(x2-x1)/Nexp;
end

x3 = dxdt0 * t+(x6+x0)/2;

if pf>=p4
    x4 = dxdts4 * t+(x6+x0)/2;
    x5 = dxdts4 * t+(x6+x0)/2;
    xright = x4;
elseif pf<p4
    x4 = dxdtfirst4 * t+(x6+x0)/2;
    x5  = dxdtlast4 * t+(x6+x0)/2;
    xright(:,1) = x4+(1:Nexp)*(x5-x4)/Nexp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% exact solution per cell                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rhoexact   = zeros(N,1);
uexact     = zeros(N,1);
pexact     = zeros(N,1);
betaexact  = zeros(N,1);
alphaexact = zeros(N,1);
eexact     = zeros(N,1);

% left initial state
for j=1:floor((x1+L/2)/dx+1/2)
    rhoexact(j)   = rho1;
    uexact(j)     = u1;
    pexact(j)     = p1;
    betaexact(j)  = beta1;
    alphaexact(j) = alpha1;
    eexact(j)     = p1/((gamma1-1)*rho1);
end

rhosol   = [rho1;rho1];
usol     = [u1;u1];
psol     = [p1;p1];
betasol  = [beta1;beta1];
alphasol = [alpha1;alpha1];
esol     = [e1;e1];

% left wave state
if pf<p1
    for j=ceil((x1+L/2)/dx+1/2):ceil((x2+L/2)/dx+1/2)
        xexp          = (j-1/2)*dx;
        uexact(j)     = 2/(gammaL+1)*((xexp-(x6+x0)/2)/t+a1)+(gammaL-1)/(gammaL+1)*u1;
        aexpleft      = -(gammaL-1)/(gammaL+1)*((xexp-(x6+x0)/2)/t-u1)+2/(gammaL+1)*a1;
        rhoexact(j)   = rho1*(aexpleft/a1)^(2/(gammaL-1));
        pexact(j)     = p1*(aexpleft/a1)^((2*gammaL)/(gammaL-1));
        betaexact(j)  = beta1;
        alphaexact(j) = alpha1;
        eexact(j)     = e1;
    end
        usol(3:(Nexp+2),1)     = 2/(gammaL+1)*((xleft-(x6+x0)/2)/t+a1)+(gammaL-1)/(gammaL+1)*u1;
        asol(1:Nexp,1)         = -(gammaL-1)/(gammaL+1)*((xleft-(x6+x0)/2)/t-u1)+2/(gammaL+1)*a1;
        psol(3:(Nexp+2),1)     = p1*(asol(1:Nexp)/a1).^((2*gammaL)/(gammaL-1));
        rhosol(3:(Nexp+2),1)   = rho1*(asol(1:Nexp)/a1).^(2/(gammaL-1));
        betasol(3:(Nexp+2),1)  = beta1*ones(1,Nexp);
        alphasol(3:(Nexp+2),1) = alpha1*ones(1,Nexp);
        esol(3:(Nexp+2),1)     = psol(3:(Nexp+2),1)./(((1./(alphasol(3:(Nexp+2),1)/gamma1+(1-alphasol(3:(Nexp+2),1))/gamma2))-1).*rhosol(3:(Nexp+2),1));
else
    rhosol(3,1)   = rho2;
    usol(3,1)     = uf;
    psol(3,1)     = pf;
    betasol(3,1)  = beta1;
    alphasol(3,1) = alpha1;
    esol(3,1)     = e2;
end

%  final state 2
for j=ceil((x2+L/2)/dx+1/2):ceil((x3+L/2)/dx+1/2)
    rhoexact(j)   = rho2;
    uexact(j)     = uf;
    pexact(j)     = pf;
    betaexact(j)  = beta1;
    alphaexact(j) = alpha1;
    eexact(j)     = e2;
end

rhosol((length(rhosol)+1):(length(rhosol)+2),1)       = [rho2;rho2];
usol((length(usol)+1):(length(usol)+2),1)             = [uf;uf];
psol((length(psol)+1):(length(psol)+2),1)             = [pf;pf];
betasol((length(betasol)+1):(length(betasol)+2),1)    = [beta1;beta1];
alphasol((length(alphasol)+1):(length(alphasol)+2),1) = [alpha1;alpha1];
esol((length(esol)+1):(length(esol)+2),1)             = [e2;e2];

% final state 3
for j=ceil((x3+L/2)/dx+1/2):floor((x4+L/2)/dx+1/2)
    rhoexact(j)   = rho3;
    uexact(j)     = uf;
    pexact(j)     = pf;
    betaexact(j)  = beta4;
    alphaexact(j) = alpha4;
    eexact(j)     = e3;
end

rhosol((length(rhosol)+1):(length(rhosol)+2),1)       = [rho3;rho3];
usol((length(usol)+1):(length(usol)+2),1)             = [uf;uf];
psol((length(psol)+1):(length(psol)+2),1)             = [pf;pf];
betasol((length(betasol)+1):(length(betasol)+2),1)    = [beta4;beta4];
alphasol((length(alphasol)+1):(length(alphasol)+2),1) = [alpha4;alpha4];
esol((length(esol)+1):(length(esol)+2),1)             = [e3;e3];

if pf<p4
    for j=ceil((x4+L/2)/dx+1/2):floor((x5+L/2)/dx+1/2)
        xexp          = (j-1/2)*dx;
        uexact(j)     = 2/(gammaR+1)*((xexp-(x6+x0)/2)/t-a4)+(gammaR-1)/(gammaR+1)*u4;
        aexpright     = (gammaR-1)/(gammaR+1)*((xexp-(x6+x0)/2)/t-u4)+2/(gammaR+1)*a4;
        rhoexact(j)   = rho4*(aexpright/a4)^(2/(gammaR-1));
        pexact(j)     = p4*(aexpright/a4)^((2*gammaR)/(gammaR-1));
        betaexact(j)  = beta4;
        alphaexact(j) = alpha4;
        eexact(j)     = pexact(j)/(((1/(alphaexact(j)/gamma1+(1-alphaexact(j))/gamma2))-1)*rhoexact(j));
    end

    usol((length(usol)+1):(length(usol)+Nexp),1)             = 2/(gammaR+1)*((xright-(x6+x0)/2)/t-a4)+(gammaR-1)/(gammaR+1)*u4;
    asol(1:Nexp,1)                                           = (gammaR-1)/(gammaR+1)*((xright-(x6+x0)/2)/t-u4)+2/(gammaR+1)*a4;
    psol((length(psol)+1):(length(psol)+Nexp),1)             = p4*(asol(1:Nexp)/a4).^((2*gammaR)/(gammaR-1));
    rhosol((length(rhosol)+1):(length(rhosol)+Nexp),1)       = rho4*(asol(1:Nexp)/a4).^(2/(gammaR-1));
    betasol((length(betasol)+1):(length(betasol)+Nexp),1)    = beta4*ones(Nexp,1);
    alphasol((length(alphasol)+1):(length(alphasol)+Nexp),1) = alpha4*ones(Nexp,1);
%     esol((length(esol)+1):(length(esol)+Nexp),1)     = psol((length(esol)+1):(length(esol)+Nexp),1)./(((1./(esol((length(esol)+1):(length(esol)+Nexp),1)/gamma1+(1-alphasol((length(esol)+1):(length(esol)+Nexp),1))/gamma2))-1)*rhosol((length(esol)+1):(length(esol)+Nexp),1));
else
    rhosol((length(rhosol)+1),1)     = rho3;
    usol((length(usol)+1),1)         = uf;
    psol((length(psol)+1),1)         = pf;
    betasol((length(betasol)+1),1)   = beta4;
    alphasol((length(alphasol)+1),1) = alpha4;
%     esol((length(esol)+1),1)     = e3;
end

for j=ceil((x5+L/2)/dx+1/2):N
    rhoexact(j)   = rho4;
    uexact(j)     = u4;
    pexact(j)     = p4;
    betaexact(j)  = beta4;
    alphaexact(j) = alpha4;
    eexact(j)     = e4;
end

rhosol((length(rhosol)+1):(length(rhosol)+2),1)       = [rho4;rho4];
usol((length(usol)+1):(length(usol)+2),1)             = [u4;u4];
psol((length(psol)+1):(length(psol)+2),1)             = [p4;p4];
betasol((length(betasol)+1):(length(betasol)+2),1)    = [beta4;beta4];
alphasol((length(alphasol)+1):(length(alphasol)+2),1) = [alpha4;alpha4];
% esol((length(esol)+1):(length(esol)+2),1)             = [e4;e4];

xplot = [x0;x1;x1;xleft;x3;x3;xright;x5;x5;x6];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Errors                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L1(1) = 1/length(rho)*sum(abs(rhoexact-rho'));
L1(2) = 1/length(u)*sum(abs(uexact-u'));
L1(3) = 1/length(p)*sum(abs(pexact-p'));
L1(4) = 1/length(beta)*sum(abs(betaexact-beta'));
L1(5) = 1/length(alpha)*sum(abs(alphaexact-alpha'));
L1(6) = 1/length(e)*sum(abs(eexact-e'));

L2(1) = sqrt(1/length(rho)*sum((rhoexact-rho').^2));
L2(2) = sqrt(1/length(u)*sum((uexact-u').^2));
L2(3) = sqrt(1/length(p)*sum((pexact-p').^2));
L2(4) = sqrt(1/length(beta)*sum((betaexact-beta').^2));
L2(5) = sqrt(1/length(alpha)*sum((alphaexact-alpha').^2));
L2(6) = sqrt(1/length(e)*sum((eexact-e').^2));

Linf(1) = max(abs(rhoexact-rho'));
Linf(2) = max(abs(uexact-u'));
Linf(3) = max(abs(pexact-p'));
Linf(4) = max(abs(betaexact-beta'));
Linf(5) = max(abs(alphaexact-alpha'));
Linf(6) = max(abs(eexact-e'));

% L(6) is the error for the internal energy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create figures                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xplot = (xplot-1/2)/2;
% density
figure
hold on
plot(xplot,rhosol,'b','linewidth',1)
plot(x,rho,'o')%,'MarkerSize',7,'linewidth',1)
set(gca,'fontsize',16)
grid on
xlim([-L/2 L/2])
Y = ylim;
title('density')
xlabel('x')
ylabel('\rho')
legend('exact','approx')
text(0.5,Y(2)*19/20,strcat('L_1 = ',num2str(L1(1),'%e')),'fontsize',14)
text(0.5,Y(2)*18/20,strcat('L_2 = ',num2str(L2(1),'%e')),'fontsize',14)
text(0.5,Y(2)*17/20,strcat('L_\infty = ',num2str(Linf(1),'%e')),'fontsize',14)

% velocity
figure
hold on
plot(xplot,usol,'b','linewidth',1)
plot(x,u,'o','MarkerSize',7,'linewidth',1)
set(gca,'fontsize',14)
grid on 
xlim([-L/2 L/2])
Y = ylim;
title('velocity')
xlabel('x')
ylabel('u')
legend('exact','approx')
text(0.5,Y(2)*19/20,strcat('L_1 = ',num2str(L1(2),'%e')),'fontsize',14)
text(0.5,Y(2)*18/20,strcat('L_2 = ',num2str(L2(2),'%e')),'fontsize',14)
text(0.5,Y(2)*17/20,strcat('L_\infty = ',num2str(Linf(2),'%e')),'fontsize',14)

% pressure
figure
hold on
plot(xplot,psol,'b','linewidth',1)
plot(x,p,'o','MarkerSize',7,'linewidth',1)
set(gca,'fontsize',16)
grid on
xlim([-L/2 L/2])
Y = ylim;
title('pressure')
xlabel('x')
ylabel('p')
legend('exact','approx')
text(0.5,Y(2)*19/20,strcat('L_1 = ',num2str(L1(3),'%e')),'fontsize',14)
text(0.5,Y(2)*18/20,strcat('L_2 = ',num2str(L2(3),'%e')),'fontsize',14)
text(0.5,Y(2)*17/20,strcat('L_\infty = ',num2str(Linf(3),'%e')),'fontsize',14)

% fractions
figure
subplot(1,2,1)
hold on
plot(xplot,alphasol,'b','linewidth',1)
plot(x,alpha,'o','MarkerSize',7,'linewidth',1)
set(gca,'fontsize',14)
grid on
if mean(u)>=0
    xlim([0 L/2])
else
    xlim([-L/2 0])
end
ylim([-0.2 1.2])
title('volume fraction')
xlabel('x')
ylabel('\alpha')
legend('exact','approx')
% text(0.5,Y(2)*19/20,strcat('L_1 = ',num2str(L1(4),'%e')),'fontsize',14)
% text(0.5,Y(2)*18/20,strcat('L_2 = ',num2str(L2(4),'%e')),'fontsize',14)
% text(0.5,Y(2)*17/20,strcat('L_\infty = ',num2str(Linf(4),'%e')),'fontsize',14)
subplot(1,2,2)
hold on
plot(xplot,betasol,'b','linewidth',1)
plot(x,beta,'o','MarkerSize',7,'linewidth',1)
set(gca,'fontsize',14)
grid on
if mean(u)>=0
    xlim([0 L/2])
else
    xlim([-L/2 0])
end
ylim([-0.2 1.2])
title('mass fraction')
xlabel('x')
ylabel('\beta')
legend('exact','approx')
% text(0.5,Y(2)*19/20,strcat('L_1 = ',num2str(L1(5),'%e')),'fontsize',14)
% text(0.5,Y(2)*18/20,strcat('L_2 = ',num2str(L2(5),'%e')),'fontsize',14)
% text(0.5,Y(2)*17/20,strcat('L_\infty = ',num2str(Linf(5),'%e')),'fontsize',14)

% % Internal energy
% figure
% hold on
% plot(xplot,esol,'r','linewidth',1)
% plot(x,e,'o','MarkerSize',7,'linewidth',1)
% set(gca,'fontsize',14)
% grid on
% xlim([0 1])
% Y = ylim;  % ylim([0 1.2])
% title('internal energy')
% xlabel('x')
% ylabel('e')
% legend('exact','approx')
% text(0.3,Y(2)*19/20,strcat('L_1 = ',num2str(L1(6),'%e')),'fontsize',14)
% text(0.3,Y(2)*18/20,strcat('L_2 = ',num2str(L2(6),'%e')),'fontsize',14)
% text(0.3,Y(2)*17/20,strcat('L_\infty = ',num2str(Linf(6),'%e')),'fontsize',14)