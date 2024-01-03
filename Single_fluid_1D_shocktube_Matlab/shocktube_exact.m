clear %all
close all
clc

%Initial conditions
% % Case 1 Double shockwave
% T_final = 0.12;
% 
% rho1 = 1.0;
% u1   = 5.0;
% p1   = 2.0;
% 
% rho4 = 1.0;
% u4   = 0;
% p4   = 0.1;
% 
% gamma = 1.4;

% Case 2a Sod problem
T_final = 0.25;

rho1 = 1.0;
u1   = 0.0;
p1   = 1.0;

rho4 = 0.125;
u4   = 0.0;
p4   = 0.1;

gamma = 1.4;

% % Case 2b reverse Sod problem
% T_final = 0.25;
% 
% rho4 = 1.0;
% u4   = 0.0;
% p4   = 1.0;
% 
% rho1 = 0.125;
% u1   = 0.0;
% p1   = 0.1;
% 
% gamma = 1.4;

% % Case 3 Double expansion wave
% T_final = 0.15;
% 
% rho1 = 1.0;
% u1   = -2.0;
% p1   = 0.4;
% 
% rho4 = 1.0;
% u4   = 2.0;
% p4   = 0.4;
% 
% gamma = 1.4;

%initial speed of sound
a1   = sqrt(gamma*p1/rho1);
a4   = sqrt(gamma*p4/rho4);

ptest = 1;
ptestold = 0;

%SHOCK OR EXPANSION
if u1 > u4
    if p1 > p4
        while abs(ptest-ptestold)>1e-4
            ptestold = ptest;            
            ptest = p4 + rho4*a4*sqrt(1+(gamma+1)/(2*gamma)*(ptestold-p4)/p4)*(u1-u4);
        end
        if ptest > p1
            left = 1;  disp('left wave is a expansion wave')
            right = 1; disp('right wave is a shockwave')
        elseif ptest < p1
            left = 1;  disp('left wave is a shockwave')
            right = 1; disp('right wave is an expansion wave')
        else
            left = 1;  disp('left wave is a shockwave')
            right = 0;
        end
    
    elseif p1 < p4
        while abs(ptest-ptestold)>1e-4
            ptestold = ptest;            
            ptest = p1+ rho1*a1*sqrt(1+(gamma+1)/(2*gamma)*(ptestold-p1)/p1)*(u1-u4);
        end
        if ptest > p4
            left = 1;  disp('left wave is a shockwave')
            right = 1; disp('right wave is a shockwave')
        elseif ptest < p4
            left = 1;  disp('left wave is a shockwave')
            right = 1; disp('right wave is an expansion wave')
        else
            left = 1;  disp('left wave is a shockwave')
            right = 0;
        end
            
    else
        left = 1;  disp('left wave is a shockwave')
        right = 1; disp('right wave is a shockwave')
    end
    
elseif u1 < u4
    if p1 > p4
        while abs(ptest-ptestold)>1e-4
            ptestold = ptest;            
            ptest = p1+rho1*a1*(gamma-1)/(2*gamma)*(1-ptestold/p1)/(1-(ptestold/p1)^((gamma-1)/(2*gamma)))*(u1-u4);
        end
        if ptest > p4
            left = 2;  disp('left wave is an expansion wave')
            right = 1; disp('right wave is a shockwave')
        elseif ptest < p4
            left = 2;  disp('left wave is an expansion wave')
            right = 2; disp('right wave is an expansion wave')
        else
            left = 2;  disp('left wave is an expansion wave')
            right = 0;
        end
        
    elseif p1 < p4
        while abs(ptest-ptestold)>1e-4
            ptestold = ptest;            
            ptest = p4+rho4*a4*(gamma-1)/(2*gamma)*(1-ptestold/p4)/(1-(ptestold/p4)^((gamma-1)/(2*gamma)))*(u1-u4);
        end
        if ptest > p1
            left = 2;  disp('left wave is an expansion wave')
            right = 1; disp('right wave is a shockwave')
        elseif ptest < p1
            left = 2;  disp('left wave is an expansion wave')
            right = 2; disp('right wave is an expansion wave')
        else
            left = 2;  disp('left wave is an expansion wave')
            right = 0;
        end
                
    else
        left = 2;  disp('left wave is an expansion wave')
        right = 2; disp('right wave is an expansion wave')
    end
    
else
    if p1 > p4
        left = 2;  disp('left wave is an expansion wave')
        right = 1; disp('right wave is a shockwave')
    elseif p1 < p4
        left = 1;  disp('left wave is a shockwave')
        right = 2; disp('right wave is an expansion wave')
    else
        left = 0;
        right = 0;
    end
end



%initial guess final state pressure
pf01 = (((gamma-1)/2*(u1-u4)+a1+a4)/(a1*p1^((2*gamma)/(gamma-1))+a4*p4^((2*gamma)/(gamma-1))))^((2*gamma)/(gamma-1));
pf02 = (rho1*a1*p4+rho4*a4*p1-rho1*a1*rho4*a4*(u4-u1))/(rho1*a1+rho4*a4);

i = 0;
pfnew = 10;
pfold = 0;

while abs(pfnew-pfold)>1e-4
i = i + 1;

if i==1
pf = pf02;
elseif i==2
pf = pf01;
else
pf = pfnew;
end

if left == 1 %shock
m1 = rho1*a1*sqrt(1+(gamma+1)/(2*gamma)*(pf-p1)/p1);
elseif left == 2 % expansion
m1 = rho1*a1*(gamma-1)/(2*gamma)*(1-pf/p1)/(1-(pf/p1)^((gamma-1)/(gamma*2)));
else
m1 = 0;
end

if right == 1 %shock
m4 = rho4*a4*sqrt(1+(gamma+1)/(2*gamma)*(pf-p4)/p4);
elseif right == 2 %expansion
m4 = rho4*a4*(gamma-1)/(2*gamma)*(1-pf/p4)/(1-(pf/p4)^((gamma-1)/(gamma*2)));
else
m4 = 0;
end

%pressure of final state
G = (m1*p4+m4*p1-m1*m4*(u4-u1))/(m1+m4) - pf;

if i>1
pfnew = pf - (G*(pf-pfold))/(G-Gold);
end

%save data
Gold = G;
pfold = pf;
end
pf = real(pfnew);

%velocity of final state
if left == 1 %shock
m1 = rho1*a1*sqrt(1+(gamma+1)/(2*gamma)*(pf-p1)/p1);
elseif left == 2 %expansion
m1 = rho1*a1*(gamma-1)/(2*gamma)*(1-pf/p1)/(1-(pf/p1)^((gamma-1)/(gamma*2)));
else
m1 = 0;
end
    
if right == 1 % shock 
m4 = rho4*a4*sqrt(1+(gamma+1)/(2*gamma)*(pf-p4)/p4);
elseif right == 2 %expansion
m4 = rho4*a4*(gamma-1)/(2*gamma)*(1-pf/p4)/(1-(pf/p4)^((gamma-1)/(gamma*2)));    
else
m4 = 0;
end

uf = (m1*u1+m4*u4-(p4-p1))/(m1+m4);


%density
if left == 1 %shock
rho2 = rho1*(1+(gamma+1)/(gamma-1)*pf/p1)/((gamma+1)/(gamma-1)+pf/p1);
elseif left == 2 %expansion
rho2 = rho1*(pf/p1)^(1/gamma);
else
rho2 = rho1;
end

if right ==1 %shock
rho3 = rho4*(1+(gamma+1)/(gamma-1)*pf/p4)/((gamma+1)/(gamma-1)+pf/p4);
elseif right == 2 %expansion
rho3 = rho4*(pf/p4)^(1/gamma);
else
rho3 = rho4;
end

%speed of sound
a2 = sqrt(gamma*pf/rho2);
a3 = sqrt(gamma*pf/rho3);


%direction characteristics
if left == 1 %shock
dxdts1 = u1 - a1*sqrt(1+(gamma+1)/(2*gamma)*(pf-p1)/p1);
elseif left == 2 %expansion
dxdtfirst1 = u1 - a1;    
dxdtlast1  = uf - a2;
else
end

if right == 1 %shock
dxdts4 = u4 + a4*sqrt(1+(gamma+1)/(2*gamma)*(pf-p4)/p4);
elseif right == 2 %expansion
dxdtfirst4 = u4 + a4;
dxdtlast4  = uf + a3;
else
end

%contact discontinuity
dxdt0 = uf;



tlevel = 1;

dt = T_final/50;
t = dt;

xsleft = [];
xsright = [];
xexpleft = [];
xexpright = [];
xcd = [];

while (t <= T_final)

if left == 1
xsleft = dxdts1 * t;
elseif left ==2
xfirstleft = dxdtfirst1 * t;
xlastleft  = dxdtlast1 * t;
else
end

if right == 1
xsright = dxdts4 * t;
elseif right == 2
xfirstright = dxdtfirst4 * t;
xlastright  = dxdtlast4 * t;
else
end

xcd    = dxdt0 * t;


N = 100;
j = 0:N;

if left == 2
dxleft = (xlastleft - xfirstleft)/N;
xexpleft = xfirstleft+j*dxleft;
uexpleft = zeros(1,N+1); aexpleft = zeros(1,N+1);
for i=1:N+1
%Method of characteristics
uexpleft(i)  = 2/(gamma+1)*(xexpleft(i)/t+a1)+(gamma-1)/(gamma+1)*u1;
aexpleft(i)  = -(gamma-1)/(gamma+1)*(xexpleft(i)/t-u1)+2/(gamma+1)*a1;
end

%pressure and density inside expansion
pexpleft  = p1*(aexpleft/a1).^((2*gamma)/(gamma-1));
rhoexpleft  = rho1*(aexpleft/a1).^(2/(gamma-1));
end

if right == 2
dxright = (xfirstright - xlastright)/N;
xexpright = xlastright+j*dxright;
uexpright = zeros(1,N+1); aexpright = zeros(1,N+1);
for i=1:N+1
%Method of characteristics
uexpright(i) = 2/(gamma+1)*(xexpright(i)/t-a4)+(gamma-1)/(gamma+1)*u4;
aexpright(i) = (gamma-1)/(gamma+1)*(xexpright(i)/t-u4)+2/(gamma+1)*a4;
end

%pressure and density inside expansion
pexpright = p4*(aexpright/a4).^((2*gamma)/(gamma-1));
rhoexpright = rho4*(aexpright/a4).^(2/(gamma-1));
end

x = [-0.5,xsleft,xsleft,xexpleft,xcd,xcd,xsright,xsright,xexpright,0.5];

if left == 1
    rholeft = [rho1, rho2];
    uleft   = [u1,uf];
    pleft   = [p1,pf];
elseif left == 2
    rholeft = rhoexpleft;
    uleft   = uexpleft;
    pleft   = pexpleft;
else
    rholeft = [];
    uleft   = [];
    pleft   = [];    
end
if right == 1
    rhoright = [rho3, rho4];
    uright   = [uf,u4];
    pright   = [pf,p4];    
elseif right == 2
    rhoright = rhoexpright;
    uright   = uexpright;
    pright   = pexpright;
else
    rhoright = [];
    uright   = [];
    pright   = [];
end

rhosolution = [rho1,rholeft,rho2,rho3,rhoright,rho4];
usolution   = [u1,uleft,uf,uf,uright,u4];
psolution   = [p1,pleft,pf,pf,pright,p4];

Colors = [0.6350, 0.0780, 0.1840
          0.0000, 0.4470, 0.7410
          0.4660, 0.6740, 0.1880];

subplot(2,2,1)
plot(x,rhosolution,'-','linewidth',3,'Color',Colors(1,:))
ylabel('density')

subplot(2,2,2)
plot(x,usolution,'-','linewidth',3,'Color',Colors(1,:))
ylabel('velocity')

subplot(2,2,3)
plot(x,psolution,'-','linewidth',3,'Color',Colors(1,:))
ylabel('pressure')

subplot(2,2,4)
hold on
axis([-0.5 0.5 0 T_final])
if left == 1
plot([0,xsleft],[0,t],'Color',Colors(1,:))
elseif left == 2
xfanleft = zeros(1,5);
for i = 0:4
xfanleft(i+1) = xfirstleft+i/4*(xlastleft-xfirstleft);
plot([0,xfanleft(i+1)],[0,t],'Color',Colors(2,:))
end
else
end

if right == 1
plot([0,xsright],[0,t],'Color',Colors(1,:))
elseif right == 2
xfanright = zeros(1,5);
for i = 0:4
xfanright(i+1) = xlastright+i/4*(xfirstright-xlastright);
plot([0,xfanright(i+1)],[0,t],'Color',Colors(2,:))
end
else
end

plot([0,xcd],[0,t],':','Color',Colors(3,:))
hold off

F(tlevel) = getframe;
tlevel = tlevel + 1;

t = t + dt;

end
