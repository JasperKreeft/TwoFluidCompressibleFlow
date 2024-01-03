clear all
close all
clc

u = 1.5;
uL = 1.5;
du = -1e-3;
rho = 1.2;
rhoL = 1.2;
p = 1;
pL = 1;
alpha = 0.4;
alphaL = 0.4;
beta = 0.8;
betaL = 0.8;
gamma1 = 1.8;
gamma2 = 1.2;
eL = alphaL*pL/(rhoL*(gamma1-1))+(1-alphaL)*pL/(rhoL*(gamma2-1));

A = [rhoL uL pL alphaL betaL 0];
tau = -eps;

while tau<0 %for i= 1:round((0.5-uL)/du)

u = u+du;

e1 = alpha*p/(beta*rho*(gamma1-1));
e1_rho1 = -alpha^2*p/(beta^2*rho^2*(gamma1-1));
e1_p = alpha/(beta*rho*(gamma1-1));
e2 = (1-alpha)*p/((1-beta)*rho*(gamma2-1));
e2_rho2 = -(1-alpha)^2*p/((1-beta)^2*rho^2*(gamma2-1));
e2_p = (1-alpha)/((1-beta)*rho*(gamma2-1));

alpha_u = (-beta*(1-beta)*rhoL*uL/u^2*(beta/alpha*e1_rho1*e2_p-(1-beta)/(1-alpha)*e2_rho2*e1_p)/(beta*e1_p+(1-beta)*e2_p)+alpha*(uL-u)+(u-uL-pL/(rhoL*uL))*(beta*e1_p)/(beta*e1_p+(1-beta)*e2_p)+alpha*pL/(rhoL*uL))/...
(-beta*(1-beta)*rhoL*uL/u*(beta/alpha^2*e1_rho1*e2_p+(1-beta)/(1-alpha)^2*e2_rho2*e1_p)/(beta*e1_p+(1-beta)*e2_p)+u*(uL-u)+pL*u/(rhoL*uL));

alpha = alpha+alpha_u*du;

beta = betaL;

rho = rhoL*uL/u;

e = eL+1/2*(uL^2+u^2)-uL*u+pL/(rhoL*uL)*(uL-u);

p = rho*e/(alpha/(gamma1-1)+(1-alpha)/(gamma2-1));

tau = rhoL*uL*(u-uL)+p-pL;

A = [A; rho u p alpha beta tau];

end

plot(A(:,2),A(:,1))
hold on
grid on
plot(A(:,2),A(:,3))
plot(A(:,2),A(:,4))
plot(A(:,2),A(:,5))
plot(A(:,2),A(:,6))

legend('\rho','p','\alpha','\beta','\tau')

t0 = find(A(:,6)>-1e-3 & A(:,6)<1e-3)

plot([A(t0(end),2) A(t0(end),2)],[min(min(A)) max(max(A))],'y')

mu = 1.5e-3;

dxdu(1) = 0;
dxdu = [dxdu;mu./A(2:end,6)];
dx=dxdu*du;

x(1)=0;
for i=1:length(dxdu)-1
x(i+1)=x(i)+(dxdu(i)+dxdu(i+1))*du;
end
figure(2)
hold on
grid on
plot(x(1:t0(end)),A(1:t0(end),1))
plot(x(1:t0(end)),A(1:t0(end),2))
plot(x(1:t0(end)),A(1:t0(end),3))
plot(x(1:t0(end)),A(1:t0(end),4))
plot(x(1:t0(end)),A(1:t0(end),5))
plot(x(1:t0(end)),A(1:t0(end),6))