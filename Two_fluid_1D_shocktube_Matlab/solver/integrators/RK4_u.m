function [dp,drho,dalpha]=RK4_u(rho,p,alpha,du,lr)

global gamma1 gamma2
global pi1 pi2

a      = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
b      = (gamma1*(p+pi1)-gamma2*(p+pi2))/(gamma1*(p+pi1)/alpha+gamma2*(p+pi2)/(1-alpha));
dp     = lr*du*rho*a;
drho   = lr*du*rho/a;
dalpha = lr*du*b/a;