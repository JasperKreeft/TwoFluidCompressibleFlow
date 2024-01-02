function [sint]=Sint(rho,u,p,beta,alpha,u_end,gamma1,gamma2,lr)

% lr is a sign function for left running wave (-) or right running wave (+)

% NRK4 = 10;
% sint = 0;
% du   = (u_end - u)/NRK4;
% 
% for k=1:NRK4
%     [dp_k1,drho_k1,dalpha_k1] = RK4_p(rho,p,alpha,0,0,0,gamma1,gamma2,lr*du);
%     [dp_k2,drho_k2,dalpha_k2] = RK4_p(rho,p,alpha,drho_k1,dp_k1,dalpha_k1,gamma1,gamma2,lr*du);
%     [dp_k3,drho_k3,dalpha_k3] = RK4_p(rho,p,alpha,drho_k2,dp_k2,dalpha_k2,gamma1,gamma2,lr*du);
%     [dp_k4,drho_k4,dalpha_k4] = RK4_p(rho,p,alpha,drho_k3,dp_k3,dalpha_k3,gamma1,gamma2,lr*du);
%     p     = p + (dp_k1 + 2*dp_k2 + 2*dp_k3 + dp_k4)/6;
%     rho   = rho + (drho_k1 + 2*drho_k2 + 2*drho_k3 + drho_k4)/6;
%     alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
%     u     = u+du;
%     a     = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
%     bS    = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
%     sint  = sint + (bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u)*du;
% end


du   = (u_end - u);

a  = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
bS = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
s0 = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;

[dp_k1,drho_k1,dalpha_k1] = RK4_p(rho,p,alpha,        0,      0,          0,gamma1,gamma2,lr*du/3);
[dp_k2,drho_k2,dalpha_k2] = RK4_p(rho,p,alpha,  drho_k1,  dp_k1,  dalpha_k1,gamma1,gamma2,lr*du/3);
[dp_k3,drho_k3,dalpha_k3] = RK4_p(rho,p,alpha,  drho_k2,  dp_k2,  dalpha_k2,gamma1,gamma2,lr*du/3);
[dp_k4,drho_k4,dalpha_k4] = RK4_p(rho,p,alpha,2*drho_k3,2*dp_k3,2*dalpha_k3,gamma1,gamma2,lr*du/3);
p     = p + (dp_k1 + 2*dp_k2 + 2*dp_k3 + dp_k4)/6;
rho   = rho + (drho_k1 + 2*drho_k2 + 2*drho_k3 + drho_k4)/6;
alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
u     = u+du/3;
a     = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
bS    = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
s1    = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;

[dp_k1,drho_k1,dalpha_k1] = RK4_p(rho,p,alpha,        0,      0,          0,gamma1,gamma2,lr*du/3);
[dp_k2,drho_k2,dalpha_k2] = RK4_p(rho,p,alpha,  drho_k1,  dp_k1,  dalpha_k1,gamma1,gamma2,lr*du/3);
[dp_k3,drho_k3,dalpha_k3] = RK4_p(rho,p,alpha,  drho_k2,  dp_k2,  dalpha_k2,gamma1,gamma2,lr*du/3);
[dp_k4,drho_k4,dalpha_k4] = RK4_p(rho,p,alpha,2*drho_k3,2*dp_k3,2*dalpha_k3,gamma1,gamma2,lr*du/3);
p     = p + (dp_k1 + 2*dp_k2 + 2*dp_k3 + dp_k4)/6;
rho   = rho + (drho_k1 + 2*drho_k2 + 2*drho_k3 + drho_k4)/6;
alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
u     = u+du/3;
a     = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
bS    = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
s2    = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;

[dp_k1,drho_k1,dalpha_k1] = RK4_p(rho,p,alpha,        0,      0,          0,gamma1,gamma2,lr*du/3);
[dp_k2,drho_k2,dalpha_k2] = RK4_p(rho,p,alpha,  drho_k1,  dp_k1,  dalpha_k1,gamma1,gamma2,lr*du/3);
[dp_k3,drho_k3,dalpha_k3] = RK4_p(rho,p,alpha,  drho_k2,  dp_k2,  dalpha_k2,gamma1,gamma2,lr*du/3);
[dp_k4,drho_k4,dalpha_k4] = RK4_p(rho,p,alpha,2*drho_k3,2*dp_k3,2*dalpha_k3,gamma1,gamma2,lr*du/3);
p     = p + (dp_k1 + 2*dp_k2 + 2*dp_k3 + dp_k4)/6;
rho   = rho + (drho_k1 + 2*drho_k2 + 2*drho_k3 + drho_k4)/6;
alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
u     = u+du/3;
a     = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
bS    = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
s3    = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;

sint = du/8*(s0+3*s1+3*s2+s3);
% Simpson's Rule (third order)