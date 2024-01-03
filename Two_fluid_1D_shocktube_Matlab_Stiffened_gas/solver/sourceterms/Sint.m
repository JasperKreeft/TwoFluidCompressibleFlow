function [sint]=Sint(rho,u,p,beta,alpha,u_end,gamma1,gamma2,pi1,pi2,lr)

% lr is a sign function for left running wave (-) or right running wave (+)

du   = (u_end - u);

a     = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
bS    = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
s0 = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;

[dp_k1,drho_k1,dalpha_k1] = RK4_p(rho,p,alpha,        0,      0,          0,gamma1,gamma2,pi1,pi2,lr*du/3);
[dp_k2,drho_k2,dalpha_k2] = RK4_p(rho,p,alpha,  drho_k1,  dp_k1,  dalpha_k1,gamma1,gamma2,pi1,pi2,lr*du/3);
[dp_k3,drho_k3,dalpha_k3] = RK4_p(rho,p,alpha,  drho_k2,  dp_k2,  dalpha_k2,gamma1,gamma2,pi1,pi2,lr*du/3);
[dp_k4,drho_k4,dalpha_k4] = RK4_p(rho,p,alpha,2*drho_k3,2*dp_k3,2*dalpha_k3,gamma1,gamma2,pi1,pi2,lr*du/3);
p     = p + (dp_k1 + 2*dp_k2 + 2*dp_k3 + dp_k4)/6;
rho   = rho + (drho_k1 + 2*drho_k2 + 2*drho_k3 + drho_k4)/6;
alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
u     = u+du/3;
a     = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
bS    = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
s1    = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;

[dp_k1,drho_k1,dalpha_k1] = RK4_p(rho,p,alpha,        0,      0,          0,gamma1,gamma2,pi1,pi2,lr*du/3);
[dp_k2,drho_k2,dalpha_k2] = RK4_p(rho,p,alpha,  drho_k1,  dp_k1,  dalpha_k1,gamma1,gamma2,pi1,pi2,lr*du/3);
[dp_k3,drho_k3,dalpha_k3] = RK4_p(rho,p,alpha,  drho_k2,  dp_k2,  dalpha_k2,gamma1,gamma2,pi1,pi2,lr*du/3);
[dp_k4,drho_k4,dalpha_k4] = RK4_p(rho,p,alpha,2*drho_k3,2*dp_k3,2*dalpha_k3,gamma1,gamma2,pi1,pi2,lr*du/3);
p     = p + (dp_k1 + 2*dp_k2 + 2*dp_k3 + dp_k4)/6;
rho   = rho + (drho_k1 + 2*drho_k2 + 2*drho_k3 + drho_k4)/6;
alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
u     = u+du/3;
a     = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
bS    = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
s2    = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;

[dp_k1,drho_k1,dalpha_k1] = RK4_p(rho,p,alpha,        0,      0,          0,gamma1,gamma2,pi1,pi2,lr*du/3);
[dp_k2,drho_k2,dalpha_k2] = RK4_p(rho,p,alpha,  drho_k1,  dp_k1,  dalpha_k1,gamma1,gamma2,pi1,pi2,lr*du/3);
[dp_k3,drho_k3,dalpha_k3] = RK4_p(rho,p,alpha,  drho_k2,  dp_k2,  dalpha_k2,gamma1,gamma2,pi1,pi2,lr*du/3);
[dp_k4,drho_k4,dalpha_k4] = RK4_p(rho,p,alpha,2*drho_k3,2*dp_k3,2*dalpha_k3,gamma1,gamma2,pi1,pi2,lr*du/3);
p     = p + (dp_k1 + 2*dp_k2 + 2*dp_k3 + dp_k4)/6;
rho   = rho + (drho_k1 + 2*drho_k2 + 2*drho_k3 + drho_k4)/6;
alpha = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
u     = u+du/3;
a     = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
bS    = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
s3    = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;

sint = du/8*(s0+3*s1+3*s2+s3);
% Simpson's Rule (third order)


% sint = sint;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% du   = (u_end - u);
% a    = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
% bS   = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
% sint = (bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u)*du;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simpson linear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 
% du   = (u_end - u)/3;
% 
% a  = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
% bS = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
% s0 = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;
% 
% u     = u+du;
% p     = p+lr*rho*a*du;
% rho   = rho+lr*rho/a*du;
% alpha = alpha+lr*bS/a*du;
% a     = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
% bS    = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
% s1    = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;
% 
% u     = u+du;
% p     = p+lr*rho*a*du;
% rho   = rho+lr*rho/a*du;
% alpha = alpha+lr*bS/a*du;
% a     = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
% bS    = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
% s2    = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;
% 
% u     = u+du;
% p     = p+lr*rho*a*du;
% rho   = rho+lr*rho/a*du;
% alpha = alpha+lr*bS/a*du;
% a     = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
% bS    = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
% s3    = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;
% 
% 
% sint = (3*du)/8*(s0+3*s1+3*s2+s3);
% % Simpson's Rule (third order)
% 
% 
% sint = sint;