function [dp_new,drho_new,dalpha_new]=RK4_p(rho,p,alpha,drho_old,dp_old,dalpha_old,gamma1,gamma2,du)

p          = p+dp_old/2;
alpha      = alpha+dalpha_old/2;
rho        = rho+drho_old/2;
a          = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
b          = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
dp_new     = du*rho*a;
drho_new   = du*rho/a;
dalpha_new = du*b/a;