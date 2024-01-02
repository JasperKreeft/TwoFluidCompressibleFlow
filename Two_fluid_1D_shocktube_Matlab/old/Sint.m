function sint=Sint(rho,u,p,beta,alpha,u_end,lr)

global gamma1 gamma2
global pi1 pi2

% lr is a sign function for left running wave (-) or right running wave (+)

s = [0 0 0 0];

du   = (u_end - u);

a  = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
bS = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
s(1) = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;

for k=1:3
[dp_k1,drho_k1,dalpha_k1] = RK4_u(rho          ,p        ,alpha            ,du/3,lr);
[dp_k2,drho_k2,dalpha_k2] = RK4_u(rho+drho_k1/2,p+dp_k1/2,alpha+dalpha_k1/2,du/3,lr);
[dp_k3,drho_k3,dalpha_k3] = RK4_u(rho+drho_k2/2,p+dp_k2/2,alpha+dalpha_k2/2,du/3,lr);
[dp_k4,drho_k4,dalpha_k4] = RK4_u(rho+drho_k3  ,p+dp_k3  ,alpha+dalpha_k3  ,du/3,lr);
p      = p     + (dp_k1     + 2*dp_k2     + 2*dp_k3     + dp_k4    )/6;
rho    = rho   + (drho_k1   + 2*drho_k2   + 2*drho_k3   + drho_k4  )/6;
alpha  = alpha + (dalpha_k1 + 2*dalpha_k2 + 2*dalpha_k3 + dalpha_k4)/6;
u      = u+du/3;
a      = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
bS     = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
s(k+1) = bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u;
end

sint = du/8*(s(1)+3*s(2)+3*s(3)+s(4));
% Simpson's Rule (third order)