%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
% Written by:    ir. Jasper Kreeft  (2007)     %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sint]=Sint_lin(rho,u,p,beta,alpha,u_end,gamma1,gamma2,pi1,pi2,lr)

% lr is a sign function for left running wave (-) or right running wave (+)

du   = (u_end - u);
a    = sqrt(1/(alpha/(gamma1*(p+pi1))+(1-alpha)/(gamma2*(p+pi2)))/rho);
bS   = alpha*(1-alpha)*(gamma1*(p+pi1)-gamma2*(p+pi2))/((1-alpha)*gamma1*(p+pi1)+alpha*gamma2*(p+pi2));
% a    = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
% bS   = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
sint = (bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u)*du;