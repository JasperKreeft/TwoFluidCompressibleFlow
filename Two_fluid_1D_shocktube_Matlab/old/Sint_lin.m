%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
% Written by:    ir. Jasper Kreeft  (2007)     %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sint]=Sint_lin(rho,u,p,beta,alpha,u_end,lr)

global gamma1 gamma2

% lr is a sign function for left running wave (-) or right running wave (+)

du   = (u_end - u);
a    = sqrt(1/(alpha/gamma1+(1-alpha)/gamma2)*p/rho);
bS   = alpha*(1-alpha)*(gamma1-gamma2)/((1-alpha)*gamma1+alpha*gamma2);
sint = (bS*p*(1+lr*u/a)+lr*(alpha-beta)*rho*a*u)*du;