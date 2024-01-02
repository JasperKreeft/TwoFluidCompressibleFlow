function [e,a] = EOS(rho,p,gamma)


e = p/(gamma-1)/rho;
a = sqrt(gamma*p/rho);

