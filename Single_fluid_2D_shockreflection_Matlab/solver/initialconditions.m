function [W,gamma]=initialconditions(NX,NY)

gamma = 1.4;

W = zeros(NY,NX,4);

W(:,:,1) = 1;
W(:,:,2) = 2.9;
W(:,:,3) = 0;
W(:,:,4) = 1/gamma;