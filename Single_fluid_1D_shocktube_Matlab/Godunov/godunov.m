clear all
%close all
clc

T_final = 0.195;
t = 0;

godunovinitials

while t < T_final
t = t + dt
    
godunovmethod;
fluxstukgodunov;

for k=2:N+1
q(:,k) = q(:,k) + dt/dx*(f(:,k-1)-f(:,k));
end
%boundary conditions
%ONLY SOFT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
q(:,1) = q(:,2);
q(:,N+2) = q(:,N+1);

plot(w(3,:),'-o')
pause(0.01)
end

% figure(2)
