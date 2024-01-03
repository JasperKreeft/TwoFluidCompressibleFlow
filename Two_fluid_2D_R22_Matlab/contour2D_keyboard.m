close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Fortran file                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W = WRK3;

x     = W(:,:,1);
y     = W(:,:,2);
rho   = W(:,:,3);
u     = W(:,:,4);
v     = W(:,:,5);
p     = W(:,:,6);
beta  = W(:,:,7);
alpha = W(:,:,8);


NX = 200;
NY = 100;

% Velocity magnitude
W(:,:,9)=sqrt(W(:,:,4).^2+W(:,:,5).^2);
figure(2)
contourf(W(:,:,1),W(:,:,2),W(:,:,9),50)
% % Vorticity
% for i=2:NY-1
% for j=2:NX-1
% OOO(i,j,10)=(OOO(i,j+1,5)-OOO(i,j-1,5))/(OOO(i,j+1,1)-OOO(i,j-1,1))-(OOO(i+1,j,4)-OOO(i-1,j,4))/(OOO(i+1,j,2)-OOO(i-1,j,2));
% end
% end
% 
% figure(2)
% contour(OOO(:,:,1),OOO(:,:,2),OOO(:,:,10),8)

% figure('Position',[1,1024,1280,1024])

% figure(1)
% subplot(3,2,1)
% pcolor(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,3))
% shading interp
% subplot(3,2,2)
% pcolor(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,4))
% shading interp
% subplot(3,2,3)
% pcolor(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,5))
% shading interp
% subplot(3,2,4)
% pcolor(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,6))
% shading interp
% subplot(3,2,5)
% pcolor(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,7))
% shading interp
% subplot(3,2,6)
% pcolor(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,8))
% shading interp

% figure(1)
% subplot(3,2,1)
% contour(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,3),[1.42:0.05:2.0 4.4:0.2:7])
% subplot(3,2,2)
% contour(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,4),50)
% subplot(3,2,3)
% contour(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,5),50)
% subplot(3,2,4)
% contour(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,6),50)
% subplot(3,2,5)
% contour(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,7),10)
% subplot(3,2,6)
% contour(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,8),10)

% figure(1)
% subplot(3,2,1)
% contourf(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,3),100)
% subplot(3,2,2)
% contourf(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,4),50)
% subplot(3,2,3)
% contourf(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,5),50)
% subplot(3,2,4)
% contourf(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,6),50)
% subplot(3,2,5)
% contourf(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,7),10)
% subplot(3,2,6)
% contourf(OOO(:,201:600,1),OOO(:,201:600,2),OOO(:,201:600,8),10)

% figure(1)
% subplot(3,2,1)
% contour(OOO(:,:,1),OOO(:,:,2),OOO(:,:,3),50)
% subplot(3,2,2)
% contour(OOO(:,:,1),OOO(:,:,2),OOO(:,:,4),50)
% subplot(3,2,3)
% contour(OOO(:,:,1),OOO(:,:,2),OOO(:,:,5),50)
% subplot(3,2,4)
% contour(OOO(:,:,1),OOO(:,:,2),OOO(:,:,6),50)
% subplot(3,2,5)
% contour(OOO(:,:,1),OOO(:,:,2),OOO(:,:,7),10)
% subplot(3,2,6)
% contour(OOO(:,:,1),OOO(:,:,2),OOO(:,:,8),10)

figure(1)
subplot(3,2,1)
contourf(W(:,:,1),W(:,:,2),W(:,:,3),50)
subplot(3,2,2)
contourf(W(:,:,1),W(:,:,2),W(:,:,4),50)
subplot(3,2,3)
contourf(W(:,:,1),W(:,:,2),W(:,:,5),50)
subplot(3,2,4)
contourf(W(:,:,1),W(:,:,2),W(:,:,6),50)
subplot(3,2,5)
contourf(W(:,:,1),W(:,:,2),W(:,:,7),10)
subplot(3,2,6)
contourf(W(:,:,1),W(:,:,2),W(:,:,8),10)

subplot(3,2,1)
title('\rho')
colorbar; axis equal; axis tight
subplot(3,2,2)
title('u')
colorbar; axis equal; axis tight
subplot(3,2,3)
title('v')
colorbar; axis equal; axis tight
subplot(3,2,4)
title('p')
colorbar; axis equal; axis tight
subplot(3,2,5)
title('\beta')
colorbar; axis equal; axis tight
subplot(3,2,6)
title('\alpha')
colorbar; axis equal; axis tight


% figure(3)
% contourf(OOO(:,:,6),20)
% hold on
% surf(OOO(:,:,6))