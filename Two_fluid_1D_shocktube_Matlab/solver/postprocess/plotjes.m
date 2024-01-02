kleur = '.-r';

figure(10)
subplot(2,2,1)
plot(x,W(1,:,end),kleur,'linewidth',1)
title('mixture density')
grid on
hold on
xlim([-.25 .25])
subplot(2,2,2)
plot(x,W(2,:,end),kleur,'linewidth',1)%*sqrt(10^5)
grid on
hold on
title('velocity')
xlim([-.25 .25])
subplot(2,2,3)
plot(x,W(3,:,end),kleur,'linewidth',1)%*10^5
grid on
hold on
title('pressure')
xlim([-.25 .25])
subplot(2,2,4)
plot(x,W(5,:,end),kleur,'linewidth',1)
title('epoxy volume fraction')
grid on
hold on
xlim([-.25 .25])
% axis([0 1 0.46 0.62])
% set(gca,'ytick',0.46:0.02:0.62)