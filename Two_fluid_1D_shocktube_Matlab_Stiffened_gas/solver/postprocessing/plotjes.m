kleur = 'ob';
set(0, 'DefaultLineMarkerSize', 6)

figure%(10)
% subplot(2,2,1)
plot(x,W(1,:),kleur,'linewidth',1)
set(gca,'box','off')
% title('mixture density')
% grid on
hold on
xlim([0 1])
figure% subplot(2,2,2)
plot(x,W(2,:)*sqrt(1e5),kleur,'linewidth',1)%
set(gca,'box','off')
% grid on
hold on
% title('velocity')
figure% subplot(2,2,3)
plot(x,W(3,:)*1e5,kleur,'linewidth',1)%
set(gca,'box','off')
% grid on
hold on
% title('pressure')
xlim([0 1])
figure% subplot(2,2,4)
plot(x,W(5,:),kleur,'linewidth',1)
set(gca,'box','off')
% title('epoxy volume fraction')
% grid on
hold on
% axis([0 1 0.46 0.62])
% set(gca,'ytick',0.46:0.02:0.62)


figure
plot(x,W(4,:),kleur,'linewidth',1)
set(gca,'box','off')
hold on

rE = W(5,:).*W(3,:)*1e5/(gamma1-1)+(1-W(5,:)).*(W(3,:)*1e5+gamma2*pi2*1e5)/(gamma2-1)+0.5*W(1,:).*W(2,:).^2*1e5;
rE = sum(rE)

figure
plot(x,rE,kleur,'linewidth',1)
set(gca,'box','off')
hold on

set(0, 'DefaultLineMarkerSize', 9)

rE0_l = W(5,1).*W(3,1)*1e5/(gamma1-1)+(1-W(5,1)).*(W(3,1)*1e5+gamma2*pi2*1e5)/(gamma2-1)+0.5*W(1,1)*W(2,1).^2*1e5;
rE0_r = W(5,end).*W(3,end)*1e5/(gamma1-1)+(1-W(5,end)).*(W(3,end)*1e5+gamma2*pi2*1e5)/(gamma2-1)+0.5*W(1,end)*W(2,end).^2*1e5;
n = length(x);
rE0 = (rE0_l+rE0_r)/2*n