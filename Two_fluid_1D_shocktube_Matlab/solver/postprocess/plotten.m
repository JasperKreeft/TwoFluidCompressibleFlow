if ani>1
    W = Cons2Prim(Q);
end

if ani==2
    figure(1)
    plot(x,W(1,:),'+-b','linewidth',2)
    grid on
    title('density')
    xlim([-L/2 L/2])
elseif ani==3
    figure(1)
    subplot(2,4,1:2)
    plot(x,W(1,:),'.-','linewidth',1)
    title('density')
    grid on
    xlim([-L/2 L/2])
    if max(W(1,:))-min(W(1,:))<1e-2; xyl = [-L/2 L/2 mean(W(1,:))+1e-2*[-1 1]]; axis(xyl); end
    subplot(2,4,3:4)
    plot(x,W(2,:),'.-','linewidth',1)
    grid on
    title('velocity')
    xlim([-L/2 L/2])
    if max(W(2,:))-min(W(2,:))<1e-2; xyl = [-L/2 L/2 mean(W(2,:))+1e-2*[-1 1]]; axis(xyl); end
    subplot(2,4,5:6)
    plot(x,W(3,:),'.-','linewidth',1)
    grid on
    title('pressure')
    xlim([-L/2 L/2])
    if max(W(3,:))-min(W(3,:))<1e-2; xyl = [-L/2 L/2 mean(W(3,:))+1e-2*[-1 1]]; axis(xyl); end
    subplot(2,4,7)
    plot(x,W(4,:),'.-','linewidth',1)
    title('mass fraction')
    grid on
    axis([-L/2 L/2 -0.2 1.2])
    subplot(2,4,8)
    plot(x,W(5,:),'.-','linewidth',1)
    title('volume fraction')
    grid on
    axis([-L/2 L/2 -0.2 1.2])
elseif ani==4
    z = zeros(2,N);
    z(1,:) = W(1,:);
    z(2,:) = W(1,:);
    figure(1)
    pcolor(x,y,z)
    shading interp
    axis equal
    axis off
    colorbar('horiz')
    colormap(hot)
end


if movie
    set(gcf,'Visible','off')
    frame = getframe(gcf,[0 0 1920 1080]);
    writeVideo(writerObj,frame);
end

pause(0.01)
