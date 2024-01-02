if ani==2
    figure(1)
    plot(x_cellcenter,W(1,:),'+-b','linewidth',2)
    grid on
    title('density')
    pause(0.1)
elseif ani==3
    figure(1)
    subplot(2,2,1)
    plot(x_cellcenter,W(1,:),'.-','linewidth',1)
    title('density')
    grid on
    xlim([0 1])
    if     nr==1; ylim([0 1.1]);
    elseif nr==2; ylim([0 1100]);
    elseif nr==3; ylim([0 1.1]);
    elseif nr==4; ylim([0 11]);
    elseif nr==5; ylim([0 1.1]);
    end
    subplot(2,2,2)
    plot(x_cellcenter,W(2,:),'.-','linewidth',1)
    grid on
    title('velocity')
    xlim([0 1])
    if     nr==1; ylim([0.9 1.1]);
    elseif nr==2; ylim([0.9 1.1]);
    elseif nr==3; ylim([0 1.1]);
    elseif nr==4; ylim([0 2]);
    elseif nr==5; ylim([-2.5 2.5]);
    end
    subplot(2,2,3)
    plot(x_cellcenter,W(3,:),'.-','linewidth',1)
    grid on
    title('pressure')
    xlim([0 1])
    if     nr==1; ylim([0.9 1.1]);
    elseif nr==2; ylim([0.9 1.1]);
    elseif nr==3; ylim([0 1.1]);
    elseif nr==4; ylim([0 11]);
    elseif nr==5; ylim([0 0.5]);
    end
    subplot(2,2,4)
    e = W(3,:)./((gamma-1)*W(1,:));
    plot(x_cellcenter,e,'.-','linewidth',1)
    title('internal energy')
    grid on
    xlim([0 1])
    if     nr==1; ylim([0 30]);
    elseif nr==2; ylim([0 3]);
    elseif nr==3; ylim([1.5 3]);
    elseif nr==4; ylim([0 5]);
    elseif nr==5; ylim([0 1.1]);
    end
    pause(0.05)
elseif ani==4
    xp = [ x_cellcenter ; x_cellcenter ];
    yp = y'*ones(1,N);
    zp = zeros(2,N);
    zp(1,:) = W(1,:);
    zp(2,:) = W(1,:);
    figure(1)
    pcolor(xp,yp,zp)
    shading interp
    axis([ 0 L y])
    axis equal
    axis tight
    axis off
    cbar = colorbar('horiz');
    cbar.Label.String = 'Density';
    set(gca,'clim',[0 round(max(W(1,:)))]);
    colormap jet % cool
    title(['Time = ' num2str(t(n),'%4.3f') ' sec'])
    pause(0.01)
end