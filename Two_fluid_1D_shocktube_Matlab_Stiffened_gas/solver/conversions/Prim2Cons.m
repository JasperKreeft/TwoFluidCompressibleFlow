function q = Prim2Cons(w,gamma1,gamma2,pi1,pi2)

% Convert from primitive (w) to conservative (q) variables
q(1,:) = w(1,:);
q(2,:) = w(1,:).*w(2,:);
q(3,:) = (w(3,:)+gamma1*pi1).*w(5,:)/(gamma1-1)+(w(3,:)+gamma2*pi2).*(1-w(5,:))/(gamma2-1)+0.5*w(1,:).*w(2,:).^2;
q(4,:) = w(1,:).*w(4,:);
q(5,:) = (w(3,:)+gamma1*pi1).*w(5,:)/(gamma1-1)+0.5*w(1,:).*w(2,:).^2.*w(4,:);