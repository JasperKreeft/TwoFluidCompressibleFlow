function w = Cons2Prim(q,gamma1,gamma2,pi1,pi2)

% global t N

% Convert from conservative (q) to primative variables
w(1,:) = q(1,:);
w(2,:) = q(2,:)./q(1,:);
w(4,:) = q(4,:)./q(1,:);

% w(3,:) = (gamma1-gamma2)*(q(5,:)-((q(4,:).*q(2,:).^2)./(2*q(1,:).^2)))+(gamma2-1)*(q(3,:)-(q(2,:).^2./(2*q(1,:))));
% w(5,:) = ((gamma1-1)*(q(5,:)-q(4,:).*q(2,:).^2./(2*q(1,:).^2)))./((gamma1-gamma2)*(q(5,:)-q(4,:).*q(2,:).^2./(2*q(1,:).^2))+(gamma2-1)*(q(3,:)-q(2,:).^2./(2*q(1,:))));

for i=1:length(w(1,:))

    Q1 = q(3,i)-q(2,i).^2./(2*q(1,i));
    Q2 = q(5,i)-q(4,i).*q(2,i).^2./(2*q(1,i).^2);

% a=(gamma2*pi2-gamma1*pi1);
% b=(gamma1-1)*Q2+(gamma2-1)*(Q1-Q2)+(gamma1*pi1-gamma2*pi2);
% c=-(gamma1-1)*Q2;
% 
% w(5,i) = (-b+sqrt(b^2-4*a*c))/(2*a);

    R = roots([(gamma2*pi2-gamma1*pi1) (gamma1-1)*Q2+(gamma2-1)*(Q1-Q2)+(gamma1*pi1-gamma2*pi2) -(gamma1-1)*Q2]);

    w(5,i) = max((R>=0).*(R<=1+10*eps).*R);
    
    w(3,i) = (gamma1-gamma2)*Q2+(gamma2-1)*Q1-w(5,i)*gamma1*pi1-(1-w(5,i))*gamma2*pi2;

% if abs(w(5,i))==0 && i<104
%     keyboard
% end

end


