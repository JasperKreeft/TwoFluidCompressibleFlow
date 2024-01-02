function w = Cons2Prim(q,gamma)

% Convert from conservative (q) to primative variables
w(1,:) = q(1,:);
w(2,:) = q(2,:)./q(1,:);
w(3,:) = (gamma-1)*(q(3,:)-(q(2,:).^2)./(2*q(1,:)));