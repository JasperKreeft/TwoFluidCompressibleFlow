function ri = Cons2Riem(q,gamma)

p      = (gamma-1)*(q(3,:)-(q(2,:).^2)./(2*q(1,:)));

% Convert from conservative (q) to Riemann invariant variables
ri(1,:) = q(2,:)./q(1,:);
ri(2,:) = sqrt(gamma*p./q(1,:));
ri(3,:) = log(p./(q(1,:).^gamma));