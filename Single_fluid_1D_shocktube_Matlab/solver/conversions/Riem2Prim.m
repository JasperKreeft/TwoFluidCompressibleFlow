function w = Riem2Prim(ri,gamma)

% Convert from conservative (q) to primative variables
w(1,:) = (1/gamma*ri(2,:).^2.*exp(-ri(3,:))).^(1/(gamma-1));
w(2,:) = ri(1,:);
w(3,:) = (ri(2,:).^2/gamma.*exp(-ri(3,:)/gamma)).^(gamma/(gamma-1));