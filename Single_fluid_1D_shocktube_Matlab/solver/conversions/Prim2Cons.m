function q = Prim2Cons(w,gamma)


%Convert from primary (w) to conservation (q) variables
q(1,:) = w(1,:);
q(2,:) = w(1,:).*w(2,:);
q(3,:) = w(3,:)/(gamma-1)+1/2*w(1,:).*w(2,:).^2;