function q = Prim2Cons(w,gamma)

%Convert from primary (w) to conservation (q)
q(:,:,1) = w(:,:,1);
q(:,:,2) = w(:,:,1).*w(:,:,2);
q(:,:,3) = w(:,:,1).*w(:,:,3);
q(:,:,4) = w(:,:,4)/(gamma-1)+1/2*w(:,:,1).*(w(:,:,2).^2+w(:,:,3).^2);