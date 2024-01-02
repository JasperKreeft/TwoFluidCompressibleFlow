function w = Cons2Prim(q,gamma)


%Convert from conservation (Q) to primary (W) variables
w(:,:,1) = q(:,:,1);
w(:,:,2) = q(:,:,2)./q(:,:,1);
w(:,:,3) = q(:,:,3)./q(:,:,1);
w(:,:,4) = (gamma-1)*(q(:,:,4)-(q(:,:,2).^2+q(:,:,3).^2)./(2*q(:,:,1)));