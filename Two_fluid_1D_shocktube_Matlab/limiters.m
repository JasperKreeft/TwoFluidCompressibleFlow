function Wf = limiters(N,w,limiter)

eps = 1e-12;

dWL=zeros(5,N);
dWR=zeros(5,N);

if limiter == 0
    % nothing
    
else

    for j=2:N-1
        for k=1:5
            r = (w(k,j+1)-w(k,j)+eps)/(w(k,j)-w(k,j-1)+eps);
            
            % Minmod limiter
            if limiter == 1
                phiL = max(0,min(1,r));
                phiR = max(0,min(1,1/r));
            
            % Koren limiter
            elseif limiter == 2
                phiL = max(0,min(2*r,min(1/3+2/3*r,2)));
                phiR = max(0,min(2/r,min(1/3+2/3/r,2)));
            end

            dWL(k,j) = 1/2*phiL*(w(k,j)-w(k,j-1));
            dWR(k,j) = -1/2*r*phiR*(w(k,j)-w(k,j-1));
        end
    end

    %Boundary conditions
    dWR(:,1) =  1/2*w(:,1)-1/2*w(:,2);
    dWL(:,1) = -1/2*w(:,1)+1/2*w(:,2);
    dWR(:,N) =  1/2*w(:,N-1)-1/2*w(:,N);
    dWL(:,N) = -1/2*w(:,N-1)+1/2*w(:,N);

end

Wf(:,2*(1:N)-1) = w + dWR;
Wf(:,2*(1:N))   = w + dWL;