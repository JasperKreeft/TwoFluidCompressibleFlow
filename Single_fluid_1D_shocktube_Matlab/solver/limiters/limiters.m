function Wf = limiters(N,w,limiter)

eps = 1e-15;

dWL=zeros(3,N);
dWR=zeros(3,N);

if limiter == 0
    % nothing
    
else
    for j=2:N-1
        for k=1:3
            r = (w(k,j+1)-w(k,j)+eps)/(w(k,j)-w(k,j-1)+eps);

            % Minmod limiter
            if limiter == 1
                phiL = max(0,min(1,r));
                phiR = max(0,min(1,1/r));
                
            % Koren limiter
            elseif limiter == 2
                phiL = max(0,min(2*r,min(1/3+2/3*r,2)));
                phiR = max(0,min(2/r,min(1/3+2/3/r,2)));
                
            % Van Leer limiter
            elseif limiter == 3
                phiL = (r+abs(r))/(1+abs(r));
                phiR = (1/r+abs(1/r))/(1+abs(1/r));
                
            % Van Albada limiter
            elseif limiter == 4
                phiL = (r^2+r)/(r^2+1);
                phiR = ((1/r)^2+1/r)/((1/r)^2+1);

            % Superbee limiter
            elseif limiter == 5
                phiL = max([0,min(2*r,1),min(r,2)]);
                phiR = max([0,min(2/r,1),min(1/r,2)]);
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