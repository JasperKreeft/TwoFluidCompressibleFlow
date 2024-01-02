function Qf = limitersCons(N,q,limiter)

eps = 1e-15;

dQL=zeros(3,N);
dQR=zeros(3,N);

if limiter == 0
    % nothing
    
else
    for j=2:N-1
        for k=1:3
            r = (q(k,j+1)-q(k,j)+eps)/(q(k,j)-q(k,j-1)+eps);

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
            
            dQL(k,j) = 1/2*phiL*(q(k,j)-q(k,j-1));
            dQR(k,j) = -1/2*r*phiR*(q(k,j)-q(k,j-1));
            
        end
    end

    %Boundary conditions
    dQR(:,1) =  1/2*q(:,1)-1/2*q(:,2);
    dQL(:,1) = -1/2*q(:,1)+1/2*q(:,2);
    dQR(:,N) =  1/2*q(:,N-1)-1/2*q(:,N);
    dQL(:,N) = -1/2*q(:,N-1)+1/2*q(:,N);
    
end

Qf(:,2*(1:N)-1) = q + dQR;
Qf(:,2*(1:N))   = q + dQL;