function RIf = limitersRiem(N,ri,limiter)

eps = 1e-15;

dRIL=zeros(3,N);
dRIR=zeros(3,N);

if limiter == 0
    % nothing
    
else
    for j=2:N-1
        for k=1:3
            r = (ri(k,j+1)-ri(k,j)+eps)/(ri(k,j)-ri(k,j-1)+eps);

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
            
            dRIL(k,j) = 1/2*phiL*(ri(k,j)-ri(k,j-1));
            dRIR(k,j) = -1/2*r*phiR*(ri(k,j)-ri(k,j-1));
            
        end
    end

    %Boundary conditions
    dRIR(:,1) = 0; % 1/2*ri(:,1)-1/2*ri(:,2);
    dRIL(:,1) = -1/2*ri(:,1)+1/2*ri(:,2);
    dRIR(:,N) =  1/2*ri(:,N-1)-1/2*ri(:,N);
    dRIL(:,N) = 0; %-1/2*ri(:,N-1)+1/2*ri(:,N);
    
end

RIf(:,2*(1:N)-1) = ri + dRIR;
RIf(:,2*(1:N))   = ri + dRIL;