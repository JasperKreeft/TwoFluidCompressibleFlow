function Wf = limitersPrim(N,w,limiter)

epss = 1e-15;

dWL=zeros(3,N);
dWR=zeros(3,N);

if limiter == 0
    % nothing
    
else
    dw = diff(w,1,2);
    
    %Boundary conditions
    dWL(:,1) =  dw(:,1)/2;
    dWL(:,N) =  0; %  dw(:,N-1)/2; % 
    dWR(:,1) =  0; % -dw(:,1)/2;   % 
    dWR(:,N) = -dw(:,N-1)/2;

    % Domain
    rL = (dw(:,2:N-1)+epss)./(dw(:,1:N-2)+epss);
    rR = 1./rL;

        % Minmod limiter
        if limiter == 1
            phiL = max(0,min(1,rL));
            phiR = max(0,min(1,rR));

        % Koren limiter
        elseif limiter == 2
            phiL = max(0,min(2*rL,min(1/3+2/3*rL,2)));
            phiR = max(0,min(2*rR,min(1/3+2/3*rR,2)));

        % Van Leer limiter
        elseif limiter == 3
            phiL = (rL+abs(rL))./(1+abs(rL));
            phiR = (rR+abs(rR))./(1+abs(rR));

        % Van Albada limiter
        elseif limiter == 4
            phiL = (rL.^2+rL)./(rL.^2+1);
            phiR = (rR.^2+rR)./(rR.^2+1);

        % Superbee limiter
        elseif limiter == 5
            phiL = max([0,min(2*rL,1),min(rL,2)]);
            phiR = max([0,min(2*rR,1),min(rR,2)]);
        end

        dWL(:,2:N-1) =  1/2*phiL.*dw(:,1:N-2);
        dWR(:,2:N-1) = -1/2*phiR.*dw(:,2:N-1);
    
end

Wf(:,2*(1:N)-1) = w + dWR;
Wf(:,2*(1:N))   = w + dWL;