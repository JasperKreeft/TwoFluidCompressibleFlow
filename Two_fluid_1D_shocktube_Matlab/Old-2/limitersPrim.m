%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                              %
% Written by: Jasper Kreeft  (2007)            %
% Updated by: Jasper Kreeft  (2015)            %
%                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Wf = limitersPrim(w,limiter)

global N

epss = 1e-15;

dWL=zeros(5,N);
dWR=zeros(5,N);

if limiter == 0
    % nothing
    
else
    
    dw = diff(w,1,2);

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


%         % Limited downwind limiter (anti-diffusive)
%         phiL_LD = max(0,min(2*rL([4 5],:),2));
%         phiR_LD = max(0,min(2*rR([4 5],:),2));
%         
%         dWL([4 5],2:N-1) =  1/2*phiL_LD.*dw([4 5],1:N-2);
%         dWR([4 5],2:N-1) = -1/2*phiR_LD.*dw([4 5],2:N-1);

    %Boundary conditions
    dWL(:,1) = -dw(:,1)/2;
    dWL(:,N) = dw(:,N-1)/2; 
    dWR(:,1) = -dw(:,1)/2;
    dWR(:,N) = dw(:,N-1)/2;
    
%     dWR(:,1) =  1/2*w(:,1)-1/2*w(:,2);
%     dWL(:,1) = -1/2*w(:,1)+1/2*w(:,2);
%     dWR(:,N) =  1/2*w(:,N-1)-1/2*w(:,N);
%     dWL(:,N) = -1/2*w(:,N-1)+1/2*w(:,N);


end

Wf(:,2*(1:N)-1) = w + dWR;
Wf(:,2*(1:N))   = w + dWL;
