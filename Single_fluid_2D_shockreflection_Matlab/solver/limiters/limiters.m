function [Wfx,Wfy]  = limiters(NX,NY,w,limiter,eps)

eps = 1e-12;

Wfx = zeros(NY,2*NX,4);
Wfy = zeros(2*NY,NX,4);

dWLx = zeros(NY,NX,4);
dWRx = zeros(NY,NX,4);
dWBy = zeros(NY,NX,4);
dWAy = zeros(NY,NX,4);

if limiter == 0 % No limiter
	% nothing

else
    for i=2:NY-1
        for j=2:NX-1
            for k=1:4
                if j==1
                    dWLx(i,j,k) = 1/2*(w(i,j+1,k)-w(i,j,k));
                    dWRx(i,j,k) = 0;
                elseif j==NX
                    dWLx(i,j,k) = 0;
                    dWRx(i,j,k) = 1/2*(w(i,j-1,k)-w(i,j,k));
                else
                    r = (w(i,j+1,k)-w(i,j,k)+eps)/(w(i,j,k)-w(i,j-1,k)+eps);
                    if limiter == 1 % Minmod
                        phiL = max(0,min(1,r));
                        phiR = max(0,min(1,1/r));
                    elseif limiter == 2 % Koren
                        phiL = max(0,min(2*r,min(1/3+2/3*r,2)));
                        phiR = max(0,min(2/r,min(1/3+2/3/r,2)));
                    elseif limiter == 3 % Van Leer                        
                        phiL = (r+abs(r))/(1+r);
                        phiR = (1/r+abs(1/r))/(1+1/r);
                    elseif limiter == 4 % Van Albada
                        phiL = (r^2+r)/(r^2+1);
                        phiR = ((1/r)^2+1/r)/((1/r)^2+1);
                    elseif limiter == 5 % HQUICK                        
                        phiL = 2*(r+abs(r))/(r+3);
                        phiR = 2*((1/r)+abs(1/r))/((1/r)+3);
                    end 
                    dWLx(i,j,k) = 1/2*phiL*(w(i,j,k)-w(i,j-1,k));
                    dWRx(i,j,k) = -1/2*r*phiR*(w(i,j,k)-w(i,j-1,k));
                end
                if i==1
                    dWAy(i,j,k) = 1/2*(w(i+1,j,k)-w(i,j,k));
                    dWBy(i,j,k) = 0;
                elseif i==NY
                    dWAy(i,j,k) = 0;
                    dWBy(i,j,k) = 1/2*(w(i-1,j,k)-w(i,j,k));
                else
                    r = (w(i+1,j,k)-w(i,j,k)+eps)/(w(i,j,k)-w(i-1,j,k)+eps);
                    if limiter == 1 % Minmod
                        phiA = max(0,min(1,r));
                        phiB = max(0,min(1,1/r));
                    elseif limiter == 2 % Koren
                        phiA = max(0,min(2*r,min(1/3+2/3*r,2)));
                        phiB = max(0,min(2/r,min(1/3+2/3/r,2)));
                    elseif limiter == 3 % Van Leer
                        phiA = (r+abs(r))/(1+r);
                        phiB = (1/r+abs(1/r))/(1+1/r);
                    elseif limiter == 4 % Van Albada                        
                        phiA = (r^2+r)/(r^2+1);
                        phiB = ((1/r)^2+1/r)/((1/r)^2+1);  
                    elseif limiter == 5 % HQUICK
                        phiA = 2*(r+abs(r))/(r+3);
                        phiB = 2*((1/r)+abs(1/r))/((1/r)+3);
                    end
                    dWAy(i,j,k) = 1/2*phiA*(w(i,j,k)-w(i,j-1,k));
                    dWBy(i,j,k) = -1/2*r*phiB*(w(i,j,k)-w(i,j-1,k));
                end
            end
        end
    end
end


for k=1:4
    Wfx(:,2*(1:NX)-1,k) = w(:,:,k)+dWRx(:,:,k);
    Wfx(:,2*(1:NX),k)   = w(:,:,k)+dWLx(:,:,k);
    Wfy(2*(1:NY)-1,:,k) = w(:,:,k)+dWAy(:,:,k);
    Wfy(2*(1:NY),:,k)   = w(:,:,k)+dWBy(:,:,k);
end