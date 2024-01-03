function [W,Wfx,Wfy]=limiterPrim(dx,dy,Q)

%--------------------------------------------------%
%                                                  %
% Limiters, used to avoid non-monotonicity         %
%                                                  %
%                                                  %
% Primitive variant                                %
%--------------------------------------------------%

global NX NY

lim = 2;

epss = rand(1)*100*eps;

Wfx  = zeros(NY,2*NX,8);
Wfy  = zeros(2*NY,NX,8);
dWLx = zeros(NY,NX,8);
dWRx = zeros(NY,NX,8);
dWAy = zeros(NY,NX,8);
dWBy = zeros(NY,NX,8);

W = Cons2Prim(Q);

% x-direction
% First-Order

if lim == 0
    % nothing
else
    dw =diff(W,1,2);
    dWLx(:,1,3:8)  = dw(:,1,3:8)/2;
    dWLx(:,NX,3:8) = 0;
    dWRx(:,1,3:8)  = 0;
    dWRx(:,NX,3:8) = -dw(:,NX-1,3:8)/2;
    rL = (dw(:,2:NX-1,3:8)+epss)./(dw(:,1:NX-2,3:8)+epss);
    rR = 1./rL;
    if lim==1
        % Minmod limiter
        phiL = (rL>=0).*(rL+(rL>=1).*(1-rL));
        phiR = (rR>=0).*(rR+(rR>=1).*(1-rR));
    elseif lim==2
        % Koren limiter
        phiL = (rL>=0).*( 2*rL + (1/3-4/3*rL).*(rL>=0.25) + (5/3-2/3*rL).*(rL>=2.5) );
        phiR = (rR>=0).*( 2*rR + (1/3-4/3*rR).*(rR>=0.25) + (5/3-2/3*rR).*(rR>=2.5) );
    end
    dWLx(:,2:NX-1,3:8) =  1/2*phiL.*dw(:,1:NX-2,3:8);
    dWRx(:,2:NX-1,3:8) = -1/2*phiR.*dw(:,2:NX-1,3:8);
end    
    
%y-direction
if lim==0
    % nothing
else
    dw =diff(W,1,1);
    dWAy(1,:,3:8)  = dw(1,:,3:8)/2;
    dWAy(NY,:,3:8) = 0;
    dWBy(1,:,3:8)  = 0;
    dWBy(NY,:,3:8) = -dw(NY-1,:,3:8)/2;
    rA = (dw(2:NY-1,:,3:8)+epss)./(dw(1:NY-2,:,3:8)+epss);
    rB = 1./rA;
    if lim==1
        % Minmod limiter
        phiA = (rA>=0).*(rA+(rA>=1).*(1-rA));
        phiB = (rB>=0).*(rB+(rB>=1).*(1-rB));
    elseif lim==2
        % Koren limiter
        phiA = (rA>=0).*( 2*rA + (1/3-4/3*rA).*(rA>=0.25) + (5/3-2/3*rA).*(rA>=2.5) );
        phiB = (rB>=0).*( 2*rB + (1/3-4/3*rB).*(rB>=0.25) + (5/3-2/3*rB).*(rB>=2.5) );
    end
    dWAy(2:NY-1,:,3:8) =  1/2*phiA.*dw(1:NY-2,:,3:8);
    dWBy(2:NY-1,:,3:8) = -1/2*phiB.*dw(2:NY-1,:,3:8);
end
    
    
Wfx(:,2*(1:NX)-1,1)   = W(:,:,1)-0.5*dx;
Wfx(:,2*(1:NX)-1,2)   = W(:,:,2);
Wfx(:,2*(1:NX),1)     = W(:,:,1)+0.5*dx;
Wfx(:,2*(1:NX),1)     = W(:,:,2);
Wfy(2*(1:NY)-1,:,1)   = W(:,:,1);
Wfy(2*(1:NY)-1,:,2)   = W(:,:,2)-0.5*dy;
Wfy(2*(1:NY),:,1)     = W(:,:,1);
Wfy(2*(1:NY),:,2)     = W(:,:,2)+0.5*dy;
for k = 1:6
    Wfx(:,2*(1:NX)-1,k+2) = W(:,:,k+2) + dWRx(:,:,k);
    Wfx(:,2*(1:NX),k+2)   = W(:,:,k+2) + dWLx(:,:,k);
    Wfy(2*(1:NY)-1,:,k+2) = W(:,:,k+2) + dWBy(:,:,k);
    Wfy(2*(1:NY),:,k+2)   = W(:,:,k+2) + dWAy(:,:,k);
end
