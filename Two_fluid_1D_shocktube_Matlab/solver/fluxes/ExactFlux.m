function [F,sL,sR,Lambda] = ExactFlux(Wf,Lambda)

% Exact Riemann solver for two-fluid mixture is based on the work of Saurel

global N
global gamma1 gamma2
global pi1 pi2

epss = 1e-12;

F  = zeros(5,N+1);
sL = zeros(1,N); % values use cell-center index instead of flux index
sR = zeros(1,N); % values use cell-center index instead of flux index

% i = j - 1/2
for i  = 2:N

    rho1   = Wf(1,2*i-2);
    u1     = Wf(2,2*i-2);
    p1     = Wf(3,2*i-2);
    beta1  = Wf(4,2*i-2);
    alpha1 = Wf(5,2*i-2);
    a1     = sqrt(1/(alpha1/(gamma1*(p1+pi1))+(1-alpha1)/(gamma2*(p1+pi2)))/rho1);

    rho4   = Wf(1,2*i-1);
    u4     = Wf(2,2*i-1);
    p4     = Wf(3,2*i-1);
    beta4  = Wf(4,2*i-1);
    alpha4 = Wf(5,2*i-1);
    a4     = sqrt(1/(alpha4/(gamma1*(p4+pi1))+(1-alpha4)/(gamma2*(p4+pi2)))/rho4);
    
    % Check for negative densities and pressures
    if rho1<=0
        disp('rho1 negative')
        rho1 = epss;
    end

    if p1<=0
        disp('p1 negative')
        p1 = epss;
    end

    if rho4<=0
        disp('rho4 negative')
        rho4 = epss;
    end

    if p4<=0
        disp('p4 negative')
        p4 = epss;
    end

    % Equal state in both cells
    if abs(rho1-rho4)+abs(u1-u4)+abs(p1-p4)+abs(beta1-beta4)+abs(alpha1-alpha4)<1e-10
    
       F(:,i) = [rho1*u1;
                 rho1*u1^2+p1;
                 alpha1*gamma1*(p1+pi1)*u1/(gamma1-1)+(1-alpha1)*gamma2*(p1+pi2)*u1/(gamma2-1)+1/2*rho1*u1^3;
                 rho1*u1*beta1;
                 alpha1*gamma1*(p1+pi1)*u1/(gamma1-1)+1/2*rho1*u1^3*beta1];

        Lambda  = max(abs([u1 u1-a1 u4+a4 Lambda]));

    % non-moving interface
    elseif abs(p1-p4)+abs(u1)+abs(u4)<1e-10

        F(2,i) = p1;

    % Two-fluid problem
    else

        alpha1 = alpha1+epss;
        alpha4 = alpha4+epss;
        beta1 = beta1+epss;
        beta4 = beta4+epss;

        rho1_1 = (beta1)/(alpha1)*rho1;
        rho1_2 = (1-beta1)/(1-alpha1)*rho1;
%         [~,~,c1] = eoscmix(p1,rho1,rho1_1,alpha1);
        c1 = sqrt(1/(alpha1/(gamma1*(p1+pi1))+(1-alpha1)/(gamma2*(p1+pi2)))/rho1);
        z1 = rho1*c1;

        rho4_1 = (beta4)/(alpha4)*rho4;
        rho4_2 = (1-beta4)/(1-alpha4)*rho4;
%         [~,~,c4] = eoscmix(p4,rho4,rho4_1,alpha4);
        c4 = sqrt(1/(alpha4/(gamma1*(p4+pi1))+(1-alpha4)/(gamma2*(p4+pi2)))/rho4);
        z4 = rho4*c4;

        p23 = (z1*p4+z4*p1+z4*z1*(u1-u4))/(z1+z4);
        if p23<0
            p23=(p4+p1)/2;
        end

        iq=0;

        while 1

            %
            % Left state
            %
            if p23>p1
            % Shockwave left
                rho2_1 = rhohug(gamma1,pi1,p23,p1,rho1_1);
                rho2_2 = rhohug(gamma2,pi2,p23,p1,rho1_2);
                m1_1 = debitchoc(gamma1,pi1,p23,p1,rho1_1);
                m1_2 = debitchoc(gamma2,pi2,p23,p1,rho1_2);
                m1   = sqrt(1/(beta1/m1_1^2+(1-beta1)/m1_2^2));
                u2   = u1-(p23-p1)/m1;
                sl   = u1-m1/rho1;
            else
            % expansion left
                rho2_1 = rhoisent(gamma1,pi1,p23,p1,rho1_1);
                rho2_2 = rhoisent(gamma2,pi2,p23,p1,rho1_2);
                somme  = gausslg(p1,p23,rho1_1,rho1_2,beta1,(1-beta1));
                u2     = u1-somme;
                sl     = u1-c1;
            end


        %
        %       Right state
        %
            if p23>p4
            % Shockwave right
                rho3_1 = rhohug(gamma1,pi1,p23,p4,rho4_1);
                rho3_2 = rhohug(gamma2,pi2,p23,p4,rho4_2);
                m4_1 = debitchoc(gamma1,pi1,p23,p4,rho4_1);
                m4_2 = debitchoc(gamma2,pi2,p23,p4,rho4_2);
                m4   = sqrt(1/(beta1/m4_1^2+(1-beta4)/m4_2^2));

                m4_1_2=rho4_1*((gamma1+1)*(p23+pi1)+(gamma1-1)*(p4+pi1))/2;
                m4_2_2=rho4_2*((gamma2+1)*(p23+pi2)+(gamma2-1)*(p4+pi2))/2;

                m4   = sqrt(1/(beta1/m4_1_2+(1-beta4)/m4_2_2));

                u3   = u4+(p23-p4)/m4;
                sr   = u4+m4/rho4;
            else
            % expansion
                rho3_1 = rhoisent(gamma1,pi1,p23,p4,rho4_1);
                rho3_2 = rhoisent(gamma2,pi2,p23,p4,rho4_2);
                somme  = gausslg(p4,p23,rho4_1,rho4_2,beta1,(1-beta4));
                u3     = u4+somme;
                sr     = u4+c4;
            end


            f=u3-u2;

           if abs(f)>1e-5
               if iq==0
                  derivf=1e-3;
               else
                  derivf=(f-fp)/dp;
               end
               iq=1;
               p23new = p23-f/derivf;
               if p23new<1e-9
                   p23new=1e-9;
               end
               dp=p23new-p23;
               p23=p23new;
               fp=f;

           else
               break
           end

        end % while

        u23 = u3;

        rho2   = 1/(beta1/rho2_1+(1-beta1)/rho2_2);
        alpha2 = rho2*beta1/rho2_1;
        rho3   = 1/(beta1/rho3_1+(1-beta4)/rho3_2);
        alpha3 = rho3*beta1/rho3_1;
        [~,~,c2] = eoscmix(p23,rho2,rho2_1,alpha2);
        [~,~,c3] = eoscmix(p23,rho3,rho3_1,alpha3);

%--------------------------------------------------------------------------
%     Flux
%--------------------------------------------------------------------------

        sstarl=u23-c2;
        sstarr=u23+c3;

        lambda0 = u23;        
        if p23>p1
            lambda1 = sl;
            lambda2 = sl;
        elseif p23<=p1
            lambda1 = sl;
            lambda2 = sstarl;
        end

        if ~exist('lambda1','var')
            keyboard
        end

        if p23>p4
            lambda4 = sr;
            lambda3 = sr;
        elseif p23<=p4
            lambda3 = sstarr;
            lambda4 = sr;
        end
        
%         keyboard
        if lambda1>=0
            rhof   = rho1;
            uf     = u1;
            pf     = p1;
            betaf  = beta1;
            alphaf = alpha1;
        elseif lambda1<0 && lambda2>=0
            keyboard
            uf   = (gamma-1)/(gamma+1)*u1+2/(gamma+1)*a1;
            af   = u23;
            pf   = p1*(af/a1)^((2*gamma)/(gamma-1));
            rhof = rho1*(af/a1)^(2/(gamma-1));
        elseif lambda2<0 && lambda0>=0
            rhof   = rho2;
            uf     = u23;
            pf     = p23;
            betaf  = beta1;
            alphaf = alpha2;
        elseif lambda0<0 && lambda3>=0
            rhof   = rho3;
            uf     = u23;
            pf     = p23;
            betaf  = beta4;
            alphaf = alpha3;
        elseif lambda3<0 && lambda4>=0 %??????????????
            keyboard
            uf   = (gamma-1)/(gamma+1)*u4-2/(gamma+1)*a4;
            af   = -u23;
            pf   = p4*(af/a4)^(2*gamma/(gamma-1));
            rhof = rho4*(af/a4)^(2/(gamma-1));
        elseif lambda4<0
            rhof   = rho4;
            uf     = u4;
            pf     = p4;
            betaf  = beta4;
            alphaf = alpha4;
        end

        F(1,i) = rhof*uf;
        F(2,i) = rhof*uf^2+pf;
        F(3,i) = alphaf*gamma1*(pf+pi1)*uf/(gamma1-1)+(1-alphaf)*gamma2*(pf+pi2)*uf/(gamma2-1)+1/2*rhof*uf^3;
        F(4,i) = betaf*rhof*uf;
        F(5,i) = alphaf*gamma1*(pf+pi1)*uf/(gamma1-1)+1/2*rhof*uf^3*betaf;

        Lambda  = max(abs([lambda0 lambda1 lambda2 lambda3 lambda4 Lambda]));


    end %if L==R
    
end %for i loop

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c1,c2,cmix] = eoscmix(p,rho,rho1,alpha1)

    global gamma1 gamma2
    global pi1 pi2

    rho2=(rho-alpha1*rho1)/(1-alpha1+eps);
    c1=sqrt(gamma1*(p+pi1)/rho1);
    c2=sqrt(gamma2*(p+pi2)/rho2);
    cmix=sqrt(1/(rho*(alpha1/(rho1*c1^2)+(1-alpha1)/(rho2*c2^2))));

end

function rho = rhohug(gamma,pii,p,p0,rho0)

%     v=((gamma+1)*(p0+pinf)+(gamma-1)*(p+pinf))/rho0;
%     v=v/((gamma+1)*(p+pinf)+(gamma-1)*(p0+pinf));
%     rho=1/v;

    rho = rho0 * ((gamma+1)*(p+pii)+(gamma-1)*(p0+pii)) / ((gamma+1)*(p0+pii)+(gamma-1)*(p+pii));

end

function m = debitchoc(gamma,pii,p,p0,rho0)

    m=sqrt(rho0*((gamma+1)*(p+pii)+(gamma-1)*(p0+pii))/2);

end



function int = gausslg(p0,p,rho01,rho02,y1,y2)

    global gamma1 gamma2
    global pi1 pi2

    z = [0.2386191861,0.6612093865,0.9324695142];
    w = [0.4679139346,0.3607615730,0.1713244924];
    xm =0.5*(p0+p);
    xr =0.5*(p-p0);
    int=0.;
    for i=1:3
        dx =xr*z(i);
        xmp=xm+dx;
        xmm=xm-dx;
        [~,z1] = rhoisent(gamma1,pi1,xmp,p0,rho01);
        [~,z2] = rhoisent(gamma2,pi2,xmp,p0,rho02);
        zmp=sqrt(y1/(z1+eps)^2+y2/(z2+eps)^2);

        [~,z1] = rhoisent(gamma1,pi1,xmm,p0,rho01);
        [~,z2] = rhoisent(gamma2,pi2,xmm,p0,rho02);
        zmm=sqrt(y1/(z1+eps)^2+y2/(z2+eps)^2);

        int=int+w(i)*(zmp+zmm);
    end
    int=xr*int;

end


function [rho,z]=rhoisent(gamma,pinf,p,p0,rho0)

  	rho=rho0*((p+pinf)/(p0+pinf))^(1/gamma);
% 	c=sqrt(gamma*(p+pinf)/rho);
	z=sqrt(rho*gamma*(p+pinf)); % =rho*c;

end