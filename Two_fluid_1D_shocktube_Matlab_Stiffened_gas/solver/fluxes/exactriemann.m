function [F,sL,sR,Lambda] = exactriemann(N,Wf,gamma1,gamma2,pi1,pi2,Lambda)

eps = 1e-12;
F  = zeros(5,N+1);
sL = zeros(1,N); % waarden zijn gedefinieerd voor de cel index en niet voor de flux index
sR = zeros(1,N); % waarden zijn gedefinieerd voor de cel index en niet voor de flux index

% i = j - 1/2
for i  = 2:N
% if i==101
%     keyboard
% end
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
    
    rho1_1 = beta1*rho1/alpha1;
    rho1_2 = (1-beta1)*rho1/(1-alpha1);
    rho4_1 = beta4*rho4/alpha4;
    rho4_2 = (1-beta4)*rho4/(1-alpha4);

    % Check for negative densities and pressures
    if rho1<=0
        disp('rho1 negative')
        rho1 = eps;
    end

    if p1<=0
        disp('p1 negative')
        p1 = eps;
    end

    if rho4<=0
        disp('rho4 negative')
        rho4 = eps;
    end

    if p4<=0
        disp('p4 negative')
        p4 = eps;
    end

    % Equal state in both cells
    if abs(rho1-rho4)<1e-6 && abs(u1-u4)<1e-6 && abs(p1-p4)<1e-6 && abs(beta1-beta4)<1e-6 && abs(alpha1-alpha4)<1e-6
    
       F(:,i) = [rho1*u1;
                 rho1*u1^2+p1;
                 alpha1*gamma1*(p1+pi1)*u1/(gamma1-1)+(1-alpha1)*gamma2*(p1+pi2)*u1/(gamma2-1)+1/2*rho1*u1^3;
                 rho1*u1*beta1;
                 alpha1*gamma1*(p1+pi1)*u1/(gamma1-1)+1/2*rho1*u1^3*beta1];

        sL(i-1) = 0;
        sR(i)   = 0;

        Lambda  = max(abs([u1 u1-a1 u4+a4 Lambda]));

    % non-moving interface
    elseif abs(p1-p4)<1e-6 && abs(u1)<1e-6 && abs(u4)<1e-6

        F(:,i) = [0;p1;0;0;0];

        sL(i-1) = 0;
        sR(i)   = 0;


    % Two-fluid problem
    else

        rho1=alpha1*rho1_1+(1-alpha1)*rho1_2;
        rho4=alpha4*rho4_1+(1-alpha4)*rho4_2;

        beta1=alpha1*rho1_1/rho1;
        beta4 = alpha4*rho4_1/rho4;

        [c1_1,c1_2,c1] = eoscmel(p1,rho1,rho1_1,alpha1,gamma1,gamma2,pi1,pi2);
        z1=rho1*c1;

        [c4_1,c4_2,c4] = eoscmel(p4,rho4,rho4_1,alpha4,gamma1,gamma2,pi1,pi2);
        z4=rho4*c4;

        p23=(z1*p4+z4*p1+z4*z1*(u1-u4))/(z1+z4);
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
                m1=1.d0/(beta1/m1_1^2+(1-beta1)/m1_2^2);
                m1=sqrt(m1);
                u2=u1-(p23-p1)/m1;
                sl=u1-m1/rho1;
            else
            % expansion left
                rho2_1 = rhoisent(gamma1,pi1,p23,p1,rho1_1);
                rho2_2 = rhoisent(gamma2,pi2,p23,p1,rho1_2);
                [somme] = gausslg(p1,p23,rho1_1,rho1_2,beta1,(1-beta1),gamma1,gamma2,pi1,pi2);
                u2=u1-somme;
                sl=u1-c1;
            end


        %
        %       Right state
        %
            if p23>p4
            % Shockwave right
                rho3_1 = rhohug(gamma1,pi1,p23,p4,rho4_1);
                rho3_2 = rhohug(gamma2,pi2,p23,p4,rho4_2);
                [m4_1] = debitchoc(gamma1,pi1,p23,p4,rho4_1);
                [m4_2] = debitchoc(gamma2,pi2,p23,p4,rho4_2);
                m4=1.d0/(beta1/m4_1^2+(1-beta4)/m4_2^2);
                m4=sqrt(m4);
                u3=u4+(p23-p4)/m4;
                sr=u4+m4/rho4;
            else
            % expansion
                rho3_1 = rhoisent(gamma1,pi1,p23,p4,rho4_1);
                rho3_2 = rhoisent(gamma2,pi2,p23,p4,rho4_2);
                [somme] = gausslg(p4,p23,rho4_1,rho4_2,beta1,(1-beta4),gamma1,gamma2,pi1,pi2);
                u3=u4+somme;
                sr=u4+c4;
            end


            f=u3-u2;

           if abs(f)>1.e-5
               if iq==0
                  derivf=1.e-3;
               else
                  derivf=(f-fp)/dp;
               end
               iq=1;
               p23new=p23-f/derivf;
               if p23new<1.e-9
                   p23new=1.e-9;
               end
               dp=p23new-p23;
               p23=p23new;
               fp=f;

           else
               break
           end

        end % while

        u23 = u3;

        %-----------------------------------------------------------------------
        %     Echantillonage du flux 
        %-----------------------------------------------------------------------   

        rho2   = 1/(beta1/rho2_1+(1-beta1)/rho2_2);
        alpha2 = rho2*beta1/rho2_1;
        rho3   = 1/(beta1/rho3_1+(1-beta4)/rho3_2);
        alpha3 = rho3*beta1/rho3_1;
        [c2_1,c2_2,c2] = eoscmel(p23,rho2,rho2_1,alpha2,gamma1,gamma2,pi1,pi2);
        [c3_1,c3_2,c3] = eoscmel(p23,rho3,rho3_1,alpha3,gamma1,gamma2,pi1,pi2);

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

        if p23>p4
            lambda4 = sr;
            lambda3 = sr;
        elseif p23<=p4
            lambda3 = sstarr;
            lambda4 = sr;
        end
% keyboard
        if lambda1>=0
            rhof   = rho1;
            uf     = u1;
            pf     = p1;
            betaf  = beta1;
            alphaf = alpha1;
        elseif lambda1<0 && lambda2>=0 %???????????????
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

function [c1,c2,cmel] = eoscmel(p,rho,rho1,alpha1,gamma1,gamma2,pi1,pi2)

    rho2=(rho-alpha1*rho1)/(1-alpha1);
    c1=sqrt(gamma1*(p+pi1)/rho1);
    c2=sqrt(gamma2*(p+pi2)/rho2);
    vitson=alpha1/rho1/(c1*c1)+(1-alpha1)/rho2/(c2*c2);
    vitson=1.d0/(rho*vitson);
    cmel=sqrt(vitson);

end

function [rho] = rhohug(gamma,pinf,p,p0,rho0)

    gm=gamma-1.d0;
    gp=gamma+1.d0;

    v0=1.d0/rho0;
    v=v0*(gp*(p0+pinf)+gm*(p+pinf));
    v=v/(gp*(p+pinf)+gm*(p0+pinf));
    rho=1.d0/v;

end

function [m] = debitchoc(gamma,pinf,p,p0,rho0)

    gm=gamma-1.d0;
    gp=gamma+1.d0;

    c0=sqrt(gamma*(p0+pinf)/rho0);
    m=rho0*c0*sqrt((gp*(p+pinf)/(p0+pinf)+gm)/2.d0/gamma);

end



function [somme] = gausslg(p0,p,rho01,rho02,y1,y2,gamma1,gamma2,pi1,pi2)

%!     implicit real*8(a-h,o-z)

    z = [0.2386191861d0,0.6612093865d0,0.9324695142d0];
    w = [0.4679139346d0,0.3607615730d0,0.1713244924d0];
    xm =0.5d0*(p0+p);
    xr =0.5d0*(p-p0);
    somme=0.d0;
    for i=1:3
        dx =xr*z(i);
        xmp=xm+dx;
        xmm=xm-dx;
        [rho,z1] = rhoisent(gamma1,pi1,xmp,p0,rho01);
        [rho,z2] = rhoisent(gamma2,pi2,xmp,p0,rho02);
        zmp=sqrt(y1/z1^2+y2/z2^2);

        [rho,z1] = rhoisent(gamma1,pi1,xmm,p0,rho01);
        [rho,z2] = rhoisent(gamma2,pi2,xmm,p0,rho02);
        zmm=sqrt(y1/z1^2+y2/z2^2);

        somme=somme+w(i)*(zmp+zmm);
    end
    somme=xr*somme;

end


function [rho,z]=rhoisent(gamma,pinf,p,p0,rho0)
%     calcul de la densite le long de l'isentrope et
%     de l'impedance acoustique (utile pour invariants de Riemann)

  	rho=rho0*((p+pinf)/(p0+pinf))^(1.d0/gamma);
	c=sqrt(gamma*(p+pinf)/rho);
	z=rho*c;

end