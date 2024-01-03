% IC_Bubble
clear all
close all
clc


NX = 200;
NY = 100;

W = zeros(NY,NX,8);

%-------------------------------------------%
%                                           %
% Initial conditions and settings           %
%                                           %
%-------------------------------------------%


%-------------------------------------------%
%       Generate grid                       %
%-------------------------------------------%

X = 0.1600;
Y = 2*0.0445;
R = 0.0250;

dx = X/NX;
dy = Y/NY;

for i=1:NY
    for j=1:NX
        W(i,j,1)=(j-0.5)*dx-0.080;
        W(i,j,2)=Y/2-(i-0.5)*dy;
    end
end

%AIR
gamma1 = 1.4;

%air_0
Wair0(1) = 1.4;
Wair0(2) = 0;
Wair0(3) = 0;
Wair0(4) = 1;
Wair0(5) = 1;
Wair0(6) = 1;

%air_1
Wair1(1) = 1.92691;
Wair1(2) = -0.33361;
Wair1(3) = 0;
Wair1(4) = 1.5698;
Wair1(5) = 1;
Wair1(6) = 1;

%Helium of R22
%gamma2 = 1.648  % Helium
gamma2 = 1.249;   % R22

%Helium + 28% air
WHe(1) = 0.25463;
WHe(2) = 0;
WHe(3) = 0;
WHe(4) = 1;
WHe(5) = 0;
WHe(6) = 0;

%R22
WR22(1) = 4.41540;
WR22(2) = 0;
WR22(3) = 0;
WR22(4) = 1;
WR22(5) = 0;
WR22(6) = 0;

for i=1:NY
    for j=1:NX
        xmin = min(sqrt((W(i,j,1)-0.5*dx)^2),sqrt((W(i,j,1)+0.5*dx)^2));
        xmax = max(sqrt((W(i,j,1)-0.5*dx)^2),sqrt((W(i,j,1)+0.5*dx)^2));
        ymin = min(sqrt((W(i,j,2)-0.5*dy)^2),sqrt((W(i,j,2)+0.5*dy)^2));
        ymax = max(sqrt((W(i,j,2)-0.5*dy)^2),sqrt((W(i,j,2)+0.5*dy)^2));

        rmax = sqrt(xmax^2+ymax^2);
        rmin = sqrt(xmin^2+ymin^2);

        if (rmax < R)
        % completely inside bubble

            % W(i,j,3:8)=WHe
            W(i,j,3:8)=WR22;

        elseif (rmin > R)
        %completely outside bubble

            W(i,j,3:8)=Wair0;

        else
        %mixture

            %stel
            y1 = ymin;
            x1 = sqrt(R^2-ymin^2);
            if (x1 > xmax)
                x1 = xmax;
                y1 = sqrt(R^2-xmax^2);
            end

            %stel
            y2 = ymax;
            x2 = sqrt(R^2-ymax^2);
            if (x2 < xmin)
                x2 = xmin;
                y2 = sqrt(R^2-xmin^2);
            end

            if (x1<xmax)
                if (y2<ymax)
                    opp = dx*dy-0.5*(x1-xmin)*(y2-ymin);
                elseif (y2==ymax)
                    opp = 0.5*(xmax-x1+xmax-x2)*dy;
                end
            elseif (x1==xmax)
                if (y2<ymax)
                    opp = 0.5*(ymax-y1+ymax-y2)*dx;
                elseif (y2==ymax)
                    opp = 0.5*(xmax-x2)*(ymax-y1);
                end
            end

            W(i,j,8) = opp/(dx*dy);
            W(i,j,3) = W(i,j,8)*Wair0(1)+(1-W(i,j,8))*WR22(1);
            W(i,j,4) = 0;
            W(i,j,5) = 0;
            W(i,j,6) = 1;
            W(i,j,7) = W(i,j,8)*Wair0(1)/W(i,j,3);

        end
    end
end

for i=1:NY
    for j=(3*NX/4):NX
        W(i,j,3:8)=Wair1;
    end
end

disp('Saving Matlab-readable output file: 0.mat')
save('0','W');

% contour2D