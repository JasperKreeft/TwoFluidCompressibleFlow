figure
for n=1:size(speeds,3)
for i=1:N+1

        if speeds(i,1,n) >= 0 && speeds(i,3,n)>=0
            if speeds(i,2,n)>=0 && speeds(i,5,n)>=0
               cas(i,n) = 1;
            elseif speeds(i,2,n)>=0 && speeds(i,5,n)<=0
               cas(i,n) = 2; 
            elseif speeds(i,2,n)<=0 && speeds(i,5,n)>=0
               cas(i,n) = 3;
            elseif speeds(i,2,n)<=0 && speeds(i,5,n)<=0
             cas(i,n) = 4;
            end

        elseif  speeds(i,1,n) >= 0 && speeds(i,3,n)<=0
            if speeds(i,2,n)>=0 && speeds(i,5,n)>=0
              cas(i,n) = 5;
            elseif speeds(i,2,n)>=0 && speeds(i,5,n)<=0
           cas(i,n) = 6;
            elseif speeds(i,2,n)<=0 && speeds(i,5,n)>=0
              cas(i,n) = 7;
            elseif speeds(i,2,n)<=0 && speeds(i,5,n)<=0
              cas(i,n) = 8;
            end

        elseif  speeds(i,1,n) <= 0 && speeds(i,4,n)>=0
            if speeds(i,2,n)>=0 && speeds(i,5,n)>=0
               cas(i,n) = 9;
            elseif speeds(i,2,n)>=0 && speeds(i,5,n)<=0
               cas(i,n) = 10;
            elseif speeds(i,2,n)<=0 && speeds(i,5,n)>=0
              cas(i,n) = 11;
            elseif speeds(i,2,n)<=0 && speeds(i,5,n)<=0
                cas(i,n) = 12;
            end

        elseif  speeds(i,1,n) <= 0 && speeds(i,4,n)<=0
            if speeds(i,2,n)>=0 && speeds(i,5,n)>=0
              cas(i,n) = 13;
            elseif speeds(i,2,n)>=0 && speeds(i,5,n)<=0
                cas(i,n) = 14;
            elseif speeds(i,2,n)<=0 && speeds(i,5,n)>=0
                cas(i,n) = 15;
            elseif speeds(i,2,n)<=0 && speeds(i,5,n)<=0
                cas(i,n) = 16;
            end
        end %flux choice

end

plot(cas(:,n),'x')
title(['n = ' num2str(n)])
pause
end

for i=1:N+1
plot(speeds(i,:,end),'o')
title(['i = ' num2str(i)])
ylim([-.15 .25])
grid
pause
end