mov = avifile('velocity.avi');
mov.Quality = 75;
for i=1:n
close(figure(1));
figure(1);
hold on;
plot(x(1:N),W(1,1:N,i),'.b')
plot(x(1:N),W(2,1:N,i),'.r')
plot(x(1:N),W(3,1:N,i),'.g')
plot(x(1:N),W(4,1:N,i),'.c');
plot(x(1:N),W(5,1:N,i),'.m');
% axis([0 1 0.9 1.1]);
legend('\beta','\alpha',2);%'\rho','p',
grid;
G = getframe(gca);
mov = addframe(mov,G);
end
mov = close(mov);