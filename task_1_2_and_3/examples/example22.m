%Nelson Campos
% Systems Identification 
% Example 2.2 of the book System Identification: Theory for the User (2nd Edition).

clc
close all
clear all

Ar = 2;
N0 = 10;
s = 10;
N = s*N0;

t = 1:N;

for r=(-N0/2+1):N0/2
    u = (1/sqrt(N0))*Ar*exp(2*pi*i*t*r/N0);
end

figure(1)
hold on
plot(t,u), xlabel('Tempo (t)'),ylabel('u(t)'), title('Sinal de entrada')

t=1:N;
w=linspace(0,pi,129);
for index_w=1:size(w,2)
    Un(index_w)=sum(u.*exp(-i*w(index_w)*t));
end
Un=1/sqrt(N)*Un;

figure(2)
plot(w,abs(Un).^2,'r'),xlabel('Frequência (w)'),ylabel('|U(w)|^2'), title('DEP de u(t)')

hold off
