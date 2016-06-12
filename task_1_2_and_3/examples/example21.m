% Systems Identification 
% Example 2.1 of the book System Identification: Theory for the User (2nd Edition).

clc
close all
clear all

N0 = 5;
s = 20;
N = s*N0;
A = 2;
t = [1:1:N];
w0 = 2*pi/N0;

u = A*cos(w0*t);

figure(1)
hold on
plot(t,u), xlabel('Tempo (t)'),ylabel('u(t) = A*cos(w0*t)'), title('Sinal de entrada')

t=1:N;
w=linspace(0,pi,129);
for index_w=1:size(w,2)
    Un(index_w)=sum(u.*exp(-i*w(index_w)*t));
end
Un=1/sqrt(N)*Un;

figure(2)
plot(w,abs(Un).^2,'r'),xlabel('Frequência (w)'),ylabel('|U(w)|^2'), title('DEP de u(t)')

hold off
