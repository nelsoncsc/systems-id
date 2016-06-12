%Nelson Campos
% Systems Identification 
% Example 2.3 of the book System Identification: Theory for the User (2nd Edition).

clc
close all
clear all

N0 = 5;
s = 20;
N = s*N0;
A = 4;
t = [1:1:N]; 
w0 = 2*pi/N0;

u = A*cos(w0*t);

tau=-100:1:100;
Rs=zeros(1,size(tau,2));
for indiceTau=1:size(tau,2)
    u2 = A*cos(w0*(t-tau(indiceTau)));
    Rs(indiceTau)=Rs(indiceTau)+u*u2';
end
Rs=1/N*Rs;
    
    
w=linspace(-pi,pi,129);
Phi_s=zeros(1,size(w,2));
for index_w=1:size(w,2)
    Phi_s(index_w)=Rs*exp(-i*w(index_w)*tau)';
end

Phi_s=1/sqrt(size(tau,2))*Phi_s; %lenght of tau must be the same of the time

% figure(1)
% hold on
% plot(t,u), xlabel('Tempo (t)'),ylabel('s(t) = A*cos(w0*t)'), title('Sinal de entrada')
% figure(2)
% plot(t,u2), xlabel('Tempo (t)'),ylabel('s(t) = A*cos(w0*t)'), title('Sinal de entrada')
% hold off

% Rs = sum(u.*u2)/N %N = s*N0 and u está definido no intervalo t=[0,N=sN0]
% w = linspace(-100,100,1000);
% for index_w=1: size(w,2)
%     Phi_s(index_w) = Rs*exp(-i*w(index_w)*tau);
% end

plot(w, abs(Phi_s)), xlabel('Frequência (w)'),ylabel('|\phi_s(w)|^2'), title('Espectro do sinal s(t)')

   
