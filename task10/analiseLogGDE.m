%Nelson Campos
%This task estimates the percentage of the CPU of a computed excited with some loads (rendering a video, play a music, ...) 
%Details in the Systems_ID_Presentation.pdf

close all;
clear all;
clc;

% vetor de valores de potencia
dados=load('logIdentificacaoVideo.txt');
dados2=load('logIdentificacaoMusica.txt');
dados3=load('logIdentificacaoEbizzy.txt');
dados4=load('logIdentificacaoOS.txt');

%__________________________________________________________________________________________________________________________________
%For dados
intervalo = 200:1:1000;
intervalo2 = 1:length(dados4);

Tm = 1;
t=Tm * intervalo';

u = dados(intervalo,6);
y = dados(intervalo,7);



%Definition of the regression vector
phi = zeros(2,size(u,1));
for(k = 2:size(u,1))
    phi(1, k) = -y(k-1);
    phi(2, k) = -y(k-1)*u(k-1);
end

N = size(u,1);
%Definition of theta
%theta = (((sum(phi*phi'))/N).^-1).*(sum(phi*y)/N);
theta = inv(phi*phi')*(phi*y);
%theta = lsqnonneg(phi', y);

%From the obtained predictor
y2 = phi'*theta;

figure(1)
subplot(3,1,1)
stairs(t, u,'-r');
ylabel('Entrada, u')
xlabel('Tempo, t')
subplot(3,1,2)
stairs(t, y,'-');
ylabel('Saída, y')
xlabel('Tempo, t')
subplot(3,1,3)
plot(1:N, y, 'b-', 1:N, y2, 'r'), xlabel('Tempo, t'), ylabel('Saída, y'), legend('Measured', 'Estimated'), title('The estimated function for Video')

V1 = sum((y-y2).^2)/N
%__________________________________________________________________________________________________________________________________
%For dados2
intervalo2 = 1:length(dados2);

Tm = 1;
t=Tm * intervalo2';

u = dados2(intervalo2,6);
y = dados2(intervalo2,7);

%Definition of the regression vector
phi = zeros(2,size(u,1));
for(k = 2:size(u,1))
    phi(1, k) = -y(k-1);
    phi(2, k) = -y(k-1)*u(k-1);
end

N = size(u,1);
%Definition of theta
%theta = (((sum(phi*phi'))/N).^-1).*(sum(phi*y)/N);
theta = inv(phi*phi')*(phi*y);
%theta = lsqnonneg(phi', y);

%From the obtained predictor
y2 = phi'*theta;
figure(2)
subplot(3,1,1)
stairs(t, u,'-r');
ylabel('Entrada, u')
xlabel('Tempo, t')
subplot(3,1,2)
stairs(t, y,'-');
ylabel('Saída, y')
xlabel('Tempo, t')
subplot(3,1,3)
plot(1:N, y, 'b-', 1:N, y2, 'r'), xlabel('Tempo, t'), ylabel('Saída, y'), legend('Measured', 'Estimated'), title('The estimated function for Music')

V2 = sum((y-y2).^2)/N

%__________________________________________________________________________________________________________________________________
%For dados3
intervalo3 = 1:length(dados3);

Tm = 1;
t=Tm * intervalo3';

u = dados3(intervalo3,6);
y = dados3(intervalo3,7);

%Definition of the regression vector
phi = zeros(2,size(u,1));
for(k = 2:size(u,1))
    phi(1, k) = -y(k-1);
    phi(2, k) = -y(k-1)*u(k-1);
end

N = size(u,1);
%Definition of theta
%theta = (((sum(phi*phi'))/N).^-1).*(sum(phi*y)/N);
theta = inv(phi*phi')*(phi*y);
%theta = lsqnonneg(phi', y);

%From the obtained predictor
y2 = phi'*theta;
figure(3)
subplot(3,1,1)
stairs(t, u,'-r');
ylabel('Entrada, u')
xlabel('Tempo, t')
subplot(3,1,2)
stairs(t, y,'-');
ylabel('Saída, y')
xlabel('Tempo, t')
subplot(3,1,3)
plot(1:N, y, 'b-', 1:N, y2, 'r'), xlabel('Tempo, t'), ylabel('Saída, y'), legend('Measured', 'Estimated'), title('The estimated function for Ebizzy')

V3 = sum((y-y2).^2)/N
%__________________________________________________________________________________________________________________________________
%For dados4
intervalo4 = 1:length(dados4);

Tm = 1;
t=Tm * intervalo4';

u = dados4(intervalo4,6);
y = dados4(intervalo4,7);

%Definition of the regression vector
phi = zeros(2,size(u,1));
for(k = 3:size(u,1))
    phi(1, k) = -y(k-1);
    phi(2, k) = u(k-1);
    phi(3, k) = y(k-2)*u(k-1);
    phi(4, k) = y(k-1)*u(k-2);
end

N = size(u,1);
%Definition of theta
%theta = (((sum(phi*phi'))/N).^-1).*(sum(phi*y)/N);
theta = inv(phi*phi')*(phi*y);
%theta = lsqnonneg(phi', y);

%From the obtained predictor
y2 = phi'*theta;
figure(4)
subplot(3,1,1)
stairs(t, u,'-r');
ylabel('Entrada, u')
xlabel('Tempo, t')
subplot(3,1,2)
stairs(t, y,'-');
ylabel('Saída, y')
xlabel('Tempo, t')
subplot(3,1,3)
plot(1:N, y, 'b-', 1:N, y2, 'r'), xlabel('Tempo, t'), ylabel('Saída, y'), legend('Measured', 'Estimated'), title('The estimated function for OS')

V4 = sum((y-y2).^2)/N


