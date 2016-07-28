%Nelson Campos
% Systems Identification 
% Example 6.1 of the book System Identification: Theory for the User (2nd Edition) 
% and question 7C.1 of the 1st edition of the same book

clear all;
close all;

N = 200;
u = idinput(N, 'prbs'); %define a N-point PRBS input;
t = 1:N;
NFFT=2^nextpow2(N);
Fa=1;

figure(1), plot(1:100, u(1:100)), title('INPUT #1')
hold on

e = randn(size(u)); %white noise with variance 1

y = zeros(size(u)); %The output signal
for k = 3 : size(u,1)
    y(k) = 1.5*y(k-1)-0.7*y(k-2)+u(k-1)+0.5*u(k-2)+e(k);
end

figure(2), plot(1:100, y(1:100)), title('OUTPUT #1')
hold off

%The frequency range
w = linspace(0, 4e0, NFFT);
epsilon = linspace(-pi,pi,1000);

%Computation for w

%Un is the Fourier Transform of u(t)
for index_w=1:size(w,2)
    Un(index_w)=sum(u.*exp(-i*w(index_w)*t'));
end
Un=1/sqrt(N)*Un;

%Yn is the Fourier Transform of y(t)
for index_w=1:size(w,2)
    Yn(index_w)=sum(y.*exp(-i*w(index_w)*t'));
end
Yn=1/sqrt(N)*Yn;

Gn_empirical = Yn./Un; %Gn_empirical is the Empirical Transfer Function

%Computation for epsilon

%Un is the Fourier Transform of u(t)
for index_e=1:size(epsilon,2)
    Un(index_e)=sum(u.*exp(-i*epsilon(index_e)*t'));
end
Un=1/sqrt(N)*Un;

%This part was used to create the function function_Hc
%Hd=tf([0 1 0.5],[1 -1.5 0.7],1)
    %Hc=d2c(Hd)
    
%Yn is the Fourier Transform of y(t)
for index_e=1:size(epsilon,2)
    Yn(index_e)=sum(y.*exp(-i*epsilon(index_e)*t'));
end
Yn=1/sqrt(N)*Yn;

vector_gamma = [10 50 200]; %To select the gamma value
type = 3; %To select the window type
for index_gamma = 1:size(vector_gamma,2)
    Un2 = abs(Un).^2;
    %To calculate the Estimated Transfer Function
    gamma = vector_gamma(index_gamma);
    num = zeros(1, size(epsilon,2));
    den = zeros(1, size(epsilon,2));
    Gn = zeros(1, size(w,2));
    for index_w = 1 : size(w,2)
        for index_e = 1 : size(epsilon,2)
            theta = epsilon(index_e)-w(index_w); %theta = epsilon - wo
            Window = windowFunction(gamma, theta, 1);
            num(index_e) = Window*Un2(index_e)*Yn(index_e)/Un(index_e);
            den(index_e) = Window*Un2(index_e);
        end
        Gn(index_w) = sum(num)/sum(den);
    end
    
    figure(3)
    subplot(2,2,1)
    loglog(w,abs(Gn_empirical),'b',w,abs(function_Hc(1i*w)),'k','Linewidth',3), title('Empirical Transfer Function');
    subplot(2,2,index_gamma+1)
    
    loglog(w,abs(Gn),'b',w,abs(function_Hc(1i*w)),'k','Linewidth',3), title('Estimated Transfer Function');
    %hold off
end

%Now using the expression (6.51)
num2 = zeros(1, size(w,2));
den2 = zeros(1, size(w,2));
Gn2 = zeros(1, size(w,2)); %The estimated Transfer Function

Ryu = xcorr(y, u);
Ruu = xcorr(u,u);
Phi_yu = fft(Ryu);
Phi_u = fft(Ruu);
Gn2 = Phi_yu./Phi_u;

figure(4), semilogx(w,20*log10(abs(Gn2(1:NFFT)))), title('Estimated Transfer Function with equation (6.51)');

%Task #5 starts here
x = zeros(1, N);
%definition of x(t)
for k = 3: N
    x(k) = u(k-2); 
end

%definition of zeta(t)
zeta = zeros(4, N);
for k = 3: N
    zeta(1, k) = -x(k-1);
    zeta(2,k) = -x(k-2);
    zeta(3,k) = u(k-1);
    zeta(4,k) = u(k-2);
end

%definition of phi
phi = zeros(4,N);
for k = 3:N
    phi(1,k) = -y(k-1);
    phi(2,k) = -y(k-2);
    phi(3,k) = u(k-1);
    phi(4,k) = u(k-2);
end

%definition of theta_IV
theta_IV = ((sum(zeta*phi')/N).^-1)*(sum(zeta*y)/N);
y_IV = phi'*theta_IV'+e;
figure(5),
plot(y_IV), title('y(t) retrieved from \theta_{IV}')

%c2: M and N being the LS of A and B respectively
theta_c2 = ((sum(phi*phi')).^-1)*(sum(phi*y));

xc2 = zeros(1, N);
%definition of x(t)
for k = 3: N
    xc2(k) = 0.7801*xc2(k-1)+0.7629*xc2(k-2)-48.6379*u(k)+98.5326*u(k-1); 
end

%definition of zeta(t)
zeta_c2 = zeros(4, N);
for k = 3: N
    zeta_c2(1, k) = -xc2(k-1);
    zeta_c2(2,k) = -xc2(k-2);
    zeta_c2(3,k) = u(k-1);
    zeta_c2(4,k) = u(k-2);
end

theta_c2_IV = ((sum(zeta_c2*phi')/N).^-1)*(sum(zeta_c2*y)/N);
yc2_IV = phi'*theta_c2_IV'+e;
figure(6),
plot(yc2_IV), title('y(t) retrieved from \theta_{IV} in case 2')

%c3: y(t) = B(q)*u(t)/A(q)+e(t), K(q) = 1/A(q), M(q) = A(q), N(q) = B(q) 
xc3 = zeros(1, N);
%definition of x(t)
for k = 3: N
    xc3(k) = y(k)-e(k); 
end

%definition of zeta(t)
zeta_c3 = zeros(4, N);
for k = 3: N
    zeta_c3(1, k) = 1.5*zeta_c3(1,k-1)-0.7*zeta_c3(1, k-2)-xc3(k-1);
    zeta_c3(2,k) = 1.5*zeta_c3(2,k-1)-0.7*zeta_c3(2, k-2)-xc3(k-2);
    zeta_c3(3,k) = 1.5*zeta_c3(3,k-1)-0.7*zeta_c3(3, k-2)+u(k-1);
    zeta_c3(4,k) = 1.5*zeta_c3(4,k-1)-0.7*zeta_c3(4, k-2)+u(k-2);
end

theta_c3_IV = ((sum(zeta_c3*phi')/N).^-1)*(sum(zeta_c3*y)/N);
yc3_IV = phi'*theta_c3_IV'+e;
figure(7),
plot(yc3_IV), title('y(t) retrieved from \theta_{IV} in case 3')

