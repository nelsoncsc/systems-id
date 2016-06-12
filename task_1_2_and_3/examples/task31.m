%Nelson Campos
% Systems Identification 
% This code reproduces the Figures 3.15, 3.16 and 3.17 of the book Modeling of Dynamic Systems

close all
clear all

t = 0:200; %The time interval

e = randn(1,size(t,2)); %randon variable with normal distribution N~(0,1)

%This part generates the signal of the Figure 3.15(a)
d1a = -0.9;
c0a = 1;

wa = zeros(1,size(t,2));
for k=2:size(t,2)
    wa(k) = -d1a*wa(k-1)+c0a*e(k);
end

%This part generates the signal of the Figure 3.15(b)
d1b = 0.9;
c0b = 1;
wb = zeros(1,size(t,2));
for k=2:size(t,2)
    wb(k) = -d1b*wb(k-1)+c0b*e(k);
end

%This part generates the signal of the Figure 3.15(c)
d1c = -0.5;
d2c = 0.7;
c0c = 1;
c1c = 0.5;
wc = zeros(1,size(t,2));
for k=3:size(t,2)
    wc(k) = -d1c*wc(k-1)-d2c*wc(k-2)+c0c*e(k)+c1c*e(k-1); 
end

%This part builds the Probability Density Funtion of e(t) defined as
%P(e(t)) = 0, with probability 0.98
%P(e(t)) = sqrt(50), with probability 0.01
%P(e(t)) = -sqrt(50), with probability 0.01

r = rand(1,size(t,2)); %Notice that r is generated with rand (values between 0 and 1)
P = [0.98, 0.99, 1];
X = [0 sqrt(50) -sqrt(50)];

count = [0 0 0];
ed = zeros(1, size(t,2));
for k=1:size(r,2)
    counter = 1;
    while(r(k) > P(counter))
        counter = counter+1;
    end
    ed(k) = X(counter); %ed is the new random variable
    count(counter) = count(counter)+1;
end

average = count/size(r,2);
status1 = ['P(e(t))=0 = '  ,num2str(average(1))];
status2= ['P(e(t)) = sqrt(50)) = ', num2str(average(2))];
status3 = ['P(e(t)) = -sqrt(50)) = '  ,num2str(average(3))];
disp(status1)
disp(status2)
disp(status3)
disp('Coded by Nelson Campos')

%This part generates the signal of the Figure 3.15(d)
d1d = -0.5;
d2d = 0.7;
c0d = 1;
c1d = 0.5;
wd = zeros(1,size(t,2));
for k=3:size(t,2)
    wd(k) = -d1d*wd(k-1)-d2d*wd(k-2)+c0d*ed(k)+c1d*ed(k-1); 
end

figure(1),
subplot(2,2,1)
plot(t,wa), title('3.15 (a)')
subplot(2,2,2)
plot(t,wb), title('3.15 (b)')
subplot(2,2,3)
plot(t,wc), title('3.15 (c)')
subplot(2,2,4)
plot(t,wd), title('3.15 (d)')
hold on

%This part makes the covariance function of wa(t)
tau = -20:20;
R_wa = zeros(1, size(tau,2));
for index_tau = 1: size(tau,2);
    for k=1:size(t,2)
        if(k-tau(index_tau) >=1 && k-tau(index_tau) <= size(t,2))
            wa2(k) = wa(k-tau(index_tau));
        else
            wa2(k) = 0;    
        end
    end
    R_wa(index_tau) = R_wa(index_tau)+wa*wa2'; 
end
R_wa = R_wa/size(t,2);

%This part makes the covariance function of wb(t)
tau = -20:20;
R_wb = zeros(1, size(tau,2));
for index_tau = 1: size(tau,2);
    for k=1:size(t,2)
        if(k-tau(index_tau) >=1 && k-tau(index_tau) <= size(t,2))
            wb2(k) = wb(k-tau(index_tau));
        else
            wb2(k) = 0;    
        end
    end
    R_wb(index_tau) = R_wb(index_tau)+wb*wb2'; 
end
R_wb = R_wb/size(t,2);

%This part makes the covariance function of wc(t)
tau = -20:20;
R_wc = zeros(1, size(tau,2));
for index_tau = 1: size(tau,2);
    for k=1:size(t,2)
        if(k-tau(index_tau) >=1 && k-tau(index_tau) <= size(t,2))
            wc2(k) = wc(k-tau(index_tau));
        else
            wc2(k) = 0;    
        end
    end
    R_wc(index_tau) = R_wc(index_tau)+wc*wc2'; 
end
R_wc = R_wc/size(t,2);

%This part makes the covariance function of wd(t)
tau = -20:20;
R_wd = zeros(1, size(tau,2));
for index_tau = 1: size(tau,2);
    for k=1:size(t,2)
        if(k-tau(index_tau) >=1 && k-tau(index_tau) <= size(t,2))
            wd2(k) = wd(k-tau(index_tau));
        else
            wd2(k) = 0;    
        end
    end
    R_wd(index_tau) = R_wd(index_tau)+wd*wd2'; 
end
R_wd = R_wd/size(t,2);

figure(2),
subplot(2,2,1)
plot(tau,R_wa), title('3.16 (a)')
subplot(2,2,2)
plot(tau,R_wb), title('3.16 (b)')
subplot(2,2,3)
plot(tau,R_wc), title('3.16 (c)')
subplot(2,2,4)
plot(tau,R_wd), title('3.16 (d)')
hold on

%This part calculates the spectra of the signals
w = logspace(-2,1,1000);
Phi_wa = zeros(1,size(w,2));
for index_w=1:size(w,2)
    Phi_wa(index_w)=R_wa*exp(-i*w(index_w)*tau)';
end
Phi_wa = Phi_wa/sqrt(size(tau,2));

Phi_wb = zeros(1,size(w,2));
for index_w=1:size(w,2)
    Phi_wb(index_w)=R_wb*exp(-i*w(index_w)*tau)';
end
Phi_wb = Phi_wb/sqrt(size(tau,2));

Phi_wc = zeros(1,size(w,2));
for index_w=1:size(w,2)
    Phi_wc(index_w)=R_wc*exp(-i*w(index_w)*tau)';
end
Phi_wc = Phi_wc/sqrt(size(tau,2));

Phi_wd = zeros(1,size(w,2));
for index_w=1:size(w,2)
    Phi_wd(index_w)=R_wd*exp(-i*w(index_w)*tau)';
end
Phi_wd = Phi_wd/sqrt(size(tau,2));

figure(3),
subplot(2,2,1)
semilogx(w,20*log10(abs(Phi_wa))), title('3.17 (a)')
subplot(2,2,2)
semilogx(w,20*log10(abs(Phi_wb))), title('3.17 (b)')
subplot(2,2,3)
semilogx(w,20*log10(abs(Phi_wc))), title('3.17 (c)')
subplot(2,2,4)
semilogx(w,20*log10(abs(Phi_wd))), title('3.17 (d)')
hold off

disp('The spectra is periodic with period 2*pi/T, where T=1 as defined here')
