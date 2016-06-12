%Nelson Campos
% Systems Identification
% This task estimates the parameters of the RLC circuit. 
% The simulink file is attached in the same directory of the task7.m
%Data acquisiton

%Continuous Data
timeC = data_c.time;
dataC = [data_c.signals.values];
inputC = dataC(:,1);
outputC = dataC(:,2);
figure(1), plot(timeC, inputC, 'r', timeC, outputC, 'g'), xlabel('time'), ylabel('y'), title('Continuous data'), legend('Input', 'Output')

%Discrete Data
timeD = data_d.time;
dataD = [data_d.signals.values];
inputD = dataD(:,1);
outputD = dataD(:,2);
figure(2), plot(timeD, inputD, 'r', timeD, outputD, 'g'), xlabel('time'), ylabel('y'), title('Discrete data'), legend('Input', 'Output')

N = size(inputD,1); %number of discrete points

%From the non-linear model, the predicted output is defined as bellow:
%y(t) = -a1*y(t-1) - a2*y(t-2) + b1*u(t-1) + b2*u(t-2)

u = inputD;
y = outputD;
%Definition of the non-linear regressor
phi = zeros(2,N);
for(k = 3:N)
    phi(1, k) = -y(k-1);
    phi(2, k) = -y(k-2);
    phi(3,k) = u(k-1);
    phi(4,k) = u(k-2);
end

%Definition of theta
%theta = (((sum(phi*phi'))/N).^-1).*(sum(phi*y)/N);
theta = inv(phi*phi')*(phi*y);
%theta=theta';

%From the non-linear predictor
y2 = phi'*theta;
figure(3), plot(1:N, outputD, 's', 1:N, y2, 'r'), xlabel('x'), ylabel('y'), legend('Measured', 'Estimated'), title('The estimated function')

num = [0 theta(3) theta(4)];
den = [1 theta(1) theta(2)];
sys = tf(num, den, 50e-3);
sysc = d2c(sys);
s = tf('s');
%From the model, the zero could be neglected and sysc_approx = (0.9998/(s^2+s+1));
%Where the ideal model is H(s) = (R/LC)/(s^2+s*R/L+1/LC)
sysc_approx = (0.9998/(s^2+s+1));
%Then, R ~= L ~= C = 1

