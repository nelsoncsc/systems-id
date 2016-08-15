%Nelson Campos
% Systems Identification
% This task estimates the parameters of a Solar Heating System. Although the system is non-linear, it is linear on the parameters.
% The data were collected using the http://arohatgi.info/WebPlotDigitizer/
% References: TOOLS FOR SEMIPHYSICAL MODELLING

clear all
close all

%The data acquisition

%The number of collected data
N = 4948;

%The solar intensity is collected from the intensity.csv file
dataI = csvread('intensity.csv');
dataI = dataI(1:N,:);
I = dataI(:,2);

%The pump velocity is collect from pump2.csv file
dataPump = csvread('pump2.csv');

%The matrix dimensions must agree
dataPump = dataPump(1:N,:);
u = dataPump(:,2);

%The output data
dataOut = csvread('y.csv');
y = dataOut(:,2);

%From the non-linear model, the predicted output is defined as bellow:
%y(t) = theta1*y(t-1) + theta2*u(t-1)*I(t-2)

%Definition of the non-linear regressor
phi = zeros(2,N);
for(k = 3:N)
    phi(1, k) = y(k-1);
    phi(2, k) = u(k-1)*I(k-2);
end

%Definition of theta
%theta = (((sum(phi*phi'))/N).^-1).*(sum(phi*y)/N);
theta = inv(phi*phi')*(phi*y);
%From the non-linear predictor
y2 = phi'*theta;
plot(1:N, y, 'b', 1:N, y2, 'rx'), xlabel('x'), ylabel('Storage Temperature'), legend('Measured', 'Estimated'), title('A non-linear model')
