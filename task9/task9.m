%November 30, 2015
%Nelson Campos 
%nelson.campos@ee.ufcg.edu.br
%This file computes the Levy's Method for Complex-Curve Fitting
%Reference: Complex-Curve Fitting

k = 1:14;
w = [0.0 0.1 0.2 0.5 0.7 1.0 2.0 4.0 7.0 10.0 20.0 40.0 70.0 100.0]';

for k = 1: 14
    F(k) = (1+i*w(k))/(1+0.1*i*w(k)-0.01*w(k)^2); %The collected data
    F1(k) = 0.05591*i*w(k)*(1+0.0028*i*w(k))/(1+0.398*i*w(k)); %delta(Vfd)/delta(Id)
    F2(k) = 0.4928*i*w(k)*(1+0.0028*i*w(k))/((1+5.06*i*w(k))*(1+0.0049*i*w(k))); %delta(Ifd)/delta(Id)
    magnitude(k) = abs(F(k));
    phase(k) = 180*angle(F(k))/pi;
    R(k) = real(F1(k));
    I(k) = imag(F1(k));
end

%Setting initial conditions
lambda_0 = 0;
lambda_2 = 0;
lambda_4 = 0;
S0 = 0;
S2 = 0;
S4 = 0;
T1 = 0;
T3 = 0;
U2 = 0;
U4 = 0;

%Computing the matrix elements
lambda_0 = sum(w.^0);
lambda_2 = sum(w.^2);
lambda_4 = sum(w.^4);
S0 = sum((w.^0).*R');
S2 = sum((w.^2).*R');
S4 = sum((w.^4).*R');
T1 = sum(w.*I');
T3 = sum((w.^3).*I');
U2 = sum((w.^2).*(I'.^2+R'.^2));
U4 = sum((w.^4).*(I'.^2+R'.^2));  

%Defining the matrix
M = [lambda_0 0 -lambda_2 T1 S2; 
     0 lambda_2 0 -S2 T3; 
     lambda_2 0 -lambda_4 T3 S4;
     T1 -S2 -T3 U2 0; 
     S2 T3 -S4 0 U4];

C = [S0 T1 S2 0 U2]';

%Obtaining the estimated coefficents
N = inv(M)*C;

for k=1:14
    G(k) = (N(1)+N(2)*i*w(k)-N(3)*w(k)^2)/((1+N(4)*i*w(k)-N(5)*w(k)^2));
end

figure(1), plot(w, abs(F1), 'k-', w, abs(G), 'ro'), xlabel('Frequency (rad/s)'), ylabel('Magnitude(w)'), title('Complex Curve Fitting: Levy´s Method')

%For the second function
for k = 1: 14
    F(k) = (1+i*w(k))/(1+0.1*i*w(k)-0.01*w(k)^2); %The collected data
    F1(k) = 0.05591*i*w(k)*(1+0.0028*i*w(k))/(1+0.398*i*w(k)); %delta(Vfd)/delta(Id)
    F2(k) = 0.4928*i*w(k)*(1+0.0028*i*w(k))/((1+5.06*i*w(k))*(1+0.0049*i*w(k))); %delta(Ifd)/delta(Id)
    magnitude(k) = abs(F(k));
    phase(k) = 180*angle(F(k))/pi;
    R(k) = real(F2(k));
    I(k) = imag(F2(k));
end

%Setting initial conditions
lambda_0 = 0;
lambda_2 = 0;
lambda_4 = 0;
S0 = 0;
S2 = 0;
S4 = 0;
T1 = 0;
T3 = 0;
U2 = 0;
U4 = 0;

%Computing the matrix elements
lambda_0 = sum(w.^0);
lambda_2 = sum(w.^2);
lambda_4 = sum(w.^4);
S0 = sum((w.^0).*R');
S2 = sum((w.^2).*R');
S4 = sum((w.^4).*R');
T1 = sum(w.*I');
T3 = sum((w.^3).*I');
U2 = sum((w.^2).*(I'.^2+R'.^2));
U4 = sum((w.^4).*(I'.^2+R'.^2));  

%Defining the matrix
M = [lambda_0 0 -lambda_2 T1 S2; 
     0 lambda_2 0 -S2 T3; 
     lambda_2 0 -lambda_4 T3 S4;
     T1 -S2 -T3 U2 0; 
     S2 T3 -S4 0 U4];

C = [S0 T1 S2 0 U2]';

%Obtaining the estimated coefficents
N = inv(M)*C;

for k=1:14
    G(k) = (N(1)+N(2)*i*w(k)-N(3)*w(k)^2)/((1+N(4)*i*w(k)-N(5)*w(k)^2));
end

figure(2), plot(w, abs(F2), 'k-', w, abs(G), 'ro'), xlabel('Frequency (rad/s)'), ylabel('Magnitude(w)'), title('Complex Curve Fitting: Levy´s Method')

