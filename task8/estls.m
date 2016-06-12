%Nelson Campos
% Systems Identification
% This task estimates the parameters of an induction machine using non-linear least squares estimation.
% References: Nonlinear parameter estimation of steady-state induction machine models

clear all
close all
%
global N V p w medidas
%
Rr=3.84;
Xr=6.789;
Xs=1.658;
Rs=1.93;
Xm=38.7;
%
N=25;
V=220;
f=60;
w=2*pi*f;
p=2;
In=5.8/50;
Pn=1.5e3/50;
Tn=8.0/50;
s=[0.02:(1-0.02)/N:1.0]';
for i=1:length(s),
%
A=Rs*(1+Xr/Xm)+(1+Xs/Xm)*Rr/s(i);
B=Xr+Xs*(1+Xr/Xm)-Rs*Rr/Xm/s(i);
C=1+Xr/Xm;
D=Rr/Xm/s(i);
Im(i,1)=V*sqrt((C*C+D*D)/(A*A+B*B));
Zeq(i,1) = sqrt(-1)*Xm*(Rr/s(i)+sqrt(-1)*Xr)/(Rr/s(i)+sqrt(-1)*(Xr+Xm)) + Rs+sqrt(-1)*Xs;
theta = angle(Zeq(i,1));
fpm(i,1) = cos(theta);
Pam(i,1) = 3*V*V/abs(Zeq(i,1))*fpm(i,1);
Prm(i,1) = 3*V*V/abs(Zeq(i,1))*sin(theta);
%
Im(i,1)=V*sqrt((C*C+D*D)/(A*A+B*B))+In*randn();
Pm(i,1)=3*V*V*(A*C-B*D)/(A*A+B*B)+Pn*randn();
Tm(i,1)=3*V*V*(p/w)*(Rr/s(i))/(A*A+B*B)+Tn*randn();
end
%
medidas=[s Im];
%
x0 = (1/1)*[3  6  2 1 30]';
lb = [0  0  0 0 0]';
ub = [10 10 5 5 50]';
%
options=optimset('Display','iter','TolCon',1e-12,'TolFun',1e-12,'TolX',1e-12);
x=lsqnonlin(@fcnls,x0,lb,ub,options);
msg = ['ideal - Rr = ' num2str(Rr) ', Xr = ' num2str(Xr) ', Xs = ' num2str(Xs) ', Rs = ' num2str(Rs) ', Xm = ' num2str(Xm)];
disp(msg)
%
Rr=x(1);
Xr=x(2);
Xs=x(3);
Rs=x(4);
Xm=x(5);
%
msg = ['final - Rr = ' num2str(Rr) ', Xr = ' num2str(Xr) ', Xs = ' num2str(Xs) ', Rs = ' num2str(Rs) ', Xm = ' num2str(Xm)];
disp(msg)
%
for i=1:length(s),
A=Rs*(1+Xr/Xm)+(1+Xs/Xm)*Rr/s(i);
B=Xr+Xs*(1+Xr/Xm)-Rs*Rr/Xm/s(i);
C=1+Xr/Xm;
D=Rr/Xm/s(i);
Ie(i,1)=V*sqrt((C*C+D*D)/(A*A+B*B));
Zeq(i,1) = sqrt(-1)*Xm*(Rr/s(i)+sqrt(-1)*Xr)/(Rr/s(i)+sqrt(-1)*(Xr+Xm)) + Rs+sqrt(-1)*Xs;
theta = angle(Zeq(i,1));
fpe(i,1) = cos(theta);
Pae(i,1) = 3*V*V/abs(Zeq(i,1))*fpe(i,1);
Pre(i,1) = 3*V*V/abs(Zeq(i,1))*sin(theta);


%
Ie(i,1)=V*sqrt((C*C+D*D)/(A*A+B*B));
Pe(i,1)=3*V*V*(A*C-B*D)/(A*A+B*B);
Te(i,1)=3*V*V*(p/w)*(Rr/s(i))/(A*A+B*B);
end
%
figure(1)
subplot(6,1,1)
plot(s,Im,'ro',s,Ie,'k-'), grid, title('I \times s')
subplot(6,1,2)
plot(s,Pm,'ro',s,Pe,'k-'), grid, title('P \times s')
subplot(6,1,3)
plot(s,Tm,'ro',s,Te,'k-'), grid, title('T \times s')
subplot(6,1,4)
plot(s,fpm,'ro',s,fpe,'k-'), grid, title('fp \times s')
subplot(6,1,5)
plot(s,Pam,'ro',s,Pae,'k-'), grid, title('Pa \times s')
subplot(6,1,6)
plot(s,Prm,'ro',s,Pre,'k-'), grid, title('Q \times s')

