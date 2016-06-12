function f = fcnls(x)
global N V p w medidas
%
s=medidas(:,1);
Im=medidas(:,2);
%
Rr=x(1);
Xr=x(2);
Xs=x(3);
Rs=x(4);
Xm=x(5);
%
for i=1:length(s),
    A=Rs*(1+Xr/Xm)+(1+Xs/Xm)*Rr/s(i);
    B=Xr+Xs*(1+Xr/Xm)-Rs*Rr/Xm/s(i);
    C=1+Xr/Xm;
    D=Rr/Xm/s(i);
    Ic(i,1)=V*sqrt((C*C+D*D)/(A*A+B*B));
end
%
f = Im-Ic;
