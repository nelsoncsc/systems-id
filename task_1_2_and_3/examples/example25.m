% Systems Identification 
% Example 2.5 of the book System Identification: Theory for the User (2nd Edition).

%clear all
%close all

t = 1:1000;
e = randn(1,1000); %white noise with variance 1 and mean zero

figure(1)
plot(t,e)
hold on

v = zeros(1, size(t,2));
for index_t=3:size(t,2)%v(t)-1.5v(t-1)+0.7v(t-2)= e(t)+0.5e(t-1)
    v(index_t) = 1.5*v(index_t-1)-0.7*v(index_t-2)+e(index_t)+0.5*e(index_t-1);
end

figure(2)
plot(t,v)
hold off

w=linspace(0,pi,size(t,2));
for index_w=1:size(w,2)
    En(index_w)=sum(e.*exp(-i*w(index_w)*t));
end
En = En/sqrt(size(t,2)); %E(s)is the Fourier Transform of the e(t)

for index_w=1:size(w,2)
    Vn(index_w)=sum(v.*exp(-i*w(index_w)*t));
end
Vn = Vn/sqrt(size(t,2)); %V(s) is the Fourier Transform of the v(t)

Hn = Vn./En; %H(s)=V(s)/E(s) is the Transfer Function
Phi_v = var(e)*abs(Hn).^2; % Phi_v(w) = var(e)*|H(w)|^2 is the spectra of v(t)
plot(w,Phi_v)

