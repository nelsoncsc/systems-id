%Nelson Campos
% Systems Identification 
% Example 3.1 of the book System Identification: Theory for the User (2nd Edition).

clc
close all
clear all

t = 0:1:4;

e = [1 2 3 4 5];%randn(1, size(t,2)); %definition of e(t)
c = 2;
v = e;

%definition of v(t)
for k=2:size(t,2)
    v(k) = e(k)+c*e(k-1);
end

e2 = zeros(1, size(t,2)); %e2 is e(t) retrieved from v(t)
%This part makes the convolution
tau = 0:size(t,2);
for index_tau = 1: size(t,2);
    for k=1:size(t,2)
        if(k-tau(index_tau) >=1 && k-tau(index_tau) <= size(t,2))
            v2(k) = v(k-tau(index_tau));
        else
            v2(k) = 0;    
        end
    end
    e2 = e2 + (-c)^tau(index_tau)*v2; %e(t) is retrieved here
end

%Then the one step ahead predictor is v(t|t-1) = v(t)-e(t)
v2 = v-e;
