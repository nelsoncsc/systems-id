function [ Window ] = Windows(gamma, theta, type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
switch type
    case 1 %Parzen
        Window = 2/(pi*gamma^3)*(2+cos(theta))*(sin(gamma*theta/4)/sin(theta/2))^4;
    case 2 %Barlett
        Window = 1/(2*gamma*pi)*((sin(gamma*theta/2))/(sin(theta/2)))^2;
    case 3 %Hamming
        H1 = 1/(4*pi)*(sin((gamma+0.5)*theta))/(sin(theta/2));
        H2 = 1/(8*pi)*(sin((gamma+0.5)*(theta-(pi/gamma))))/(sin((theta-(pi/gamma))/2));
        H3 = 1/(8*pi)*(sin((gamma+0.5)*(theta+(pi/gamma))))/(sin((theta+(pi/gamma))/2));
        Window = H1+H2+H3;
    otherwise
        Window = 1;
end

end

