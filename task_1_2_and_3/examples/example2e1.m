%Nelson Campos
% Systems Identification 
% Example 2E.1 of the book System Identification: Theory for the User (2nd Edition).

Phi_v = [1 1; 1.64 1.6]; %The spectra of the process v(t)
coeff= Phi_v(:,2); %coefficients of the R(q) = (1+c1*q^-1)/(1+a1*q^-1)

R = tf([1 0.5*coeff(1)], [1 0.5*coeff(2)],1,'Variable', 'q^-1')
