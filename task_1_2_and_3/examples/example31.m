%Nelson Campos
% Systems Identification 
% Example 3.1 of the book System Identification: Theory for the User (2nd Edition).

t = 1:5;
e = [1 2 3 4 5]; %randn(1, size(t,2))
a = 0.5;

e2 = zeros(1, size(t,2)); %e(t-k)
tau = 0:size(t,2);
v = zeros(1,size(t,2)); %v(t)
for index_tau = 1: 5;
    for k=1:size(t,2)
        if(k-tau(index_tau) >=1&& k-tau(index_tau) <= size(t,2))
            e2(k) = e(k-tau(index_tau));
        else
            e2(k) = 0;    
        end
    end
    v = v + a^tau(index_tau)*e2; %e(t) is retrieved here
end

%Then the one step ahead predictor is v(t|t-1) = v(t)-e(t)
v2 = v-e;
