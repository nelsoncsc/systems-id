function [ z ] = function_Hc( q )
    z =  (0.254*q + 1.82)./(q.^2 + 0.3567*q + 0.2426);
end
