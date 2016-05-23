%Omid55
function [ fitness ] = ObjectiveFunction( x,sp )

alpha = 0.001;
fitness = 0;
M = length(x);
for i = 1 : M
    fitness = fitness + sum(sp(i,:) .* ( alpha - (x(i)~=x) ) - log( 1 + exp( alpha - (x(i)~=x) ) ));
end

end

