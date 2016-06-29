%Omid55
function [ fitness ] = ObjectiveFunction( x,sp,lx )

alpha = 0.001;
fitness = 0;
M = length(x);
% %calculating the P matrix
% P = zeros(M,M);
% for i=1:M
%     for j=1:M
%         if i == j
%             continue;
%         end
%         
%         P(i,j) = 1 / (1 + exp(abs(x(i) - x(j)) - alpha));
%         if sp(i,j) > 0
%             fitness = fitness + log(P(i,j));
%         else
%             fitness = fitness + log(1 - P(i,j));
%         end
%     end
% end

for i = 1 : M
    fitness = fitness + sum(sp(i,:) .* ( alpha - abs(x(i) - x) ) - log( 1 + exp( alpha - abs(x(i) - x) ) ));
end

if nargin > 2
    const = 0.01;
    sigma = 1;
    notFarPart = 0;
    for i=1:M
        notFarPart = notFarPart - ((x(i) - lx(i)) ^ 2 / (2*sigma^2)) + const;
    end
    fitness = fitness + notFarPart;
end

end

