%Omid55
function [ sp ] = BAnet( N,m,m0 )

sp = sparse(ones(m0) - diag(ones(m0,1)));   % clique creator
probablity = repmat(1:m0,1,m0);

n = m0;
while n < N
    lastProbability = probablity;
    n = n + 1;
    selected = zeros(m,1);
    for i = 1 : m
        index = randi(length(probablity));
        selected(i) = probablity(index);
        probablity(probablity == probablity(index)) = [];
    end
    for i=1:length(selected)        
        sp(selected(i),n) = 1;
        sp(n,selected(i)) = 1;    % for an undirected network
    end
    probablity = [lastProbability selected' ones(1,length(selected))*n];
end

end

