%Omid55
function [ sp ] = BAnet( N,m,m0 )

sp = sparse(ones(m0) - diag(ones(m0,1)));   % clique creator
probablity = repmat(1:m0,1,m0);

n = m0;
while n < N
    n = n + 1;
    selected = [];
    while length(selected) < m
        if(isempty(probablity))
            selected = 1;
        else
            prob = probablity;
            for i=1:length(selected)
                inds = find(prob == selected(i));
                prob(inds) = [];
            end
            index = randi(length(prob));
            selected = [selected; prob(index)];
            selected = unique(selected);
        end
    end
    for i=1:length(selected)        
        sp(selected(i),n) = 1;
        sp(n,selected(i)) = 1;
    end
    probablity = [probablity selected' ones(1,length(selected))*n];
end

% sp = max(sp,sp');

end

