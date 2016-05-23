%Omid55
%Forest Fire Network Creator Code
%gives N for number nodes and PF for Forward burning probability and PB for Backward burning probability.
function [ sp ] = ForestFireCreator( N,PF,PB )

initCliqueSize = 3;
sp = sparse(ones(initCliqueSize) - diag(ones(initCliqueSize,1)));   % clique creator

for i=initCliqueSize+1:N
    ambassador = 1 + floor( rand(1) * (i-1) );
    sp(i,ambassador) = 1;
    sp(ambassador,i) = 1;
    nodes = getRecursivelyNodes(ambassador,sp,PF,PB);
    nodes(nodes == i)=[];
    sp(i,nodes) = 1;
    sp(nodes,i) = 1;
end

sp = max(sp,sp');

end

