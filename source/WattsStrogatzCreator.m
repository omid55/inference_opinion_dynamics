%Omid55
function [ sp ] = WattsStrogatzCreator( N,K,P )

probability = 1 - P;
sp = CreateRegularLattice(N,K);
for i = 1 : N
    adj = find(sp(i,:) == 1);
    available = setdiff(1:N,[adj i]);
    pr = rand(length(adj),1);
    selected = find(pr > probability);
    index = 1 + floor(rand(length(selected),1)*length(available));
    sp(i,selected) = 0;
    sp(selected,i) = 0;
    sp(i,available(index)) = 1;
    sp(available(index),i) = 1;
end

end

