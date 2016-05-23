%Omid55
% Get n distinct items from range 
function [ result ] = GetDistinctItems( list,n )

result = zeros(n,1);
for i=1:n
    index = randi(length(list));
    result(i) = list(index);
    list(index) = [];
end

end

