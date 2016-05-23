%Omid55
%Give count from the best measures
function [ result ] = GiveGoods( measures,count )

sorted = sort(measures,'descend');
threshold = sorted(count);
result = find(measures > threshold);
indices = find(measures == threshold);
if size(measures,1) > 1     % column major
    result = [result; indices(1:count-length(result))];   
else                                % row major
    result = [result indices(1:count-length(result))];   
end

end

