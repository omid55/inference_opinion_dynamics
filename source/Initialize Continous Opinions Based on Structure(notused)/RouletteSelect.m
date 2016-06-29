%Omid55
function [ selectedIndex ] = RouletteSelect( probabilities )

r = rand(1,1) * sum(probabilities);
s = 0;
j = 0;

while s <= r
    j = j + 1;
    s = s + probabilities(j);
end

selectedIndex = j;


end

