%Omid55
function [ value ] = StayInBound( value,rangeBegin,rangeEnd )

if value > rangeEnd
    value = rangeEnd;
end

if value < rangeBegin
    value = rangeBegin;
end
    
end

