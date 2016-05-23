%Omid55
function [ nodes ] = getRecursivelyNodes( ambassador,sp,PF,PB )

fAdj = OuterNodes(sp,ambassador);
bAdj = InnerNodes(sp,ambassador);
nodes = [];

for j=1:length(fAdj)
    if  fAdj(j)~=ambassador && isempty(find(nodes==fAdj(j), 1)) && rand(1)>(1-PF)
        nodes = [nodes fAdj(j) getRecursivelyNodes(fAdj(j),sp,PF*0.1,PB*0.1)];  %nodes = [nodes fAdj(j) getRecursivelyNodes(fAdj(j),sp,PF*(1-PF),PB*(1-PB))];
    end
end

for j=1:length(bAdj)
    if  bAdj(j)~=ambassador && isempty(find(nodes==bAdj(j), 1)) && rand(1)>(1-PB)
        nodes = [nodes bAdj(j) getRecursivelyNodes(bAdj(j),sp,PF*0.1,PB*0.1)];  %nodes = [nodes bAdj(j) getRecursivelyNodes(bAdj(j),sp,PF*(1-PF),PB*(1-PB))];
    end
end

%%%nodes(find(nodes == ambassador))=[];
%nodes = unique(nodes);

end

