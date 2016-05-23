%Omid55
function [ nodes ] = OuterNodes( sp,nodeIndex )

% nodes = [];
% for i=1:size(sp,1)
%     if(sp(nodeIndex,i) ~= 0)
%         nodes = [nodes i];
%     end
% end

adjSp = sp(nodeIndex,:);
nodes = find(adjSp>0);

end

