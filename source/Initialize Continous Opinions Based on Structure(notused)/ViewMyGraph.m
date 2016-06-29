%Omid55
%View Graph based on bio graph design
function [  ] = ViewMyGraph( sp,lables )

N = size(sp,1);
labs=cell(N,1);
for ag=1:N
    if nargin > 1
        labs{ag} = [num2str(ag) ': ' num2str(lables(ag))];
    else
       labs{ag} = num2str(ag); 
    end
end
bg = biograph(sp,labs);
view(bg);
pause;

% close all biograph windows
child_handles = allchild(0);
names = get(child_handles,'Name');
k = find(strncmp('Biograph Viewer', names, 15));
close(child_handles(k))

end

