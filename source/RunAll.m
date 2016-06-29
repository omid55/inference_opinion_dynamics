%Omid55
function [  ] = RunAll(  )

% addpath(genpath('/home/omid55/Tools/cvx'));
% warning off
% addpath(genpath('~/Tools/matlab_bgl'));
% addpath('Initialize Discrete Opinions Based on Structure');

%% All and All
% for net_type = 1 : 4
%     for soc_type = 1 : 5
%         maincvx(net_type,soc_type);
%     end
% end

% fig = figure;
% subplot(2,2,1);
% maincvx(1,1);
% subplot(2,2,2);
% maincvx(2,1); 
% subplot(2,2,3);
% maincvx(3,1);
% subplot(2,2,4);
% maincvx(4,1);
% saveas(fig,'Figure1.fig');
% print('Figure1.png', '-dpng', '-r800');
% close all;


fig = figure;
for c = 1 : 2
    subplot(1,2,c);
    maincvx(c+4,1);
end
saveas(fig,'Figure2.fig');
print('Figure2.png', '-dpng', '-r800');

end