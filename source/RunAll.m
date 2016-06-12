%Omid55
function [  ] = RunAll(  )


%% All and All
% for net_type = 1 : 4
%     for soc_type = 1 : 5
%         maincvx(net_type,soc_type);
%     end
% end


fig = figure;
subplot(2,2,1);
maincvx(1,1);
subplot(2,2,2);
maincvx(2,1);
subplot(2,2,3);
maincvx(3,1);
subplot(2,2,4);
maincvx(4,1);
saveas(fig,'Figure1.fig');
print('Figure1.png', '-dpng', '-r400');
close all;

fig = figure;
subplot(2,3,1);
maincvx(1,1);
subplot(2,3,2);
maincvx(1,2);
subplot(2,3,3);
maincvx(1,3);
subplot(2,3,4);
maincvx(1,4);
subplot(2,3,5);
maincvx(1,5);
saveas(fig,'Figure2.fig');
print('Figure2.png', '-dpng', '-r400');

end