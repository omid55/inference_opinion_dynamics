%% Model Networks
clc;
close all;
clear;

fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 8 5]);
set(gcf, 'PaperSize',[8.5 11]);
my_cnt = 1;
my_ch = 'a';
ppath = 'RESULTS/';
lw = 1;

for num = 1 : 4
    subplot(2,2,num);
    load([ppath num2str(num) '_1/FinalData.mat']);
    
    hold on;
    times = (1:size(randD,2)) * N;
    myScale = 1/(N + InformedAgentsSize);
    errorbar(times,myScale*mean(degsD),myScale*std(degsD)/sqrt(size(degsD,1)),'b','LineWidth',lw);
    errorbar(times,myScale*mean(betwD),myScale*std(betwD)/sqrt(size(betwD,1)),'Color',[1 0.69 0.39],'LineWidth',lw);
    errorbar(times,myScale*mean(closD),myScale*std(closD)/sqrt(size(closD,1)),'c','LineWidth',lw);
    errorbar(times,myScale*mean(randD),myScale*std(randD)/sqrt(size(randD,1)),'r','LineWidth',lw);
    errorbar(times,myScale*mean(rankD),myScale*std(rankD)/sqrt(size(rankD,1)),'k','LineWidth',lw);
    errorbar(times,myScale*mean(anotherD),myScale*std(anotherD)/sqrt(size(anotherD,1)),'m','LineWidth',lw);
    errorbar(times,myScale*mean(mineD),myScale*std(mineD)/sqrt(size(mineD,1)),'g','LineWidth',lw);
    if(num == 1)
        leg = legend('Degree','Betweenness','Closeness','EPN','PageRank','Yamagishi','HESP','Location','Southeast');
        set(leg,'FontSize',7);
    end
    t = title(['\bf(' char(my_ch) ')']);
    set(t, 'FontSize', 10);
    my_ch = my_ch + 1;
    axis([0 10500 0 1.1]);
end

[~,h3] = suplabel('\bfRatio of Followers','y');
set(h3,'FontSize',12);
[~,h3] = suplabel('\bfIteration #','x');
set(h3,'FontSize',12);

print('Figure2.png', '-dpng', '-r200');
print('Figure2_large.png', '-dpng', '-r800');
saveas(fig,'Figure2.fig');


% -----------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------
% -----------------------------------------------------------------------------------------------


%% Real Networks
clc;
close all;
clear;

fig = figure;
set(gcf, 'PaperPosition',[0.25 2.5 12 9]);
set(gcf, 'PaperSize',[8.5 11]);
my_cnt = 1;
my_ch = 'a';
ppath = 'RESULTS/';
lw = 1;

for num = 5 : 10
    subplot(2,3,num-4);
    load([ppath num2str(num) '_1/FinalData.mat']);
    
    hold on;
    times = (1:size(randD,2)) * N;
    myScale = 1/(N + InformedAgentsSize);
    errorbar(times,myScale*mean(degsD),myScale*std(degsD)/sqrt(size(degsD,1)),'b','LineWidth',lw);
    errorbar(times,myScale*mean(betwD),myScale*std(betwD)/sqrt(size(betwD,1)),'Color',[1 0.69 0.39],'LineWidth',lw);
    errorbar(times,myScale*mean(closD),myScale*std(closD)/sqrt(size(closD,1)),'c','LineWidth',lw);
    errorbar(times,myScale*mean(randD),myScale*std(randD)/sqrt(size(randD,1)),'r','LineWidth',lw);
    errorbar(times,myScale*mean(rankD),myScale*std(rankD)/sqrt(size(rankD,1)),'k','LineWidth',lw);
    errorbar(times,myScale*mean(anotherD),myScale*std(anotherD)/sqrt(size(anotherD,1)),'m','LineWidth',lw);
    errorbar(times,myScale*mean(mineD),myScale*std(mineD)/sqrt(size(mineD,1)),'g','LineWidth',lw);
    if(num == 5)
        leg = legend('Degree','Betweenness','Closeness','EPN','PageRank','Yamagishi','HESP','Location','Southeast');
        set(leg,'FontSize',10);
    end
    t = title(['\bf(' char(my_ch) ')']);
    set(t, 'FontSize', 10);
    my_ch = my_ch + 1;
end

[~,h3] = suplabel('\bfRatio of Followers','y');
set(h3,'FontSize',15);
[~,h3] = suplabel('\bfIteration #','x');
set(h3,'FontSize',15);

print('Figure3.png', '-dpng', '-r200');
print('Figure3_large.png', '-dpng', '-r800');
saveas(fig,'Figure3.fig');

