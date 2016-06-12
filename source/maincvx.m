% Omid55
% Learn Parameter tetta by the opinion data with cvx
%   With Bayesian Learning
function [  ] = maincvx( net_type,soc_type )

%% net_type:
%1- BA
%2- ER
%3- WS
%4- FF

%% soc_type:
%1- high degree apporach
%2- random approach
%3- high closeness apporach
%4- high betweenness apporach
%5- low degree apporach


%% Clear everything
%clear;
%---for subplot---%close all;
clc;


%% ----------------===== All Parameters =====----------------
% ---------------------=====================================---------------------
N = 100;                                  % number of nodes
% N = 20;                                  % << SIMPLE >>
% net_type = 1;
% soc_type = 1;
MaxDatasetIteration = 10;        % iteration count
K = 5;                                       % number of opinions
times = 10;
gamma = 0.7;
%desiredOp = 1;
MaxSimulationIteration = 50;
InformedAgentsSize = ceil(N/10);
MaxOptimizationIteration = 0;       % 0 means just one single iteration
InformedAgentsStrength = 1000;   % for betweenness make than 10 times bigger
PowerfullStrength = InformedAgentsStrength/2;
DE_config.N = 100;
DE_config.MaxIteration = 20;
DE_config.Coef = 0.7;
DE_config.Beta = 0.8;
DE_config.Pr = 0.6;
DE_config.Nv = 1;   % number of difference vectors
DE_config.withFigure = 0;    % DE optimization without any figure
% ---------------------=====================================---------------------


%% Generating dataset
[data,net,initOp,soc_pow] = GenerateDataset(N,MaxDatasetIteration,K-1,PowerfullStrength,DE_config,net_type,soc_type);
desiredOp = K;

% %% Loading the dataset
% load('OpinionDataset');

% dd = [1 2; 1 3;2 3; 2 4;3 4];
% N = max(max(dd(:,1)),max(dd(:,2)));
% net = sparse(dd(:,1),dd(:,2),1,N,N);
% net = max(net,net');
% initOp = [1 1 0 0];
% data = [1 2 0;4 3 1;3 4 1;1 5 1];
% %data = [1 2 0;4 3 0;2 4 0];

% dd = [1 2; 1 3;2 3];
% N = max(max(dd(:,1)),max(dd(:,2)));
% net = sparse(dd(:,1),dd(:,2),1,N,N);
% net = max(net,net');
% initOp = [0 0 1];
% data = [1 1 1;2 2 1];

% dd = [1 2; 1 3;2 3;2 4;2 5;4 5];
% N = max(max(dd(:,1)),max(dd(:,2)));
% net = sparse(dd(:,1),dd(:,2),1,N,N);
% net = max(net,net');
% initOp = [1 1 0 0 0];
% data = [1 1 0;4 2 1;5 3 1];


%% Learning the parameter by cvx optimizer
loglikelihoods = [];
data = sortrows(data,2);  % sorting data through time steps
ahat = full(sum(net)');    %ones(N,1);
ahat = ahat / max(ahat);
%Eps = 10^-5;
cnt = 1;
%figure;

degz = full(sum(net));

while 1
    cvx_begin
    %cvx_precision best
    variable a(N,1)
    
    op = initOp;
    newOp = op;
    LL = 0;
    for i = 1 : size(data,1)
        v = data(i,1);
        t = data(i,2);
        k = data(i,3);
        
        lastOpV = op(v);
        newOp(v) = k;
        
        if i > 1 && t ~= data(i-1,2)
            op = newOp;        % updating the opinion of node v
        end
        
        adjV = find(net(v,:) == 1);
        compliant = adjV(op(adjV) == k);
        su1 = 0;
        if ~isempty(compliant)
            su1 = sum( exp(ahat(compliant)) .* a(compliant) )  /  sum( exp(ahat(compliant)) );
        end
        others = setdiff(adjV,compliant);
        if isempty(others)
            continue;
        end
        su2 = sum(exp(a(adjV)));
        h = adjV(degz(adjV)==max(degz(adjV)));
        prior = ahat(h(1)) - log(sum(exp(ahat(adjV))));
        
        if lastOpV == k
            self = (K*(1-gamma)+gamma) / K;
        else
            self = gamma / K;
        end
        
        LL = LL + (prior + log(self) + su1 - log(su2));
    end
    
    maximize LL
%    % we do not need this constraint
%     subject to
%         a > 0;
    cvx_end
    
    save('CvxFinishedData');
    
    loglikelihoods = [loglikelihoods cvx_optval];
    %save('LoglikelihoodsData','loglikelihoods');
    
    %     diff = norm(ahat - a)
    %     if diff < Eps
    %         ALPHA = ahat
    %         break;
    %     end
    
    SP = exp(a);
    [r,p] = corr(soc_pow,SP)
    %[r,p] = corr(soc_pow,SP,'type','Spearman')
    
    ahat = a;
    
    if cnt > MaxOptimizationIteration || (cnt > 1  &&  cvx_optval < loglikelihoods(end))
        break;
    end
    
    cnt = cnt + 1;
end

if length(loglikelihoods) > 1
    figure;
    plot(1:length(loglikelihoods),loglikelihoods);
    title('\bfLogLikelihood');
end

ranks = HITS(net);
[hitsR,hitsP] = corr(ranks,SP,'type','Spearman')


%% ----------------------------------------------------------------------------------
%% ----------------------------------------------------------------------------------
%% Voter simulation with informed agents

% % % % by proposed
% % % newSP = zeros(N,1);
% % % for i = 1 : N
% % %     adjI = find(net(i,:) == 1);
% % %     tmp = 0;
% % %     for j = 1 : length(adjI)
% % %         tmp = tmp + mean(SP(net(adjI(j),:) == 1));
% % %     end
% % %     newSP(i) = SP(i) + tmp / length(adjI);
% % % end
% % % SP = newSP;
% % % sortedSPs = sort(SP,'descend');
% % % targets = [];
% % % c = 1;
% % % while length(targets) < InformedAgentsSize
% % %     goods = find(SP == sortedSPs(c));
% % %     targets = unique([targets (goods(op(goods) ~= desiredOp))]);
% % %     c = c + 1;
% % % end
% % % targets = targets(1:InformedAgentsSize);   % for removing excess targets
% % % by proposed
% % newSP = zeros(N,1);
% % for i = 1 : N
% %     newSP(i) = SP(i) + mean(SP(net(i,:) == 1));
% % end
% % SP = newSP;
% % sortedSPs = sort(SP,'descend');
% % targets = [];
% % c = 1;
% % while length(targets) < InformedAgentsSize
% %     goods = find(SP == sortedSPs(c));
% %     targets = unique([targets (goods(op(goods) ~= desiredOp))]);
% %     c = c + 1;
% % end
% % targets = targets(1:InformedAgentsSize);   % for removing excess targets
% by proposed
sortedSPs = sort(SP,'descend');     % sortedSPs = sort(SP,'descend');  % <<<<<<<<<<<  << CHECK HERE >>  >>>>>>>>>>>
targets = [];
c = 1;
while c<N && length(targets) < InformedAgentsSize
    goods = find(SP == sortedSPs(c));    % goods = find(SP == sortedSPs(c));  % <<<<<<<<<<<  << CHECK HERE >>  >>>>>>>>>>>
    targets = unique([targets (goods(op(goods) ~= desiredOp))]);
    c = c + 1;
end
if length(targets) < InformedAgentsSize
    targets = GiveGoods(SP,InformedAgentsSize);
end
targets = targets(1:InformedAgentsSize);   % for removing excess targets


% by high degs
degs = sum(net);
sortedDegs = sort(degs,'descend');
highDegs = [];
c = 1;
while c<N && length(highDegs) < InformedAgentsSize
    goods = find(degs == sortedDegs(c));
    highDegs = unique([highDegs (goods(op(goods) ~= desiredOp))]);
    c = c + 1;
end
if length(highDegs) < InformedAgentsSize
    highDegs = GiveGoods(degs,InformedAgentsSize);
end
highDegs = highDegs(1:InformedAgentsSize);   % for removing excess targets
% % % highDegs = GiveGoods(degs,InformedAgentsSize);


%by betweenness
betw = betweenness_centrality(net)';
sorted = sort(betw,'descend');
highBetweenness = [];
c = 1;
while c<N && length(highBetweenness) < InformedAgentsSize
    goods = find(betw == sorted(c));
    highBetweenness = unique([highBetweenness (goods(op(goods) ~= desiredOp))]);
    c = c + 1;
end
if length(highBetweenness) < InformedAgentsSize
    highBetweenness = GiveGoods(betw,InformedAgentsSize);
end
highBetweenness = highBetweenness(1:InformedAgentsSize);   % for removing excess targets
% % % highBetweenness = GiveGoods(betw,InformedAgentsSize);


%by closeness
asp = all_shortest_paths(net);
clos = (N-1) ./ sum(asp);
sorted = sort(clos,'descend');
highCloseness = [];
c = 1;
while c<N && length(highCloseness) < InformedAgentsSize
    goods = find(clos == sorted(c));
    highCloseness = unique([highCloseness (goods(op(goods) ~= desiredOp))]);
    c = c + 1;
end
if length(highCloseness) < InformedAgentsSize
    highCloseness = GiveGoods(clos,InformedAgentsSize);
end
highCloseness = highCloseness(1:InformedAgentsSize);   % for removing excess targets
% % % highCloseness = GiveGoods(clos,InformedAgentsSize);


% % by random
% randTar = randi(N,[1,InformedAgentsSize]);    % << CHECK HERE >>  REMOVE THIS AND UNCOMMENT BELLOW

% by my influence maximization's metric:
EPN1_measure = zeros(N,1);
indegs = sum(net,1);
outdegs = sum(net,2);
eps = 0.0001;
for i = 1 : N
    f1 = ((outdegs(i) + eps) / (indegs(i) + eps));
    f2 = (mean(indegs(net(i,:) == 1)) / max(indegs));
    EPN1_measure(i) = f1 * f2;
end
[~,EP1] = sort(EPN1_measure,'descend');
EP1_targets = EP1(1:InformedAgentsSize)';
randTar = EP1_targets;


%% Influence maximization simulations
randD = zeros(times,MaxSimulationIteration+1);
degsD = zeros(times,MaxSimulationIteration+1);
mineD = zeros(times,MaxSimulationIteration+1);
closD = zeros(times,MaxSimulationIteration+1);
betwD = zeros(times,MaxSimulationIteration+1);

for time = 1 : times
    time
    
    %     % by random
    %     randTar = randi(N,[1,InformedAgentsSize]);
    
    resRand = VoterSimulationMethod(op,net,randTar,MaxSimulationIteration,K,desiredOp,soc_pow,InformedAgentsStrength);
    resDeg = VoterSimulationMethod(op,net,highDegs,MaxSimulationIteration,K,desiredOp,soc_pow,InformedAgentsStrength);
    resBetweenness = VoterSimulationMethod(op,net,highBetweenness,MaxSimulationIteration,K,desiredOp,soc_pow,InformedAgentsStrength);
    resCloseness = VoterSimulationMethod(op,net,highCloseness,MaxSimulationIteration,K,desiredOp,soc_pow,InformedAgentsStrength);
    resMine = VoterSimulationMethod(op,net,targets,MaxSimulationIteration,K,desiredOp,soc_pow,InformedAgentsStrength);
    
    randD(time,:) = resRand.desOpPop;
    degsD(time,:) = resDeg.desOpPop;
    betwD(time,:) = resBetweenness.desOpPop;
    closD(time,:) = resCloseness.desOpPop;
    mineD(time,:) = resMine.desOpPop;
end


%% Result Figure
%---for subplot---%fig = figure;
%---for subplot---%set(gcf, 'PaperPosition',[0.25 2.5 5 3.5]);
hold on;
times = (1:size(randD,2)) * N;
myScale = 10;    % << CHECK HERE >> << SET IT TO 1 >>
errorbar(times,myScale*mean(randD),myScale*std(randD)/sqrt(size(randD,1)),'r');
errorbar(times,myScale*mean(degsD),myScale*std(degsD)/sqrt(size(degsD,1)),'b');
errorbar(times,myScale*mean(betwD),myScale*std(betwD)/sqrt(size(betwD,1)),'Color',[1 0.69 0.39]);
errorbar(times,myScale*mean(closD),myScale*std(closD)/sqrt(size(closD,1)),'c');
errorbar(times,myScale*mean(mineD),myScale*std(mineD)/sqrt(size(mineD,1)),'g');
leg = legend('Random','High Degree','High Betweenness','High Closeness','Proposed','Location','Southeast');
set(leg,'FontSize',7);
%axis([0,times(end)+1000,myScale*(min(mean(randD))-20),myScale*(N+InformedAgentsSize+20)]);
%axis('tight');
xlabel('\bfTimestep');
ylabel('\bfDesired Opinion Population');
%set(gca,'XScale','log');
%set(gca,'YScale','log');
% % % saveas(fig,'ErrorBar.fig');
% % % print -loos -dtiff ErrorBar.tiff;

%     figure;
%     plot(1:length(resMine.desOpPop),resMine.desOpPop,'-o',1:length(resRand.desOpPop),resRand.desOpPop,'-x',1:length(resDeg.desOpPop),resDeg.desOpPop,'-p');
%     legend('proposed','rand','degree');


%% Save all of data
path = 'NewResults';
mkdir(path);
f1 = num2str(net_type);
f2 = num2str(soc_type);
newpath = [path '/' f1 '_' f2];
mkdir(newpath);
%---for subplot---%name = [newpath '/ErrorBar'];
%---for subplot---%saveas(fig,[name '.fig']);
%---for subplot---%eval(['print -loos -dtiff ' name '.tiff']);
save([newpath '/FinalData']);


end