%Omid55
%Initialize Opinions based on latent space model
function [ net,initOp ] = InitializeOpinions( N,K,DE_config,netType )

%% Social Network
averageDegree = 12;
switch netType
    case 1
    %BA
    net = BAnet(N,averageDegree/2,averageDegree/2);
    name = 'BA';

    case 2
    %ER
    erdosProbability = averageDegree/N;
    net = erdos_reyni(N,erdosProbability);
    name = 'ER';
    
    case 3
    %WS
    net = WattsStrogatzCreator(N,averageDegree,0.1);
    name = 'WS';
    
    case 4
    %FF
    net = ForestFireCreator(N,0.25,0.25);
    name = 'FF';
end
disp([name '('  num2str(N) ')']);


%% Init Opinions
% -- init by random --
initOp = randi(K,1,N);

% % -- init by optimization --
% initOp = ExpectOpinions(net,K,DE_config);
%save('InitOpinions','initOp','net');


% JUST A SAMPLE
% % My Simple Network
% ed = [1 2;1 3;1 4;1 5;1 7;5 6;4 6;6 7;7 8;9 8;9 13;9 10;9 11;9 12;11 12;10 11;10 12];
% N = max(max(ed));
% sp = sparse(ed(:,1),ed(:,2),1,N,N);

% Simulation Part
%x = ExpectOpinions(sp);
%ViewMyGraph(sp,x);
% JUST A SAMPLE

end

