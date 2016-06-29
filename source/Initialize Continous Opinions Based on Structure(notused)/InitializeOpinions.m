%Omid55
%Initialize Opinions based on latent space model
function [ net,initOp ] = InitializeOpinions( N )

tic;
%% Social Network
% BA
averageDegree = 12;
A = BAnet(N,averageDegree/2,averageDegree/2);
net = sparse(A);
initOp = ExpectOpinions(net);
save('InitOpinions','initOp','net');
toc
% % My Simple Network
% ed = [1 2;1 3;1 4;1 5;1 7;5 6;4 6;6 7;7 8;9 8;9 13;9 10;9 11;9 12;11 12;10 11;10 12];
% N = max(max(ed));
% sp = sparse(ed(:,1),ed(:,2),1,N,N);


%% Simulation Part
%x = ExpectOpinions(sp);
%ViewMyGraph(sp,x);

% hist(initOp);


end

