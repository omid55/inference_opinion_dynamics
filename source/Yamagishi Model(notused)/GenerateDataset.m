%Omid55
function [  ] = GenerateDataset(  )

%% Parameters
N = 100;                        % number of nodes
J = 5;                            % number of attributes of each person
MaxIteration = 1;     % iteration count
K = 10;                             % number of opinions
tetta = [0.1; 0.3; 0.5; 0.01; 0.09];   % its size is J
x = rand(N,J);
op = randi(K,[N,1]);
initOp = op;
originalTetta = tetta;


%% Social Network
% BA
averageDegree = 12;
net = BAnet(N,averageDegree/2,averageDegree/2);


%% Simulation Body
%Voter Model with social power
soc_pow = exp(x*tetta);
data = zeros(MaxIteration,3);
cnt = 1;
for it = 1 : MaxIteration
   for v = 1 : N
        adj = find(net(v,:) == 1);
        prob = zeros(K,1);
        for k = 1 : K
            compliant = adj(op(adj) == k);
            prob(k) = sum(soc_pow(compliant)) / sum(soc_pow(adj));
        end
        newOpinion = RouletteSelect(prob);
        op(v) = newOpinion;
        data(cnt,:) = [v cnt op(v)];
        cnt = cnt + 1;
   end
end

save('OpinionDataset.mat','data','x','net','initOp','originalTetta','N');


end

