%Omid55
function [ data,net,initOp,soc_pow ] = GenerateDataset( N,MaxIteration,K,PowerfullStrength,DE_config,netType,soc_type )

%% Initializing the opinions
[net, op] = InitializeOpinions(N,K,DE_config,netType);
initOp = op;

% update number of nodes
N = size(net, 1);


%% Calculating the social powers
soc_pow = CalculateSocialPowers(net,PowerfullStrength,soc_type)';


%% Simulation Body
%Voter Model with social power
data = zeros(MaxIteration,3);
cnt = 1;
gamma = 0.7;
xs = op;
%newOp = zeros(size(op));
for it = 1 : MaxIteration
    %v = randi(N);
    for v = 1 : N
        adj = find(net(v,:) == 1, 1);
        if isempty(adj)
            continue;
        end
        prob = zeros(K,1);
        for k = 1 : K
            compliant = adj(op(adj) == k);
            
%             if op(v) == k
%                 self = (K*(1-gamma)+gamma) / K; 
%             else
%                 self = gamma / K;
%             end
%             prob(k) = self * sum(soc_pow(compliant)) / sum(soc_pow(adj));

            prob(k) = sum(soc_pow(compliant)) / sum(soc_pow(adj));
        end
        newOpinion = RouletteSelect(prob);
        op(v) = newOpinion;   %%%newOp(v) = newOpinion;
        data(cnt,:) = [v cnt newOpinion];   %%%data(cnt,:) = [v it newOpinion];
        cnt = cnt + 1;
    end
    %%%op = newOp;
    
    %if mod(it,N) == 0
        xs = [xs; op];
    %end
end
xs = [xs; op];

% figure;
% title('Opinion Changing');
% xlabel('Iterations');
% ylabel('Opinion');
% for i=1:N
%     plot(xs(:,i));
%     hold on;
% end

%save('OpinionDataset','data','net','initOp','soc_pow');


end

