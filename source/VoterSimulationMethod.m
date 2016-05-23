%Omid55
% Voter model simulation method with social power
function [ result ] = VoterSimulationMethod( op,net,targets,MaxIteration,K,desiredOp,soc_pow,InformedAgentsStrength )

N = size(net,1);

%% Network Edge Adding
for i=1:length(targets)
    net(targets(i),N+i) = 1;
    net(N+i,targets(i)) = 1;
end
op(N+1:N+length(targets)) = desiredOp;   % agents' opinion

soc_pow(N+1:N+length(targets)) = InformedAgentsStrength;

xs = op(1:N);
newOp = zeros(1,N);
desOpPop = zeros(1,MaxIteration+1);   %desOpPop == desired opinion population
desOpPop(1) = sum(op == desiredOp);
prob = zeros(K,1);
for it = 1 : MaxIteration
    for v = 1 : N        
        adj = find(net(v,:) == 1);
        for k = 1 : K
            compliant = adj(op(adj) == k);
            prob(k) = sum(soc_pow(compliant)) / sum(soc_pow(adj));
        end
        newOpinion = RouletteSelect(prob);
        newOp(v) = newOpinion;
    end
    
    xs = [xs; newOp];
    desOpPop(it+1) = sum(op == desiredOp);

    op(1:N) = newOp;
%     ViewMyGraph(net,op);
end
xs = [xs; newOp];

% figure;
% title('Opinion Changing');
% xlabel('Iterations');
% ylabel('Opinion');
% for i=1:N
%     plot(xs(:,i));
%     hold on;
% end

result.lastOp = op;
result.desOpPop = desOpPop;

end
