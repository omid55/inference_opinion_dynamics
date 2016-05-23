%Omid55
% Voter model simulation method with social power
function [ result ] = VoterSimulationMethod( op,net,targets,MaxIteration,K,desiredOp,soc_pow )

N = size(net,1);
op(targets) = desiredOp;   % desired opinion

xs = op(1:N);
newOp = zeros(1,N);
desOpPop = zeros(1,MaxIteration+1);   %desOpPop == desired opinion population
desOpPop(1) = sum(op == desiredOp);   %meanOps(1) = mean(op);
prob = zeros(K,1);
for it = 1 : MaxIteration
    for v = 1 : N
        if isempty(find(targets == v, 1))    % we suppose that initial persons will never change their minds
            
            adj = find(net(v,:) == 1);
            for k = 1 : K
                compliant = adj(op(adj) == k);
                prob(k) = sum(soc_pow(compliant)) / sum(soc_pow(adj));
            end
            newOpinion = RouletteSelect(prob);
            %newOpinion = find(prob == max(prob));   % << CHECK HERE>>
            newOp(v) = newOpinion;
            
        end
    end
    
    xs = [xs; newOp];
    desOpPop(it+1) = sum(op == desiredOp);

    op = newOp;
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

