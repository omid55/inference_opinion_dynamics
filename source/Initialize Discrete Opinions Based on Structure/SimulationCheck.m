%Omid55
function [  ] = SimulationCheck(  )

%% Clear Everything
clc;
close all;

load('OpinionsData.mat');
x = initOp;
N = size(net,1);
%x = 2 * rand(1,N) - 1;


%% Simulation Body
MaxIterations = 5000;
Randomness = 0.01;
Mu = 0.5;
% U = 0.5;
rangeBegin = -1;
rangeEnd = 1;
xs = [];
for it = 1 : MaxIterations
    %for i = 1 : N
    i = randi(N);
    adj = find(net(i,:)==1);
    if isempty(adj)
        continue;
    end
    
  %% Deffuant Expression with Bounded Confidence
    j = adj(randi(length(adj)));
    U = rand(1) * Randomness + 1 - abs(x(i));
    if abs(x(i) - x(j)) < U
        x(i) = StayInBound(x(i) + Mu * (x(j) - x(i)),rangeBegin,rangeEnd);
        x(j) = StayInBound(x(j) + Mu * (x(i) - x(j)),rangeBegin,rangeEnd);
    end
   
% % %    %% Hegselmann and Krause Bounded Confidence Model
% % %     op = 0;
% % %     leng = 0;
% % %     for k = 1 : length(adj)
% % %         j = adj(k);
% % %         U = rand(1) * Randomness + 1 - abs(x(i));
% % %         if sp(i,j) > 0 && abs(x(i) - x(j)) <= U
% % %             op = op + (x(j) - x(i));
% % %             leng = leng + 1;
% % %         end
% % %     end
% % %     if leng == 0
% % %         continue;
% % %     end
% % %     op = op / leng;
% % %     op = x(i) + Mu * op;
% % %     x(i) = StayInBound(op,rangeBegin,rangeEnd);

    %end
    if mod(it,N) == 0
        xs = [xs; x];
    end
end
xs = [xs; x];

title('Opinion Changing');
xlabel('Iterations');
ylabel('Opinion');
for i=1:N
    plot(xs(:,i));
    hold on;
end

end

