%Omid55
%Convex Estimation
function [  ] = Main(  )

%% Clear Everything
clc;
close all;


%% Network Creation
% BA
N = 30;
averageDegree = 12;
y = BAnet(N,averageDegree/2,averageDegree/2);

% %network_file = 'out.petster-friendships-hamster.txt';
% %network_file = 'USairport500-net.txt';
% network_file = 'out.ego-facebook.txt';
% a = importdata(network_file);
% N = max(max(a(:,1)),max(a(:,2)));
% net = sparse(a(:,1),a(:,2),1,N,N);


%% Optimization with cvx
alpha = 0.001;
cvx_begin
    variables x(1,N)

    fitness = 0;
    for i = 1 : N
        fitness = fitness + sum(  y(i,:) .* ( alpha - (x(i) - x) ) - log( 1 + exp( alpha - (x(i) - x) ) )  );
    end
     
%     fitness = 0;
%     for i = 1 : N
%         for j = 1 : N
%             if i ~= j
%                 fitness = fitness + y(i,j) * ( alpha - (x(i)-x(j)) ) - log( 1 + exp( alpha - (x(i)-x(j)) ) );
%             end
%         end
%     end

    maximize fitness
cvx_end

hist(x);

end
