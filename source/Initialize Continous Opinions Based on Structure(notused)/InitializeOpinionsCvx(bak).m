%Omid55
%Convex Estimation
function [  ] = Main(  )

%% Clear Everything
clc;
close all;


%% Network Creation
% BA
N = 100;
averageDegree = 12;
y = BAnet(N,averageDegree/2,averageDegree/2);


% %network_file = 'out.petster-friendships-hamster.txt';
% %network_file = 'USairport500-net.txt';
% network_file = 'out.ego-facebook.txt';
% a = importdata(network_file);
% N = max(max(a(:,1)),max(a(:,2)));
% net = sparse(a(:,1),a(:,2),1,N,N);


%% Optimization with cvx
D = zeros(N,N);
alpha = 1;
cvx_begin
    variables d(N*(N-1)/2,1)

    fitness = 0;
    for i = 1 : N
        for j = 1 : N
            if i ~= j
                p = min(i,j);
                q = max(i,j);
                ij = (p - 1) * N - (p - 1) * p / 2 + (q - p);
                fitness = fitness + y(i,j) * (alpha - d(ij)) - log(1 + exp(alpha - d(ij)));
            end
        end
    end

    maximize fitness

    subject to
        d >= 0;
cvx_end

for i = 1 : N
    for j = 1 : N
        if i ~= j
            p = min(i,j);
            q = max(i,j);
            ij = (p - 1) * N - (p - 1) * p / 2 + (q - p);
            D(i,j) = d(ij);
        end
    end
end

LL = sum(sum( y * (repmat(alpha,[N,N]) - D) - log(1 + exp(repmat(alpha,[N,N]) - D)) ))
save('Data.mat');

figure;
histfit(D(:));
title('D');

end
