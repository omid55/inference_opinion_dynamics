%Omid
%HITS ranking method
function [ ranks ] = HITS( net )

%% HITS Ranking
epsilon = 10^-10;
N = size(net,1);
MaxIts = 10^4;

h = ones(N,1);    % initialization of hubs array
a = ones(N,1);    % initialization of authorities array

delta = 99999;
t = 0;
while delta > epsilon && t < MaxIts

    t = t + 1
    lastA = a;
    
    disp('O Operation ...');
    %O Operation
    for i=1:N
        h(i) = sum(a( net(i,:) == 1 ));
    end

    disp('I Operation ...');
    %I Operation
    for i=1:N
        a(i) = sum(h( net(:,i) == 1 ));
    end
    
    disp('Normalization ...');
    % Normalization
    a = a / max(a);   % normalization with max weight
    %a = a / norm(a);    % normalization with norm 2
    
    delta = norm(a - lastA);

end

ranks = a;

end

