%Omid55
function [ sp ] = CalculateSocialPowers( net,PowerfullStrength,soc_type )    % IMPORTANT: just one time calculation and not more

N = size(net,1);

switch(soc_type)
    case 1
    % high degree apporach
    degs = sum(net);
    sp = degs;

    case 2
    % random approach
    sp = rand(1,N);
    
    case 3
    % high closeness apporach
    asp = all_shortest_paths(net);
    sp = (N-1) ./ sum(asp);

    case 4
    % high betweenness apporach
    sp = betweenness_centrality(net)';

    case 5
    % low degree apporach
    degs = sum(net);
    sp = max(degs) - degs + 10;
    sp = sp * 3;
end

if soc_type ~= 5
    % for normalization
    sp = max(sum(net)) * sp / max(sp);
end


sorted = sort(sp,'descend');
len = ceil(N/10);
threshold = sorted(len);
highSPs = find(sp > threshold);
indices = find(sp == threshold);
highSPs = [highSPs indices(1:len-length(highSPs))];        
list = 1 : N;
list(highSPs) = [];
selected = list(randperm(length(list),ceil(N/20)));
%selected = randi(N,1,ceil(N/20));

% sorted = sort(degs,'descend');
% threshold = sorted(ceil(N/20));
% highDegs = find(degs > threshold);
% indices = find(degs == threshold);
% selected = [highDegs indices(1:ceil(N/20)-length(highDegs))];

sp(selected) = PowerfullStrength;
% sp(highSPs(1:length(selected)/2)) = PowerfullStrength;
% sp(selected(1:length(selected)/2)) = PowerfullStrength;


%sp = exp(-degs);
%sp = 1./degs;   
%sp = exp(degs);
%sp = degs;

% %% --== Gaussian noise added ==--
% stand = 100;
% N = length(degs);
% noise = stand * randn(1,N);
% sp = degs(1:N) + noise;
% sp = sp - min(sp);
% meanDeg = mean(degs(1:N));
% sp = meanDeg * sp / mean(sp);
% DegsWithSpCorr = corr(full(degs)',full(sp)')

% % Matlab's noise
% N = length(degs);
% sp = awgn(degs,10,'measured');
% s
% p = sp - min(sp);
% meanDeg = mean(degs(1:N));
% sp = meanDeg * sp / mean(sp);
end

