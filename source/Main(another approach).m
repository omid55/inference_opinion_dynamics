%Omid55
%Learn Parameter tetta by the opinion data with fsolve
function [  ] = Main(  )

% Loading the dataset
load('OpinionDataset');
N = size(net,1);

% dd = [1 2; 1 3;2 3; 2 4;3 4];
% N = max(max(dd(:,1)),max(dd(:,2)));
% net = sparse(dd(:,1),dd(:,2),1,N,N);
% net = max(net,net');
% initOp = [1 1 0 0];
% data = [1 2 0;4 3 1;3 4 1;1 5 1];


%% Learning the parameter by EM algorithm
data = sortrows(data,2);  % sorting data through time steps

alpha0 = zeros(N,1);
alphaBar = zeros(N,1);
%options = optimoptions('fsolve','Display','iter'); % Option to display output
options = 0;
diff = 1;
cnt = 0;
while diff > 0.01
    [a,fval] = fsolve(@myfunc ,alpha0,options,data,net,initOp,alphaBar); % Call solver
    diff = norm(alphaBar - a)
    alphaBar = a;
        
    cnt = cnt + 1;
    %if mod(cnt,5) == 0
        SP = exp(a)
        degs = sum(net)';
        [r,p] = corr(exp(degs),SP)
        ooooooooooooooooooooo = 0;
    %end
end

end