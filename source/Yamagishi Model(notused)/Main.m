%Omid55
%Learn Parameter t by the opinion data
%TESTING
function [  ] = Main(  )

%% Loading the dataset
load('OpinionDataset');


%% Learning the parameter by EM algorithm
J = 5;
tHat = zeros(J,1);
loglikelihood = 0;
while 1
    
    cvx_begin
        %variables t(J,1) x(N,J)
        variable t(J,1)

        op = initOp;
        for i = 1 : size(data,1)
            v = data(i,1);
            k = data(i,3);
            op(v) = k;                                     % updating the opinion of node v
            adj = find(net(v,:) == 1);
            compliant = adj(op(adj) == k);

            v1 = 0;
            su1 = sum(exp(x(compliant,:) * tHat));
            for j = 1 : length(compliant)
                u = compliant(j);
                v1 = v1 + (exp(x(u,:)*tHat)/su1) * x(u,:) * t;
            end

            v2 = sum(exp(x(adj,:)*t));
            loglikelihood = loglikelihood + (v1 - log(v2));
        end

        maximize loglikelihood
        subject to
        t > 0
    cvx_end
    
    tHat = t

end

Tetta = t


end

