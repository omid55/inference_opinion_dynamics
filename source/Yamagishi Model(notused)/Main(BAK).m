%Omid55
%Learn Parameter tetta by the opinion data
function [  ] = Main(  )

%% Loading the dataset
load('OpinionDataset');


%% Learning the parameter by EM algorithm
J = 5;
Eps1 = 10^-10;
Eps2 = 10 ^ -3;
tetta = zeros(J,1);
tettaHat = tetta;
%lastNr = Inf;
it = 0;
while 1
    
    it = it + 1;
    gradient = zeros(1,J);
    op = initOp;
    for i = 1 : size(data,1)
        v = data(i,1);
        k = data(i,3);
        op(v) = k;                                     % updating the opinion of node v
        adj = find(net(v,:) == 1);
        compliant = adj(op(adj) == k);
        
        v1 = 0;
        su1 = sum(exp(x(compliant,:) * tettaHat));
        for j = 1 : length(compliant)
            u = compliant(j);
            v1 = v1 + (exp(x(u,:)*tettaHat)/su1) * x(u,:);
        end
        
        v2 = 0;
        su2 = sum(exp(x(adj,:) * tetta));
        for j = 1 : length(adj)
            u = adj(j);
            v2 = v2 + (exp(x(u,:)*tetta)/su2) * x(u,:);
        end
        
        gradient = gradient + (v1 - v2);
    end
    
    nr = norm(gradient)
%     if nr > lastNr
%         oooo = 0;
%     end
%     lastNr = nr;
    if nr < Eps1
        break;
    end
    
    hessian = zeros(J,J);
    for i = 1 : size(data,1)
        v = data(i,1);
        adj = find(net(v,:) == 1);
        su = sum(exp(x(adj,:) * tetta));
        v1 = 0;
        for j = 1 : length(adj)
            u = adj(j);
            v1 = v1 + (exp(x(u,:)*tetta)/su) * x(u,:)' * x(u,:);
        end
        v2 = 0;
        for j = 1 : length(adj)
            u = adj(j);
            v2 = v2 + (exp(x(u,:)*tetta)/su) * x(u,:);
        end
        
        hessian = hessian - (v1 - v2 * v2');
    end
    
    deltaTetta = - hessian \ gradient';
    tetta = tetta + deltaTetta;
    %norm(deltaTetta)
    if abs(tettaHat - tetta) < Eps1
        TETTA = tetta
        originalTetta
        break;
    end
    if norm(deltaTetta) < 0.01
        tettaHat = tetta;
        tetta = zeros(J,1);
        continue;
    end
    
end

end

