%Omid55
%My function to be solved
function [ F ] = myfunc( alpha,data,net,initOp,alphaBar )

% F = MyCFunction(alpha,data,net,initOp,alphaBar);
% oooo = 0;

    N = length(alpha);
    F = zeros(N,1);
    for u = 1 : N
        item = exp(alpha(u));
        
        adjU = find(net(u,:) == 1);
        ss = 0;
        for i = 1 : length(adjU)
            vPrime = adjU(i);
            cvPrime = length(find(data(:,1) == vPrime));
            adjvPrime = net(vPrime,:) == 1;
            
            s = sum(exp(alpha(adjvPrime)));
%             s = 0;
%             for j = 1 : length(adjvPrime)
%                 uPrime = adjvPrime(j);
%                 s = s + exp(alpha(uPrime));
%             end

            ss = ss + cvPrime/s;
        end
        
        item = item * ss;
        
        sss = 0;
        op = initOp;
        newOp = op;
        for i = 1 : size(data,1)
            v = data(i,1);
            t = data(i,2);
            k = data(i,3);
            newOp(v) = k;
            
            if i > 1 && t ~= data(i-1,2)
                op = newOp;        % updating the opinion of node v
            end
            
            if net(v,u) ~= 0 && op(u) == k   % if u and v were compliant
                adjV = find(net(v,:) == 1);
                compliantV = adjV(op(adjV) == k);
                if isempty(compliantV)
                    continue;
                end
                
                s = sum(exp(alphaBar(compliantV)));
%                 s = 0;
%                 for j = 1 : length(compliantV)
%                     uu = compliantV(j);
%                     s = s + exp(alphaBar(uu));
%                 end
                
                sss = sss + exp(alphaBar(u))/s;
            end
        end
        
        item = item - sss;
        
        F(u) = item;
    end
    
end

