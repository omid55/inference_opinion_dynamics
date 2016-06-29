%Omid55
%Arguments :
%   sp is the sparse network
%   lx is the last expected opinions
%   x is the expected opinions
function [ x ] = ExpectOpinions( sp,lx ) 

M = size(sp,1);

%% Nonconvex Optimization with DE
disp('DE');
N = 100;
MaxIteration = 20;
Eps = 10 ^ -10;
Beta = 0.8;
Pr = 0.6;
Nv = 1;   % number of difference vectors
% DE/rand/1/binomial
InitPopulation = 2 * rand(N,M) - 1;
if ~exist('lx','var')
    [x,fitness] = DE_Optimizer(Beta,Pr,Nv,Eps,MaxIteration,InitPopulation,sp);
else
    [x,fitness] = DE_Optimizer(Beta,Pr,Nv,Eps,MaxIteration,InitPopulation,sp,lx);
end


end