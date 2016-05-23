%Omid55
%Arguments :
%   sp is the sparse network
%   lx is the last expected opinions
%   x is the expected opinions
function [ x ] = ExpectOpinions( sp,K,config ) 

M = size(sp,1);

%% Nonconvex Optimization with DE
disp('DE');
N = config.N;
MaxIteration = config.MaxIteration;
Coef = config.Coef;
Beta = config.Beta;
Pr = config.Pr;
Nv = config.Nv;   % number of difference vectors
withFigure = config.withFigure;
% DE/rand/1/binomial
InitPopulation = randi(K,N,M);
%InitPopulation = Scale * InitPopulation / max(InitPopulation);    % because we want to prevent people to be extremist then we make them have smaller opinions with a Scale
[x,~] = DE_Optimizer(Beta,Pr,Nv,Coef,MaxIteration,InitPopulation,sp,K,withFigure);


end