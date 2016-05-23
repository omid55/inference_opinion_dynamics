%Omid55
% DE Optimizer
function [ BestChromosome,BestFitness ] = DE_Optimizer( Beta,Pr,Nv,Eps,MaxIteration,InitPopulation , sp,lx )

%% Differential Evolution (DE)
population = InitPopulation;

if ~exist('lx','var')
    fitnesses = CalculateFitnesses(population,sp);
else
    fitnesses = CalculateFitnesses(population,sp,lx);
end
it = 0;
bests = [];
means = [];
while it < MaxIteration
    lastDe = population;

    it = it + 1
    bests = [bests; max(fitnesses)];
    means = [means; mean(fitnesses)];
    plot(1:it,bests,1:it,means);
    legend('Max of Fitnesses','Mean of Fitnesses');
    xlabel('Generations');
    ylabel('Fitnesses');
    title('Differential Evolution');

    if ~exist('lx','var')
        [population] = CreateTrialVectorAndCrossOver(population,fitnesses,Beta,Pr,Nv,sp);
        fitnesses = CalculateFitnesses(population,sp);
    else
        [population] = CreateTrialVectorAndCrossOver(population,fitnesses,Beta,Pr,Nv,sp,lx);
        fitnesses = CalculateFitnesses(population,sp,lx);
    end
    
    if norm(population - lastDe) < Eps
        break;
    end
end
best = find(fitnesses == max(fitnesses));
bestIndex = best(1);
BestChromosome = population(bestIndex,:);
BestFitness = fitnesses(bestIndex);

end

