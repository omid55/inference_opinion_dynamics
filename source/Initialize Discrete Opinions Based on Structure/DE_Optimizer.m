%Omid55
% DE Optimizer
function [ BestChromosome,BestFitness ] = DE_Optimizer( Beta,Pr,Nv,Coef,MaxIteration,InitPopulation,sp,K,withFigure )

%% Differential Evolution (DE)
population = InitPopulation;

fitnesses = CalculateFitnesses(population,sp);
it = 0;
bests = [];
means = [];
if withFigure == 1
    figure;
end
while it < MaxIteration
    tic;

    it = it + 1;
    if withFigure == 1
        bests = [bests; min(fitnesses)];
        means = [means; mean(fitnesses)];
        figure(1);
        plot(0:it,bests,0:it,means);
        legend('Best Cost','Mean Cost');
        xlabel('Generations');
        ylabel('Cost');
        title('Differential Evolution Algorithm');
    end
    
    [population,fitnesses] = CreateTrialVectorAndCrossOver(population,fitnesses,Beta,Pr,Nv,sp,K);
    
    time = toc;
    disp(['DE: iteration ' num2str(it) ' done in ' num2str(time) 's.']);
    
    diff = length(find(fitnesses <= mean(fitnesses)));
    if diff >= Coef * size(population,1)
        break;
    end
end
best = find(fitnesses == max(fitnesses));
bestIndex = best(1);
BestChromosome = population(bestIndex,:);
BestFitness = fitnesses(bestIndex);

end

