%Omid55
function [ fitnesses ] = CalculateFitnesses( population,sp )

%% Fitness Caculation
fitnesses = zeros(size(population,1),1);
parfor i=1:size(population,1)
    fitnesses(i) = ObjectiveFunction(population(i,:),sp);
end

end

