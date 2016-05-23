%Omid55
function [ fitnesses ] = CalculateFitnesses( population,sp,lx )

%% Fitness Caculation
fitnesses = zeros(size(population,1),1);
for i=1:size(population,1)
    if nargin == 2
        fitnesses(i) = ObjectiveFunction(population(i,:),sp);
    else
        fitnesses(i) = ObjectiveFunction(population(i,:),sp,lx);
    end
end

end

