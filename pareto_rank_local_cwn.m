
% pareto_rank.m
% July 31, 2014
% Written by Hieu V. Dang

function f = pareto_rank_local_cwn(population,Nvar,Nfun)

Npop = size(population,3);
front = 1;
F(front).f = [];
individual = [];

for i=1:Npop
    individual(i).n = 0; %number of individuals that dominate this individual
    individual(i).p = [];   % Individuals whichs this individual dominate
    
    for j=1:Npop
        dom_less = 0;
        dom_equal = 0;
        dom_more = 0;
        
        for k=1:Nfun
            if (population(Nvar+k,1,i) < population(Nvar+k,1,j))
                dom_less = dom_less +1;
            elseif (population(Nvar+k,1,i)==population(Nvar+k,1,j))
                dom_equal = dom_equal +1;
            else
                dom_more = dom_more +1;
            end
        end
        if dom_less ==0 && dom_equal ~= Nfun
            individual(i).n = individual(i).n + 1;
        elseif dom_more == 0 && dom_equal ~= Nfun
            individual(i).p = [individual(i).p j];
        end
    end
    
    if individual(i).n == 0
        population(Nvar+Nfun+1,1,i) = 1;
        F(front).f = [F(front).f i];
    end
end

% Find the subsequent fronts
while ~isempty(F(front).f)
   Q = [];
   for i = 1 : length(F(front).f)
       if ~isempty(individual(F(front).f(i)).p)
        	for j = 1 : length(individual(F(front).f(i)).p)
            	individual(individual(F(front).f(i)).p(j)).n = ...
                	individual(individual(F(front).f(i)).p(j)).n - 1;
        	   	if individual(individual(F(front).f(i)).p(j)).n == 0
               		population(Nfun + Nvar+ 1,1,individual(F(front).f(i)).p(j)) = ...
                        front + 1;
                    Q = [Q individual(F(front).f(i)).p(j)];
                end
            end
       end
   end
   front =  front + 1;
   F(front).f = Q;
end


f = population;


        