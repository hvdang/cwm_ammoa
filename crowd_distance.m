% crowd_distance.m
% May 14, 2014
% Written by Hieu V. Dang

function f = crowd_distance(population,Nvar,Nfun)

%Npop = size(population,1);

K = Nvar+Nfun;


max_rank = max(population(K+1,1,:));

previous_index=0;
for i=1:max_rank
    % find the index for current rank
    %temp_pop = [];
    current_index = max(find(population(K+1,1,:)==i));
    
    temp_pop = population(:,:,previous_index+1:current_index);
    %temp_pop1 = population(previous_index+1:current_index,:);
    
    for j=1:Nfun
        % sorting based on objective
        [~,index_obj] = sort(temp_pop(Nvar+j,1,:));
        
        sorted_based_objective = [];
        for k=1:length(index_obj)
            sorted_based_objective(:,:,k) = temp_pop(:,:,index_obj(k));
        end
        
        fmax = sorted_based_objective(Nvar+j,1,length(index_obj));
        fmin = sorted_based_objective(Nvar+j,1,1);
        
        temp_pop(K+1+j,1,index_obj(length(index_obj))) = Inf;
        temp_pop(K+1+j,1,index_obj(1)) = Inf;
        
        for k=2:length(index_obj)-1
            next_obj = sorted_based_objective(Nvar+j,1,k+1);
            previous_obj = sorted_based_objective(Nvar+j,1,k-1);
            if (fmax - fmin == 0)
                temp_pop(K+1+j,1,index_obj(k)) = Inf;
            else
                temp_pop(K+1+j,1,index_obj(k)) = 100*(next_obj - previous_obj)/(fmax - fmin);
            end
        end
    end
    distance = zeros(1,1,current_index);
    for j = 1:Nfun
        distance(1,1,:) = distance(1,1,:) + temp_pop(K+1+j,1,:);
    end
    temp_pop(K+3,1,:) = distance;
    y = temp_pop;%(K+3,1,:);
    z(:,:,previous_index+1:current_index) = y;
end
f = z();
    
    