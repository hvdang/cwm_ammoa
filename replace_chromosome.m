% replace_chromosome.m
% Feburary 20, 2013

% Written by  Hieu V. Dang

function f = replace_chromosome(intermediate_chromosome,V,M,pop)

% function f  = replace_chromosome(intermediate_chromosome,pro,pop)
% This function replaces the chromosomes based on rank and crowding
% distance. Initially until the population size is reached each front is
% added one by one until addition of a complete front which results in
% exceeding the population size. At this point the chromosomes in that
% front is added subsequently to the population based on crowding distance.

N = size(intermediate_chromosome,3);

% Get the index for the population sort based on the rank
[temp,index] = sort(intermediate_chromosome(M+V+1,1,:));
clear temp;

% Now sort the individuals based on the index
for i=1:N
    sorted_chromosome(:,:,i) = intermediate_chromosome(:,:,index(i));
end
 
% Find the maximum rank in the current population
max_rank = max(intermediate_chromosome(M+V+1,1,:));

% Start adding each front based on rank and crowding distance until te
% whole population is filled
previous_index = 0;
for i=1:max_rank
    % Get the index for current rank
    current_index = max(find(sorted_chromosome(M+V+1,1,:)==i));
    % Check to see if the population is filled if all the individuals with
    % rank i is added to the population
    if current_index > pop
        % If so find the number of individuals with current rank i
        remaining = pop - previous_index;
        % Get information about the individuals in the current rank i
        temp_pop = sorted_chromosome(:,:,previous_index+1:current_index);
        % Sort the individuals with rank i in the descending order based on
        % the crowding distance
        [temp_sort,temp_sort_index] = sort(temp_pop(M+V+3,1,:),'descend');
        clear temp_sort;
        % Start filling individuals into the population in descending order
        % until the population is filled
        for j=1:remaining
            f(:,:,previous_index + j) = temp_pop(:,:,temp_sort_index(j));
        end
        return;
    elseif current_index < pop
        % Add all the individuals with rank i into the population
        f(:,:,previous_index + 1:current_index) = sorted_chromosome(:,:,previous_index+1:current_index);
    else
        % Add all the individuals with rank i into the population
        f(:,:,previous_index+1:current_index)= sorted_chromosome(:,:,previous_index+1:current_index);
        return
    end
    % Get the index for the last added individual
    previous_index = current_index;
end

% sort the chromosome according to the objective weights

