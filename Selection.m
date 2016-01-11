% Selection.m
% August 10, 2014
% Written by Hieu V. Dang

function f = Selection(population,pool_size,tour_size,Pt)

[K,~,Npop] = size(population);

% the last element contains the Reny entropy
rHr = K-1;
% Until the mating pool is filled, performed tournament selection
candidate = zeros(1,tour_size);

it=1;
while it < pool_size
%for i=1:pool_size
    %select tour_size individuals at random
    %candidate = zeros(1,tour_size);
    
    for j=1:tour_size
        candidate(j) = round(Npop*rand(1));
        if candidate(j)==0
            candidate(j)=1;
        end
        if j>1
            % make sure that same candidate is not choosen
            if ~isempty(find(candidate(1:j-1) == candidate(j)))
                candidate(j) = round(Npop*rand(1));
                if candidate(j) == 0
                    candidate(j)=1;
                end
            end
        end
    end
    
    % Collect information about the selected candidate
    RRSE = zeros(1,tour_size);

    for j=1:tour_size
        RRSE(j) = population(rHr,1,candidate(j));
    end
    % sorting RRSE ascending 
    [~,indxR] = sort(RRSE);
    
    Z = rand(1);
    if  Z <Pt
        f(:,:,it) = population(:,:,indxR(1));
        it = it+1;
    elseif Z < Pt*(1-Pt) 
        f(:,:,it) = population(:,:,indxR(2));
        it = it+1;
    elseif Z < Pt*((1-Pt)^2)
        f(:,:,it) = population(:,:,indxR(3));
        it = it+1;
    else
        continue
    end

  
end


    
            
    
