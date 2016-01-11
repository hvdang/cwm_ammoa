% local_search_stage2.m for JSSPA problem

% this function refines each chromosome in the population to find the
% better neighbours. 

%Written by Hieu V. Dang
% Aoril 1, 2015

function [improved_population, gdAvg] = local_search_stage2(population,Nvar,Nfun,Q,Nc,Var_min,Var_max,Nnb)

K = Nvar + Nfun;

Npop = size(population,3);

%Nnb = 100;
inter_pop = zeros(K,Nc,Nnb+1);

% tabu params
nTL = Npop;
TL = zeros(K,Nc,nTL); 
% Tabu list setting
inTL=1;
%dom_rate = zeros(1,Npop);
%lamda = ones(1,Nfun);

for i=1:Npop
    
    % Generate the random normalized weight vector
    Zn = rand(1);
    lamda = zeros(1,Nfun);
    lamda(1) = 1-Zn^(1/(Nfun-1));
    temp=0;
    for j=2:Nfun-1
        temp = temp + lamda(j-1);
        lamda(j) = (1-temp)*(1-Zn/(1/(Nfun-j)));
    end
    lamda(Nfun) = 1 - temp - lamda(Nfun-1);
    
    % generate Nnb neighbours of population(i,:).
            
    chrome_neighbors = randomPointsInCircle_cwn(Nnb,population(1:Nvar,:,i),0.1,Var_min,Var_max);
    
    % evaluate objectives of Nnb neighbours and their information criterion
    for j=1:Nnb
        chrome_neighbors(Nvar+1:K,1,j) = evaluate_objective(chrome_neighbors(1:Nvar,:,j),Nvar,Q,Nc);
    end
    inter_pop(1:K,:,1:Nnb) = chrome_neighbors;
    inter_pop(1:K,:,Nnb+1) = population(1:K,:,i);
    
      
    inter_pop = pareto_rank_local_cwn(inter_pop,Nvar,Nfun);
   
    
    % calculate the dominance param
    neighbors_rank = inter_pop(K+1,1,1:Nnb);
    db_indx = find(neighbors_rank < inter_pop(K+1,1,Nnb+1));
    
    if ~isempty(db_indx)
        dom_rate(i) = length(db_indx);
        % select the best neighbor
        for j=1:length(db_indx)
            sum_obj(j) = lamda*inter_pop(Nvar+1:K,1,db_indx(j));
        end
        [~,min_sumObj] = min(sum_obj);        
        
        better_neighbor = inter_pop(:,:,min_sumObj);
    else
        dom_rate(i) = 0;
        better_neighbor = zeros(K,Nc);%population(i,1:K);
    end
    
    % Store the better solutions
    % check if better_neighbor is in the TL
    flag_tl = 0;
    for k=1:nTL
        if isequal(better_neighbor(1:Nvar,:),TL(1:Nvar,:,k))==1 
            flag_tl = flag_tl +1;
        end
    end
    
    % update Tabu
    if (flag_tl ==0) && (inTL < nTL)
        inTL = inTL +1;
        TL(1:K,:,inTL) = better_neighbor(1:K,:);
    elseif (inTL >= nTL) && (flag_tl ==0)
        inTL = 1;
        TL(1:K,:,inTL) = better_neighbor(1:K,:);
    end
end

    
% Calculate the dominance measure
 gdAvg = sum(dom_rate)/(Npop*Nnb)  ;

% Sort the TL
w_objs = zeros(1,nTL);
for i=1:nTL
    w_objs(i) = sum(TL(Nvar+1:K,1,i));
end
indx = find(w_objs~=0);
%[~,tl_index] = sort(w_objs(indx));
%length(indx)

TL_sort = TL(:,:,indx); 

improved_population = TL_sort;