
% RRSE.m
% Evaluate the Renyi relaitive summation entropy

% August 10, 2014
% Written by Hieu V. Dang

function [f] = RRSE(population,Nvar,Nfun)

Npop = size(population,3);
Nc = floor(Npop/3);
K = Nvar + Nfun;

alpha1 = 2:1:20;
for i=1:Npop
    
    % Select Nc individuals at random to compute the winning probability
    win_times = zeros(1,Nfun);
    candidate = zeros(1,Nc);
    prob = zeros(1,Nfun);
    
    for j=1:Nc
        % select an individual at random
        candidate(j) = round(Npop*rand(1));
        % the array has to start from 1
        if candidate(j) ==0
            candidate(j) = 1;
        end
        % Check that the same candidate is not chosen
        if j>1
            if ~isempty(find(candidate(1:j-1)==candidate(j)))
                candidate(j) = round(Npop*rand(1));
                if candidate(j) ==0
                    candidate(j) =1;
                end
            end
        end
    end
    
    % Compute the winning probability
    for j=1:Nc
        for k=1:Nfun
            temp = population(Nvar+k,1,i) - population(Nvar+k,1,candidate(j));
            if  temp < 0
                win_times(k) = win_times(k) + 1;
            end
        end
    end
    
    
    % winning probabilities of k objectives of chromosome i
    for k=1:Nfun
        prob(k) = win_times(k)/Nc;
     end
    
    % Generalized Renyi entropy for chromosome i
    temp1 = 0;

   for a = 1:length(alpha1)
    for k=1:Nfun
        if prob(k) ==0
            temp1 = temp1 +prob(k);
        else
            temp1 = temp1 + (prob(k))^(-alpha1(a));
        end
    end
    if temp1 == 0
        Rrse(a) = 100;
    else
        Rrse(a) = (log2(temp1))/(alpha1(a)-1);
    end
   end
    Rrse1 = sum(Rrse)/length(alpha1);
    population(K+2,1,i) = Rrse1;
end
f = population;
%% Find the Renyi distance


    
        
        
            
        
    
    

        