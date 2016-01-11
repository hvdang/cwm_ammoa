
% tabu_search_stage1.m
% August 4, 2014,
% Written by Hieu V. Dang

function [improved_population] = tabu_search_stage1(cluster,Nvar,Nfun,Q,Nc,min_range,max_range)

Nclus = length(cluster);


K = Nvar + Nfun;
improved_temp = cell(1,Nclus);

if rand(1) <0.4% check probability of local search at stage 1

for jj = 1:Nclus

population = cluster{jj};


Npop= size(population,3);
n = 10;%  floor(Npop);
maxIters = 5*n;
nTL = n; % length of tabu list
TL = zeros(K,Nc,nTL); % Tabu list

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

    %lamda = [1 1];


       
    % select the best offspring as the initial population for local search
    rEntropy = zeros(1,Npop);
    for i=1:Npop
        rEntropy(i) = sum(lamda.*population(Nvar+1:K,1,i)');
        %rEntropy(i) = population(K+2,1,i);
    end
    [~,minEntropy] = min(rEntropy);
    %clear temp;
    
    initial_offspring = population(1:K,:,minEntropy);
    x = initial_offspring;
    
    init_w_objs = sum(lamda.*x(Nvar+1:K,1,1)');
    better_w_objs = init_w_objs;
    
    %Tabu list setting
    inTL = 1;
    TL(1:K,:,inTL) = x(1:K,:,:);
    
    
    %%%%%Main loop%%%%%%%
    %S = zeros(n,K);
    fv = zeros(1,n);
    itrc = 0;
    while(itrc < maxIters)
        
        S = initialize_variables(n,Nfun,Nvar,Q,Nc,min_range,max_range);
        for i=1:n
            fv(i) = sum(lamda.*S(Nvar+1:K,1,i)');
        end
        
        % Search for better neighbors
        for i=1:n
            % check S is in the TL
            flag_tl = 0;
            for j=1:nTL
                if isequal(S(1:Nvar,:,i),TL(1:Nvar,:,i))==1
                    flag_tl = flag_tl +1;
                end
            end
            
            b_flag =0;
            if(fv(i) < better_w_objs) && (flag_tl==0)
                xnew = S(1:K,:,i);
                fnew = fv(i);
                b_flag=1;
            else
                xnew = x(1:K,:,:);
                fnew = better_w_objs;
            end
            
            % update Tabulist
            if (inTL < nTL) && (b_flag==1)
                inTL = inTL +1;
               TL(1:K,:,inTL) = xnew;
            elseif (inTL >= nTL) && (b_flag==1) 
                inTL = 1;
                TL(1:K,:,inTL) = xnew;
            end
        end
        x(1:K,:,:) = xnew(1:K,:,:);
        better_w_objs = fnew;
        itrc = itrc +1;
    end
%end

% Sort the TL
w_objs = zeros(1,nTL);
for i=1:nTL
    w_objs(i) = sum(lamda.*TL(Nvar+1:K,1,i)');
end
indx = find(w_objs~=0);
[~,tl_index] = sort(w_objs(indx));

TL_sort = TL(:,:,tl_index);
%improved_population = TL_sorted(1:n,1:K);
if length(tl_index) < n
    improved_temp{jj} = TL_sort;        
else
    improved_temp{jj} = TL_sort(:,:,1:n);
end

end   %end For JJ

improved_pop = [];
k=0;
for i=1:Nclus-1
    ln = size(improved_temp{i},3);
    for ii = 1:ln
        improved_pop(:,:,k+ii) = improved_temp{i}(:,:,ii);
    end
    k = k+1;
end

improved_population = improved_pop; 


else
    improved_population = [];
end  % end if
     
