% initialize_variables.m

% initialize the population for CWN_AMMOA

% Written by Hieu V. Dang


% chromosome = initialize_variables(pop,MO,V,tmin,tmax,Pmin,Pmax,min_deci,max_deci);

function chromosome = initialize_variables(pop,MO,V,Q,Nc,min_range,max_range)

K = MO + V; % array consists of V decision variables and MO objectives

%% Initialize each chromosome
%chromosome = cell(1,pop);
chromosome = zeros(K,Nc,pop);
%t = zeros(Q,1);
%p = zeros(Q,Nc);
%gama = zeros(Q,Nc);
for j=1:pop
    for i=1:V
        %temp = tmin + (tmax-tmin)*rand();
        %t(i)= temp;
        for k=1:Nc
            %t(i,k) = temp; 
            %p(i,k) = Pmin + (Pmax-Pmin)*rand();
            %gama(i,k) = min_deci + (max_deci - min_deci)*rand();
            chromosome(i,k,j) = min_range(i) + (max_range(i) - min_range(i))*rand(1);
        end
    end
    
    %chromosome(1:Q,1,j) = t;
    %chromosome(Q+1:2*Q,:,j) = p;
    %chromosome(2*Q+1:V,:,j) = gama;
    %F = evaluate_objective(chromosome(:,:,j),V,Q,Nc);
    
    chromosome(V+1:K,1,j) = evaluate_objective(chromosome(:,:,j),V,Q,Nc);
end



