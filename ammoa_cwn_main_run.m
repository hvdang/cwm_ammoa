% cwn_ammoa_main_run.m
% this program is for multiobjective spectrum sensing and power allocation
% based on AMOMA

% April 1, 2015
% Written by Hieu V. Dang

close all
clear all
clc

tic

%% Parameters for CWN
Q=6; % 3 SUs
Np=1; % 1 PU
Nc = 8; % 2 channels
Tf = 2e-2;% s, time frame for sensing    
tmin = 1e-4; %min of sensing time 
tmax = 5e-3; %max of sensing time
Pmindb = -30; %dBm
Pmin = 10^(Pmindb/10); %mw
Pmaxdb = -3; % dbm
Pmax = 10^(Pmaxdb/10); %mw
thresh_min = -30; %dbm
min_deci = 10^(thresh_min/10);
thresh_max = -3; %dbm
max_deci = 10^(thresh_max/10);
fq = 100; %hz

%% Parameters for memetic algorithms
% population size
Npop = 100;
% number of generation
maxIters =400;
% number of objectives
Nfun = 2;
% number of decision variables
Nvar = 3*Q;
K = Nvar+Nfun;
pool = Npop; %round(pop/2); % size of mating pool
tour = 10; % tournament size.

% distribution indices for crossover and mutation operators 
mu = 10;
mum = 20;

% The min and max values for each decision variable
for i=1:Q
    min_range(i) = tmin ;
    max_range(i) = tmax ;
end
for i=Q+1:2*Q
    min_range(i) = Pmin ;
    max_range(i) = Pmax ;
end
for i=2*Q+1:Nvar
    min_range(i) = min_deci;
    max_range(i) = max_deci;
end


%% Initialize population
init_population = initialize_variables(Npop,Nfun,Nvar,Q,Nc,min_range,max_range);
population = init_population;

%% Find the Pareto ranks and crowding distance
population = pareto_rank(population,Nvar,Nfun);
population = crowd_distance(population,Nvar,Nfun);

%% Compute the RRSE
population = RRSE(population,Nvar,Nfun);

%% Start the evolutionary process
population_array = cell(1,maxIters);
Pt=1; 

IQ = zeros(1,maxIters);
IQ0 = IQ;
t0=5;

dqAverage = [];

%% Main loop
for i=1:maxIters
    
    % select parents for reproduction
    parent_population = Selection(population,pool,tour,Pt);
    
    % Perform crossover and mutation
    offspring_population = realcrossover(parent_population,Nfun,Nvar,Q,Nc,mu,min_range,max_range);
    offspring_population = realmutation(offspring_population,Nfun,Nvar,Q,Nc,mum,min_range,max_range);
    N_offs = size(offspring_population,3);
    
    % Combine and update populations
    inter0_population(1:K,:,1:Npop) = population(1:K,:,:);
    inter0_population(1:K,:,Npop+1:Npop+N_offs) = offspring_population(1:K,:,:);
    inter0_population = pareto_rank(inter0_population,Nvar,Nfun);
    inter0_population = crowd_distance(inter0_population,Nvar,Nfun);
    population0 = replace_chromosome(inter0_population,Nvar,Nfun,Npop);
    population0 = RRSE(population0,Nvar,Nfun);
    Npop0 = size(population0,3);
    
    % Doing clustering for local search
     [CLuster,IQ(i),RE] = clustering(population0,K,10,1000,0.000000001);
     
     if i > t0
        temp=0;
        for j=i-t0:i-1
            temp = temp + IQ(j)/(2*(i-j));
        end
        IQ0(i) = temp/(t0);
     end
    
     if IQ(i) >= IQ0(i)
         % Perform Tabu local search stage 1
         improved_population = tabu_search_stage1(CLuster,Nvar,Nfun,Q,Nc,min_range,max_range);
         if ~isempty(improved_population)
             N_impr = size(improved_population,3);
         else
             N_impr=0;
         end
         
         % Intermediate population
         inter_population(1:K,:,1:Npop0) = population0(1:K,:,:);
         if N_impr > 0
             inter_population(1:K,:,Npop0+1:Npop0+N_impr) = improved_population;
         end
         Ninter = size(inter_population,3);
         
         % Pareto ranks & Crowding distance
         inter_population = pareto_rank(inter_population,Nvar,Nfun);
         inter_population = crowd_distance(inter_population,Nvar,Nfun);
     
         % Chromosomes replacement 
         population1 = replace_chromosome(inter_population,Nvar,Nfun,Npop);
         
         % local search stage 2 
         %impr_pop1 = tabu_search_stage2(population1,Nvar,Nfun,Q,Nc,min_range,max_range); 
         [impr_pop1,dqAverage0] = local_search_stage2(population1,Nvar,Nfun,Q,Nc,min_range,max_range,200);
         dqAverage = [dqAverage dqAverage0];
         if ~isempty(impr_pop1)
             Nimpr1 = size(impr_pop1,3);
         else
             Nimpr1=0;
         end
         
         % inter pop 
         inter_pop1(1:K,:,1:Npop) = population1(1:K,:,:); 
         inter_pop1(1:K,:,Npop+1:Npop+Nimpr1) = impr_pop1(1:K,:,:);
         
         % Pareto ranks and crowdistance
         inter_pop1 = pareto_rank(inter_pop1,Nvar,Nfun);
         inter_pop1 = crowd_distance(inter_pop1,Nvar,Nfun);
         
         % Replacement
         population = replace_chromosome(inter_pop1,Nvar,Nfun,Npop);
         
         % compute RRSE values 
         population = RRSE(population,Nvar,Nfun);
     else
         Pt=Pt-0.05;
         if Pt < 0.1
             Pt = 1;
         end
     end
     population_array{i} = population;
    
    
end

toc

%% Visualize pareto ranks
population = population_array{44};
for i=1:Npop
    f1_init(i) = -init_population(Nvar+1,1,i);
    f2_init(i) = init_population(Nvar+2,1,i);
end

for i=1:Npop
    f1(i) = -population(Nvar+1,1,i);
    f2(i) = population(Nvar+2,1,i);
end
% Visualize the pareto rank
figure(1);
for i=1:Npop
plot([f2(i),f2_init(i)],[f1(i),f1_init(i)],':k'); hold on
end
%hold on
plot(f2,f1,'or','LineWidth',2); 
hold on
plot(f2_init,f1_init,'*b','LineWidth',2);
hold off
xlabel('Averaged Interference, I');
ylabel('Averaged throughput, R');
axis tight;

figure(2);
for i=1:Npop
plot([f1(i),f1_init(i)],[f2(i),f2_init(i)],':k'); hold on
end
%hold on
plot(f1,f2,'or','LineWidth',2);
hold on
plot(f1_init,f2_init,'*b','LineWidth',2);
hold off
ylabel('Averaged Interference, I');
xlabel('Averaged throughput, R');
axis tight;


% visualize dqAverage
ld= length(dqAverage);
figure(3);
plot(1:ld,dqAverage);
axis tight

%% save files
%save initial population
savefile1 ='Results/init_pop100.mat';
save(savefile1,'init_population');

% save the array of population
savefile2 = 'Results/pop100_400iters.mat';
save(savefile2,'population_array');


savefile3 = 'Results/dqAverage_400iters.mat';
save(savefile3,'dqAverage');

