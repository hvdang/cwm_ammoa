% demo_results_visualization.m

% April 4, 2015

% Written by Hieu V. Dang


close all 
clear all
clc

Q=10;
Nc=12;
V = 3*Q;

%% visualize initial population of 100
load 'results\init_pop100.mat';
pop = size(chromosome_init,3);

% Visualize pareto ranks
for i=1:pop
    f1_init(i) = -chromosome_init(V+1,1,i);
    f2_init(i) = chromosome_init(V+2,1,i);
end

figure(1);
plot(f1_init,f2_init,'*b','LineWidth',2);
ylabel('Averaged Interference, I');
xlabel('Averaged Throughput, R');
axis tight;

%% Visualize obtained population
load 'results\pops100_400iters.mat';
population = population_array{400};
Npop = size(population,3);

for i=1:Npop
    f1(i) = - population(V+1,1,i);
    f2(i) = population(V+2,1,i);
end

figure(2);
plot(f1,f2,'or','LineWidth',2);
ylabel('Averaged Interference, I');
xlabel('Averaged Throughput, R');
axis tight;

%% obtained variables from the best selected chromosome

best_chromosome = population(:,:,73);

sense_time = best_chromosome(1:Q,1);
deci_threshold = best_chromosome(Q+1:2*Q,:);
power = best_chromosome(2*Q+1:V,:);

% Plot the sensing time
q = 1:Q;
figure(3);
%[t_best1,t_best2] = hist(bsol(1:Q,1));
%plot(t_best2,t_best1,'-o');
bar(q,sense_time,0.3);
ylabel('Sensing time (Second)');
xlabel('Secondary user, i');
axis tight

%% plot decision threshold
figure(4)
plot(1:Nc,deci_threshold(10,:),'-r*','LineWidth',2);
xlabel('Channel, k');
ylabel('Decision threshold');
axis tight;


%% PLOT power 
figure(5)
plot(1:Nc,power(10,:),'-bo','LineWidth',2);
xlabel('Channel, k');
ylabel('Power, mw');
axis tight;

