% evaluate_objective.m

% Calculate the average throughput objective and average interference
% objective for AMMOA run

% Written by  Hieu V. Dang


function f = evaluate_objective(chromosome,V,Q,Nc)

%% Objective function 1
fq = 100;
Tf = 1e-2;% s, time frame for sensing 

tau = chromosome(1:Q,1);
power = chromosome(Q+1:2*Q,:);
gamma = chromosome(2*Q+1:V,:);

% weights
lamda = ones(1,Q)/Q;

% variance parameters
var_noise = 2e-4; % mw noise variance of channel
var_pt = 0.05; % 1W power of PU on all channels

% Compute intermediate Mu and Var parameters'
MUqk0 = var_noise;
MUqk1 = var_pt + var_noise;
SIGqk0_2 = var_noise^2;
SIGqk1_2 = (var_pt + var_noise)^2;

SIGqk0 = sqrt(SIGqk0_2);
SIGqk1 = sqrt(SIGqk1_2);

% Compute probabilities of false alarms and missed detections
prob_fa = zeros(Q,Nc);
prob_de = zeros(Q,Nc);
prob_miss = zeros(Q,Nc);

for i=1:Q
    for j=1:Nc
        temp1 = sqrt(tau(i)*fq)*(gamma(i,j) - MUqk0)/SIGqk0;
        prob_fa(i,j) = qfunc(temp1);
        temp2 = sqrt(tau(i)*fq)*(gamma(i,j) - MUqk1)/SIGqk1;
        prob_de(i,j) = qfunc(temp2);
        prob_miss(i,j) = 1 - prob_de(i,j);
    end
end
clear temp1 temp2;


%% Compute objective 1 (throughput)

[Hqq,Hrq] = channel_tf(Q,Nc);

%calculate achievable rate
rqk = zeros(Q,Nc);
for i=1:Q
    for j=1:Nc
        temp=0;
        for i1=1:Q
            if i1~=i
                temp = temp + (Hrq(i,j)^2)*power(i,j);
            end
        end
        rqk(i,j) = log2(1+(power(i,j)*Hqq(i,j)^2)/(var_noise+temp));
    end
end
clear temp;

% The throughputs
Rq = zeros(1,Q);
for i=1:Q
    %for j=1:K
        temp = 0;
        for j1=1:Nc
            temp = temp + (1-prob_fa(i,j1))*rqk(i,j1);
            %temp = temp + rqk(i,j1);
        end
        Rq(i) = (1-tau(i)/Tf)*temp;
    %end
end
clear temp;
f(1) = -lamda*Rq';

%% Compute Objective 2, Interference
Sq = zeros(1,Q);
for i=1:Q
    temp=0;
    for k=1:Nc
        temp = temp + prob_miss(i,k)*power(i,k);
    end
    Sq(i) = temp;
end
f(2) = lamda*Sq';


