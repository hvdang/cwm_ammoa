% Clustering
% May 24, 2014
% Written by Hieu V. Dang

function [cluster,IQ,U] = clustering(population,K,Cnum,maxIter,thrE)

m = 2;
winSize = 0;

[Npop] = size(population,3);

for i=1:Npop
RenyE(i) = population(K+2,1,i);
end



% Initialize randomly the memberships array U
U = rand(Npop,Cnum);

% FLICM
[U,c,iter] = FLICM(RenyE,U,m,Cnum,winSize,maxIter,thrE);

[~,clus] = max(U,[],2);

cluster = cell(1,Cnum);

for i=1:Cnum
    q = find(clus == i);
    cluster{i} = population(:,:,q);
end

% Calculate the clustering quality index
for i=1:Cnum
    pop = cluster{i};
    temp = 0;
    l = size(pop,3);
    for j=1:l
        temp = temp + (pop(K+2,1,j)-c(i))^2;
    end
    inDist(i) = temp/l;
end
IH = sum(inDist);

for i=1:Cnum
    temp = 0;
    for j=1:Cnum
        temp = temp + (c(i) - c(j))^2;
    end
    cenDist(i) = temp/Cnum;
end
IS = sum(cenDist);

IQ = IH;