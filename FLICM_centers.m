% FLICM_centers.m
% May 12, 2014
% Written by Hieu V. Dang

function c = FLICM_centers(data,U,cnum,m)

Len = length(data);
c = zeros(1,cnum);

for k=1:cnum
    sSum=0;
    for i=1:Len
        sSum = sSum + U(i,k)^m;
        c(k) = c(k) + (U(i,k)^m)*data(i);
    end
    c(k) = c(k)/sSum;
end
