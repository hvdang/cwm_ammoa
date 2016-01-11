% FLICM.m
% May 12, 2014
% Written by Hieu V. Dang

function [Uout,c,iter] = FLICM(data,Uinit,m,cnum,winsize,maxIter,thrE)

Len = length(data); %size(data,1); 
sStep = (winsize-1)/2;
U = zeros(size(Uinit));
%Uold = zeros(size(U));
Uout = zeros(size(U));

U = Uinit;
Uold = U;
c = FLICM_centers(data,U,cnum,m);

itrs = 0;
dmax = 10;

while(dmax > thrE && itrs <= maxIter)
    % Calculation of the new array U
    for i=1:Len
        for k=1:cnum
            vec = data(i);
            sSum = abs(vec-c(k))^2;  % eucludean metric
            
            % Compute G
           % for ii=-sStep:sStep
           %     x = i + ii;
           %     dist = abs(x-i);
            %    if (x>0 && x<=Len && ii~=0)
            %        vec = data(x);
            %        sSum = sSum + (1/(1+dist))*((1-Uold(x,k)^m)*abs(vec-c(k)));
            %    end
            %end
            d(k) = sSum;
        end
        
        % Compute U
        for k=1:cnum
            dd = d(k);
            sSum = 0;
            for ii=1:cnum
                sSum = sSum + (dd/d(ii))^(2/(m-1));
                U(i,k) = 1/sSum;
            end
        end
    end
    
    % Next center calculation
    c = FLICM_centers(data,U,cnum,m);
    
    dmax = -1;
    % copy new array U to old one
    for i=1:Len
        for k=1:cnum
            if(dmax < abs(Uold(i,k) - U(i,k)))
                dmax = abs(Uold(i,k) - U(i,k));
            end
            Uold(i,k) = U(i,k);
            Uout(i,k) = U(i,k);
        end
    end
    itrs = itrs + 1;
end
iter = itrs;
    
            

