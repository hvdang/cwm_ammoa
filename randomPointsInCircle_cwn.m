% randomPointsInCircle.m
% This function to generate a random numbers in a point circles 

% Written by Hieu V Dang


function X_nb = randomPointsInCircle_cwn(Nnb, c, r,Var_min,Var_max)
%n: number of generating samples
% c: current set of variables
% r: radius = 0.05
c = c';

d = size(c,2);
n = size(c,1);
X_nb = zeros(d,n,Nnb);

for t = 1:Nnb

rad = rand(n,1).*r;
X = ones(n,d).*repmat(rad,1,d);
X1 = X;
phi = zeros(n,d-1);
phi(:,1:d-2) = rand(n,d-2).*pi;
phi(:,d-1) = rand(n,1).*2.*pi;
for j = 1:n
    for k = 1:d
        if k > 1
            for l = 1:k-1
                X(j,k) = X(j,k).*sin(phi(j,l));
            end;
        end;
        if k < d
            X(j,k) = X(j,k).*cos(phi(j,k));
        end;
    end
end;
C = c;%repmat(c,n,1)
%X1 = repmat(c,n,1) + X;
X1 = C + X;
% check the constraints
for i=1:n
    for j=1:d
        if(X1(i,j) < Var_min(j))
            X1(i,j) = C(i,j) + abs(X(i,j));
        end
        if (X1(i,j) > Var_max(j))
            X1(i,j) = max(C(i,j) - X(i,j),X(i,j)-C(i,j));
        end
    end
end

X_nb(:,:,t) = X1';

end