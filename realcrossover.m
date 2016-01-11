% realcrossover.m
% Feb 23, 2013

% Written by Hieu V. Dang


function f = realcrossover(parent_chromosome,M,V,Q,Nc,mu,min_range,max_range)

N = size(parent_chromosome,3);
%clear kk1 kk2;
%child_1 = zeros(N,M+V,kk2);
%child_2 = zeros(N,M+V,kk2); child = zeros(N,M+V,kk2);


y=1;n=1; ncross = 0;  %ncross_p=0; %child = [];child_1 = []; child_2 = [];
for i=1:N
    rnd = rand(1);
   if y >= N
       y = y-N+1;
       n = n-N+1;
   end
           
               
   %Check whether the cross-over to be performed
    if (rnd < 0.8)
        
        % Loop over no of variables
        
        for j=1:V
            for k=1:Nc
            %select two parents
            par1 = parent_chromosome(j,k,y);
            par2 = parent_chromosome(j,k,y+1);
                                  
            yL = min_range(j);
            yU = max_range(j);
                               
            rnd = rand;
            % Check whether variable is selected or not
            if (rnd <= 0.8)
                
                ncross = ncross + 1;
                % Variable selected
                if(abs(par1-par2) > 0.000001)
                    if(par2 > par1)
                        y2 = par2;
                        y1 = par1;
                    else
                        y2 = par1;
                        y1 = par2;
                    end
                    if((y1-yL) > (yU - y2))
                        beta = 1 + abs(2*(yU-y2)/(y2-y1));
                    else
                        beta = 1 + abs(2*(y1-yL)/(y2-y1));
                    end
                    
                    % Find alpha
                    expp = mu + 1;
                    beta = 1/beta;
                   
                    alpha = 2 - beta^expp;
                    
                    rnd = rand;
                    if(rnd <= 1/alpha)
                        alpha = alpha*rnd;
                        expp = 1/(mu+1);
                        betaq = alpha^expp;
                    else
                        alpha = alpha*rnd;
                        alpha = 1/(2-alpha);
                        expp = 1/(mu+1);
                        betaq = alpha^expp;
                    end
                    
                    % Generating two children
                    child1 = 0.5*((y1+y2) - betaq*(y2-y1));
                    child2 = 0.5*((y1+y2) + betaq*(y2-y1));
                else
                    betaq = 1;
                    y1 = par1;
                    y2 = par2;
                    
                    % Generating two children
                    child1 = 0.5*((y1+y2) - betaq*(y2-y1));
                    child2 = 0.5*((y1+y2) + betaq*(y2-y1));
                end
                
                if(child1 < yL) child1 = yL; end
                if(child1 > yU) child1 = yU; end
                if(child2 < yL) child2 = yL; end
                if(child2 > yU) child2 = yU; end
            else
                % Copying the children to parents
                child1 = par1;
                child2 = par2;
            end
            child_1(j,k,i) = child1;
            child_2(j,k,i) = child2;
            
            end
        end
        
        % Evaluate the objective function for the offspring
        child_1(V+1:M+V,1,i) = evaluate_objective(child_1(1:V,:,i),V,Q,Nc);
        child_2(V+1:M+V,1,i) = evaluate_objective(child_1(1:V,:,i),V,Q,Nc);
        
        child(1:M+V,:,n) = child_1(1:M+V,:,i);
        child(1:M+V,:,n+1) = child_2(1:M+V,:,i);
        %n = n+2;
    else
        child(1:M+V,:,n) = parent_chromosome(1:M+V,:,y);
        child(1:M+V,:,n+1) = parent_chromosome(1:M+V,:,y+1);
    end
    n = n+2;
    y = y+2;
  
end

 f = child;       
        
            
                
                    
                    
                    