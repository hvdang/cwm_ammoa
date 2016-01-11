% realmutation.m
% Feb 23, 2013

% Written by Hieu V. Dang

function f = realmutation(chromosome,M,V,Q,Nc,mum,min_range,max_range)

pop = size(chromosome,3);

for i=1:pop
    rnd = rand;
    % For each chromosome find whether to do mutation or not
   if (rnd < 0.01)
    for j=1:V
        
       for k=1:Nc
            
        rnd = rand;
        % For each variable find whether to do mutation or not
        y = chromosome(j,k,i);
        if(rnd < 0.5)
            yL = min_range(j);
            yU = max_range(j);
            
            if(y > yL)
                % Calculate delta
                if((y-yL) < (yU-y))
                    delta = (y-yL)/(yU-yL);
                else
                    delta = (yU - y)/(yU-yL);
                end
                
                rnd = rand;
                
                indi = 1/(mum+1);
                
                if(rnd <= 0.5)
                    xy = 1 - delta;
                    val = 2*rnd + (1-2*rnd)*(xy^(mum+1));
                    deltaq = val^indi - 1;
                else
                    xy = 1-delta;
                    val = 2*(1-rnd) + 2*(rnd - 0.5)*(xy^(mum+1));
                    deltaq = 1 - val^indi;
                end
                
                % Change the value for the parent
                y = y + deltaq*(yU - yL);
                
                if(y < yL) y = yL; end
                if(y > yU) y = yU; end
                
                child(j,k,i) = y;
            else
                xy = rand(1);
                child(j,k,i) = xy*(yU - yL) + yL;
            end
        else
            child(j,k,i) = y;
        end
      end
    end
    chromosome(1:V,:,i) = child(1:V,:,i);
    chromosome(V+1:M+V,1,i) = evaluate_objective(chromosome(1:V,:,i),V,Q,Nc); 
    end
    
   %end
end
f = chromosome;


                
                    
    