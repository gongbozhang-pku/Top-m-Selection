function [PCS,EOC]=EAm(k,n0,T,sigma,num,m,truemu)

PCS = zeros(num,T);
EOC = zeros(num,T);

tic

parfor t = 1:num
    
    X0 = normrnd(truemu .* ones(n0,k), sqrt(sigma) .* ones(n0,k));
    estmean = mean(X0);
    estvar = var(X0);
    
    N = n0 * ones(1,k);
    
    pv = (N ./ estvar) .^ (-1);
    pm = pv .* (N .* estmean ./ estvar);
    
    Nrv = normrnd(0, 1, 1, T); 
    
    [~,rb] = sort(truemu);
    rb = sort(rb(1:m));
    

    for budget = 1:T
        
        [~,id4] = sort(pm);
        id4 = sort(id4(1:m));
        
        PCS(t,budget) = prod(id4==rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(truemu(id4))-sum(truemu(rb)))/num;
        
        q = fix(budget/k);
        res = budget-k*q;
        if res>0
            id2 = res;
        else
            id2 = k;
        end
        
        mm = estmean(id2);
        x = truemu(id2) + (sigma(id2)) .^ (1/2) .* Nrv(budget);
        estmean(id2) = (estmean(id2) .* N(id2) + x) ./ (N(id2) + 1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        N(id2) = N(id2) + 1;
        
        pv(id2) = (N(id2)./estvar(id2))^(-1);
        pm(id2) = pv(id2).*(N(id2).*estmean(id2)./estvar(id2));
    end
end
toc

PCS = sum(PCS);
EOC = sum(EOC);
end