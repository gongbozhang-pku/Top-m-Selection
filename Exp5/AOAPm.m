function [PCS,EOC] = AOAPm(k,n0,T,sigma,num,m,truemu)

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
        
        [~,id1] = sort(pm);
        id1 = sort(id1(1:m));
        
        PCS(t,budget) = prod(id1==rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(truemu(id1))-sum(truemu(rb)))/num;
        
        cm = setdiff((1:k),id1);
        
        V = zeros(1,k);
        
        for i=1:k
            nv = pv;
            M = N;
            M(i) = N(i)+1;
            nv(i) = (M(i)/estvar(i))^(-1);
            a = 1:m;
            b = 1:(k-m);
            V1 = (pm(id1(a))'-pm(cm(b))).^2./(nv(id1(a))'+nv(cm(b)));
            V(i) = min(min(V1));
        end
        [~,id2]=max(V);
            
        mm = estmean(id2);
        x = truemu(id2) + (sigma(id2)) .^ (1/2) .* Nrv(budget);
        estmean(id2) = (estmean(id2) .* N(id2) + x) ./ (N(id2) + 1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        N(id2) = N(id2) + 1;
        
        pv(id2) = (N(id2)./estvar(id2)).^(-1);
        pm(id2) = pv(id2).*(N(id2).*estmean(id2)./estvar(id2));
    end
end
toc

PCS = sum(PCS);
EOC = sum(EOC);
end