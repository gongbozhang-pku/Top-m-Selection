function [PCS,EOC]=OCBAmsa(k,n0,T,sigma,num,m,truemu)

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
    wt = ones(1,k)/k;
    
    for budget = 1:T
        [~,id1] = sort(pm);
        id1 = sort(id1(1:m));
        
        PCS(t,budget) = prod(id1==rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(truemu(id1))-sum(truemu(rb)))/num;
        
        [pmm,id4] = sort(estmean);
        c = (estvar(id4(m+1))*pmm(m)+estvar(id4(m))*pmm(m+1))/(estvar(id4(m))+estvar(id4(m+1)));
        
        w = zeros(1,k);
        delta = estmean-c;
        Omega = 1:k;
        AW1 = (estvar(Omega).*delta(k).^2)./(estvar(k).*delta(Omega).^2);
        w(k) = 1/sum(AW1);
        w(Omega) = AW1*w(k);
        
        df = (k*n0+budget)*w-(k*n0+budget-1)*wt;
        
        [~,id2] = max(df);
        
        mm = estmean(id2);
        x = truemu(id2) + (sigma(id2)) .^ (1/2) .* Nrv(budget);
        estmean(id2) = (estmean(id2) .* N(id2) + x) ./ (N(id2) + 1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        N(id2) = N(id2) + 1;
        
        pv(id2) = (N(id2)./estvar(id2)).^(-1);
        pm(id2) = pv(id2).*(N(id2).*estmean(id2)./estvar(id2));
        wt = N./sum(N);
    end
    
end
toc

PCS = sum(PCS);
EOC = sum(EOC);
end