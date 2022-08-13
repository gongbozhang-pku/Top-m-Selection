function [PCS,EOC] = OCBAmjia(k,n0,T,sigma,num,m,truemu)
 
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
        
        [~,s] = min(estmean,[],2);
        [B,id4] = sort(estmean);
        
        a = min(B(m)-B(m-1),B(m+1)-B(m));
        b = min(B(m+1)-B(m),B(m+2)-B(m+1));
        if a >= b
            w1 = zeros(1,k);
            Omega1 = setdiff((1:k),id4(m));
            alpha1 = zeros(1,k);
            alpha1(Omega1) = (estmean(Omega1)-B(m)).^2;
            AW11 = (estvar(Omega1)*alpha1(s))./(estvar(s)*alpha1(Omega1));
            AW21 = (estvar(id4(m))./estvar(Omega1)).*AW11.^2;
            w1(s) = 1/((sum(AW21))^(1/2)+sum(AW11));
            w1(Omega1) = AW11*w1(s);
            w1(id4(m)) = w1(s)*(sum(AW21))^(1/2);
            w = w1;
        else
            w2 = zeros(1,k);
            Omega2 = setdiff((1:k),id4(m+1));
            alpha2 = zeros(1,k);
            alpha2(Omega2) = (estmean(Omega2)-B(m+1)).^2;
            AW12 = (estvar(Omega2)*alpha2(s))./(estvar(s)*alpha2(Omega2));
            AW22 = (estvar(id4(m+1))./estvar(Omega2)).*AW12.^2;
            w2(s) = 1/((sum(AW22))^(1/2)+sum(AW12));
            w2(Omega2) = AW12*w2(s);
            w2(id4(m+1)) = w2(s)*(sum(AW22))^(1/2);
            w = w2;
        end
        
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