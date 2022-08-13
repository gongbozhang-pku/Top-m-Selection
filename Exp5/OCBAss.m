function [PCS,EOC]=OCBAss(k,n0,T,sigma,num,m,truemu)

PCS = zeros(num,T);
EOC = zeros(num,T);

tic

parfor t=1:num
    
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
        
        [~,id3] = sort(estmean);
        id3 = sort(id3(1:m));
        
        Sn = setdiff((1:k),id3);
        u0 = sum(N(id3).^2./estvar(id3));
        un = sum(N(Sn).^2./estvar(Sn));
        
        i = 1:m;
        j = 1:(k-m);
        I = (estmean(id3(i))'-estmean(Sn(j))).^2./(estvar(id3(i))'./wt(id3(i))'+estvar(Sn(j))./wt(Sn(j)));
        [a,b] = find(I==min(min(I)));
        
        
        if u0<un
            id2 = id3(a);
        else
            id2 = Sn(b);
        end
        
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