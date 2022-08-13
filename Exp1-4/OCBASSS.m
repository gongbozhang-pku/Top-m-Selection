function [PCS,EOC]=OCBASSS(k,n0,T,mu0,sigma0,v,num,m)

PCS = zeros(num,T);
EOC = zeros(num,T);

tic

parfor t=1:num
    
    mu = normrnd(mu0.* ones(n0,k),sqrt(sigma0).* ones(n0,k));
    
    X0 = normrnd(mu.* ones(n0,k), sqrt(v).* ones(n0,k));
    
    estmean = mean(X0);
    estvar = var(X0);
    %estvar = var(X0,1);
    
    N = n0*ones(1,k);
    pv = (1./sigma0+N./estvar).^(-1);
    pm = pv.*(mu0./sigma0+N.*estmean./estvar);
    
    Nrv = normrnd(0,1,1,T);
    
    [~,rb] = sort(mu,'descend');
    %[~,rb] = sort(mu);
    rb = sort(rb(1:m));
    wt = ones(1,k)/k;
    
    for budget=1:T
        [~,id1] = sort(pm,'descend');
        %[~,id1]=sort(estmean,'descend');
        %[~,id1]=sort(estmean);
        %[~,id1]=sort(pm);
        PCS(t,budget) = prod(id1==rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(mu(rb))-sum(mu(id1)))/num;
        
        [~,id3] = sort(estmean,'descend');
        id3 = id3(1:m);
        
        Sn = setdiff((1:k),id3);
        
        i = 1:m;
        j = 1:(k-m);
        I = (estmean(id3(i))'-estmean(Sn(j))).^2./(estvar(id3(i))'./wt(id3(i))'+estvar(Sn(j))./wt(Sn(j)));
        [a,b] = find(I==min(min(I)));
        
        h = rand(1);
        if h<0.5
            id2 = id3(a);
        else
            id2 = Sn(b);
        end
        
        mm = estmean(id2);
        x = mu(id2)+(v(id2)).^(1/2).*Nrv(budget);
        estmean(id2) = (estmean(id2) .* N(id2) + x) ./ (N(id2) + 1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        %estvar(id2) = (N(id2)./(N(id2)+1)).*(estvar(id2)+(mm-x).^2./(N(id2)+1));
        N(id2) = N(id2)+1;
        wt = N./sum(N);
        
        pv(id2) = (1./sigma0(id2)+N(id2)./estvar(id2)).^(-1);
        pm(id2) = pv(id2).*(mu0(id2)./sigma0(id2)+N(id2).*estmean(id2)./estvar(id2));
    end
end
toc

PCS = sum(PCS);
EOC = sum(EOC);
end