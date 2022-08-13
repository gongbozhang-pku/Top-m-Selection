function [PCS,EOC]=AOAPm(k,n0,T,mu0,sigma0,v,num,m)

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
    
    for budget=1:T
        
        [~,id1] = sort(pm,'descend');
        %[~,id1] = sort(pm);
        id1 = sort(id1(1:m));
        
        PCS(t,budget) = prod(id1==rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(mu(rb))-sum(mu(id1)))/num;
        
        cm = setdiff((1:k),id1);
        
        V = zeros(1,k);
        for i = 1:k
            nv = pv;
            M = N;
            M(i) = N(i)+1;
            nv(i) = (1/sigma0(i)+M(i)/estvar(i))^(-1);
            a = 1:m;
            b = 1:(k-m);
            V1 = (pm(id1(a))'-pm(cm(b))).^2./(nv(id1(a))'+nv(cm(b)));
            V(i) = min(min(V1));
        end
        [~,id2] = max(V);
        
        mm = estmean(id2);
        x = mu(id2)+(v(id2)).^(1/2).*Nrv(budget);
        estmean(id2) = (estmean(id2).*N(id2)+x)./(N(id2)+1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        %estvar(id2) = (N(id2)./(N(id2)+1)).*(estvar(id2)+(mm-x).^2./(N(id2)+1));
        N(id2) = N(id2)+1;
        
        pv(id2) = (1./sigma0(id2)+N(id2)./estvar(id2)).^(-1);
        pm(id2) = pv(id2).*(mu0(id2)./sigma0(id2)+N(id2).*estmean(id2)./estvar(id2));
    end
end
PCS = sum(PCS);
EOC = sum(EOC);
end