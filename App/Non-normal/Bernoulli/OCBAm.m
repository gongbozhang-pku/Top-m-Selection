function [PCS,EOC] = OCBAm(k,n0,T,num,m)

PCS = zeros(num,T);
EOC = zeros(num,T);

tic

parfor t = 1:num
    
    alpha0 = unifrnd(1,20,1,k);
    beta0 = unifrnd(1,20,1,k);
    theta  =  betarnd(alpha0,beta0);
    
    X0 = binornd(ones(n0,k),theta.*ones(n0,k));
    alpha = alpha0+sum(X0);
    beta = beta0+n0-sum(X0);
    pm = alpha./(alpha+beta);
    estmean = mean(X0);
    estvar = var(X0);
    N = n0*ones(1,k);
    
    [~,rb] = sort(theta,'descend');
    rb = sort(rb(1:m));
    wt = ones(1,k)/k;
    
    for budget = 1:T
        [~,id9] = sort(pm,'descend');
        id9 = sort(id9(1:m));
        
        PCS(t,budget) = prod(id9 == rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(theta(rb))-sum(theta(id9)))/num;
        
        [pmm,id4] = sort(estmean,'descend');
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
        x = binornd(1,theta(id2));
        alpha(id2) = alpha(id2)+x;
        beta(id2) = beta(id2)+1-x;
        pm(id2) = alpha(id2)/(alpha(id2)+beta(id2));

        estmean(id2) = (estmean(id2).*N(id2)+x)./(N(id2)+1); %得到replication的alternative的估计分布的新均值
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        N(id2) = N(id2)+1;
        wt = N./sum(N);
    end
end

toc

PCS = sum(PCS);
EOC = sum(EOC);
end