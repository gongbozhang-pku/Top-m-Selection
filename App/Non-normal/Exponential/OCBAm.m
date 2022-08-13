function [PCS,EOC] = OCBAm(k,n0,T,num,m)

PCS = zeros(num,T);
EOC = zeros(num,T);

tic

parfor t = 1:num
    
    alpha0 = unifrnd(2,10,1,k);
    beta0 = unifrnd(1,2,1,k);
    lambda  =  gamrnd(alpha0,beta0);
    truemu = 1./lambda;
    
    X0  =  exprnd(repmat(truemu,n0,1),[n0 k]);
    estmean = mean(X0);
    estvar = var(X0);
    N = n0*ones(1,k);
    
    alpha = alpha0+n0;
    beta = beta0./(1+beta0.*n0.*estmean);
    pm = alpha.*beta;
    
    [~,rb] = sort(truemu,'descend');
    rb = sort(rb(1:m));
    wt = ones(1,k)/k;
    
    for budget = 1:T
        [~,id1] = sort(1./pm,'descend');
        id1 = sort(id1(1:m));
        PCS(t,budget) = prod(id1 == rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(truemu(rb))-sum(truemu(id1)))/num;
        
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
        x = exprnd(truemu(id2));
        estmean(id2) = (estmean(id2).*N(id2)+x)./(N(id2)+1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        N(id2) = N(id2)+1;
        
        alpha(id2) = alpha(id2)+1;
        beta(id2) = beta0(id2)/(1+beta0(id2)*N(id2)*estmean(id2));
        pm(id2) = alpha(id2)*beta(id2);
        wt = N./sum(N);
    end
end

toc

PCS = sum(PCS);
EOC = sum(EOC);
end