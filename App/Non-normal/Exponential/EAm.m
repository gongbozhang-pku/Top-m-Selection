function [PCS,EOC] = EAm(k,n0,T,num,m)

PCS = zeros(num,T);
EOC = zeros(num,T);

tic

parfor t = 1:num
    
    alpha0 = unifrnd(2,10,1,k);
    beta0 = unifrnd(1,2,1,k);
    lambda = gamrnd(alpha0,beta0);
    truemu = 1./lambda;
    
    X0 = exprnd(repmat(truemu,n0,1),[n0 k]);
    estmean = mean(X0);
    N = n0*ones(1,k);
    
    alpha = alpha0+n0;
    beta = beta0./(1+beta0.*n0.*estmean);
    pm = alpha.*beta;
    
    [~,rb] = sort(truemu,'descend');
    rb = sort(rb(1:m));
    
    for budget = 1:T
        [~,id4] = sort(1./pm,'descend');
        id4 = sort(id4(1:m));
        
        PCS(t,budget) = prod(id4 == rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(truemu(rb))-sum(truemu(id4)))/num;
        
        q = fix(budget/k);
        res = budget-k*q;
        if res>0
            id2 = res;
        else
            id2 = k;
        end
        
        x = exprnd(truemu(id2));
        estmean(id2) = (estmean(id2).*N(id2)+x)./(N(id2)+1);
        N(id2) = N(id2)+1;
        
        alpha(id2) = alpha(id2)+1;
        beta(id2) = beta0(id2)/(1+beta0(id2)*N(id2)*estmean(id2));
        pm(id2) = alpha(id2)*beta(id2);
    end
end
toc

PCS  =  sum(PCS);
EOC  =  sum(EOC);
end