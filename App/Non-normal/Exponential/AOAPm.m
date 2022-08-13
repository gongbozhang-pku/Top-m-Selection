function [PCS,EOC] = AOAPm(k,n0,T,num,m)

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
    pv = alpha.*beta.^2;
    
    [~,rb] = sort(truemu,'descend');
    rb = sort(rb(1:m));
    
    for budget = 1:T
        [~,id1] = sort(1./pm,'descend');
        id1 = sort(id1(1:m));
        
        PCS(t,budget) = prod(id1 == rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(truemu(rb))-sum(truemu(id1)))/num;
        
        cm = setdiff((1:k),id1);
        
        V  =  zeros(1,k);
        
        for i = 1:k
            nv = pv;
            nv(i) = pm(i)*(1-pm(i))/(alpha(i)+beta(i)+2);
            a  =  1:m;
            b  =  1:(k-m);
            V1  =  (pm(id1(a))'-pm(cm(b))+10^(-5)).^2./(nv(id1(a))'+nv(cm(b)));
            V(i) = min(min(V1));
        end
        [~,id2] = max(V);
        
        x = exprnd(truemu(id2));
        estmean(id2) = (estmean(id2).*N(id2)+x)./(N(id2)+1);
        N(id2) = N(id2)+1;
        
        alpha(id2) = alpha(id2)+1;
        beta(id2) = beta0(id2)./(1+beta0(id2).*N(id2).*estmean(id2));
        pm(id2) = alpha(id2).*beta(id2);
        pv(id2) = alpha(id2).*beta(id2).^2;
    end
end

toc

PCS  =  sum(PCS);
EOC  =  sum(EOC);
end