function [PCS,EOC] = AOAPm(k,n0,T,num,m)

PCS = zeros(num,T);
EOC = zeros(num,T);

tic

parfor t = 1:num
    
    alpha0 = unifrnd(1,20,1,k);
    beta0 = unifrnd(1,20,1,k);
    theta = betarnd(alpha0,beta0);
    
    X0 = binornd(ones(n0,k),theta.*ones(n0,k));
    alpha = alpha0+sum(X0);
    beta = beta0+n0-sum(X0);
    pm = alpha./(alpha+beta);
    pv = pm.*(1-pm)./(alpha+beta+1);
    
    [~,rb] = sort(theta,'descend');
    rb = sort(rb(1:m));
    
    for budget = 1:T
        [~,id1] = sort(pm,'descend');
        id1 = sort(id1(1:m));
        PCS(t,budget) = prod(id1 == rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(theta(rb))-sum(theta(id1)))/num;
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
        
        x = binornd(1,theta(id2));
        alpha(id2) = alpha(id2)+x;
        beta(id2) = beta(id2)+1-x;
        pm(id2) = alpha(id2)./(alpha(id2)+beta(id2));
        pv(id2) = pm(id2).*(1-pm(id2))./(alpha(id2)+beta(id2)+1);
    end
end

toc

PCS  =  sum(PCS);
EOC  =  sum(EOC);
end