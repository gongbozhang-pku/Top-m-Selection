function [PCS,EOC] = EAm(k,n0,T,num,m)

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
    
    [~,rb] = sort(theta,'descend');
    rb = sort(rb(1:m));

    for budget = 1:T
        [~,id9] = sort(pm,'descend');
        id9 = sort(id9(1:m));
        
        PCS(t,budget) = prod(id9 == rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(theta(rb))-sum(theta(id9)))/num;
        
        q = fix(budget/k);
        res = budget-k*q;
        if res>0
            id2 = res;
        else
            id2 = k;
        end
        
        x = binornd(1,theta(id2));
        alpha(id2) = alpha(id2)+x;
        beta(id2) = beta(id2)+1-x;
        pm(id2) = alpha(id2)/(alpha(id2)+beta(id2));
    end
end

toc

PCS = sum(PCS);
EOC = sum(EOC);
end