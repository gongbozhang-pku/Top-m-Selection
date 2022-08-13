function [PCS,EOC] = OCBAmjia(k,n0,T,num,m)

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
    estvar = var(X0);
    N = n0*ones(1,k);
    
    alpha = alpha0+n0;
    beta = beta0./(1+beta0.*n0.*estmean);
    pm = alpha.*beta;
    
    [~,rb] = sort(truemu,'descend');
    rb = sort(rb(1:m));
    wt = ones(1,k)/k;
    
    for budget = 1:T
        [~,id3] = sort(1./pm,'descend');
        id3 = sort(id3(1:m));
        
        PCS(t,budget) = prod(id3 == rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(truemu(rb))-sum(truemu(id3)))/num;
        
        [B,id4] = sort(estmean,'descend');
        [~,s] = min(estmean,[],2);
        
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
            w  =  w1;
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
        x = exprnd(truemu(id2));
        estmean(id2) = (estmean(id2).*N(id2)+x)./(N(id2)+1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        N(id2) = N(id2)+1;
        
        alpha(id2) = alpha(id2)+1;
        beta(id2) = beta0(id2)./(1+beta0(id2).*N(id2).*estmean(id2));
        pm(id2) = alpha(id2)*beta(id2);
        wt = N./sum(N);
    end
end

toc

PCS = sum(PCS);
EOC = sum(EOC);
end