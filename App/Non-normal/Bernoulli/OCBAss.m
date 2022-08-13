function [PCS,EOC] = OCBAss(k,n0,T,num,m)

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
        
        [~,id3] = sort(estmean,'descend');
        id3 = sort(id3(1:m));
        
        Sn  =  setdiff((1:k),id3);
        u0  =  sum(N(id3).^2./estvar(id3));
        un  =  sum(N(Sn).^2./estvar(Sn));
        
        i = 1:m;
        j = 1:(k-m);
        I  =  (estmean(id3(i))'-estmean(Sn(j))).^2./(estvar(id3(i))'./wt(id3(i))'+estvar(Sn(j))./wt(Sn(j)));
        [a,b] = find(I == min(min(I)));
        
        if u0<un
            id2 = id3(a);
        else
            id2 = Sn(b);
        end
        
        mm = estmean(id2);
        x = binornd(1,theta(id2));
        alpha(id2) = alpha(id2)+x;
        beta(id2) = beta(id2)+1-x;
        pm(id2) = alpha(id2)/(alpha(id2)+beta(id2));

        estmean(id2) = (estmean(id2).*N(id2)+x)./(N(id2)+1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        N(id2) = N(id2)+1;
        wt = N./sum(N);
    end
end

toc

PCS = sum(PCS);
EOC = sum(EOC);
end