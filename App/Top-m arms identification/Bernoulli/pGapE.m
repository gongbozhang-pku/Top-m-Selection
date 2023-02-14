function [PCS,EOC]=pGapE(k,n0,T,num,m,c)

PCS=zeros(num,T);
EOC=zeros(num,T);

parpool('local',24);

parfor_progress(num/10);

tic

parfor t=1:num
    
    alpha0=unifrnd(1,10,1,k);
    beta0=unifrnd(1,10,1,k);
    
    theta = betarnd(alpha0,beta0);
    
    X0 = binornd(ones(n0,k),theta.*ones(n0,k));
    
    estmean=mean(X0);
    
    alpha=alpha0+sum(X0);
    beta=beta0+n0-sum(X0);
    
    N=n0*ones(1,k);
    
    pm=alpha./(alpha+beta);
    
    alphaa = alpha;
    betaa = beta;
    pmm = pm;
    estesstmean = estmean;
    NN = N;
    
    [~,rb]=sort(theta,'descend');
    rb=sort(rb(1:m));
    
    for budget=1:T
        
        [~,id4]=sort(pm,'descend');
        id4=sort(id4(1:m));
        
        PCS(t,budget)=prod(id4==rb) / num;
        EOC(t,budget)=EOC(t,budget)+(sum(theta(rb))-sum(theta(id4)))/num;
        
        alpha = alphaa;
        beta = betaa;
        
        pm = pmm;
        
        estmean = estesstmean;
        
        N = NN;
        
        for draw = 1:budget
            
            [~,id1]=sort(estmean,'descend');
            Delta1 = estmean(id1(1:m))-estmean(id1(m+1));
            Delta2 = estmean(id1(m))-estmean(id1((m+1):k));
            Delta = [Delta1 Delta2];
            
            H = sum(1./(Delta).^2);
            
            [~,id2] = max(-Delta+sqrt(c*budget/H./N(id1)));
            
            id2 = id1(id2);
            
            x=binornd(1,theta(id2));
            alpha(id2)=alpha(id2)+x;
            beta(id2)=beta(id2)+1-x;
            pm(id2)=alpha(id2)/(alpha(id2)+beta(id2));
            
            estmean(id2)=(estmean(id2).*N(id2)+x)./(N(id2)+1);
            
            N(id2)=N(id2)+1;
        end
    end
    
    if mod(t,10) == 0
        parfor_progress;
    end
end

parfor_progress(0);

poolobj = gcp('nocreate');
delete(poolobj);

toc
PCS = sum(PCS);
EOC = sum(EOC);
end