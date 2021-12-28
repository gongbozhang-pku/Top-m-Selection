function [PCS,EOC]=EAm(k,n0,T,mu0,sigma0,v,num,m)

PCS=zeros(1,T);
EOC=zeros(1,T);
mu=zeros(1,k);
X0=zeros(n0,k);

tic
for t=1:num
    for i=1:k
        mu(i)=normrnd(mu0(i),(sigma0(i))^(1/2));
    end
    for i=1:n0
        X0(i,:)=normrnd(mu,v.^(1/2));
    end
    estmean=mean(X0);
    estvar=var(X0);
    N=n0*ones(1,k);
    pv=(1./sigma0+N./estvar).^(-1);
    pm=pv.*(mu0./sigma0+N.*estmean./estvar);
    Nrv=normrnd(0,1,1,T);
    [~,rb]=sort(mu,'descend');
    %[~,rb]=sort(mu);
    rb=sort(rb(1:m));
    
    for budget=1:T
        [~,id4]=sort(pm,'descend');
        %[~,id4]=sort(pm);
        id4=sort(id4(1:m));
        if id4==rb
            PCS(budget)=PCS(budget)+1/num;
        end
        EOC(budget)=EOC(budget)+(sum(mu(rb))-sum(mu(id4)))/num;
        
        q=fix(budget/k);
        res=budget-k*q;
        if res>0
            id2=res;
        else
            id2=k;
        end
        mm=estmean(id2);
        x=mu(id2)+(v(id2)).^(1/2).*Nrv(budget);
        estmean(id2)=(estmean(id2).*N(id2)+x)./(N(id2)+1);
        estvar(id2)=(N(id2)./(N(id2)+1)).*(estvar(id2)+(mm-x).^2./(N(id2)+1));
        N(id2)=N(id2)+1;
        
        pv(id2)=(1./sigma0(id2)+N(id2)./estvar(id2))^(-1);
        pm(id2)=pv(id2).*(mu0(id2)./sigma0(id2)+N(id2).*estmean(id2)./estvar(id2));
    end
end
toc
end
