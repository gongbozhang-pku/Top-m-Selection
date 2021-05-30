function [PCS,EOC]=OCBAss(k,n0,T,mu0,sigma0,v,num,m)

tic

PCS=zeros(1,T);
EOC=zeros(1,T);
mu=zeros(1,k);
X0=zeros(n0,k);

for t=1:num
    for i=1:k
        mu(i)=normrnd(mu0(i),(sigma0(i))^(1/2));
    end
    for i=1:n0
        X0(i,:)=normrnd(mu,v.^(1/2));
    end
    estmean=mean(X0);
    estvar=var(X0,1);
    N=n0*ones(1,k);
    %     pv=(1./sigma0+N./estvar).^(-1);
    %     pm=pv.*(mu0./sigma0+N.*estmean./estvar);
    Nrv=normrnd(0,1,1,T);
    [~,rb]=sort(mu,'descend');
    %[~,rb]=sort(mu);
    rb=sort(rb(1:m));
    wt=ones(1,k)/k;
    
    for budget=1:T
        [~,id1]=sort(estmean,'descend');
        %[~,id1]=sort(estmean);
        %[~,id1]=sort(pm,'descend');
        %[~,id1]=sort(pm);
        id1=sort(id1(1:m));
        
        if id1==rb
            PCS(budget)=PCS(budget)+1/num;
        end
        EOC(budget)=EOC(budget)+(sum(mu(rb))-sum(mu(id1)))/num;
        
        [~,id3]=sort(estmean,'descend');
        id3=sort(id3(1:m));
        
        Sn = setdiff((1:k),id3);
        u0 = sum(N(id3).^2./estvar(id3));
        un = sum(N(Sn).^2./estvar(Sn));
        
        for i=1:m
            for j=1:(k-m)
                I(i,j)=(mu(id3(i))-mu(Sn(j)))^2/(estvar(id3(i))/wt(id3(i))+estvar(Sn(j))/wt(Sn(j)));
            end
        end
        [a,b]=find(I==min(min(I)));
        
        if u0<un
            id2=id3(a);
        else
            id2=Sn(b);
        end
        
        mm=estmean(id2);
        x=mu(id2)+(v(id2)).^(1/2).*Nrv(budget);
        estmean(id2)=(estmean(id2).*N(id2)+x)./(N(id2)+1);
        estvar(id2)=(N(id2)./(N(id2)+1)).*(estvar(id2)+(mm-x).^2./(N(id2)+1));
        N(id2)=N(id2)+1;
        wt=N./sum(N);
        
        %         pv(id2)=(1./sigma0(id2)+N(id2)./estvar(id2))^(-1);
        %         pm(id2)=pv(id2).*(mu0(id2)./sigma0(id2)+N(id2).*estmean(id2)./estvar(id2));
    end
end
toc
end