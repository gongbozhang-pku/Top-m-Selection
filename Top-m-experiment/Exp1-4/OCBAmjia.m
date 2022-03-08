function [PCS,EOC]=OCBAmjia(k,n0,T,mu0,sigma0,v,num,m)

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
    estvar=var(X0,1);
    %estvar=var(X0);
    N=n0*ones(1,k);
    %     pv=(1./sigma0+N./estvar).^(-1);
    %     pm=pv.*(mu0./sigma0+N.*estmean./estvar);
    Nrv=normrnd(0,1,1,T);
    [~,rb]=sort(mu,'descend');
    %[~,rb]=sort(truemu);
    rb=sort(rb(1:m));
    wt=ones(1,k)/k;
    
    for budget=1:T
        [~,id1]=sort(estmean,'descend');
        %[~,id1]=sort(estmean);
        %[~,id1]=sort(pm,'descend');
        %[~,id1]=sort(pm);
        id3=sort(id1(1:m));
        if id3==rb
            PCS(budget)=PCS(budget)+1/num;
        end
        EOC(budget)=EOC(budget)+(sum(mu(rb))-sum(mu(id3)))/num;
        
        [~,s]=min(estmean,[],2);
        [B,id4]=sort(estmean,'descend');
        
        a = min(B(m)-B(m-1),B(m+1)-B(m));
        b = min(B(m+1)-B(m),B(m+2)-B(m+1));
        if a >= b
            w1 = zeros(1,k);
            Omega1=setdiff((1:k),id4(m));
            alpha1=zeros(1,k);
            alpha1(Omega1)=(estmean(Omega1)-B(m)).^2;
            AW11=(estvar(Omega1)*alpha1(s))./(estvar(s)*alpha1(Omega1));
            AW21=(estvar(id4(m))./estvar(Omega1)).*AW11.^2;
            w1(s)=1/((sum(AW21))^(1/2)+sum(AW11));
            w1(Omega1)=AW11*w1(s);
            w1(id4(m))=w1(s)*(sum(AW21))^(1/2);
            w = w1;
        else
            w2 = zeros(1,k);
            Omega2=setdiff((1:k),id4(m+1));
            alpha2=zeros(1,k);
            alpha2(Omega2)=(estmean(Omega2)-B(m+1)).^2;
            AW12=(estvar(Omega2)*alpha2(s))./(estvar(s)*alpha2(Omega2));
            AW22=(estvar(id4(m+1))./estvar(Omega2)).*AW12.^2;
            w2(s)=1/((sum(AW22))^(1/2)+sum(AW12));
            w2(Omega2)=AW12*w2(s);
            w2(id4(m+1))=w2(s)*(sum(AW22))^(1/2);
            w = w2;
        end
        
        df=(k*n0+budget)*w-(k*n0+budget-1)*wt;
        %         for i=1:k
        %             qujian (i) = sum (w(1:i));
        %         end
        %         suijishu = rand(1);
        %         id2 = length(qujian(suijishu > qujian)) + 1;
        
        [~,id2]=max(df);
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