function [PCS,EOC]=OCBAmsa(k,n0,T,num,m,mu,PlanN1,PlanN2,PlanN3,PlanN4,TotalPlan)

PCS=zeros(1,T);
EOC=zeros(1,T);
X0=zeros(n0,k);

tic

for t=1:num
    
    for i=1:81
        Route = zeros(8,5);
        Route((1:2),:) = PlanN1(TotalPlan(i,(1:2)),:);
        Route((3:4),:) = PlanN2(TotalPlan(i,(3:4)),:);
        Route((5:6),:) = PlanN3(TotalPlan(i,(5:6)),:);
        Route((7:8),:) = PlanN4(TotalPlan(i,(7:8)),:);
        for j=1:n0
            X0(j,i) = evacuation(Ntpmu,Ntpsigma,Ncpmu,Ncpsigma,Route);
        end
    end
    
    estmean=mean(X0);
    estvar=var(X0,1);
    %estvar=var(X0);
    N=n0*ones(1,k);
    
    [~,rb]=sort(mu);
    rb=sort(rb(1:m));
    wt=ones(1,k)/k;
    
    for budget=1:T
        [~,id1]=sort(estmean);
        id1=sort(id1(1:m));
        if id1==rb
            PCS(budget)=PCS(budget)+1/num;
        end
        EOC(budget)=EOC(budget)+(sum(mu(id1))-sum(mu(rb)))/num;
        
        [pmm,id4]=sort(estmean);
        c=(estvar(id4(m+1))*pmm(m)+estvar(id4(m))*pmm(m+1))/(estvar(id4(m))+estvar(id4(m+1)));
        
        w = zeros(1,k);
        delta = estmean-c;
        Omega=1:k;
        AW1=(estvar(Omega).*delta(k).^2)./(estvar(k).*delta(Omega).^2);
        w(k)=1/sum(AW1);
        w(Omega)=AW1*w(k);
        
        df=(k*n0+budget)*w-(k*n0+budget-1)*wt;
        
        [~,id2]=max(df);
        
        mm=estmean(id2);
        Route = zeros(8,5);
        Route((1:2),:) = PlanN1(TotalPlan(id2,(1:2)),:);
        Route((3:4),:) = PlanN2(TotalPlan(id2,(3:4)),:);
        Route((5:6),:) = PlanN3(TotalPlan(id2,(5:6)),:);
        Route((7:8),:) = PlanN4(TotalPlan(id2,(7:8)),:);
        x=evacuation(Ntpmu,Ntpsigma,Ncpmu,Ncpsigma,Route);
        estmean(id2)=(estmean(id2).*N(id2)+x)./(N(id2)+1);
        estvar(id2)=(N(id2)./(N(id2)+1)).*(estvar(id2)+(mm-x).^2./(N(id2)+1));
        N(id2)=N(id2)+1;
        wt=N./sum(N);
    end
end
toc
end