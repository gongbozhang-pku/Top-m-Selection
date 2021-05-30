function [PCS,EOC]=EAm(k,n0,T,num,m,mu,PlanN1,PlanN2,PlanN3,PlanN4,TotalPlan)

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
    pv=(N./estvar).^(-1);
    pm=pv.*(N.*estmean./estvar);
    
    [~,rb]=sort(mu);
    rb=sort(rb(1:m));
    
    for budget=1:T
        [~,id4]=sort(pm);
        id4=sort(id4(1:m));
        if id4==rb
            PCS(budget)=PCS(budget)+1/num;
        end
        EOC(budget)=EOC(budget)+(sum(mu(id4))-sum(mu(rb)))/num;
        
        q=fix(budget/k);
        res=budget-k*q;
        if res>0
            id2=res;
        else
            id2=k;
        end
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
        
        pv(id2)=(N(id2)./estvar(id2))^(-1);
        pm(id2)=pv(id2).*(N(id2).*estmean(id2)./estvar(id2));
    end
end
toc
end