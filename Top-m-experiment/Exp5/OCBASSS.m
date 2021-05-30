function [PCS,EOC]=OCBASSS(k,n0,T,num,m,mu,PlanN1,PlanN2,PlanN3,PlanN4,TotalPlan)

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
        
        [~,id3]=sort(estmean);
        id3=id3(1:m);
        
        Sn = setdiff((1:k),id3);
        
        for i=1:m
            for j=1:(k-m)
                I(i,j)=(estmean(id3(i))-estmean(Sn(j)))^2/(estvar(id3(i))/wt(id3(i))+estvar(Sn(j))/wt(Sn(j)));
            end
        end
        
        [a,b]=find(I==min(min(I)));
        
        h = rand(1);
        if h<0.5
            id2 = id3(a);
        else
            id2 = Sn(b);
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
    end
end
toc
end