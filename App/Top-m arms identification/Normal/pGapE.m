function [PCS,EOC]=pGapE(k,n0,T,mu0,sigma0,sigma,num,m,c)

PCS=zeros(num,T);
EOC=zeros(num,T);

parfor_progress(num);

tic

parfor t=1:num

    truemu = normrnd(mu0, sqrt(sigma0));

    X0 = normrnd(truemu .* ones(n0,k), sqrt(sigma) .* ones(n0,k));

    estmean=mean(X0);
    estvar=var(X0);

    N=n0*ones(1,k);

    pv=(1./sigma0+N./estvar).^(-1);
    pm=pv.*(mu0./sigma0+N.*estmean./estvar);

    pvv = pv;
    pmm = pm;

    estestmean = estmean;
    estestvar = estvar;

    NN = N;

    TT = (1+T)*T/2;

    Nrv=normrnd(0,1,1,TT);
    [~,rb]=sort(truemu,'descend');
    rb=sort(rb(1:m));

    for budget=1:T

        [~,id4]=sort(pm,'descend');
        id4=sort(id4(1:m));

        PCS(t,budget)=prod(id4==rb) / num;
        EOC(t,budget)=EOC(t,budget)+(sum(truemu(rb))-sum(truemu(id4)))/num;

        pv = pvv;
        pm = pmm;

        estmean = estestmean;
        estvar = estestvar;

        N = NN;

        for draw = 1:budget

            [~,id1]=sort(estmean,'descend');
            Delta1 = estmean(id1(1:m))-estmean(id1(m+1));
            Delta2 = estmean(id1(m))-estmean(id1((m+1):k));
            Delta = [Delta1 Delta2];

            H = sum(1./(Delta).^2);

            [~,id2] = max(-Delta+sqrt(c*budget/H./N(id1)));

            id2 = id1(id2);

            mm=estmean(id2);
            x=truemu(id2)+(sigma(id2)).^(1/2).*Nrv((budget-1)*budget/2+draw);
            estmean(id2)=(estmean(id2).*N(id2)+x)./(N(id2)+1);
            estvar(id2)=((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
            N(id2)=N(id2)+1;

            pv(id2) = (1 ./ sigma0(id2) + N(id2) ./ estvar(id2)) .^ (-1);
            pm(id2) = pv(id2) .* (mu0(id2) ./ sigma0(id2) + N(id2) .* estmean(id2) ./ estvar(id2));
        end
    end

    parfor_progress;

end
parfor_progress(0);
toc
PCS = sum(PCS);
EOC = sum(EOC);
end