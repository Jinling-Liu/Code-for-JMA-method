function result=SimulationNPB(Zn)
%calculate the nonparametric bayesian estimation of offspring mean M based Zn

[n,d]=size(Zn);
kk=3;N=1000;M=zeros(d,d);
%=========obtain all possible values of zlk=========================
for k=1:n-1
    for i=1:d
        ZUHE_zij{k,i} = MBP_zuhe_zi(Zn(k+1,i),d);
        ZUHE_zij_f=ZUHE_zij{k,i};
        [m,~]=size(ZUHE_zij_f);
        for j=1:d
            ZUHE_ij_s=[];
            for t=1:m
                ZUHE_ij = MBP_zueh_super_fast_group2(Zn(k,i),ZUHE_zij_f(t,j),kk);
                ZUHE_ij_s=[ZUHE_ij_s;ZUHE_ij];
            end
            ZUHE_ij_S{k,i,j}=ZUHE_ij_s;
        end
    end
end

for i=1:d
    for j=1:d
        alpha_diri{i,j} = 0.5*ones(1,kk+1);
        ZUHE_S{i,j}=zeros(1,kk+1);
    end
end


t=1;
while t<N
    
    %==========update p===========
    for i=1:d
        for j=1:d
            alpha_diri_c = alpha_diri{i,j}+sum(ZUHE_S{i,j});
            p{i,j} = SimulatinSampleDirichlet(alpha_diri_c);
        end
    end

    %=========update zlk======================
    for k=1:n-1
        for i=1:d
            for j=1:d
                [n_zuhe_ij_s,~]=size(ZUHE_ij_S{k,i,j});
                ZUHE_ij_s=ZUHE_ij_S{k,i,j};pij=p{i,j};
                weights=(factorial(Zn(k,i))./prod(factorial(ZUHE_ij_s),2)).*prod(repmat(pij,n_zuhe_ij_s,1).^(ZUHE_ij_s),2);
                weights=weights/sum(weights);
                label_ZUHE_ij_S=randsample([1:n_zuhe_ij_s]', 1, true, weights);
                choose_ZUHE_ij_S=ZUHE_ij_S{k,i,j};
                ZUHE_ij_SS{i,j}=choose_ZUHE_ij_S(label_ZUHE_ij_S,:);
            end
        end
        for i=1:d
            for j=1:d
                ZUHE_S{i,j}=[ZUHE_S{i,j};ZUHE_ij_SS{i,j}];
            end
        end
    end
    t=t+1;

    if t>900
        for i=1:d
            for j=1:d
                M0(i,j)=p{i,j}*[0:3]';
            end
        end
        M=M+M0/100;
    end

end
result=M;
end







