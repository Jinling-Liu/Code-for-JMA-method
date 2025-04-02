function [V_ij_lse,V_ij_wlse]=SimulationLseWlseV(Zn)
%calculate the unweighted LSE and weighted LSE of offspring covariance V based on Zn

%unweighted LSE
[n,d] = size(Zn);n = n-1;
X = Zn(1:n,:);Y = Zn(2:n+1,:);X=[ones(n,1) X];
for ii=1:d
    beta0=((X'*X)^-1)*(X'*Y(:,ii));
    U(:,ii) = Y(:,ii)-X*beta0;
end
for i=1:d
    for j=1:d
        U_ij{i,j}=U(:,i).*U(:,j);
    end
end
V1=[];
for label_i=1:d
    for label_j=1:d
        beta0 = (X'*X)^(-1)*(X'*U_ij{label_i,label_j});
        V1=[V1;beta0(2:end)'];
    end
end

%we only estimate the nonzero part, diagonal element matrices V_1,...,V_d of V 
V_ij_lse=[];
for i=1:d
V_ij_lse=[V_ij_lse;reshape(V1(:,i),d,d)];
end

%================================================================================
%weighted LSE
Xm = Zn(1:n,:)./repmat(sqrt(sum(Zn(1:n,:),2)),1,d);
Xw = Zn(1:n,:)./repmat((sum(Zn(1:n,:),2)),1,d);
Yw = Zn(2:n+1,:)./repmat(sqrt(sum(Zn(1:n,:),2)),1,d);

for ii=1:d
    beta0=((Xm'*Xm)^(-1))*(Xm'*Yw(:,ii));
    Uw(:,ii) = Yw(:,ii)-Xm*beta0;
end
for i=1:d
    for j=1:d
        U_wij{i,j}=(Uw(:,i).*Uw(:,j));
    end
end
V2=[];Xw=[ones(n,1) Xw];
for label_i=1:d
    for label_j=1:d
        beta0 = (Xw'*Xw)^(-1)*(Xw'*U_wij{label_i,label_j});
        V2 = [V2;beta0(2:end)'];
    end
end


%we only estimate the nonzero part, diagonal element matrices V_1,...,V_d of V 
V_ij_wlse = [];
for i = 1:d
    V_ij_wlse = [V_ij_wlse;reshape(V2(:,i),d,d)];
end

end