function [M_Lse,M_wLse]=SimulationLseWlseM(Zn)
%calculate the unweighted LSE and weighted LSE of offspring mean M based on Zn

%unweighted LSE
[n,d]=size(Zn);n=n-1;
X=Zn(1:n,:);Y=Zn(2:n+1,:);

for ii=1:d
model = fitlm(X,Y(:,ii));
beta0 = model.Coefficients.Estimate;
M_Lse(:,ii)=beta0(2:end);
end

%weighted LSE
Xw=Zn(1:n,:)./repmat(sqrt(sum(Zn(1:n,:),2)),1,d);
Yw=Zn(2:n+1,:)./repmat(sqrt(sum(Zn(1:n,:),2)),1,d);

for ii=1:d
model = fitlm(Xw,Yw(:,ii));
beta0 = model.Coefficients.Estimate;
M_wLse(:,ii)=beta0(2:end);
end
end