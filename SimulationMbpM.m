function [M_hat_ma, CI_ma]=SimulationMbpM(Zn) 
%calculate the model averaging estimation of offspring mean M based on Zn
%input: Zn<-observations of d-type branching process from generations 0 to n-1

[n,d]=size(Zn);n=n-1;alpha=0.05;
pi_panety=[n*[1:d]'./(n-[1:d]'-2)];

X=Zn(1:n,:)./repmat(sqrt(sum(Zn(1:n,:),2)),1,d);
Y=Zn(2:n+1,:)./repmat(sqrt(sum(Zn(1:n,:),2)),1,d);

for j=1:d
    model = fitlm(X,Y(:,j));
    hat_eta = model.Coefficients.Estimate;
    hat_eta = hat_eta(1);
    Y(:,j) = Y(:,j)-hat_eta;
end


for type_number=1:d
for i=1:d
    for j=1:d
        Pi_i=[eye(i) zeros(i,d-i)]';Pi_j=[eye(j) zeros(j,d-j)]';
        Xi=X*Pi_i;Xj=X*Pi_j;
        Pi=Xi*(Xi'*Xi)^(-1)*Xi';Pj=Xj*(Xj'*Xj)^(-1)*Xj';
        Di=diag(1./diag((eye(n)-diag(diag(Pi))))); Dj=diag(1./diag((eye(n)-diag(diag(Pj)))));
        Phi_Matrix(i,j)=Y(:,type_number)'*(eye(n)-Pi)*Di*Dj*(eye(n)-Pj)*Y(:,type_number);
    end
end
Phi_main=Phi_Matrix;
w0=ones(d,1)/d;

%determine the parameter phi_n (the penalty term)
best_para_phi = 1;
best_obj_value = inf;
for para_phi = 0:0.001:0.1
    warning('off', 'all');
    options = optimset('Display', 'off');
    [~, fval] = quadprog(Phi_main,para_phi*pi_panety,[],[],ones(d,d),ones(d,1),zeros(d,1),ones(d,1),w0,options);
    if fval < best_obj_value
        best_obj_value = fval;
        best_para_phi = para_phi;
    end
end

%calculate the model averaging weights
warning('off', 'all');
options = optimset('Display', 'off');
w(:,type_number) = quadprog(Phi_main,best_para_phi*pi_panety,[],[],ones(d,d),ones(d,1),zeros(d,1),ones(d,1),w0,options);
end

%calculate the model averaging estimator of M
CI_ma=[];
for j = 1:d
    for i = 1:d
        Pi_i = [eye(i) zeros(i,d-i)]';
        Xi = X*Pi_i;
        beta0 = (Xi'*Xi)\(Xi'*Y(:,j));
        r_ma = Y(:,j)-Xi*beta0;
        sigma2_ma = (r_ma' * r_ma) / (size(X,1) - size(X,2));
        V_list{i} = sigma2_ma * inv(Xi' * Xi);
        beta(:,i) = [beta0;zeros(d-i,1)];
    end
   for label_v=1:d
      V_ma = w(label_v,j)^2*V_list{label_v};
   end

   se_ma = sqrt(diag(V_ma));
   z_value = norminv(1 - alpha/2);
    
   M_hat_ma(:,j)=sum(repmat(w(:,j)',d,1).*beta,2);
   CI_ma = [CI_ma;[M_hat_ma(:,j) - z_value * se_ma, M_hat_ma(:,j) + z_value * se_ma]];
end
end























