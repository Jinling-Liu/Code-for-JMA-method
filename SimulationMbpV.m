function V_ij_ma=SimulationMbpV(Zn) 
%calculate the model averaging estimation of offspring covariance V based on Zn
%input: Zn<-observations of d-type branching process from generations 0 to n-1

[n,d]=size(Zn);n=n-1;
pi_panety=[n*[1:d]'./(n-[1:d]'-2)];

Xm=Zn(1:n,:)./repmat(sqrt(sum(Zn(1:n,:),2)),1,d);
X=Zn(1:n,:)./repmat((sum(Zn(1:n,:),2)),1,d);
Y=Zn(2:n+1,:)./repmat(sqrt(sum(Zn(1:n,:),2)),1,d);
Xm=[ones(n,1) Xm];
for ii=1:d
    beta0=(Xm'*Xm)^(-1)*(Xm'*Y(:,ii));
    U(:,ii) = Y(:,ii)-Xm*beta0;
end

for i = 1:d
    for j = 1:d
        Uw{i,j} =U(:,i).*U(:,j);
    end
end

for type_number1 = 1:d
for type_number2 = 1:d   
for i = 1:d
    for j = 1:d
        Pi_i = [eye(i) zeros(i,d-i)]';Pi_j=[eye(j) zeros(j,d-j)]';
        Xi = X*Pi_i;Xj = X*Pi_j;
        Pi = Xi*(Xi'*Xi)^(-1)*Xi';Pj = Xj*(Xj'*Xj)^(-1)*Xj';
        Di = diag(1./diag((eye(n)-diag(diag(Pi))))); Dj = diag(1./diag((eye(n)-diag(diag(Pj)))));
        Phi_Matrix(i,j) = Uw{type_number1,type_number2}'*(eye(n)-Pi)*Di*Dj*(eye(n)-Pj)*Uw{type_number1,type_number2};
    end
end
Phi_main = Phi_Matrix;
w0 = ones(d,1)/d;

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
w{type_number1,type_number2} = quadprog(Phi_main,best_para_phi*pi_panety,[],[],ones(d,d),ones(d,1),zeros(d,1),ones(d,1),w0,options);
end
end

%calculate the model averaging estimator of M
V_hat_ma=[];
for type_number1 = 1:d
    for type_number2=1:d
        for i = 1:d
            Pi_i = [eye(i) zeros(i,d-i)]';
            Xi = X*Pi_i;
            beta0 = (Xi'*Xi)^(-1)*(Xi'*Uw{type_number1,type_number2});
            beta(:,i) = [beta0;zeros(d-i,1)];
        end
        V_hat_ma=[V_hat_ma sum(repmat(w{type_number1,type_number2}',d,1).*beta,2)];
    end
end
V_hat_ma=V_hat_ma';
V_ij_ma=[];
for i=1:d
V_ij_ma=[V_ij_ma;reshape(V_hat_ma(:,i),d,d)];
end
end

