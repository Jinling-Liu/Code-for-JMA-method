clear
clc
M=[0.75 0.3 0.5;0 0.5 0;0.4 0 0];
d=length(M);
%===============n=20===================
n=30;
Z0=5;Z0=Z0*ones(1,d);N=100;
for simulation_time=1:N
    tic
    simulation_time1=simulation_time
    rng(simulation_time);
    Zn = generate_mbp(Z0,M,n);
    tic;
    result{simulation_time}=SimulationNPB(Zn);
    T_NPB(simulation_time)=toc;
    tic;
    M_hat_ma{simulation_time}=SimulationMbpM(Zn);
    T_MA(simulation_time)=toc;
end
mean_M_bayesian=zeros(d,d);mean_M_ma=zeros(d,d);
for j=1:N
    M_bayesian = [result{j}];
    M_ma = [M_hat_ma{j}];

    mean_M_bayesian = mean_M_bayesian+M_bayesian/N;
    mean_M_ma = mean_M_ma+M_ma/N;
end
for j=1:N
    mse_M_bayesian(j) = norm(result{j}-M)^2;
    mse_M_ma(j) = norm(M_hat_ma{j}-M)^2;

    var_M_bayesian(j) = norm(result{j}-mean_M_bayesian)^2;
    var_M_ma(j) = norm(M_hat_ma{j}-mean_M_ma)^2;
end
Result_30=[mean(mse_M_bayesian) mean(mse_M_ma) mean(var_M_bayesian) mean(var_M_ma) mean(T_NPB) mean(T_MA)]';

%===============n=25===================
% n=25;
% Z0=10;Z0=Z0*ones(1,d);N=100;
% for simulation_time=1:N
%     tic
%     simulation_time2=simulation_time
%     rng(simulation_time);
%     Zn = generate_mbp(Z0,M,n);
%     tic;
%     result{simulation_time}=SimulationNPB(Zn);
%     T_NPB(simulation_time)=toc;
%     tic;
%     M_hat_ma{simulation_time}=SimulationMbpM(Zn);
%     T_MA(simulation_time)=toc;
% end
% mean_M_bayesian=zeros(d,d);mean_M_ma=zeros(d,d);
% for j=1:N
%     M_bayesian = [result{j}];
%     M_ma = [M_hat_ma{j}];
% 
%     mean_M_bayesian = mean_M_bayesian+M_bayesian/N;
%     mean_M_ma = mean_M_ma+M_ma/N;
% end
% for j=1:N
%     mse_M_bayesian(j) = norm(result{j}-M)^2;
%     mse_M_ma(j) = norm(M_hat_ma{j}-M)^2;
% 
%     var_M_bayesian(j) = norm(result{j}-mean_M_bayesian)^2;
%     var_M_ma(j) = norm(M_hat_ma{j}-mean_M_ma)^2;
% end
% Result_25=[mean(mse_M_bayesian) mean(mse_M_ma) mean(var_M_bayesian) mean(var_M_ma) mean(T_NPB) mean(T_MA)]';
% 
% %===============n=30===================
% n=30;
% Z0=10;Z0=Z0*ones(1,d);N=100;
% for simulation_time=1:N
%     tic
%     simulation_time3=simulation_time
%     rng(simulation_time);
%     Zn = generate_mbp(Z0,M,n);
%     tic;
%     result{simulation_time}=SimulationNPB(Zn);
%     T_NPB(simulation_time)=toc;
%     tic;
%     M_hat_ma{simulation_time}=SimulationMbpM(Zn);
%     T_MA(simulation_time)=toc;
% end
% mean_M_bayesian=zeros(d,d);mean_M_ma=zeros(d,d);
% for j=1:N
%     M_bayesian = [result{j}];
%     M_ma = [M_hat_ma{j}];
% 
%     mean_M_bayesian = mean_M_bayesian+M_bayesian/N;
%     mean_M_ma = mean_M_ma+M_ma/N;
% end
% for j=1:N
%     mse_M_bayesian(j) = norm(result{j}-M)^2;
%     mse_M_ma(j) = norm(M_hat_ma{j}-M)^2;
% 
%     var_M_bayesian(j) = norm(result{j}-mean_M_bayesian)^2;
%     var_M_ma(j) = norm(M_hat_ma{j}-mean_M_ma)^2;
% end
% Result_30=[mean(mse_M_bayesian) mean(mse_M_ma) mean(var_M_bayesian) mean(var_M_ma) mean(T_NPB) mean(T_MA)]';

% [Result_20 Result_25 Result_30]
% Result_bayesian2=[Result_20 Result_25 Result_30];



save('Result_30.mat','Result_30')





