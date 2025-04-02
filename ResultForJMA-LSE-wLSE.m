clear
clc
d=5;
M=[0.3 0.4 0 0.5 0;0.5 0 0.5 0 0.6;0 0.3 0.6 0 0 ;0.4 0 0 0 0 ;0 0.3 0 0 0]';
V=[];
for i = 1:d
    V = [V;diag(M(:,i))];
end
Z0_value = 2000;Z0 = Z0_value*ones(d,1)';
n_size = [300 400 500];N = 100;

for i_n=1:3
    n=n_size(i_n);
    for simulation_time = 1:N
        [i_n simulation_time]
        rng(simulation_time); %set random seed
        Zn = generate_mbp(Z0,M,n); %generate simulation data

        M_hat_ma{i_n,simulation_time} = SimulationMbpM(Zn);
        V_hat_ma{i_n,simulation_time} = SimulationMbpV(Zn) ;

        [M_Lse{i_n,simulation_time},M_wLse{i_n,simulation_time}] = SimulationLseWlseM(Zn);
        [V_ij_lse{i_n,simulation_time},V_ij_wlse{i_n,simulation_time}] = SimulationLseWlseV(Zn);
    end
end

%=================calculate MSE and VAR=======================================
mean_M_hat_ma = zeros(d*3,d);mean_M_Lse = zeros(d*3,d);mean_M_wLse = zeros(d*3,d);
mean_V_hat_ma = zeros(d^2*3,d);mean_V_Lse = zeros(d^2*3,d);mean_V_wLse = zeros(d^2*3,d);

for j=1:N
    M_hat_ma0 = [M_hat_ma{1,j};M_hat_ma{2,j};M_hat_ma{3,j}];
    V_hat_ma0 = [V_hat_ma{1,j};V_hat_ma{2,j};V_hat_ma{3,j}];
    M_Lse0 = [M_Lse{1,j};M_Lse{2,j};M_Lse{3,j}];
    M_wLse0 = [M_wLse{1,j};M_wLse{2,j};M_wLse{3,j}];
    V_Lse0 = [V_ij_lse{1,j};V_ij_lse{2,j};V_ij_lse{3,j}];
    V_wLse0 = [V_ij_wlse{1,j};V_ij_wlse{2,j};V_ij_wlse{3,j}];

    mean_M_hat_ma = mean_M_hat_ma+M_hat_ma0/N;
    mean_V_hat_ma = mean_V_hat_ma+V_hat_ma0/N;
    mean_M_Lse = mean_M_Lse+M_Lse0/N;
    mean_M_wLse = mean_M_wLse+M_wLse0/N;
    mean_V_Lse = mean_V_Lse+V_Lse0/N;
    mean_V_wLse = mean_V_wLse+V_wLse0/N;
end



for i = 1:3
    for j=1:N
        mse_M_ma(j,i) = norm(M_hat_ma{i,j}-M)^2;
        mse_V_ma(j,i) = norm(V_hat_ma{i,j}-V)^2;
        mse_M_Lse(j,i) = norm(M_Lse{i,j}-M)^2;
        mse_M_wLse(j,i) = norm(M_wLse{i,j}-M)^2;
        mse_V_Lse(j,i) = norm(V_ij_lse{i,j}-V)^2;
        mse_V_wLse(j,i) = norm(V_ij_wlse{i,j}-V)^2;

        var_M_ma(j,i) = norm(M_hat_ma{i,j}-mean_M_hat_ma(d*(i-1)+1:i*d,:))^2;
        var_V_ma(j,i) = norm(V_hat_ma{i,j}-mean_V_hat_ma(d^2*(i-1)+1:i*d^2,:))^2;
        var_M_Lse(j,i) = norm(M_Lse{i,j}-mean_M_Lse(d*(i-1)+1:i*d,:))^2;
        var_M_wLse(j,i) = norm(M_wLse{i,j}-mean_M_wLse(d*(i-1)+1:i*d,:))^2;
        var_V_Lse(j,i) = norm(V_ij_lse{i,j}-mean_V_Lse(d^2*(i-1)+1:i*d^2,:))^2;
        var_V_wLse(j,i) = norm(V_ij_wlse{i,j}-mean_V_wLse(d^2*(i-1)+1:i*d^2,:))^2;
    end
end

mse_var = [mean(mse_M_ma);mean(mse_V_ma);mean(mse_M_Lse);mean(mse_M_wLse);mean(mse_V_Lse);mean(mse_V_wLse);...
    mean(var_M_ma);mean(var_V_ma);mean(var_M_Lse);mean(var_M_wLse);mean(var_V_Lse);mean(var_V_wLse)]

for i_m=1:N
    M_p1=M_hat_ma{1,i_m};
    M_p2=M_hat_ma{2,i_m};
    M_p3=M_hat_ma{3,i_m};
    M_plot1(i_m,:)=reshape(M_p1,1,d^2);
    M_plot2(i_m,:)=reshape(M_p2,1,d^2);
    M_plot3(i_m,:)=reshape(M_p3,1,d^2);
end

for i_m=1:N
    V_p1=V_hat_ma{1,i_m};
    V_p2=V_hat_ma{2,i_m};
    V_p3=V_hat_ma{3,i_m};
    V_plot1(i_m,:)=[diag(V_p1(1:5,:))' diag(V_p1(6:10,:))' diag(V_p1(11:15,:))' diag(V_p1(16:20,:))' diag(V_p1(21:25,:))'];
    V_plot2(i_m,:)=[diag(V_p2(1:5,:))' diag(V_p2(6:10,:))' diag(V_p2(11:15,:))' diag(V_p2(16:20,:))' diag(V_p2(21:25,:))'];
    V_plot3(i_m,:)=[diag(V_p3(1:5,:))' diag(V_p3(6:10,:))' diag(V_p3(11:15,:))' diag(V_p3(16:20,:))' diag(V_p3(21:25,:))'];
end


%=================Subplot for M=======================================
label_plot0=[1:25];
label_plot1=reshape(label_plot0,5,5);
label_plot=reshape(label_plot1',25,1);
M_r=reshape(M,1,d^2);
for s_b=1:d^2
    subplot(d,d,s_b)
    data=M_plot3(:,label_plot(s_b));
    set(gcf, 'Position', [100, 100, 1000, 800])
    histogram(data, 'Normalization', 'pdf');
    hold on;
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');
    x_values = edges(1:end-1) + diff(edges)/2;
    density_smooth = smooth(x_values, counts, 1, 'loess');
    plot(x_values, density_smooth, 'r-', 'LineWidth', 2);
    hold on
    plot([M_r(label_plot(s_b)) M_r(label_plot(s_b))], [0 8],'g-*');
end

%=================Subplot for V=======================================

M_r=reshape(M,1,d^2);
for s_b=1:d^2
    subplot(d,d,s_b)
    data=V_plot3(:,label_plot(s_b));
    set(gcf, 'Position', [100, 100, 1000, 800])
    histogram(data, 'Normalization', 'pdf');
    hold on;
    [counts, edges] = histcounts(data, 'Normalization', 'pdf');
    x_values = edges(1:end-1) + diff(edges)/2;
    density_smooth = smooth(x_values, counts, 0.5, 'loess');
    plot(x_values, density_smooth, 'r-', 'LineWidth', 2);
    hold on
    plot([M_r(s_b) M_r(s_b)], [0 1],'g-*');
end





















