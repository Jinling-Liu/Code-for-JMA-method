function Result=var_mse(C,NN,M,d1,d2)
%calculate the MSE and VAR

for ii=1:3
    M_ma_bar=zeros(d1,d2);
    for n_mc=1:NN
        M_ma_bar=C{n_mc,ii}/NN+M_ma_bar;
    end

    for n_mc=1:NN
        mse_ma(n_mc,1)=(norm(C{n_mc,ii}-M))^(2);
        var_ma(n_mc,1)=(norm(C{n_mc,ii}-M_ma_bar))^(2);
    end

    result=[mean(mse_ma) mean(var_ma)]';

    Result(:,ii)=result;
end
end


