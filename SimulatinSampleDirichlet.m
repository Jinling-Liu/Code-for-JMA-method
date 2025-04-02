function samples = SimulatinSampleDirichlet(alpha)
    % generate num_samples samples from Dirichlet(alpha)
    % alpha is (1×K) parameter vector
    % num_samples is the sample size
    K = length(alpha);
    gamma_samples = gamrnd(repmat(alpha, 1, 1), 1, 1, K);
    samples = gamma_samples ./ sum(gamma_samples, 2);
end
