
mcN = 1; % monte carlo experiment  length

N = 10000;

b =1;
a = [1 -0.8]; % true estimates
c_hat = [0 0 0 0]';
mu = 0.03;

s=3;
e = zeros(mcN,N);
for mc_loop = 1:mcN
    x = randn(N,1);
    d = filter(b,a,x);
    r = autocorr(d);
    y = zeros(1,s)';
    y(1) = r(2);
    beta = 1;
    c = r(2);
    g = r(2);
    for k = 2:s
        beta = (1-g*g)*beta;
        g = (r(k+1) - c.'*r(2:k))/beta;
        c = [g; c - g*c(k-1:-1:1)];
        y(k) = g;
    end
    c_hat_array(:,mc_loop) = [1; -c(s:-1:1)];
end
mean(c_hat_array,2)