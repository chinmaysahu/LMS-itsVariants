%SHARF ALgorithm
mcN = 50; % monte carlo experiment  length

N = 10000;

b =1;
a = [1 0.8 0.64 0.512]; % true estimates
c_hat = [0 0 0 0]';
s = length(c_hat);
mu = 0.03;
% Montecarlo
e = zeros(mcN,N);
for mc_loop = 1:mcN
    c_hat = [0 0 0 0]';
    x = randn(N,1);
    d = filter(b,a,x);
    % LMS Iterations
    for LMS_loop = s:N-s
        u = zeros(1,s)';
        u(1) = x(LMS_loop);
        u(2:s) = -d(LMS_loop-1:-1:LMS_loop-s+1);
        e(mc_loop,LMS_loop) = d(LMS_loop) - u.'*c_hat;
        
        c_hat = c_hat + 2*mu*e(mc_loop,LMS_loop)*u;
    end
end
MSE = mean(e(:,s:N).^2);
if mcN==1
    MSE = e(:,s:N).^2;
end
plot(s:N,db(MSE))