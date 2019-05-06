% LMS Example

mcN = 50; % monte carlo experiment  length

N = 10000;

a =1;
b = [1 0.2 0 -0.8]; % true estimates

b_hat = [0 0 0 0]'; %initial b estimates
s = length(b_hat);
Rxx=eye(s);
mu = 1; % step size
norm_mu=mu/((s+1)*trace(Rxx));

% Montecarlo
e = zeros(mcN,N);
for mc_loop = 1:mcN
    x = randn(N,1);  
    d = filter(b,a,x);
%     [~,R]=corrmtx(x,10);
%     norm_mu=mu/((s+1)*trace(R));
    % LMS Iterations
    for LMS_loop = s:N
        x_e = x(LMS_loop:-1:LMS_loop-s+1);
        e(mc_loop,LMS_loop) = d(LMS_loop) - b_hat'*x_e;
        b_hat = b_hat + 2*norm_mu*x_e*e(mc_loop,LMS_loop);
    end
    b_hat
%     plot(s:N,db(e(s:N).^2))
end
MSE = mean(e(:,s:N).^2);
box on
hold on
plot(s:N,db(MSE))