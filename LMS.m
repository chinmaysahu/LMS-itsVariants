% LMS Example

mcN = 50; % monte carlo experiment  length

N = 10000;

b =1;
a = [1 -0.8]; % true estimates
% b_hat = [0 0]'; %initial b estimates L=2
% b_hat = [0 0 0 0]'; %initial b estimates L=4
b_hat = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]'; %initial b estimates L=10
s = length(b_hat);

mu = 0.01; % step size

% Montecarlo
e = zeros(mcN,N);
for mc_loop = 1:mcN
    x = randn(N,1);
    d = filter(b,a,x);
    % LMS Iterations
    for LMS_loop = s:N
        x_e = x(LMS_loop:-1:LMS_loop-s+1);
        e(mc_loop,LMS_loop) = d(LMS_loop) - b_hat'*x_e;
        b_hat = b_hat + 2*mu*x_e*e(mc_loop,LMS_loop);
    end
    b_hat
%     plot(s:N,db(e(s:N).^2))
end

MSE = mean(e(:,s:N).^2);

%  plot(s:N,db(MSE))

% d_hat=filter(b_hat,a,x);
% figure
% freqz(b,a);
% figure
% freqz(b_hat,a);
% figure
% impz(b,a);
figure
impz(b_hat,a);
