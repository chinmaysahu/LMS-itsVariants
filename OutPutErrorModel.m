%Equation Error Model
mcN = 50;
N = 10000;
a = [1 0.8 0.64 0.5120];
% a = [1 0.4 0.4*.4 .4*.4*.4];
b = 1';
c_hat = [0.1 0.1 0.1 0.1]';
s = length(c_hat);
mu = 0.005;
% Montecarlo
e = zeros(mcN,N);
for mc_loop = 1:mcN
    c_hat = [0.1 0.1 0.1 0.1]';
    x = randn(N,1);
    d = filter(b,a,x);
    beta_alpha=zeros(s,N);
    % LMS Iterations
    for LMS_loop = s+1:N-s-1
        u = zeros(1,s)';
        u(1) = x(LMS_loop);
        u(2:s) = -d(LMS_loop-1:-1:LMS_loop-s+1);
        e(mc_loop,LMS_loop) = d(LMS_loop) - u.'*c_hat;
        
        beta_alpha(1,LMS_loop) = -x(LMS_loop);
        beta_alpha(2:s,LMS_loop) = d(LMS_loop-1:-1:LMS_loop-s+1);
        
        for p=1:1:s-1
            beta_alpha(1,LMS_loop) = beta_alpha(1,LMS_loop) - c_hat(p+1)*beta_alpha(1,LMS_loop-p);
        end
        for z=2:1:s
            for p=1:1:s-1
                beta_alpha(z,LMS_loop) = beta_alpha(z,LMS_loop) - c_hat(p+1)*beta_alpha(z,LMS_loop-p);
            end
        end
        
        del = 2*e(mc_loop,LMS_loop)*beta_alpha(:,LMS_loop);
        c_hat = c_hat - 2*mu*del;
    end
   
end
MSE = mean(e(:,s:N).^2);
if mcN==1
    MSE = e(:,s:N).^2;
end
plot(s:N,db(MSE))