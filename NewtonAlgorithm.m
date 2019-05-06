% LMS Example

mcN = 50;

N = 10000;

a =1;
b = [1 0.2 0 -0.8];

% bn = [1 -0.8];

b_hat = [0 0 0 0]';
s = length(b_hat);

mu = 0.05;

% Montecarlo
e = zeros(mcN,N);
for mc_loop = 1:mcN
    x = randn(N,1);
    % x = filter(bn,1,w)/100;
    
    d = filter(b,a,x);
    
    % LMS Iterations
    
    for LMS_loop = s:N
        x_e = x(LMS_loop:-1:LMS_loop-s+1);
        e(mc_loop,LMS_loop) = d(LMS_loop) - b_hat'*x_e;
        b_hat = b_hat + 2*mu*eye(4)*x_e*e(mc_loop,LMS_loop);
    end
    
    b_hat
%     plot(s:N,db(e(s:N).^2))
 
end

MSE = mean(e(:,s:N).^2);
box on
hold on
plot(s:N,db(MSE))