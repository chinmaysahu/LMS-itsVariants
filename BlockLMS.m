% Block LMS Example

mcN =50; % monte carlo experiment  length

N = 10000;

a =1;
b = [1 0.2 0 -0.8]; % true estimates
%% Block LMS
b_hat = [0 0 0 0]';
s = length(b_hat); %Number of new samples per iteration
M =4; %Block size
mu = 0.05;
% Montecarlo
e = zeros(M,N,mcN-1);
for mc_loop = 1:mcN
    b_hat = [0 0 0 0]';
    x = randn(N,1);
    d = filter(b,a,x);
    % LMS Iterations
    %for LMS_loop = (1+M/s):N/s-M
    for BLOCK_loop = s:N/M-s
        x_e = x(BLOCK_loop*s:1:BLOCK_loop*s+M-1);
        x_b = zeros(M,length(b_hat));
        for i=1:1:M
            for j=1:1:length(b_hat)
                x_b(i,j)=x(BLOCK_loop*s+i-j);
            end
        end
        e(:,BLOCK_loop,mc_loop) = d(BLOCK_loop*s:1:BLOCK_loop*s+M-1) - x_b*b_hat;
        b_hat = b_hat + 2*mu/M*x_b.'*e(:,BLOCK_loop,mc_loop);
    end
end

MSE = mean(mean(e(:,:,:).^2,3));
if mcN==1
    MSE = mean(e(:,:,mc_loop).^2);
end
box on
hold on
plot(s:N,db(MSE(s:N)))