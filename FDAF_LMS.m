% FDAF LMS Example

mcN =50; % monte carlo experiment  length

N = 10000;

a =1;
b = [1 0.2 0 -0.8 ]; % true estimates
B_hat = fft([0 0 0 0])'; %initial b estimates
s = length(b_hat);
M=4; %block size
padding=2*s;
mu = 0.01; % step size

% Montecarlo
e = zeros(s,N/M,mcN-1);
for mc_loop = 1:mcN
    b_hat = [0 0 0 0]';
    B_hat = fft(b_hat,padding);
    x = randn(N,1);
    d = filter(b,a,x);
    for BLOCK_loop = s:N/M-s
        X_e=diag(fft(x(BLOCK_loop*s-s:1:BLOCK_loop*s+s-1)));
        y = ifft(X_e*B_hat);
        e(:,BLOCK_loop,mc_loop) = d(BLOCK_loop*s:1:BLOCK_loop*s+s-1)-y(padding-s+1:padding);
        e_padded = zeros(1,padding)';
        e_padded(padding-s+1:padding) = e(:,BLOCK_loop,mc_loop);
        PHI = X_e'*fft(e_padded,padding);
        B_hat = B_hat + mu*PHI;
    end
    b_hat = ifft(B_hat);
    b_hat_array(:,mc_loop) = b_hat(1:s)
end
mean(b_hat_array,2)
MSE = mean(e(:,:,:).^2,3);
if mcN==1
    MSE = mean(e(:,:,mc_loop).^2);
end
plot(s:N/M,db(MSE(s:N/M)))