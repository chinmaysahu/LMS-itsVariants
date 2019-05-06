%% System Identification Using Recursive Least Square (RLS) and Least Mean Square (LMS) algorithm
%% Start
clc;
clear all;
close all;
N = 1000;       % Number of samples
Bits = 2;      % For PSK modulation    
SNR = 10;       % Noise level
h = [0.9 0.2 0.5 -0.7];         % Plant impulse response
data = randi(1,N);            % Random index for input data
x = real(pskmod(data,Bits));    % Phase shit keying (PSK) modulation
r = filter(h,1,x);              % Input passed trought system(h)
d = awgn(r, SNR);               % Addition of white gaussian noise of decined SNR
%% LMS parameter
etac = 1e-3;                    % Learning rate for LMS
Wlms = zeros(size(h));                   % Initial weights of LMS
V = zeros(1,length(h));         % Input frame length of LMS
%% RLS Parameters
p=3;
lamda=1;
sigma=1;
Wrls = Wlms';                        % Initial weights of RLS
P=(sigma^-1)*eye(p+1);
U = zeros(size(Wrls));         % Input frame length of RLS
for n = 1 : N
                  
%% LMS
        V(1,2:end) = V(1,1:end-1);  % Shifting of frame window
        V(1,1) = x(n);              % Input of LMS
        
        yc = (Wlms)*V';                 % Output of LMS                           
        ec = d(n) - yc;                 % Instantaneous error of LMS 
        Wlms = Wlms +  etac * ec * V;   % Weight update rule of LMS
%% RLS
        U(2:end,1) = U(1:end-1,1);  % Shifting of frame window
        U(1,1) = x(n);              % Input of RLS
        
        y(n) = U'*Wrls;                             % Output of RLS
        alpha(n) = d(n) - y(n);                     % Instantaneous error of RLS
        g(:,n)=P*U*((lamda + U'*P*U).^-1);          % Gain vector
        P=(lamda^-1)*P - g(:,n)*U'*(lamda^-1)*P;    % RLS intermediate term
        Wrls = Wrls + alpha(n)*g(:,n);              % Weight update rule of RLS
        
    %% Normalized weight difference (NWD)
NWDc(n) = norm(Wlms-h)./norm(h);      % Normalized weight difference of LMS
NWD(n)  = norm(Wrls'-h)./norm(h);     % Normalized weight difference of RLS
end
% %% Cost function plots
figure
fsize=14; % plot text font size
plot(10*log10(NWDc),'','linewidth',4)
hold on
plot(10*log10(NWD),'r','linewidth',4)
lgh=legend(strcat('Least mean square (LMS):', int2str(SNR),' (dB)'),strcat('Recursive least square (RLS):', int2str(SNR),' (dB)'),'Location','NorthEast');
grid minor
xlabel('Iterations','FontName','Times New Roman','FontSize',fsize);
ylabel('Normalized weight difference (NWD) in (dB)','FontName','Times New Roman','FontSize',fsize);
title('Cost function (NWD vs epochs iteration)','FontName','Times New Roman','FontSize',6*fsize/5);
set(lgh,'FontName','Times New Roman','FontSize',fsize)
set(gca,'FontName','Times New Roman','FontSize',fsize)
saveas(gcf,strcat('LMS_RLS_Comparision.png'),'png')
[h;Wrls';Wlms]