clear;clc;

%% 传输参数
M = 30; % BS antennas
K = 10; % total users
n_iterations = 1e5; % 迭代次数
SNR_dB = -30:2:0; % SNR
nn = 2;
modu_num = 2^nn; % 调制阶数
J = 2; % 不活跃的用户
L = 10; % 每次每个用户发送的数据长度

%% 统计BER
errorbit_in_zf = zeros(1,length(SNR_dB));
errorbit_in_mmse = zeros(1,length(SNR_dB));
errorbit_direct_zf = zeros(1,length(SNR_dB));
errorbit_direct_mmse = zeros(1,length(SNR_dB));
errorbit_proposed_zf = zeros(1,length(SNR_dB));
errorbit_proposed_mmse = zeros(1,length(SNR_dB));

%% 仿真
parfor i = 1 : length(SNR_dB)
    snr = 10^(SNR_dB(i)/10);
    noise_var = 1/snr;
    for j = 1 : n_iterations
       %% 传输信号
        bits = randi([0 modu_num-1],K+J, L);
        xf = pskmod(bits, modu_num, pi/modu_num);
        x = xf(K, L);
        
       %% 经过信道，高斯噪声
        n = sqrt(noise_var/2)*(randn(M,L)+1i*randn(M,L)); 
        Hf = sqrt(0.5)*(randn(M,K+J)+1i*randn(M,K+J));
        H = Hf(1:M, 1:K);
        Hj = Hf(1:M, K+1:end);
        H_power = abs(H).^2;
        H_power_mean = mean(H_power,'all');
        
       %% 接收信号
        yf = Hf*xf+n;
        
       %% inapproate ZF
        w_zf = (H'*H)\H';
        xt_in_zf = w_zf*yf;
        hats_in_zf = pskdemod(xt_in_zf, modu_num, pi/modu_num);
        errorbit_in_zf(i) = errorbit_in_zf(i)+biterr(bits(1:K,:), hats_in_zf);

       %% inapproate MMSE
        Ik = eye(size(H,2));
        w_mmse = (H'*H+1/snr*Ik)\H';
        xt_in_mmse = w_mmse*yf;
        hats_in_mmse = pskdemod(xt_in_mmse, modu_num, pi/modu_num);
        errorbit_in_mmse(i) = errorbit_in_mmse(i)+biterr(bits(1:K,:), hats_in_mmse);        
        
       %% direct ZF
        w_direct_zf = (Hf'*Hf)\Hf';
        xt_direct_zf = w_direct_zf*yf;
        hats_direct_zf = pskdemod(xt_direct_zf, modu_num, pi/modu_num);
        errorbit_direct_zf(i) = errorbit_direct_zf(i)+biterr(bits ,hats_direct_zf);        
        
       %% direct MMSE
        Ik = eye(K+J);
        w_direct_mmse = (Hf'*Hf+1/snr*Ik)\Hf';
        xt_direct_mmse = w_direct_mmse*yf;
        hats_direct_mmse = pskdemod(xt_direct_mmse, modu_num, pi/modu_num);
        errorbit_direct_mmse(i) = errorbit_direct_mmse(i)+biterr(bits, hats_direct_mmse);           

       %% proposed ZF
        L21 = Hj'*H/(H'*H);
        S = Hj'*Hj - Hj'*H/(H'*H)*H'*Hj;
        inv_HfHf = zeros(K+J, K+J);
        inv_HfHf(1:K,1:K) = inv(H'*H) + L21'/(S)*L21;
        inv_HfHf(K+1:end,K+1:end) = inv(S);
        inv_HfHf(1:K,K+1:end) = -L21'/(S);
        inv_HfHf(K+1:end,1:K) = -inv(S)*L21;
        xt_proposed_zf = inv_HfHf * Hf' * yf;
        hats_proposed_zf = pskdemod(xt_proposed_zf, modu_num, pi/modu_num);
        errorbit_proposed_zf(i) = errorbit_proposed_zf(i)+biterr(bits, hats_proposed_zf);        
        
       %% proposed MMSE
        Ik = eye(K);
        U21 = (Hj'*H/(H'*H+1/snr*Ik));
        Ij = eye(J);
        V = Hj'*Hj+1/snr*Ij-Hj'*H/(H'*H+1/snr*Ik)*H'*Hj;
        inv_HfHf = zeros(K+J, K+J);
        inv_HfHf(1:K,1:K) = inv(H'*H+1/snr*Ik)+U21'/(V)*U21;
        inv_HfHf(K+1:end,K+1:end) = inv(V);
        inv_HfHf(1:K,K+1:end) = -U21'*inv(V);
        inv_HfHf(K+1:end,1:K) = -inv(V)*U21;
        xt_proposed_mmse = inv_HfHf * Hf' * yf;
        hats_proposed_mmse = pskdemod(xt_proposed_mmse, modu_num, pi/modu_num);
        errorbit_proposed_mmse(i) = errorbit_proposed_mmse(i)+biterr(bits, hats_proposed_mmse);    
    end
end

ber_in_zf  = errorbit_in_zf/(n_iterations*nn*(K)*L);
ber_in_mmse  = errorbit_in_mmse/(n_iterations*nn*(K)*L);
ber_direct_zf  = errorbit_direct_zf/(n_iterations*nn*(K+J)*L);
ber_direct_mmse  = errorbit_direct_mmse/(n_iterations*nn*(K+J)*L);
ber_proposed_zf  = errorbit_proposed_zf/(n_iterations*nn*(K+J)*L);
ber_proposed_mmse  = errorbit_proposed_mmse/(n_iterations*nn*(K+J)*L);

semilogy(SNR_dB ,ber_in_zf , 'kv-', 'Linewidth',1);
hold on;grid on;
semilogy(SNR_dB ,ber_in_mmse ,'b>-', 'Linewidth',1);
hold on;grid on;
semilogy(SNR_dB ,ber_proposed_zf,'k+-.', 'Linewidth',1);
hold on;grid on;
semilogy(SNR_dB ,ber_direct_zf , 'ko', 'Linewidth',1);
hold on;grid on;
semilogy(SNR_dB ,ber_proposed_mmse,'b--', 'Linewidth',1);
hold on; grid on;
semilogy(SNR_dB ,ber_direct_mmse ,'bh', 'Linewidth',1);
hold on;grid on;

%% 传输参数
J = 4; % 不活跃的用户

%% 统计BER
errorbit_in_zf = zeros(1,length(SNR_dB));
errorbit_in_mmse = zeros(1,length(SNR_dB));
errorbit_direct_zf = zeros(1,length(SNR_dB));
errorbit_direct_mmse = zeros(1,length(SNR_dB));
errorbit_proposed_zf = zeros(1,length(SNR_dB));
errorbit_proposed_mmse = zeros(1,length(SNR_dB));

%% 仿真
parfor i = 1 : length(SNR_dB)
    snr = 10^(SNR_dB(i)/10);
    noise_var = 1/snr;
    for j = 1 : n_iterations
       %% 传输信号
        bits = randi([0 modu_num-1],K+J, L);
        xf = pskmod(bits, modu_num, pi/modu_num);
        x = xf(K, L);
        
       %% 经过信道，高斯噪声
        n = sqrt(noise_var/2)*(randn(M,L)+1i*randn(M,L)); 
        Hf = sqrt(0.5)*(randn(M,K+J)+1i*randn(M,K+J));
        H = Hf(1:M, 1:K);
        Hj = Hf(1:M, K+1:end);
        H_power = abs(H).^2;
        H_power_mean = mean(H_power,'all');
        
       %% 接收信号
        yf = Hf*xf+n;
        
       %% inapproate ZF
        w_zf = (H'*H)\H';
        xt_in_zf = w_zf*yf;
        hats_in_zf = pskdemod(xt_in_zf, modu_num, pi/modu_num);
        errorbit_in_zf(i) = errorbit_in_zf(i)+biterr(bits(1:K,:), hats_in_zf);

       %% inapproate MMSE
       Ik = eye(size(H,2));
        w_mmse = (H'*H+1/snr*Ik)\H';
        xt_in_mmse = w_mmse*yf;
        hats_in_mmse = pskdemod(xt_in_mmse, modu_num, pi/modu_num);
        errorbit_in_mmse(i) = errorbit_in_mmse(i)+biterr(bits(1:K,:), hats_in_mmse);        
        
       %% direct ZF
        w_direct_zf = (Hf'*Hf)\Hf';
        xt_direct_zf = w_direct_zf*yf;
        hats_direct_zf = pskdemod(xt_direct_zf, modu_num, pi/modu_num);
        errorbit_direct_zf(i) = errorbit_direct_zf(i)+biterr(bits ,hats_direct_zf);        
        
       %% direct MMSE
        Ik = eye(K+J);
        w_direct_mmse = (Hf'*Hf+1/snr*Ik)\Hf';
        xt_direct_mmse = w_direct_mmse*yf;
        hats_direct_mmse = pskdemod(xt_direct_mmse, modu_num, pi/modu_num);
        errorbit_direct_mmse(i) = errorbit_direct_mmse(i)+biterr(bits, hats_direct_mmse);           

       %% proposed ZF
        L21 = Hj'*H/(H'*H);
        S = Hj'*Hj - Hj'*H/(H'*H)*H'*Hj;
        inv_HfHf = zeros(K+J, K+J);
        inv_HfHf(1:K,1:K) = inv(H'*H) + L21'/(S)*L21;
        inv_HfHf(K+1:end,K+1:end) = inv(S);
        inv_HfHf(1:K,K+1:end) = -L21'/(S);
        inv_HfHf(K+1:end,1:K) = -inv(S)*L21;
        xt_proposed_zf = inv_HfHf * Hf' * yf;
        hats_proposed_zf = pskdemod(xt_proposed_zf, modu_num, pi/modu_num);
        errorbit_proposed_zf(i) = errorbit_proposed_zf(i)+biterr(bits, hats_proposed_zf);        
        
       %% proposed MMSE
        Ik = eye(K);
        U21 = (Hj'*H*inv(H'*H+1/snr*Ik));
        Ij = eye(J);
        V = Hj'*Hj+1/snr*Ij-Hj'*H/(H'*H+1/snr*Ik)*H'*Hj;
        inv_HfHf = zeros(K+J, K+J);
        inv_HfHf(1:K,1:K) = inv(H'*H+1/snr*Ik)+U21'/(V)*U21;
        inv_HfHf(K+1:end,K+1:end) = inv(V);
        inv_HfHf(1:K,K+1:end) = -U21'*inv(V);
        inv_HfHf(K+1:end,1:K) = -inv(V)*U21;
        xt_proposed_mmse = inv_HfHf * Hf' * yf;
        hats_proposed_mmse = pskdemod(xt_proposed_mmse, modu_num, pi/modu_num);
        errorbit_proposed_mmse(i) = errorbit_proposed_mmse(i)+biterr(bits, hats_proposed_mmse);    
    end
end

ber_in_zf  = errorbit_in_zf/(n_iterations*nn*(K)*L);
ber_in_mmse  = errorbit_in_mmse/(n_iterations*nn*(K)*L);
ber_direct_zf  = errorbit_direct_zf/(n_iterations*nn*(K+J)*L);
ber_direct_mmse  = errorbit_direct_mmse/(n_iterations*nn*(K+J)*L);
ber_proposed_zf  = errorbit_proposed_zf/(n_iterations*nn*(K+J)*L);
ber_proposed_mmse  = errorbit_proposed_mmse/(n_iterations*nn*(K+J)*L);

semilogy(SNR_dB ,ber_in_zf , 'r^-', 'Linewidth',1);
hold on;grid on;
semilogy(SNR_dB ,ber_in_mmse ,'g<-', 'Linewidth',1);
hold on;grid on;
semilogy(SNR_dB ,ber_proposed_zf,'rx-.', 'Linewidth',1);
hold on;grid on;
semilogy(SNR_dB ,ber_direct_zf , 'gd', 'Linewidth',1);
hold on;grid on;
semilogy(SNR_dB ,ber_proposed_mmse,'g-.', 'Linewidth',1);
hold on;grid on;
semilogy(SNR_dB ,ber_direct_mmse ,'gh', 'Linewidth',1);

set(gcf, 'position', [200 200 400 300]);
xlabel('SNR $\overline{\gamma}$ [dB]','interpreter', 'latex');
ylabel('BER');
ylim([10^-6,1]);
legend('Inapp. ZF,\it{J=2}','Inapp. MMSE,\it{J=2}','Proposed ZF,\it{J=2}','Direct ZF,\it{J=2}',...
    'Proposed MMSE,\it{J=2}','Direct MMSE,\it{J=2}','Inapp. ZF,\it{J=4}','Inapp. MMSE,\it{J=4}',...
    'Proposed ZF,\it{J=4}','Direct ZF,\it{J=4}','Proposed MMSE,\it{J=4}','Direct MMSE,\it{J=4}');