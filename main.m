clear;clc;

%% 传输参数
M = 30; % BS antennas
K = 10; % total users
n_iterations = 1e5; % 迭代次数
SNR_dB = -30:2:0; % SNR
nn = 2;
modu_num = 2^nn; % 调制阶数
I = 2; % 不活跃的用户
L = 10; % 每次每个用户发送的数据长度

%% 统计BER
errorbit_in_zf = zeros(1,length(SNR_dB));
errorbit_in_mmse = zeros(1,length(SNR_dB));

errorbit_direct_zf = zeros(1,length(SNR_dB));
errorbit_direct_mmse = zeros(1,length(SNR_dB));

errorbit_proposed_zf = zeros(1,length(SNR_dB));
errorbit_proposed_mmse = zeros(1,length(SNR_dB));

%% Inapproate ZF
for i = 1 : length(SNR_dB)
    snr = 10^(SNR_dB(i)/10);
    noise_var = 1/snr;
    for j = 1 : n_iterations
       %% 随机活跃用户
        inactive_user = randperm(K,I);
        bits = randi([0 modu_num-1],K, L);
        x = pskmod(bits, modu_num, pi/modu_num);
        x(inactive_user,:) = 0;
        T = permutation_T(K, I, inactive_user);
        
        bits_T = T*bits;         
        bits_part = bits_T;
        bits_part = bits_part(1:K-I,:);
        xp = T*x;
        xa = xp(1:K-I, :);
        
       %% 经过信道，高斯噪声
        n = sqrt(noise_var/2)*(randn(M,L)+1i*randn(M,L)); 
        H = sqrt(0.5)*(randn(M,K)+1i*randn(M,K));
        H_power = abs(H).^2;
        H_power_mean = mean(H_power,'all');
        Hp = H * T';
        Ha = Hp(:, 1:K-I);
        Hi = Hp(:, K-I+1:end);
        y = H*x+n;
        y0 = Ha*xa+n;
        
       %% inapproate ZF
        w_zf = (H'*H)\H';
        xt_in_zf = w_zf*y;
        hats_in_zf = pskdemod(xt_in_zf, modu_num, pi/modu_num);
        hats_in_zf(inactive_user,:) = bits(inactive_user,:);
        errorbit_in_zf(i) = errorbit_in_zf(i)+biterr(bits, hats_in_zf);

       %% inapproate MMSE
        Ik = eye(size(H,2));
        w_mmse = (H'*H+1/snr*Ik)\H';
        xt_in_mmse = w_mmse*y;
        hats_in_mmse = pskdemod(xt_in_mmse, modu_num, pi/modu_num);
        hats_in_mmse(inactive_user,:) = bits(inactive_user,:);
        errorbit_in_mmse(i) = errorbit_in_mmse(i)+biterr(bits, hats_in_mmse);        
        
       %% direct ZF
        w_zf = (Ha'*Ha)\Ha';
        xt_direct_zf = w_zf*y0;
        hats_direct_zf = pskdemod(xt_direct_zf, modu_num, pi/modu_num);
        errorbit_direct_zf(i) = errorbit_direct_zf(i)+biterr(bits_part ,hats_direct_zf);        
        
       %% direct MMSE
        Ik = eye(K-I);
        w_mmse = (Ha'*Ha+1/snr*Ik)\Ha';
        xt_direct_mmse = w_mmse*y0;
        hats_direct_mmse = pskdemod(xt_direct_mmse, modu_num, pi/modu_num);
        errorbit_direct_mmse(i) = errorbit_direct_mmse(i)+biterr(bits_part, hats_direct_mmse);           

       %% proposed ZF
        inv_HpHp = T*inv(H'*H)*T';
        S = inv(inv_HpHp(K-I+1:end,K-I+1:end));
        LSL = inv_HpHp(1:K-I, K-I+1:end) * S * inv_HpHp(K-I+1:end, 1:K-I);
        inv_HaHa = inv_HpHp(1:K-I, 1:K-I) - LSL;
        xt_proposed_zf = inv_HaHa * Ha' * y0;
        hats_proposed_zf = pskdemod(xt_proposed_zf, modu_num, pi/modu_num);
        errorbit_proposed_zf(i) = errorbit_proposed_zf(i)+biterr(bits_part, hats_proposed_zf);        
        
       %% proposed MMSE
        Ik = eye(K);
        w_mmse = (H'*H+1/snr*Ik)\H';
        x0 = w_mmse*y0;
        x0_mmse = T* x0 ;
        inv_HH = T*inv(H'*H+1/snr*Ik)*T';
        V = inv(inv_HH(K-I+1:end, K-I+1:end));
        UVU = inv_HH(1:K-I, K-I+1:end)*V*inv_HH(K-I+1:end, 1:K-I);
        phi = UVU*Ha'+ inv_HH(1:K-I, K-I+1:end)*Hi';
        xt_proposed_mmse = x0_mmse(1:K-I,:) -phi*y0;
        hats_proposed_mmse = pskdemod(xt_proposed_mmse, modu_num, pi/modu_num);
        errorbit_proposed_mmse(i) = errorbit_proposed_mmse(i)+biterr(bits_part, hats_proposed_mmse);    
    end
end

ber_in_zf  = errorbit_in_zf/(n_iterations*nn*(K-I)*L);
ber_in_mmse  = errorbit_in_mmse/(n_iterations*nn*(K-I)*L);
ber_direct_zf  = errorbit_direct_zf/(n_iterations*nn*(K-I)*L);
ber_direct_mmse  = errorbit_direct_mmse/(n_iterations*nn*(K-I)*L);
ber_proposed_zf  = errorbit_proposed_zf/(n_iterations*nn*(K-I)*L);
ber_proposed_mmse  = errorbit_proposed_mmse/(n_iterations*nn*(K-I)*L);

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
I = 4; % 不活跃的用户

%% 统计BER
errorbit_in_zf = zeros(1,length(SNR_dB));
errorbit_in_mmse = zeros(1,length(SNR_dB));

errorbit_direct_zf = zeros(1,length(SNR_dB));
errorbit_direct_mmse = zeros(1,length(SNR_dB));

errorbit_proposed_zf = zeros(1,length(SNR_dB));
errorbit_proposed_mmse = zeros(1,length(SNR_dB));

%% Inapproate ZF
for i = 1 : length(SNR_dB)
    snr = 10^(SNR_dB(i)/10);
    noise_var = 1/snr;
    for j = 1 : n_iterations
       %% 随机活跃用户
        inactive_user = randperm(K,I);
        bits = randi([0 modu_num-1],K, L);
        x = pskmod(bits, modu_num, pi/modu_num);
        x(inactive_user,:) = 0;
        T = permutation_T(K, I, inactive_user);
        
        bits_T = T*bits;         
        bits_part = bits_T;
        bits_part = bits_part(1:K-I,:);
        xp = T*x;
        xa = xp(1:K-I, :);
        
       %% 经过信道，高斯噪声
        n = sqrt(noise_var/2)*(randn(M,L)+1i*randn(M,L)); 
        H = sqrt(0.5)*(randn(M,K)+1i*randn(M,K));
        H_power = abs(H).^2;
        H_power_mean = mean(H_power,'all');
        Hp = H * T';
        Ha = Hp(:, 1:K-I);
        Hi = Hp(:, K-I+1:end);
        y = H*x+n;
        y0 = Ha*xa+n;
        
       %% inapproate ZF
        w_zf = (H'*H)\H';
        xt_in_zf = w_zf*y;
        hats_in_zf = pskdemod(xt_in_zf, modu_num, pi/modu_num);
        hats_in_zf(inactive_user,:) = bits(inactive_user,:);
        errorbit_in_zf(i) = errorbit_in_zf(i)+biterr(bits, hats_in_zf);

       %% inapproate MMSE
        Ik = eye(size(H,2));
        w_mmse = (H'*H+1/snr*Ik)\H';
        xt_in_mmse = w_mmse*y;
        hats_in_mmse = pskdemod(xt_in_mmse, modu_num, pi/modu_num);
        hats_in_mmse(inactive_user,:) = bits(inactive_user,:);
        errorbit_in_mmse(i) = errorbit_in_mmse(i)+biterr(bits, hats_in_mmse);        
        
       %% direct ZF
        w_zf = (Ha'*Ha)\Ha';
        xt_direct_zf = w_zf*y0;
        hats_direct_zf = pskdemod(xt_direct_zf, modu_num, pi/modu_num);
        errorbit_direct_zf(i) = errorbit_direct_zf(i)+biterr(bits_part ,hats_direct_zf);        
        
       %% direct MMSE
        Ik = eye(K-I);
        w_mmse = (Ha'*Ha+1/snr*Ik)\Ha';
        xt_direct_mmse = w_mmse*y0;
        hats_direct_mmse = pskdemod(xt_direct_mmse, modu_num, pi/modu_num);
        errorbit_direct_mmse(i) = errorbit_direct_mmse(i)+biterr(bits_part, hats_direct_mmse);           

       %% proposed ZF
        inv_HpHp = T*inv(H'*H)*T';
        S = inv(inv_HpHp(K-I+1:end,K-I+1:end));
        LSL = inv_HpHp(1:K-I, K-I+1:end) * S * inv_HpHp(K-I+1:end, 1:K-I);
        inv_HaHa = inv_HpHp(1:K-I, 1:K-I) - LSL;
        xt_proposed_zf = inv_HaHa * Ha' * y0;
        hats_proposed_zf = pskdemod(xt_proposed_zf, modu_num, pi/modu_num);
        errorbit_proposed_zf(i) = errorbit_proposed_zf(i)+biterr(bits_part, hats_proposed_zf);        
        
       %% proposed MMSE
        Ik = eye(K);
        w_mmse = (H'*H+1/snr*Ik)\H';
        x0 = w_mmse*y0;
        x0_mmse = T* x0 ;
        inv_HH = T*inv(H'*H+1/snr*Ik)*T';
        V = inv(inv_HH(K-I+1:end, K-I+1:end));
        UVU = inv_HH(1:K-I, K-I+1:end)*V*inv_HH(K-I+1:end, 1:K-I);
        phi = UVU*Ha'+ inv_HH(1:K-I, K-I+1:end)*Hi';
        xt_proposed_mmse = x0_mmse(1:K-I,:) -phi*y0;
        hats_proposed_mmse = pskdemod(xt_proposed_mmse, modu_num, pi/modu_num);
        errorbit_proposed_mmse(i) = errorbit_proposed_mmse(i)+biterr(bits_part, hats_proposed_mmse);    
    end
end

ber_in_zf  = errorbit_in_zf/(n_iterations*nn*(K-I)*L);
ber_in_mmse  = errorbit_in_mmse/(n_iterations*nn*(K-I)*L);
ber_direct_zf  = errorbit_direct_zf/(n_iterations*nn*(K-I)*L);
ber_direct_mmse  = errorbit_direct_mmse/(n_iterations*nn*(K-I)*L);
ber_proposed_zf  = errorbit_proposed_zf/(n_iterations*nn*(K-I)*L);
ber_proposed_mmse  = errorbit_proposed_mmse/(n_iterations*nn*(K-I)*L);

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
xlabel('SNR $\overline{\gamma}$ [dB]','Interpreter','latex');
ylabel('BER');
xlim([-30,0]);
legend('Inapp. ZF,\it{I=2}','Inapp. MMSE,\it{I=2}','Proposed ZF,\it{I=2}','Direct ZF,\it{I=2}',...
    'Proposed MMSE,\it{I=2}','Direct MMSE,\it{I=2}','Inapp. ZF,\it{I=4}','Inapp. MMSE,\it{I=4}',...
    'Proposed ZF,\it{I=4}','Direct ZF,\it{I=4}','Proposed MMSE,\it{I=4}','Direct MMSE,\it{I=4}');