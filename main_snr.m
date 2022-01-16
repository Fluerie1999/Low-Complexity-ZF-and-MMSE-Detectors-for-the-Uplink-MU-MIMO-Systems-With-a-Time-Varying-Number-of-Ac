clear;clc;

%% 传输参数
M = 30; % BS antennas
K = 10; % total users
n_iterations = 1e5; % 迭代次数
SNR_dB = -30:2:10; % SNR
nn = 2;
modu_num = 2^nn; % 调制阶数
I = 2; % 不活跃的用户
L = 10; % 每次每个用户发送的数据长度

%% 统计BER
snr_in_mc_zf = zeros(1,length(SNR_dB));
snr_in_mc_mmse = zeros(1,length(SNR_dB));
snr_in_ana_zf = zeros(1,length(SNR_dB));
snr_in_ana_mmse = zeros(1,length(SNR_dB));
snr_proposed_mc_zf = zeros(1,length(SNR_dB));
snr_proposed_mc_mmse = zeros(1,length(SNR_dB));
snr_proposed_ana_zf = zeros(1,length(SNR_dB));
snr_proposed_ana_mmse = zeros(1,length(SNR_dB));

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
    
       %% inapproate MMSE
        Ik = eye(size(H,2));
        w_mmse = (H'*H+(1/snr)*Ik)\H';
        xt_in_mmse = w_mmse*y;
        
       %% Propsoed ZF
        inv_HpHp = T*inv(H'*H)*T';
        S = inv(inv_HpHp(K-I+1:end,K-I+1:end));
        LSL = inv_HpHp(1:K-I, K-I+1:end) * S * inv_HpHp(K-I+1:end, 1:K-I);
        xt_proposed_zf = (Ha'*Ha)\Ha'* y0;
        L21 = Hi'*Ha/(Ha'*Ha);
        Z = LSL*Ha'+inv_HpHp(1:K-I, K-I+1:end)*Hi';
        pho_proposed_mc_zf = Z*Z'+inv(Ha'*Ha)+snr*Z*(Ha*Ha')*Z';
        pho_proposed_mc_zf_ave = sum(diag(pho_proposed_mc_zf))/(K-I);
        
       %% proposed MMSE
        Ik = eye(K);
        w_mmse = (H'*H+(1/snr)*Ik)\H';
        x0 = w_mmse*y0;
        x0_mmse = T* x0 ;
        inv_HH = T*inv(H'*H+1/snr*Ik)*T';
        V = inv(inv_HH(K-I+1:end, K-I+1:end));
        UVU = inv_HH(1:K-I, K-I+1:end)*V*inv_HH(K-I+1:end, 1:K-I);
        phi = UVU*Ha'+ inv_HH(1:K-I, K-I+1:end)*Hi';
        xt_proposed_mmse = x0_mmse(1:K-I,:) -phi*y0;
        pho_proposed_mc_mmse = inv(Ha'*Ha+eye(K-I)/snr)+phi*(snr*(Ha*Ha')+eye(M))*phi';
        pho_proposed_mc_mmse_ave = mean(diag(pho_proposed_mc_mmse));
        
        %% 计算仿真SNR
         xt_in_mc_zf = sum(mean(abs(xt_in_zf).^2 , 2))/(K-I);
         pnoise_in_mc_zf = mean(mean(abs(w_zf*n).^2, 2));
         
         xt_in_mc_mmse = sum(mean(abs(xt_in_mmse).^2 , 2))/(K-I);
         xt_in_mmse(inactive_user,:) = 0;
         pnoise_in_mc_mmse = sum(mean(abs(x-xt_in_mmse).^2, 2))/(K-I);
         
         xt_proposed_mc_zf = mean(mean(abs(xt_proposed_zf).^2 , 2));
         pnoise_proposed_mc_zf = mean(mean(abs(inv(Ha'*Ha) *Ha'*n).^2, 2));
         
         xt_proposed_mc_mmse = mean(mean(abs(xt_proposed_mmse).^2 , 2));
         noise_proposed_mc_mmce =  T*w_mmse*n;
         noise = noise_proposed_mc_mmce(1:K-I,:)-phi*n;
         pnoise_proposed_mc_mmse = mean(mean(abs(xa-xt_proposed_mmse).^2, 2));
         
       %% 计算推导SNR
        IkI = eye(K-I);
        inv_HaHa = inv(Ha'*Ha);
        inv_HaHa_ave = mean(diag(inv_HaHa ));
        inv_HaHa_mmse = inv(Ha'*Ha+1/snr*IkI);
        inv_HaHa_mmse_ave = mean(diag(inv_HaHa_mmse));
        
        %% 实际值
        snr_in_mc_zf(i) = snr_in_mc_zf(i) +  1 / pnoise_in_mc_zf;
        snr_in_mc_mmse(i) = snr_in_mc_mmse(i) + 1 / pnoise_in_mc_mmse;
        snr_proposed_mc_zf(i) = snr_proposed_mc_zf(i) + 1 / pnoise_proposed_mc_zf;
        snr_proposed_mc_mmse(i) = snr_proposed_mc_mmse(i) + 1 / pnoise_proposed_mc_mmse;
        
        %% 理论值
        snr_in_ana_zf(i)  = snr_in_ana_zf(i) + snr/real(pho_proposed_mc_zf_ave);
        snr_in_ana_mmse(i)  = snr_in_ana_mmse(i) + snr/real(pho_proposed_mc_mmse_ave);
        snr_proposed_ana_zf(i)  = snr_proposed_ana_zf(i) + snr/inv_HaHa_ave ;
        snr_proposed_ana_mmse(i)  = snr_proposed_ana_mmse(i) + snr/inv_HaHa_mmse_ave ;
    end
end

snr_in_mc_zf = snr_in_mc_zf /n_iterations;
snr_in_mc_zf = 10*log10(snr_in_mc_zf);

snr_in_mc_mmse = snr_in_mc_mmse /n_iterations;
snr_in_mc_mmse = 10*log10(snr_in_mc_mmse);

snr_proposed_mc_zf = snr_proposed_mc_zf /n_iterations;
snr_proposed_mc_zf = 10*log10(snr_proposed_mc_zf);

snr_proposed_mc_mmse = snr_proposed_mc_mmse /n_iterations;
snr_proposed_mc_mmse = 10*log10(snr_proposed_mc_mmse);

snr_in_ana_zf = snr_in_ana_zf/n_iterations;
snr_in_ana_zf = 10*log10(snr_in_ana_zf);

snr_in_ana_mmse = snr_in_ana_mmse/n_iterations;
snr_in_ana_mmse = 10*log10(snr_in_ana_mmse );

snr_proposed_ana_zf = snr_proposed_ana_zf/n_iterations;
snr_proposed_ana_zf = 10*log10(snr_proposed_ana_zf );

snr_proposed_ana_mmse = snr_proposed_ana_mmse/n_iterations;
snr_proposed_ana_mmse = 10*log10(snr_proposed_ana_mmse );

plot(SNR_dB ,snr_proposed_mc_mmse ,'b|-', 'Linewidth',1);
hold on;grid on;
plot(SNR_dB ,snr_proposed_ana_mmse ,'b>', 'Linewidth',1);
hold on;grid on;
plot(SNR_dB ,snr_proposed_mc_zf , 'kx-', 'Linewidth',1);
hold on;grid on;
plot(SNR_dB ,snr_proposed_ana_zf , 'ko', 'Linewidth',1);
hold on;grid on;
plot(SNR_dB ,snr_in_mc_mmse ,'b^-', 'Linewidth',1);
hold on;grid on;
plot(SNR_dB ,snr_in_ana_mmse ,'bd', 'Linewidth',1);
hold on;grid on;
plot(SNR_dB ,snr_in_mc_zf  , 'k|-', 'Linewidth',1);
hold on;grid on;
plot(SNR_dB ,snr_in_ana_zf  , 'ks', 'Linewidth',1);
set(gcf, 'position', [200 200 400 300]);
xlabel('SNR $\overline{\gamma}$ [dB]', 'interpreter', 'latex'); ylabel('Post processing SNR');
legend('MMSE(MC)','MMSE(Analytical)','ZF(MC)','ZF(Analytical)',...
    'Inappropriate MMSE(MC)','Inappropriate MMSE(Analytical)','Inappropriate ZF(MC)','Inappropriate ZF(Analytical)');