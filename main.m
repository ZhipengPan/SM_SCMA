%作者：潘志鹏
%最后修订日期：2019年7月11日
%说明：本程序仿真SM-SCMA的BER性能

clear;
NT = 4;%发送天线数
M1 = NT;
b1 = log2(M1);%天线索引映射比特数
load('codebook_6users_4chips_qpsk_NR2.mat');%接收天线数的改变等效为码书维度的增加
[K,M2,V] = size(CB);
NR = K/4;
b2 = log2(M2);%SCMA码字映射信息比特数
M = M1*M2;

iter = 6;
bits_per_sym = (b1+b2)*V/K;%每个符号所占比特数
SNR = 0:2:40;
Es = 1;%叠加后SCMA码字的每个符号的平均能量为1
N0 = Es./(10.^(SNR/10));

F = zeros(K, V);
for uu = 1:V
    IND = CB(:,1,uu);
    F(IND~=0,uu) = 1;
end

mento = 1e7;%10000帧300秒，对应240000比特300秒

BER = zeros(length(SNR),1);
for snr_index = 1:length(SNR)

    ESC = 0;
    error = 0;
    tic
    for ii = 1:mento
        x = randi([0 M-1], V, 1); % log2(M)-bit symbols
        bit_x      = de2bi(x, log2(M), 'left-msb');
        
        H = sqrt(0.5)*(randn(K,V,NT)+1i*randn(K,V,NT));

        y = smscma_modu(bit_x,b1,b2,V,K,CB,H);
        noise = sqrt(0.5)*(randn(K, 1) + 1i*randn(K, 1))* sqrt(N0(snr_index));         

        ESC = ESC + mean(y.*conj(y));
        rec_y = y + noise;
        
        
        prior_LLR  = zeros(V,log2(M));

        [LLR] = SMSCMA_MPA_mex(rec_y,prior_LLR,NT,NR,CB,F,H,N0(snr_index),iter);

        bit_x_est = LLR<0 ;
        error = error + sum(sum(abs(bit_x-bit_x_est)));
        if((error>500 && ii*log2(M)*V>5e4) || ii*log2(M)*V > 1e7)
            break;
        end
        
    end
    
    disp(['Es = ',num2str(ESC/ii)]);
    BER(snr_index) = error/(ii*log2(M)*V);
    fprintf('SNR = %.1f dB, BER = %.6f\n',SNR(snr_index),BER(snr_index));
    toc
    if (BER(snr_index) < 1e-4)
        break;
    end

end
fig = figure;
semilogy(SNR(1:snr_index),BER(1:snr_index),'b-d','LineWidth',1.5,'MarkerSize',5);
legend(['N_T = ',num2str(NT),', N_R = ', num2str(NR)]);
xlabel('SNR(dB)');
ylabel('BER');
grid on;







