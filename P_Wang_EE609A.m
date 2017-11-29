% Simulate and compare the BER of a 256-PSK system and a 256-QAM system 
% with Grey coding and Eb/N0 = 0, 2, 4, 6, 8, 10 dB. (P)  
close all;
clear all;
clc;

% %{
num = 1e5; % the number of symbols
M = 256; % for MPSK and MQAM
L = sqrt(M);
N_bpsym = log2(M); % for MPSK, kbits/symbol

s_w_mpsk = randi([0 M-1],1,num); % generate the omega of the signal symbols, 0~M-1
s_s_mpsk = exp(1j*s_w_mpsk*2*pi/M); % generate the signal of symbols

ampli_mqam = -L+1:2:L-1; % MQAM amplitude
s_o_re_mqam = randi([0 L-1],1,num); % the order of symbols, e.g.0,1,2,3-->-3,-1,1,3 for 16QAM
s_o_im_mqam = randi([0 L-1],1,num); % the order of symbols, e.g.0,1,2,3-->-3,-1,1,3 for 16QAM
s_s_mqam = ampli_mqam(s_o_re_mqam+1) + 1j*ampli_mqam(s_o_im_mqam+1);

ref = 0:M-1;
map = bitxor(ref,floor(ref/2)); % map for Gray code

s_wg_mpsk = map(s_w_mpsk+1); % after gray code
s_bg_mpsk = dec2bin(s_wg_mpsk)-'0';
s_b_mpsk =  reshape(s_bg_mpsk.',1,N_bpsym*num);

s_o_re_g_mqam = map(s_o_re_mqam+1); % after gray code
s_o_im_g_mqam = map(s_o_im_mqam+1); % after gray code
s_b_re_g_mqam = dec2bin(s_o_re_g_mqam)-'0';
s_b_im_g_mqam = dec2bin(s_o_im_g_mqam)-'0';
s_b_g_mqam =  [reshape(s_b_re_g_mqam.',1,log2(L)*num),reshape(s_b_im_g_mqam.',1,log2(L)*num)];

EbN0_db = 0:2:10;
BER_mpsk = zeros(1,length(EbN0_db));
SER_mpsk = zeros(1,length(EbN0_db));

BER_mqam = zeros(1,length(EbN0_db));
SER_mqam = zeros(1,length(EbN0_db));
s_s1_re_mqam = zeros(1,length(s_s_mqam));
s_s1_im_mqam = zeros(1,length(s_s_mqam));

Eb_mpsk = 1/N_bpsym; % Eb is the average energy per bit, Es=abs(s_s)=1=log2(M)*Eb,for PSK,Es=1
N0_sd_mpsk = Eb_mpsk./10.^(EbN0_db/10); % SNR is Eb/N0, N0_sd is noise power spectral density

Es_mqam = (2*(M-1)/3); % Es=abs(s_s)_average=(2*(M-1)/3),for QAM
Eb_mqam = Es_mqam/N_bpsym; % Eb is the average energy per bit
N0_sd_mqam = Eb_mqam./10.^(EbN0_db/10); % SNR is Eb/N0, N0_sd is noise power spectral density

for m = 1:1:length(EbN0_db)
    N0_mpsk = sqrt(N0_sd_mpsk(m)/2)*(randn(1,length(s_s_mpsk))+1j*randn(1,length(s_s_mpsk))); % N0_sd/2 is two-sided power spectral density of the noise
    s_sn_mpsk = s_s_mpsk+N0_mpsk; % symbol signal with noise
    
    s_wn_mpsk = round((angle(s_sn_mpsk))/(2*pi/M)); % the angle of symbol signal with noise,-pi~pi/(2*pi/M))
    s_w1_mpsk = mod(s_wn_mpsk,M);
    SER_mpsk(m) = length(find(s_w1_mpsk-s_w_mpsk))./length(s_w1_mpsk); % SER of symbol Signal with AWGN
    
    s_w1_g_mpsk = map(s_w1_mpsk+1); % decode
    s_b1_g_mpsk = dec2bin(s_w1_g_mpsk)-'0'; % decimal number to binary number
    s_b1_mpsk =  reshape(s_b1_g_mpsk.',1,N_bpsym*num);
    BER_mpsk(m) = length(find(s_b1_mpsk-s_b_mpsk))./length(s_b1_mpsk); % BER of bit Signal with AWGN
    
    N0_mqam = sqrt(N0_sd_mqam(m)/2)*(randn(1,length(s_s_mqam))+1j*randn(1,length(s_s_mqam))); % N0_sd/2 is two-sided power spectral density of the noise
    s_sn_mqam = s_s_mqam+N0_mqam; % symbol signal with noise
    
    s_s1_re_mqam = 2*floor(real(s_sn_mqam)/2)+1;
    s_s1_re_mqam(s_s1_re_mqam>max(ampli_mqam)) = max(ampli_mqam);
    s_s1_re_mqam(s_s1_re_mqam<min(ampli_mqam)) = min(ampli_mqam);
    s_s1_im_mqam = 2*floor(imag(s_sn_mqam)/2)+1;
    s_s1_im_mqam(s_s1_im_mqam>max(ampli_mqam)) = max(ampli_mqam);
    s_s1_im_mqam(s_s1_im_mqam<min(ampli_mqam)) = min(ampli_mqam);
    s_s1_mqam = s_s1_re_mqam + 1j*s_s1_im_mqam;
    SER_mqam(m) = length(find(s_s1_mqam-s_s_mqam))./length(s_s1_mqam); % SER of symbol Signal with AWGN
    
    s_o1_re_g_mqam = map(floor((s_s1_re_mqam+L)/2)+1); 
    s_o1_im_g_mqam = map(floor((s_s1_im_mqam+L)/2)+1); 
    s_b1_re_g_mqam = dec2bin(s_o1_re_g_mqam)-'0';
    s_b1_im_g_mqam = dec2bin(s_o1_im_g_mqam)-'0';
    s_b1_g_mqam =  [reshape(s_b1_re_g_mqam.',1,log2(L)*num),reshape(s_b1_im_g_mqam.',1,log2(L)*num)];
    BER_mqam(m) = length(find(s_b1_g_mqam-s_b_g_mqam))./length(s_b1_g_mqam); % BER of bit Signal with AWGN
end

EbN0_db = linspace(0,10,100); % the number of the theory simulate points in the picture are 100
EbN0 = 10.^(EbN0_db/10);
EsN0 = N_bpsym*EbN0;
% SER_t_mpsk = 2*qfunc(sqrt(2*EsN0)*sin(pi/M));
% BER_t_mpsk = SER_t_mpsk/N_bpsym; % coherence detected MPSK
[BER_t_mpsk,SER_t_mpsk] = berawgn(EbN0_db,'psk',M,'nondiff');

% SER_t_mqam = 4*(1-1/L)*qfunc(sqrt(3*EsN0/(M-1))).*(1-(1-1/L)*qfunc(sqrt(3*EsN0/(M-1))));
% % SER_t_mqam = 4*(1-1/L)*qfunc(sqrt(2*EsN0/Es_mqam));
% % BER_t_mqam = SER_t_mqam/N_bpsym; % 
% BER_t_mqam = 2*(1-1/L)/log2(L)*qfunc(sqrt(3*log2(L)/(L^2-1)*2*EbN0));
% % SER_t_mqam = BER_t_mqam*N_bpsym; % 
[BER_t_mqam,SER_t_mqam] = berawgn(EbN0_db,'qam',M);

figure;
semilogy(EbN0_db, SER_t_mpsk,'b');
hold on;
semilogy(EbN0_db, BER_t_mpsk,'r--');
hold on;
semilogy(EbN0_db, SER_t_mqam,'g');
hold on;
semilogy(EbN0_db, BER_t_mqam,'c--');
for m = 1:1:m
    hold on;
    plot(2*(m-1),SER_mpsk(m),'b*');
    hold on;
    plot(2*(m-1),BER_mpsk(m),'ro');
    hold on;
    plot(2*(m-1),SER_mqam(m),'g*');
    hold on;
    plot(2*(m-1),BER_mqam(m),'co'); 
end
grid on;
ylabel('P');
xlabel('E_b/N_0 (dB)');
legend(['theoretical SER in ',num2str(M),'PSK'],...
    ['theoretical BER in ',num2str(M),'PSK'],...
    ['theoretical SER in ',num2str(M),'QAM'],...
    ['theoretical BER in ',num2str(M),'QAM'],...
    ['simulated SER in ',num2str(M),'PSK'],...
    ['simulated BER in ',num2str(M),'PSK'],...
    ['simulated SER in ',num2str(M),'QAM'],...
    ['simulated BER in ',num2str(M),'QAM'],...
    'Location','southwest');
% legend(['theoretical SER in ',num2str(M),'PSK'],...
%     ['theoretical BER in ',num2str(M),'PSK'],...
%     ['simulated SER in ',num2str(M),'PSK'],...
%     ['simulated BER in ',num2str(M),'PSK'],...
%     'Location','southwest');
% legend(['theoretical SER in ',num2str(M),'QAM'],...
%     ['theoretical BER in ',num2str(M),'QAM'],...
%     ['simulated SER in ',num2str(M),'QAM'],...
%     ['simulated BER in ',num2str(M),'QAM'],...
%     'Location','southwest');
% hold off;
return
%}











