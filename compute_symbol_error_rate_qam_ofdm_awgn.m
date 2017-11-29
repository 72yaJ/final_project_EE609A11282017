%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All rights reserved by Krishna Pillai, http://www.dsplog.com 
% The file may not be re-distributed without explicit authorization
% from Krishna Pillai.
% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna Pillai
% Email		: krishna@dsplog.com
% Version	: 1.0
% Date		: 01 January 2012
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computing the Symbol Error probability for a general M-QAM 
% using OFDM modulation

% function compute_symbol_error_rate_qam_ofdm_awgn()

close all; 
clear all;
clc;

figure
EsN0dB 	= 0:33; % symbol to noise ratio
M 	= [16 64 256]; % 16QAM/64QAM and 256 QAM  
color_vec1 = ['b-','m-','g-'];    
color_vec2 = ['ks-','rx-','cd-'];    
for jj= 1:length(M)
	k      = sqrt(1/((2/3)*(M(jj)-1))); 
	simSer(jj,:) = compute_symbol_error_rate(EsN0dB, M(jj));
	theorySer(jj,:) = 2*(1-1/sqrt(M(jj)))*erfc(k*sqrt((10.^(EsN0dB/10)))) ...
	              - (1-2/sqrt(M(jj)) + 1/M(jj))*(erfc(k*sqrt((10.^(EsN0dB/10))))).^2;
	semilogy(EsN0dB,theorySer(jj,:),color_vec1(jj));
	hold on;
	semilogy(EsN0dB,simSer(jj,:),   color_vec2(jj));
end
axis([0 33 10^-5 1])
grid on
legend('theory-16QAM', 'sim-16QAM', 'theory-64QAM', 'sim-64QAM', 'theory-256QAM', 'sim-256QAM');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for 16QAM/64QAM/256QAM using OFDM');
return ;

function [simSer] = compute_symbol_error_rate(EsN0dB, M)

% ofdm specifications
nFFT = 64; % fft size
nDSC = 52; % number of data subcarriers
nConstperOFDMsym = 52; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)
nOFDMsym = 10^4; % number of ofdm symbols

% modulation 
k = sqrt(1/((2/3)*(M-1))); % normalizing factor
m = [1:sqrt(M)/2]; % alphabets
alphaMqam = [-(2*m-1) 2*m-1]; 

EsN0dB_eff = EsN0dB  + 10*log10(nDSC/nFFT) + 10*log10(64/80); % accounting for the used subcarriers and cyclic prefix

for ii = 1:length(EsN0dB)

   % Transmitter
   ipMod = randsrc(1,nConstperOFDMsym*nOFDMsym,alphaMqam) + j*randsrc(1,nConstperOFDMsym*nOFDMsym,alphaMqam);
   ipMod_norm = k*reshape(ipMod,nConstperOFDMsym,nOFDMsym).'; % grouping into multiple symbolsa

   % Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
   xF = [zeros(nOFDMsym,6) ipMod_norm(:,[1:nConstperOFDMsym/2]) zeros(nOFDMsym,1) ipMod_norm(:,[nConstperOFDMsym/2+1:nConstperOFDMsym]) zeros(nOFDMsym,5)] ;
   
   % Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1 
   xt = (nFFT/sqrt(nDSC))*ifft(fftshift(xF.')).';

   % Appending cylic prefix
   xt = [xt(:,[49:64]) xt];

   % Concatenating multiple symbols to form a long vector
   xt = reshape(xt.',1,nOFDMsym*80);

   % Gaussian noise of unit variance, 0 mean
   nt = 1/sqrt(2)*[randn(1,nOFDMsym*80) + j*randn(1,nOFDMsym*80)];

   % Adding noise, the term sqrt(80/64) is to account for the wasted energy due to cyclic prefix
   yt = sqrt(80/64)*xt + 10^(-EsN0dB_eff(ii)/20)*nt;

   % Receiver
   yt = reshape(yt.',80,nOFDMsym).'; % formatting the received vector into symbols
   yt = yt(:,[17:80]); % removing cyclic prefix

   % converting to frequency domain
   yF = (sqrt(nDSC)/nFFT)*fftshift(fft(yt.')).'; 
   yMod = sqrt(64/80)*yF(:,[6+[1:nConstperOFDMsym/2] 7+[nConstperOFDMsym/2+1:nConstperOFDMsym] ]); 

   % demodulation
   y_re = real(yMod)/k;
   y_im = imag(yMod)/k;
   % rounding to the nearest alphabet
   % 0 to 2 --> 1
   % 2 to 4 --> 3
   % 4 to 6 --> 5 etc
   ipHat_re = 2*floor(y_re/2)+1;
   ipHat_re(find(ipHat_re>max(alphaMqam))) = max(alphaMqam);
   ipHat_re(find(ipHat_re<min(alphaMqam))) = min(alphaMqam);
            
   % rounding to the nearest alphabet
   % 0 to 2 --> 1
   % 2 to 4 --> 3
   % 4 to 6 --> 5 etc
   ipHat_im = 2*floor(y_im/2)+1;
   ipHat_im(find(ipHat_im>max(alphaMqam))) = max(alphaMqam);
   ipHat_im(find(ipHat_im<min(alphaMqam))) = min(alphaMqam);
    
   ipHat = ipHat_re + j*ipHat_im; 

   % converting to vector 
   ipHat_v = reshape(ipHat.',nConstperOFDMsym*nOFDMsym,1).';

   % counting the errors
   nErr(ii) = size(find(ipMod - ipHat_v ),2);

end
simSer = nErr/(nOFDMsym*nConstperOFDMsym);
end
% return;
%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creative Commons
% Attribution-Noncommercial 2.5 India
% You are free:
% to Share â€? to copy, distribute and transmit the work
% to Remix â€? to adapt the work
% Under the following conditions:
% Attribution. You must attribute the work in the manner 
% specified by the author or licensor (but not in any way 
% that suggests that they endorse you or your use of the work). 
% Noncommercial. You may not use this work for commercial purposes. 
% For any reuse or distribution, you must make clear to others the 
% license terms of this work. The best way to do this is with a 
% link to this web page.
% Any of the above conditions can be waived if you get permission 
% from the copyright holder.
% Nothing in this license impairs or restricts the author's moral rights.
% http://creativecommons.org/licenses/by-nc/2.5/in/

% Script for simulating 16-QAM transmission and reception and compare the 
% simulated and theoretical symbol error probability

% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna
% Email		: krishna@dsplog.com
% Version	: 1.0
% Date		: 9 December 2007
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% symbol error rate for 16-QAM modulation
clear
N = 2*10^3; % number of symbols
alpha16qam = [-3 -1 1 3]; % 16-QAM alphabets
Es_N0_dB = [0:20]; % multiple Es/N0 values
ipHat = zeros(1,N);
for ii = 1:length(Es_N0_dB)
    ip = randsrc(1,N,alpha16qam) + j*randsrc(1,N,alpha16qam);
    s = (1/sqrt(10))*ip; % normalization of energy to 1
    n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance

    y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise

    % demodulation
    y_re = real(y); % real part
    y_im = imag(y); % imaginary part

    ipHat_re(find(y_re< -2/sqrt(10)))           = -3;
    ipHat_re(find(y_re > 2/sqrt(10)))           =  3;
    ipHat_re(find(y_re>-2/sqrt(10) & y_re<=0))  = -1;
    ipHat_re(find(y_re>0 & y_re<=2/sqrt(10)))   =  1;

    ipHat_im(find(y_im< -2/sqrt(10)))           = -3;
    ipHat_im(find(y_im > 2/sqrt(10)))           =  3;
    ipHat_im(find(y_im>-2/sqrt(10) & y_im<=0))  = -1;
    ipHat_im(find(y_im>0 & y_im<=2/sqrt(10)))   =  1;
    ipHat = ipHat_re + j*ipHat_im;
    nErr(ii) = size(find([ip- ipHat]),2); % couting the number of errors
end

simBer = nErr/N;
theoryBer = 3/2*erfc(sqrt(0.1*(10.^(Es_N0_dB/10))));
close all
figure
semilogy(Es_N0_dB,theoryBer,'b.-','LineWidth',2);
hold on
semilogy(Es_N0_dB,simBer,'mx-','Linewidth',2);
axis([0 20 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for 16-QAM modulation') 
%}

% %{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All rights reserved by Krishna Pillai, http://www.dsplog.com 
% The file may not be re-distributed without explicit authorization
% from Krishna Pillai.
% Checked for proper operation with Octave Version 3.0.0
% Author	: Krishna Pillai
% Email		: krishna@dsplog.com
% Version	: 1.0
% Date		: 05 June 2008
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bit Error Rate for 16-QAM modulation using Gray modulation mapping

clear
N = 10^3; % number of symbols
M = 16;   % constellation size
k = log2(M); % bits per symbol

% defining the real and imaginary PAM constellation
% for 16-QAM
alphaRe = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
alphaIm = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
k_16QAM = 1/sqrt(10);

Eb_N0_dB  = 0:15; % multiple Es/N0 values
Es_N0_dB  = Eb_N0_dB + 10*log10(k);

% Mapping for binary <--> Gray code conversion
ref = 0:k-1;
map = bitxor(ref,floor(ref/2));
[tt, ind] = sort(map);                                

for ii = 1:length(Eb_N0_dB)
    
    % symbol generation
    % ------------------
    ipBit = rand(1,N*k,1)>0.5; % random 1's and 0's
    ipBitReshape = reshape(ipBit,k,N).';
    bin2DecMatrix = ones(N,1)*(2.^((k/2-1):-1:0)) ; % conversion from binary to decimal
    % real
    ipBitRe =  ipBitReshape(:,1:k/2);
    ipDecRe = sum(ipBitRe.*bin2DecMatrix,2);
    ipGrayDecRe = bitxor(ipDecRe,floor(ipDecRe/2));
    % imaginary
    ipBitIm =  ipBitReshape(:,k/2+1:k);
    ipDecIm = sum(ipBitIm.*bin2DecMatrix,2);
    ipGrayDecIm = bitxor(ipDecIm,floor(ipDecIm/2)); 
    % mapping the Gray coded symbols into constellation
    modRe = alphaRe(ipGrayDecRe+1);
    modIm = alphaIm(ipGrayDecIm+1);
    % complex constellation
    mod = modRe + 1j*modIm;
    s = k_16QAM*mod; % normalization of transmit power to one 
    
    % noise
    % -----
    n = 1/sqrt(2)*(randn(1,N) + 1j*randn(1,N)); % white guassian noise, 0dB variance 
    
    y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise

    % demodulation
    % ------------
    y_re = real(y)/k_16QAM; % real part
    y_im = imag(y)/k_16QAM; % imaginary part

    % rounding to the nearest alphabet
    ipHatRe = 2*floor(y_re/2)+1;
    ipHatRe(ipHatRe>max(alphaRe)) = max(alphaRe);
    ipHatRe(ipHatRe<min(alphaRe)) = min(alphaRe);
    ipHatIm = 2*floor(y_im/2)+1;
    ipHatIm(ipHatIm>max(alphaIm)) = max(alphaIm);
    ipHatIm(ipHatIm<min(alphaIm)) = min(alphaIm);

    % Constellation to Decimal conversion
    ipDecHatRe = ind(floor((ipHatRe+4)/2+1))-1; % LUT based
    ipDecHatIm = ind(floor((ipHatIm+4)/2+1))-1; % LUT based

    % converting to binary string
    ipBinHatRe = dec2bin(ipDecHatRe,k/2);
    ipBinHatIm = dec2bin(ipDecHatIm,k/2);

    % converting binary string to number
    ipBinHatRe = ipBinHatRe.';
    ipBinHatRe = ipBinHatRe(1:end).';
    ipBinHatRe = reshape(str2num(ipBinHatRe).',k/2,N).' ;
    
    ipBinHatIm = ipBinHatIm.';
    ipBinHatIm = ipBinHatIm(1:end).';
    ipBinHatIm = reshape(str2num(ipBinHatIm).',k/2,N).' ;

    % counting errors for real and imaginary
    nBitErr(ii) = size(find(ipBitRe- ipBinHatRe),1) + size(find(ipBitIm - ipBinHatIm),1) ;

end 
simBer = nBitErr/(N*k);
theoryBer = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(Eb_N0_dB/10))));

close all; figure
semilogy(Eb_N0_dB,theoryBer,'bs-','LineWidth',2);
hold on
semilogy(Eb_N0_dB,simBer,'mx-','LineWidth',2);
axis([0 15 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for 16-QAM modulation')
%}























