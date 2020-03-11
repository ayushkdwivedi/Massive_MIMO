% Single User Tx Beamforming

close all;
clear all;
rng('shuffle');

SNRdB = [1:1:30];
blockLength = 1000;
numBlocks = 10000;
BER = zeros(size(SNRdB));
SNR = zeros(size(SNRdB));
numAnt = 4;

for L = 1:numBlocks
    %ChNoise = zeros(1,blockLength);
    %Check for zero BER at zero noise
    ChNoise = (rand(1,blockLength) + j*randn(1,blockLength));
    bits = randi([0,1],[1,blockLength]);
    h = 1/sqrt(2)*(randn(1,numAnt) + j*randn(1,numAnt));
    
    for K = 1:length(SNRdB)
        SNR(K) = 10^(SNRdB(K)/10);
        Txbits = sqrt(SNR(K))*(2*bits-1); %Fading channel model
        tt = (h'/norm(h))*Txbits;
        Rxbits = h*tt + ChNoise; %Beamforming for single user
        
        ProcessedBits = Rxbits;
        DecodedBits = ((real(ProcessedBits))>0);
        BER(K) = BER(K) + sum(DecodedBits~=bits);
    end
end

BER = BER/(numBlocks*blockLength);
semilogy(SNRdB,BER,'b-s','linewidth',2.0);
%axis tight;
grid on;
legend('BER - Single User');
xlabel('SNR(dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR(dB) for single user bemforming')
        
    
  