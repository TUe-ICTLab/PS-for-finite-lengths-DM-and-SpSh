clear
%===================
% This script generates Figure 5 of [3]. This Figure shows the
% gap-to-capacity, i.e., eq. (13) in [3], of M=8-ASK at rate R=2.25 bit/1D
% for a family of distributions with entropies ranging from 2.25 (only
% shaping) to 3=log2(M) (only coding). These distributions characterized by
% a MB-distribution for the amplitudes and the uniform distribution for
% the signs. As the achievable rate, R_BMD which is given by eq. (63) in 
% [1] is used. R_BMD is computed using the idea rexplained in Sec. 8.3 
% of [2].
% 
% Yunus Can Gültekin, Feb. 2021
%-------------------
% [1] Böcherer, G.; Steiner, F.; Schulte, P. Bandwidth efficient and 
%     rate-matched LDPC coded modulation. IEEE Trans. Commun. 2015, 63, 4651–4665.
% [2] Cover, T.M.; Thomas, J.A. Elements of Information Theory; John Wiley & Sons: Hoboken, NJ, USA, 2006.
% [3] Gültekin, Y.C.; Fehenberger, T.; Alvarado, A.; Willems, F.M.J. Probabilistic 
%     shaping for finite blocklengths: Distribution matching and sphere shaping. Entrop. 2020, vol. 19, no. 5: 581.
%===================
R = 2.25;                % target information rate
min_SNR_lin = 2^(2*R)-1; % capacity
M = 8;                   % 8-ASK
Aset = 1:2:M-1;          % Amplitude set of 8-ASK
snr_dB = 0:20;             
snr_lin = 10.^(snr_dB/10);

% Binary reflected Gray code for 8-ASK
labeling = [0 0 0; 0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 1 1; 1 0 1; 1 0 0]; 
        
% Family of Maxwell-Boltzmann distributions
load('lambda_set.mat')
PA = zeros(length(lambda_set),M/2); % set of MB-distributed P(A)'s
for il = 1:length(lambda_set)
    K = (sum(exp(-lambda_set(il)*Aset.^2)))^-1;
    PA(il,:) = K*exp(-lambda_set(il)*Aset.^2);
end
PX = 0.5*[fliplr(PA) PA];  % set of corresponding P(X)'s assuming uniform signs
HX = sum(-PX.*log2(PX),2); % set of corresponding constellation entropies 

% Compute R_BMD for all distributions with entropies in (R, log2(M))
Rbmd = zeros(length(lambda_set),length(snr_lin));
Best_MB_snr = zeros(length(lambda_set),1);
for iMB = 1:length(lambda_set)
    Rbmd = BMDrate_MASK_AWGN (M,PX(iMB,:),labeling,snr_lin);
    [AIR, indices] = unique(Rbmd);
    Best_MB_snr(iMB) = interp1(AIR, snr_lin(indices), R);
end

% Figure 5 of [3]
figure(5);
plot(HX,10*log10(Best_MB_snr/min_SNR_lin),'blue')
xlim([R log2(M)]); ylim([0 4]); grid on; grid minor;
xlabel('H(X) [bits]'); ylabel('\Delta SNR [dB]');
%------------- END