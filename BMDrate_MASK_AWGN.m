function BMDrate = BMDrate_MASK_AWGN (M,PX,labeling,snr_lin)
%===================
% This function computes the achievable rate R_BMD given by eq. (63) in [1]
% as R_BMD = H(X) - \sum_i H(Bi|Y) where X is labeled by a binary
% string B1B2...Bm.
%-------------------
% - Channel model is expressed by Y = X + N where X takes values from an 
% M-ASK constellation, and N is white Gaussian noise.
% - The computation is done by dividing the range of Y (which is continuous
% in reality) into bins of length delta, and using the definition of
% Riemann integrability, see Theorem 8.3.1 in [2].
% - This computation technique is analogous to the computation of the
% discrete entropy of a quantized version A' of a cont. r.v. A from the
% differential entropy of A, see Sec. 8.3  in [1].
%-------------------
% [1] Böcherer, G.; Steiner, F.; Schulte, P. Bandwidth efficient and rate-matched LDPC coded modulation. IEEE Trans. Commun. 2015, 63, 4651–4665.
% [2] Cover, T.M.; Thomas, J.A. Elements of Information Theory; John Wiley & Sons: Hoboken, NJ, USA, 2006.
%===================
% M: Number of symbols in the ASK constellation
% PX: 1D symbol probabilities (row vector)
% labeling: binary labels of symbols -(M-1), -(M-3), ..., (M-1) in each row
% snr_lin: (linear) SNR values at which the BMD rate is to be computed
%===================
if abs(sum(PX)-1)<1e-10
    if isrow(PX)
    else
        PX = PX.';
    end
else
    error('Unfortunately, the 2nd Kolmogorov axiom has to be satisfied.');
end
%-------------------
nr_quant_steps = 15000; % # of steps used to quantize the support set of Y
                        % decreasing this value decreases accuracy, see
                        % the accuracy checks (***) below
%-------------------
X1D = -(M-1):2:(M-1);       % 1D constellation
m   = log2(M);              % # of address bits per symbol
Eav = sum(PX .* X1D.^2);    % 1D average energy
%-------------------
PB = zeros(m,2);            % bit-level distributions
for ibit = 1:m
    PB(ibit,0+1) = sum(PX(labeling(:,ibit)==0)) ;        % being zero
    PB(ibit,1+1) = sum(PX(labeling(:,ibit)==1)) ;        % being one
end
%-------------------
BMDrate = zeros(length(snr_lin),1);
%-------------------
for isnr = 1:length(snr_lin) 
    
    noise_var     = Eav/snr_lin(isnr);
    noise_std_dev = sqrt(noise_var);
    inv2Var       = 1 / (2*noise_var);
    invSq2piVar   = 1 / sqrt(2*pi*noise_var);
    
    % It is assumed that p(y) is negligible 5 standard deviation outside
    % the support set of X, i.e., outside (-X1D-5sigma, X1D+5sigma)
    delta = ((max(X1D)+5*noise_std_dev)-(min(X1D)-5*noise_std_dev))/nr_quant_steps;
    
    support_Y = min(X1D)-5*noise_std_dev : delta : max(X1D)+5*noise_std_dev;
    
    %----
    % Compute: R_BMD = H(X) - \sum_i H(Bi|Y) 
    %----
    %---- The minuend in R_BMD expression: Entropy H(X) = H(B1B2...Bm)
    PX_temp = PX;
    PX_temp(PX_temp==0) = [];  % Convention: 0*log(0)=0, see p.31 [1]. In MATLAB, it leads to NaN.
    H_X = sum(-PX_temp.*log2(PX_temp));
    
    %---- The subtrahend in R_BMD expression: Conditional entropies H(Bi|Y) 
    % H(Bi|Y) = H(Bi) - I(Bi;Y) = H(Bi) - h(Y) + h(Y|Bi) -> 3 terms

        %------- 1. H(Bi)'s
        H_Bi = zeros(m,1);
        for ibit = 1:m % bit-levels 1,2,...,m
            H_Bi(ibit,1) = sum(-PB(ibit,:).*log2(PB(ibit,:)));
        end

        %------- 2. h(Y)
        p_y = zeros(1,length(support_Y));
        iy = 1;
        for y = support_Y
            for ix = 1:M
                p_y(iy) = p_y(iy) + PX(ix) * invSq2piVar * exp( -(y-X1D(ix))^2 * inv2Var);
            end
            iy = iy + 1;
        end
        if abs(sum(p_y)*delta - 1) > 1e-4 % the accuracy check (***)
            error('Resolution is not enough.')
        else
        end
        p_y(p_y==0)=[]; % Convention: 0*log(0)=0, see p.31 [1].
        H_Ydelta  = - sum(delta*p_y.*log2(p_y)) - log2(delta); % eq. (8.29) in [1]
        h_Y       = H_Ydelta + log2(delta);                    % eq. (8.30) in [1], assuming delta is small enough

        %------- 3. h(Y|Bi)
        h_YgBi = zeros(m,1);
        for ibit = 1:m % bit-levels 1,2,...,m
            %---- p(y|Bi=0) & p(y|Bi=1)
            p_ygb0 = zeros(1,length(support_Y));
            p_ygb1 = zeros(1,length(support_Y));
            %---- Set of X=x for which Bi=0 & Bi=1
            X1D_b0 = X1D(labeling(:,ibit)==0);
            X1D_b1 = X1D(labeling(:,ibit)==1);
            %---- Prob. of X=x for which Bi=0 & Bi=1
            PX0 = PX(labeling(:,ibit)==0)/PB(ibit,0+1);
            PX1 = PX(labeling(:,ibit)==1)/PB(ibit,1+1);         
            %----
            iy = 1;
            for y = support_Y
                for ix = 1:M/2
                    p_ygb0(iy) = p_ygb0(iy) + PX0(ix) * invSq2piVar * exp( -(y-X1D_b0(ix))^2 * inv2Var);
                    p_ygb1(iy) = p_ygb1(iy) + PX1(ix) * invSq2piVar * exp( -(y-X1D_b1(ix))^2 * inv2Var);
                end
                iy = iy + 1;
            end
            if (abs(sum(p_ygb0)*delta - 1) > 1e-2) ||...
                    (abs(sum(p_ygb0)*delta - 1) > 1e-2) % the accuracy check (***)
                error('Resolution is not enough.')
            else
            end
            p_ygb0(p_ygb0==0)=[]; % Convention: 0*log(0)=0, see p.31 [1]. 
            p_ygb1(p_ygb1==0)=[]; % Convention: 0*log(0)=0, see p.31 [1]. 

            H_Yb0delta  = - sum(delta*p_ygb0.*log2(p_ygb0)) - log2(delta);  % eq. (8.29) in [1]
            H_Yb1delta  = - sum(delta*p_ygb1.*log2(p_ygb1)) - log2(delta);  % eq. (8.29) in [1]
            h_Ygbi0      = H_Yb0delta + log2(delta); % eq. (8.30) in [1], assuming delta is small enough
            h_Ygbi1      = H_Yb1delta + log2(delta); % eq. (8.30) in [1], assuming delta is small enough
            h_YgBi(ibit) = PB(ibit,:) * [h_Ygbi0; h_Ygbi1]; % eq. (2.10) [1]
        end
    
        %------- Combine the three terms to obtain H(Bi|Y)
        H_BgY = H_Bi - h_Y + h_YgBi;
    
    %---- R_BMD = minuend - subtrahend
    BMDrate(isnr) = max(0, H_X - sum(H_BgY)); % in case it somehow turns out to be negative
end
%---------------------- END