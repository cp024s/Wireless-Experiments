clearvars; clc;

% -------- Simulation parameters ----------------
nSym = 10^4; % Number of OFDM Symbols to transmit
EbN0dB = 0:2:20; % Bit to noise ratio
MOD_TYPE = 'MPSK'; % Modulation type - 'MPSK' or 'MQAM'
M = 64; % Choose modulation order for the chosen MOD_TYPE
N = 64; % FFT size or total number of subcarriers (used + unused) 64
Ncp = 16; % Number of symbols in the cyclic prefix

% -------- Derived Parameters --------------------
k = log2(M); % Number of bits per modulated symbol
EsN0dB = 10 * log10(k) + EbN0dB; % Convert to symbol energy to noise ratio
errors = zeros(1, length(EsN0dB)); % To store symbol errors

% Monte Carlo Simulation
for i = 1:length(EsN0dB)
    for j = 1:nSym
        % ----------------- Transmitter --------------------
        d = ceil(M * rand(1, N)); % Uniform distributed random syms from 1:M
        [X, ref] = modulation_mapper(MOD_TYPE, M, d);
        x = ifft(X, N); % IDFT
        s = add_cyclic_prefix(x, Ncp); % Add CP

        % ---------------- Channel ----------------
        r = add_awgn_noise(s, EsN0dB(i)); % Add AWGN noise r = s + n

        % ----------------- Receiver ----------------------
        y = remove_cyclic_prefix(r, Ncp, N); % Remove CP
        Y = fft(y, N); % DFT
        [~, dcap] = iqOptDetector(Y, ref); % Demapper using IQ detector

        % ---------------- Error counter ------------------
        numErrors = sum(d ~= dcap); % Count number of symbol errors
        errors(i) = errors(i) + numErrors; % Accumulate symbol errors
    end
end

simulatedSER = errors / (nSym * N);
theoreticalSER = ser_awgn(EbN0dB, MOD_TYPE, M);

% Plot theoretical curves and simulated BER points
plot(EbN0dB, log10(simulatedSER), 'ro'); hold on;
plot(EbN0dB, log10(theoreticalSER), 'r-'); grid on;
title(['Performance of ', num2str(M), '-', MOD_TYPE, ' OFDM over AWGN channel']);
xlabel('Eb/N0 (dB)'); ylabel('Symbol Error Rate');
legend('simulated', 'theoretical');

% Function definitions

function [s, ref] = mpsk_modulator(M, d)
    % Function to MPSK modulate the vector of data symbols - d
    % [s, ref] = mpsk_modulator(M, d) modulates the symbols defined by the
    % vector d using MPSK modulation, where M specifies the order of
    % M-PSK modulation and the vector d contains symbols whose values
    % in the range 1:M. The output s is the modulated output and ref
    % represents the reference constellation that can be used in demod
    ref_i = 1 / sqrt(2) * cos(((1:1:M) - 1) / M * 2 * pi);
    ref_q = 1 / sqrt(2) * sin(((1:1:M) - 1) / M * 2 * pi);
    ref = ref_i + 1i * ref_q;
    s = ref(d); % M-PSK Mapping
end

function [X, ref] = modulation_mapper(MOD_TYPE, M, d)
    % Modulation mapper for OFDM transmitter
    % MOD_TYPE - 'MPSK' or 'MQAM' modulation
    % M - modulation order, For BPSK M = 2, QPSK M = 4, 256-QAM M = 256, etc.,
    % d - data symbols to be modulated drawn from the set {1,2,...,M}
    % returns X - modulated symbols
    % ref - ideal constellation points that could be used by IQ detector
    if strcmpi(MOD_TYPE, 'MPSK')
        [X, ref] = mpsk_modulator(M, d);
        % MPSK modulation
    else
        if strcmpi(MOD_TYPE, 'MQAM')
            [X, ref] = mqam_modulator(M, d);
            % MQAM modulation
        else
            error('Invalid Modulation specified');
        end
    end
end

function y = remove_cyclic_prefix(r, Ncp, N)
    % Function to remove cyclic prefix from the received OFDM symbol r
    % r - received OFDM symbol with CP
    % Ncp - num. of samples at the beginning of r that need to be removed
    % N - number of samples in a single OFDM symbol
    % y - returns the OFDM symbol without cyclic prefix
    y = r(Ncp + 1:N + Ncp); % cut from index Ncp + 1 to N + Ncp
end

function s = add_cyclic_prefix(x, Ncp)
    % Function to add cyclic prefix to the generated OFDM symbol x
    % that is generated at the output of the IDFT block
    % x - OFDM symbol without CP (output of IDFT block)
    % Ncp - num. of samples at x's end that will be copied to its beginning
    % s - returns the cyclic prefixed OFDM symbol
    s = [x(end - Ncp + 1:end), x]; % Cyclic prefixed OFDM symbol
end

function [idealPoints, indices] = iqOptDetector(received, ref)
    % Optimum Detector for 2-dim. signals (MQAM,MPSK,MPAM) in IQ Plane
    % received - vector of form I + jQ
    % ref - reference constellation of form I + jQ
    % Note: MPAM/BPSK are one dim. modulations. The same function can be
    % applied for these modulations since quadrature is zero (Q=0).
    x = [real(received); imag(received)]'; % received vec. in cartesian form
    y = [real(ref); imag(ref)]'; % reference vec. in cartesian form
    [idealPoints, indices] = minEuclideanDistance(x, y);
end

function [idealPoints, indices] = minEuclideanDistance(x, y)
    % Function to compute the pairwise minimum Distance between two
    % vectors x and y in p-dimensional signal space and select the
    % vectors in y that provide the minimum distances.
    % x - a matrix of size mxp
    % y - a matrix of size nxp. This acts as a reference against
    % which each point in x is compared.
    % idealPoints - contain the decoded vector
    % indices - indices of the ideal points in reference matrix y
    [m, p1] = size(x);
    [n, p2] = size(y);
    if p1 ~= p2
        error('Dimension Mismatch: x and y must have the same dimension')
    end

    X = sum(x .* x, 2);
    Y = sum(y .* y, 2)';
    d = X(:, ones(1, n)) + Y(ones(1, m), :) - 2 * x * y'; % Squared Euclidean Dist.
    [~, indices] = min(d, [], 2); % Find the minimum value along DIM=2
    idealPoints = y(indices, :);
    indices = indices.';
end

function [SER] = ser_awgn(EbN0dB, MOD_TYPE, M, COHERENCE)
    % Theoretical Symbol Error Rate for various modulations over AWGN
    % EbN0dB - list of SNR per bit values
    % MOD_TYPE - 'BPSK', 'PSK', 'QAM', 'PAM', 'FSK'
    % M - Modulation level for the chosen modulation
    % - For PSK, PAM, FSK M can be any power of 2
    % - For QAM M must be an even power of 2 (square QAM only)
    % Parameter COHERENCE is only applicable for FSK modulation
    % COHERENCE = 'coherent' for coherent FSK detection
    % = 'noncoherent' for noncoherent FSK detection
    gamma_b = 10.^(EbN0dB / 10); % SNR per bit in linear scale
    gamma_s = log2(M) * gamma_b; % SNR per symbol in linear scale

    SER = zeros(size(EbN0dB));
    
    switch lower(MOD_TYPE)
        case 'bpsk'
            SER = 0.5 * erfc(sqrt(gamma_b));
        case {'psk', 'mpsk'}
            if M == 2 % for BPSK
                SER = 0.5 * erfc(sqrt(gamma_b));
            else
                if M == 4 % for QPSK
                    Q = 0.5 * erfc(sqrt(gamma_b));
                    SER = 2 * Q - Q.^2;
                else % for other higher order M-ary PSK
                    SER = erfc(sqrt(gamma_s) * sin(pi / M));
                end
            end
        case {'qam', 'mqam'}
            SER = 1 - (1 - (1 - 1 / sqrt(M)) * erfc(sqrt(3 / 2 * gamma_s / (M - 1)))).^2;
        case {'fsk', 'mfsk'}
            if strcmpi(COHERENCE, 'coherent')
                for ii = 1:length(gamma_s)
                    fun = @(q) (0.5 * erfc((-q - sqrt(2 * gamma_s(ii))) / sqrt(2))).^(M - 1) .* ...
                        1 / sqrt(2 * pi) .* exp(-q.^2 / 2);
                    SER(ii) = 1 - integral(fun, -inf, inf);
                end
            else % Default compute for noncoherent
                for jj = 1:length(gamma_s)
                    summ = 0;
                    for i = 1:M - 1
                        n = M - 1; r = i; % for nCr formula
                        summ = summ + (-1).^(i + 1) / (i + 1) * prod((n - r + 1:n) / (1:r)) * ...
                            exp(-i / (i + 1) * gamma_s(jj));
                    end
                    SER(jj) = summ; % Theoretical SER for non-coherent detection
                end
            end
        case {'pam', 'mpam'}
            SER = 2 * (1 - 1 / M) * 0.5 * erfc(sqrt(3 * gamma_s / (M^2 - 1)));
        otherwise
            display('ser_awgn.m: Invalid modulation (MOD_TYPE) selected');
    end
end

function [r, n, N0] = add_awgn_noise(s, SNRdB, L)
    % Function to add AWGN to the given signal
    % [r, n, N0] = add_awgn_noise(s, SNRdB) adds AWGN noise vector to signal
    % 's' to generate a resulting signal vector 'r' of specified SNR
    % in dB. It also returns the noise vector 'n' that is added to the
    % signal 's' and the spectral
    % [r, n, N0] = add_awgn_noise(s, SNRdB, L) adds AWGN noise vector to
    % signal 's' to generate a resulting signal vector 'r' of specified
    % SNR in dB. The parameter 'L' specifies the oversampling ratio used
    % in the system (for waveform simulation). It also returns the noise
    % vector 'n' that is added to the signal 's' and the spectral
    % density N0 of noise added
    s_temp = s;
    if iscolumn(s)
        s = s.';
    end % to return the result in the same dim as 's'
    gamma = 10^(SNRdB / 10); % SNR to linear scale
    if nargin == 2
        L = 1;
    end % if the third argument is not given, set it to 1
    if isvector(s)
        P = L * sum(abs(s).^2) / length(s); % Actual power in the vector
    else % for multi-dimensional signals like MFSK
        P = L * sum(sum(abs(s).^2)) / length(s); % if s is a matrix [MxN]
    end
    N0 = P / gamma; % Find the noise spectral density
    if(isreal(s))
        n = sqrt(N0 / 2) * randn(size(s)); % computed noise
    else
        n = sqrt(N0 / 2) * (randn(size(s)) + 1i * randn(size(s))); % computed noise
    end
    r = s + n; % received signal
    if iscolumn(s_temp)
        r = r.';
    end % return r in the original format as s
end


