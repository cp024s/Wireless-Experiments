% OFDM Channel Estimation using Least Squares Method
clc;
clear all;
close all;

%%%%% Assumptions %%%%%
Nfft = 32;          % Number of FFT subcarriers
Ng = Nfft/8;        % Number of CP symbols
Nofdm = Nfft + Ng;  % Total length of one OFDM symbol

Nps = 4;            % Spacing between the pilot subcarriers
Np = Nfft / Nps;    % Number of pilot symbols per OFDM symbol
Nm = Nfft - Np;     % Number of message symbols (unknown symbols)

% Generation of Pilot symbols (Known symbols)
Xp = 2 * (randn(1, Np) > 0) - 1;  % Pilot symbols

% Generation of message symbols
msgint = randi([0, 15], 1, Nm);   % Unknown message
Nbps = 4;                          % QAM modulation
M = 2^Nbps;                        % Number of possible QAM symbols
Es = 1;
A = sqrt(3 / (2 * (M - 1) * Es));  % Normalization factor
Data = A * qammod(msgint, M);

% Loading both pilot and message symbols onto subcarriers
X = zeros(1, Nfft);  % New OFDM symbol
l = 1;
ll = 1;

for k = 1:Nfft
    if mod(k, 4) == 1
        X(k) = Xp(l);
        l = l + 1;
    else
        X(k) = Data(ll);
        ll = ll + 1;
    end
end

% Taking IFFT and adding cyclic prefix for transmission
xtx = ifft(X, Nfft);
x_cp = [xtx(:, (Nfft - Ng + 1):Nfft), xtx];  % OFDM symbol + CP

% Generating channel
h = (1/2) * [randn + 1i * randn, randn + 1i * randn];  % Random Rayleigh channel of two taps

% Received Signal
y = conv(x_cp, h);
snr = 30;  % SNR = 30 dB for the AWGN noise
y_noise = awgn(y, snr, 'measured');  % Addition of Gaussian noise

% At the Receiver
% Removal of CP and Taking FFT for the received signal
y_minus_CP = y_noise(:, (Ng + 1):(Nfft + Ng));
Y_rx = fft(y_minus_CP, Nfft);

% Estimating the Sub-channel information at the pilot location using LS method
% Finding the pilot locations
pilot_loc = [];

for k = 1:Nfft
    if mod(k, Nps) == 1
        pilot_loc = [pilot_loc, k];
    end
end

H_LS = Y_rx(pilot_loc) ./ Xp;  % LS estimation of subchannels at pilot location

% Interpolation, to estimate the sub-channel gains in remaining subchannels
slope = (H_LS(end) - H_LS(end - 1)) / (pilot_loc(end) - pilot_loc(end - 1));
H_LS = [H_LS, H_LS(end) + slope * (Nfft - pilot_loc(end))];
pilot_loc = [pilot_loc, Nfft];
H_interpolated = interp1(pilot_loc, H_LS, [1:Nfft]);

% Performance comparison
H = fft(h, Nfft);
H_power_dB = 10 * log(abs(H .* conj(H)));
H_interpolated_power_dB = 10 * log(abs(H_interpolated .* conj(H_interpolated)));

plot(H_power_dB, 'b');
hold on;
plot(H_interpolated_power_dB, 'r:+');
legend('True Channel', 'LS Estimated Channel');
