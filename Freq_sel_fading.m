ts = 1e-5;         % Input sample period (seconds)
fd = 120;          % Maximum Doppler Shift (Hertz)
tau = [0, 1e-5, 1e-6];  % Path delays of 3 multipaths
pdb = [0, -10, -20];    % Path gains in dB

ral = rayleighchan(ts, fd, tau, pdb);  % Rayleigh Model Channel

n = 30000;         % Number of input bits
x = randi([0, 1], [n, 1]);  % Generating random input bits

M = 16;            % Modulation index of 16QAM
k = log2(M);       % Number of bits per sample in 16QAM

xd = bi2de(reshape(x, k, length(x)/k).', 'left-msb');  % Reshaping the stream of bits (x) into k bit samples for QAM modulation

modout = qammod(xd, M);  % QAM modulation

ral.StoreHistory = true;  % Enabling the Store Option for plotting the channel response

y = filter(ral, modout);  % Received signal due to transmission of x bits through the channel

figure(1)
plot(20*log10(abs(y)))  % Plotting the received power in dB
xlabel('Time (s)')
ylabel('Received Signal Power (dB)')

figure(2)
plot(ral)
