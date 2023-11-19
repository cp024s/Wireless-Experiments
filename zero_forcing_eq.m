clc;
clear all;
close all
% Generation of QAM transmit signal
a = [repmat(-3:2:3,1,4)]-[repmat(-3j,1,4) repmat(-1j,1,4)
repmat(1j,1,4) repmat(3j,1,4)];
Ak = a(randi(16,1000,1));
figure(1);
plot(real(a),imag(a),'ko');
title('16-QAM constellation');
grid;
xlabel('RE');
ylabel('IM');
% Generate Channel
h = [0.19+.56j .45-1.28j -.14-.53j -.19+.23j .33+.51j]; % LTE
channel
L1=1; L2=3;
figure(2);
stem(-L1:L2,abs(h),'r');
legend('ISI Channel');
title('Absolute values of impulse responses'); % Absolute values of
channel impulse response
% Received Signal
Rk = filter(h,1,Ak); % Received signal
figure(3);
plot(real(Rk),imag(Rk),'o');
grid on;
xlabel('RE');ylabel('IM');% Constellation with ISI
% Zeroforcing Equalization
N1=2; N2=5; % Equalizer length is N1+N2+1
H=convmtx(h.',N1+N2+1); % Convolution matrix
q_ZF=zeros(N1+N2+L1+L2+1,1); % q vector (zero forcing)
q_ZF(L1+N1+1)=1; % Puts the one in the place corresponding to
desired symbol
w_LS=((H'*H)\(H'))*q_ZF; % ZF equalizer coefficients using LS
estimate
% Received signal after equalization
Rk_eq=filter(w_LS,1,Rk); % Equalized data
% Comparisons
figure(4);
hold on;
plot(real(Rk_eq),imag(Rk_eq),'r.'); % Equalized constellation
plot(real(a),imag(a),'ko', 'MarkerSize',6,'MarkerFaceColor','k');
hold off; %Ideal constellation
grid on;
legend('Equalizeddata','Ideal constellation points');
