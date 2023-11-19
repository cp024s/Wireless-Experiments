ts=1e-6; %Input sample period(seconds)
fd=100; %Maximum Doppler Shift(Hertz)
ral=rayleighchan(ts,fd); %Rayleigh Model Channel
n=30000; %Number of input bits
x=randi(1,[n,1]); %Generating random input bits
y=filter(ral,x); %Received signal due to transmission of x bits through the channel
plot(20*log10(abs(y))) %Plotting the received power in dB
xlabel('Time (s)')
ylabel('Received Signal Power (dB)')
