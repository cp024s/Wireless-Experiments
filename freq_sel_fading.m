ts=1e-5; %Input sample period(seconds)
fd=120; %Maximum Doppler Shift(Hertz)
tau=[0,1e-5,1e-6]; %Path delays of 3 multipaths
pdb=[0,-10,-20]; %Path gains in dB
ral=rayleighchan(ts,fd,tau,pdb); %Rayleigh Model Channel
n=30000; %Number of input bits
x=randi(1,[n,1]); %Generating random input bits
y=filter(ral,x); %Received signal due to transmission of x bits through the channel
figure(1)
plot(20*log10(abs(y))) %Plotting the received power in dB
xlabel('Time (s)')
ylabel('Received Signal Power (dB)')
