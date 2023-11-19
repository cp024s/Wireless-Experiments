clc;
clear all;
close all;
%Simulate Rayleigh Channel
s=1.5; %Variance(sigma) of channel
h=s*(randn(1,10e4) + (1i*randn(1,10e4)));
[val,bin]=hist(abs(h),100);
plot(bin,val/trapz(bin,val),'r')
%Theoretical Rayleigh PDF
r=0:0.1:10;
pdf=(r/s^2).*exp(-r.^2/(2*s^2));
hold on
plot(r,pdf,'k')
title('Simulated and Theoretical PDF of Rayleigh Distribution')
xlabel('r')
ylabel('PDF f(r)')
legend('Simulated PDF','Theoretical PDF')
