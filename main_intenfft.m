Fs = 100;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 10001;             % Length of signal
t = (0:L-1)*T;        % Time vector

% S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);

% X = S + 2*randn(size(t));
% X = [];

Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
% axis([0 0.1 0 70])



%%



figure(2)
subplot(1,3,1)
plot(f,P1_85000)
subplot(1,3,2)
plot(f,P1_95000)
subplot(1,3,3)
plot(f,P1_95000-P1_85000)