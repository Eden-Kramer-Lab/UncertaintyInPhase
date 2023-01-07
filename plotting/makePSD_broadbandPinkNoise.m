Fs = 1000;
timeL = 10;

tspan = 1/Fs:1/Fs:timeL;

power  = 5.25;
fbins = 1/timeL; % our best freq bin is 1/time(secs) in Hz
realFreqs = 0:fbins:Fs/2;
gaussMean = 6;
gaussStdDev =1;
gaussFunc = exp(-(realFreqs - gaussMean).^2/(2*gaussStdDev^2));
gaussFunc = (gaussFunc./(2*sum(gaussFunc)));
gaussFunc = gaussFunc * (power/max(gaussFunc));
gaussFunc = [gaussFunc,fliplr(gaussFunc(2:end-1))];

% create signal
origNoise = randn(Fs*timeL,1)*power; % 
origNoiseFFT = fft(origNoise);
Vlo = ifft(origNoiseFFT.*gaussFunc');
[pn, pn_FT] = make_pink_noise(1,1e4,1/Fs);
pn = 7.5*pn;

plot(0:1/(timeL):1000-(1/timeL), 10*log10(abs(gaussFunc) + abs(7.5*pn_FT)),'linewidth',1.5)
hold on
plot(0:1/(timeL):1000-(1/timeL),10*log10(abs(7.5*pn_FT)),'linewidth',1.5)
ylim([-10,20])
xlim([1,50])
set(gca,'Fontsize', 14)
ylabel('Power (dB)')
xlabel('Frequency (Hz)')