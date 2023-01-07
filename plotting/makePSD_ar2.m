Fs = 1000;

coeffVal = .995;
freqs = 6;

% generate the AR(2) data   
r = coeffVal * (cos(2*pi*(freqs/Fs)) +  1i* sin(2*pi*(freqs/Fs)));
coeffs = poly([r,r']);
coeffs = -(coeffs(2:3));
sigma = 0.1;
cnt = 1;
f = [];
for g =0: 2*pi/1e4:pi
    den = 1 + coeffs(1)^2 + 2*coeffs(2) +coeffs(2)^2+ ...
            2*(coeffs(1)*coeffs(2) - coeffs(1))*cos(g) - 4*coeffs(2) * (cos(g)^2);
    f(cnt) = (sigma)/ (2*pi*(den));
cnt = cnt + 1;
end


plot(Fs*[0:2*pi/1e4:pi]/(2*pi), 10*log10(f/Fs) ,'Linewidth',1.5 )
% ylim([-10,20])
xlim([1,50])
set(gca,'Fontsize', 14)
ylabel('Power (dB)')
xlabel('Frequency (Hz)')