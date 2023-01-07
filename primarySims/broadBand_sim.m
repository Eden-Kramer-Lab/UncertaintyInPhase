function [data, snr,mk_phase,phaseBounds, hilb_phase, confLimits, poincarePhase] = ...
            broadBand_sim(Fs,timeL)

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
[pn] = make_pink_noise(1,1e4,1/Fs);
pn = 7.5*pn;

data = (Vlo+pn');
% filter data below 100 Hz:
fNQ = Fs/2;
locutoff = 100;                              
hicutoff = 120;
filtorder = 5*fix(Fs/locutoff);
MINFREQ = 0;
f=[0 locutoff/fNQ hicutoff/fNQ 1];
m=[1       1           0       0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients
data = filtfilt(filtwts,1,data);

params.tapers = [1,10,1]; %
params.Fs = Fs;
params.pad = -1;
params.fpass = [1:100];
[psds,f1]  = mtspectrumc(data,params); 
snr = sum(psds(f1>4 & f1<8))/sum(psds((f1<4 | f1 > 8) & f1 > 0)); %2.5
%%
[omega, ampEst, allQ, R,stateVec, stateCov] = ...
fit_MKModel_multSines(detrend(data),[6,60],1000,[.99,.99],[.03,.01],.001);

initParams.freqs = omega;
initParams.Fs = Fs;
initParams.ampVec = ampEst;
initParams.sigmaFreqs = allQ;
initParams.sigmaObs  =R;
[mk_phase,phaseBounds] = estimateMKphase(data,initParams,1);

%%

fNQ = Fs/2;
locutoff = 4;                               % Low freq passband = [4,7] Hz.
hicutoff = 8;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
trans  = 0.15;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients

lowAct = filtfilt(filtwts,1,data);
hilb_phase = angle(hilbert(lowAct));
[confLimits] = hilbConfLimits(data,lowAct,1000,4,.001);

%%
% defining temporal markers
timeMarker = [];
x_vals_1 = sign(detrend(data));
cnt1 = 2;
timeMarker(1) = 1;
for ii = 3:length(x_vals_1)-1
    if x_vals_1(ii) == -1
        if x_vals_1(ii+1) == 1
            timeMarker(cnt1) = ii;
            cnt1 = cnt1 +1;
        end
    end  
end

timeMarker(cnt1) = length(x_vals_1);
poincarePhase = [];
%defining phase
for ii = 2:length(timeMarker)
        amtTime = timeMarker(ii) - timeMarker(ii-1);
    poincarePhase(timeMarker(ii-1):timeMarker(ii)) = [0:1/amtTime:1]*2*pi - pi;
end

tmp = wrapTo2Pi((-pi/2) + unwrap(poincarePhase))-pi;
poincarePhase = tmp;
%%

makePlot_data_phaseEsts
makePlot_phaseDiffsAllPhaseEst