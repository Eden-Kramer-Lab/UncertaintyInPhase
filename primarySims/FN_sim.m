function [data,snr, mk_phase,phaseBounds, hilb_phase, confLimits, poincarePhase] = ...
            FN_sim()

run fitzhugh_nagumo.m
data = 20*y(1:1e4,1);% + 1*randn(length(y),1);
tspan = 1/1000:1/1000:10;
Fs = 1000;
SNR = std(20*y(:,1))/1;
%%
params.tapers = [1,10,1]; %
params.Fs = Fs;
params.pad = -1;
params.fpass = [1:100];
[psds,f1]  = mtspectrumc(data,params); 
snr = sum(psds(f1>1 & f1<3))/sum(psds((f1<1 | f1 > 3) & f1 > 0)); %

%%
[omega, ampEst, allQ, R,stateVec, stateCov] = ...
fit_MKModel_multSines(detrend(data(1:1e4)),[1.3,4],1000,[.99,.99],[5,5],.1);

initParams.freqs = omega;
initParams.Fs = Fs;
initParams.ampVec = ampEst;
initParams.sigmaFreqs = allQ;
initParams.sigmaObs  =R;
[mk_phase,phaseBounds] = estimateMKphase(data,initParams,1);
%%

fNQ = Fs/2;
locutoff = .6;                               
hicutoff = 3;
filtorder = 2*fix(Fs/locutoff);
MINFREQ = 0;
trans  = 0.25;                      % fractional width of transition zones
f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts = firls(filtorder,f,m);             % get FIR filter coefficients

lowAct = filtfilt(filtwts,1,data);
hilb_phase = angle(hilbert(lowAct));
[confLimits] = hilbConfLimits(data,lowAct,Fs,4,.001);
%%

% defining temporal markers
timeMarker = [];
% [b,a]=butter(4, 30/(Fs/2), 'low');
% tmpData = filtfilt(b,a,data);
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
    poincarePhase(timeMarker(ii-1):timeMarker(ii)) = angle(cos(2*pi*[0:1/amtTime:1]) + 1i*sin(2*pi*[0:1/amtTime:1]));
end

tmp = wrapTo2Pi(unwrap(poincarePhase)+pi/2)-pi;
poincarePhase = tmp;
%%
makePlot_data_phaseEsts
makePlot_phaseDiffsAllPhaseEst