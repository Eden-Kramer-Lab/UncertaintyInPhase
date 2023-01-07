function [data,snr, mk_phase,phaseBounds, hilb_phase, confLimits, poincarePhase] = ...
            AR_2_sim(Fs,timeL)

% ang_var2dev = @(v) sqrt(-2*log(1-v));

tspan = 1/Fs:1/Fs:timeL;

%% generate AR(2) data
coeffVal = .994;
freqs = 6;

% generate the AR(2) data   
r = coeffVal * (cos(2*pi*(freqs/Fs)) +  1i* sin(2*pi*(freqs/Fs)));
coeffs = poly([r,r']);
sigma = 0.1;

if ~isreal(coeffs)
    disp('something wrong with roots')
end

Vlo = [randn(1,1)*sigma,1];
noiseAdded = zeros(1,Fs*timeL);
for i = 2:Fs*timeL-1
    noiseAdded(i) = randn*sigma;    
    Vlo(i+1) = - coeffs(2)*Vlo(i) - coeffs(3)*Vlo(i-1) +noiseAdded(i);
end

% V = [1,1; ...
% 1/r, 1/(r')];
% V_inv = inv(V);
% u = zeros(Fs*time,2);
% for i = 1:(Fs*time)-1
%     u(i+1,:) = V_inv * Vlo(i+1:-1:i)';    
% end
% truePhase = angle(u(:,1));

Vlo = Vlo - mean(Vlo);
data = Vlo;

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

%%
params.tapers = [1,10,1]; %
params.Fs = Fs;
params.pad = -1;
params.fpass = [1:100];
[psds,f1]  = mtspectrumc(data,params); 
snr = sum(psds(f1>4 & f1<8))/sum(psds((f1<4 | f1 > 8) & f1 > 0)); %
%%
% initial parameters for the MK model
freqs = [6];
ampVec = [.99];
sigmaFreqs = [.5];
sigmaObs = .01;
window = 2000;
lowFreqBand = [4,8];

[omega, ampEst, allQ, R] = fit_MKModel_multSines(data,freqs, Fs,ampVec, sigmaFreqs,sigmaObs);

initParams.freqs = omega;
initParams.Fs = Fs;
initParams.ampVec = ampEst;
initParams.sigmaFreqs = allQ;
initParams.sigmaObs  =R;
[mk_phase,phaseBounds] = estimateMKphase(data,initParams,1);

    %%
   
% FIR being set up for theta [4,8] below:
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
hilb_phase = angle(hilbert(lowAct))';
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
%angle(cos(2*pi*[0:1/amtTime:1]) + 1i*sin(2*pi*[0:1/amtTime:1]));
end

tmp = wrapTo2Pi((-pi/2) + unwrap(poincarePhase))-pi;
poincarePhase = tmp;

%%
makePlot_data_phaseEsts
makePlot_phaseDiffsAllPhaseEst