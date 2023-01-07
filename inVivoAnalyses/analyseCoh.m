
%% State Space Model
for i = 1:2
    tic
    dataToFit = useData(1e4+1:2e4,i);
    [omega, ampEst, allQ, R, stateVec, stateCov] = fit_MKModel_multSines(zscore(dataToFit),[.5,2.4,20], Fs,[.99,.99,.99],...
                                                            [.1,.1,.1],.01);

    rhythmLoc = find(omega > .1 & omega < 1);       
    if ~isempty(rhythmLoc)

        initParams{i}.freqs = omega;
        initParams{i}.Fs = Fs;
        initParams{i}.ampVec = ampEst;
        initParams{i}.sigmaFreqs = allQ;
        initParams{i}.sigmaObs  =R;
        initParams{i}.stateVec = stateVec;
        [mk_phase,phaseBounds,~,newAllX] = estimateMKphase(zscore(useData(:,i)),initParams{i},rhythmLoc);
        initParams{i}.phaseEst = mk_phase;
        initParams{i}.phaseBounds = phaseBounds;
        initParams{i}.stateVec = newAllX;

    end
    toc
end
 % analysing the results
   
newAllX1 = initParams{1}.stateVec;
newAllX2 = initParams{2}.stateVec;
   
phaseBounds1 = initParams{1}.phaseBounds;
phaseBounds2 = initParams{2}.phaseBounds;
credWidth2 = phaseBounds2(:,2) - phaseBounds2(:,1);
credWidth1 = phaseBounds1(:,2) - phaseBounds1(:,1);
inds = find(credWidth1<prctile(credWidth1,25) & credWidth2<prctile(credWidth2,25));

figure

subplot(211)
unThrCross = [newAllX1(1,:)-1i*newAllX1(2,:)]'.*[newAllX2(1,:)+1i*newAllX2(2,:)]';
% X = transpose(newAllX1(1,:) + 1i*newAllX1(2,:));
% Y = transpose(newAllX2(1,:) + 1i*newAllX2(2,:));
% unThrCrossCorr = (abs(mean(X.*conj(Y))) ./ sqrt(mean(X.*conj(X)) * mean(Y.*conj(Y))));
unThrCrossCorr = abs(corr([newAllX1(1,:)-1i*newAllX1(2,:)]', [newAllX2(1,:)-1i*newAllX2(2,:)]'));

polarhistogram(angle(unThrCross),100)
title(['Coherence = ',num2str(round(unThrCrossCorr,2))])
set(gca,'Fontsize',14)

subplot(212)
thrCross = ([newAllX1(1,inds)-1i*newAllX1(2,inds)]'.*[newAllX2(1,inds)+1i*newAllX2(2,inds)]');
thrCrossCorr = abs(corr([newAllX1(1,inds)+1i*newAllX1(2,inds)]', [newAllX2(1,inds)+1i*newAllX2(2,inds)]'));
polarhistogram(angle(thrCross),100)
title(['Coherence = ', num2str(round(thrCrossCorr,2))])
set(gca,'Fontsize',14)

%% FIR-Hilbert
% FIR being set up for delta below:
fNQ = Fs/2;
locutoff = .25;                               
hicutoff = 1.5;
filtorder = 2*fix(Fs/locutoff);
MINFREQ = 0;
transRight  = 0.1; transLeft = 0.1;           % fractional width of transition zones
f=[MINFREQ (1-transRight/10)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+transLeft)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_delta = firls(filtorder,f,m);             % get FIR filter coefficients

lowDelta = filtfilt(filtwts_delta,1,useData);

[confLimits] = hilbConfLimits(useData(:,1), lowDelta(:,1),1000,1.75,.001);
confLimDiff(:,1) = rad2deg(confLimits(:,2) - confLimits(:,1));
[confLimits] = hilbConfLimits(useData(:,2),lowDelta(:,2),1000,1.75,.001);
confLimDiff(:,2) = rad2deg(confLimits(:,2) - confLimits(:,1));

figure
inds_hilb = find(confLimDiff(:,1)<prctile(confLimDiff(:,1),25) & confLimDiff(:,2)<prctile(confLimDiff(:,2),25));
%bit circlar here since amp is part of conf lim and coherence, so high amps
%=> high coh. but this puts no constrsint on phase.
subplot(211)
unThrCross = hilbert(lowDelta(:,1)) .* conj(hilbert(lowDelta(:,2)));
unThrCrossCorr = abs(corr(hilbert(lowDelta(:,1)) , (hilbert(lowDelta(:,2)))));
polarhistogram(angle(unThrCross),100)
title(['Coherence = ',num2str(round(unThrCrossCorr,2))])
set(gca,'Fontsize',14)

subplot(212)
thrCross = hilbert(lowDelta(inds_hilb,1)) .* conj(hilbert(lowDelta(inds_hilb,2)));
thrCrossCorr = abs(corr(hilbert(lowDelta(inds_hilb,1)) , (hilbert(lowDelta(inds_hilb,2)))));
polarhistogram(angle(thrCross),100)
title(['Coherence = ', num2str(round(thrCrossCorr,2))])
set(gca,'Fontsize',14)


% 
% [b,a]=butter(2, 5/(Fs/2), 'low');
% highDelta = filtfilt(b,a,useData);
% [b,a]=butter(2, 2.25/(Fs/2), 'high');
% highDelta = filtfilt(b,a,highDelta);
% 
% [b,a]=butter(2, 8/(Fs/2), 'high');
% beta = filtfilt(b,a,useData);
% [b,a]=butter(2, 30/(Fs/2), 'low');
% beta = filtfilt(b,a,beta);

%% Poincare section

% defining temporal markers
timeMarker = [];
% [b,a]=butter(4, 30/(Fs/2), 'low');
% tmpData = filtfilt(b,a,data);
x_vals_1 = sign(detrend(useData(:,2)));
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


tmp = wrapToPi(-(pi/2)+ conv(unwrap(poincarePhase),ones(100,1),'same')/100);
poincarePhase(2,:) = tmp;

%%
figure 

subplot(211)
unThrCross = exp(1i*phase_poincare1)' .* exp(1i*(-phase_poincare2))';
unThrCrossCorr = abs(corr(exp(1i*phase_poincare1)' , exp(1i*(phase_poincare2))'));
polarhistogram(angle(unThrCross),100)
title(['Coherence = ',num2str(round(unThrCrossCorr,2))])
set(gca,'Fontsize',14)

subplot(212)
inds_new = intersect(inds, inds_hilb);
thrCross = exp(1i*phase_poincare1(inds_new)) .* exp(1i*(-phase_poincare2(inds_new)));
thrCrossCorr = abs(corr(exp(1i*phase_poincare1(inds_new))' , exp(1i*(phase_poincare2(inds_new)))'));
polarhistogram(angle(thrCross),100)
title(['Coherence = ', num2str(round(thrCrossCorr,2))])
set(gca,'Fontsize',14)

