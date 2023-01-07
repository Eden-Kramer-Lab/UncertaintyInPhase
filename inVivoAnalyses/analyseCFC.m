
%% State Space model
num_bins = 20;

fullStateVec = initParams{2}.stateVec;
phaseBounds = initParams{2}.phaseBounds;
betaAmp = abs(zscore(fullStateVec(7,:)) + 1i*zscore(fullStateVec(8,:)));
lowDeltaAnalytic = (sum(fullStateVec([3],:),1) + 1i*sum(fullStateVec([4],:),1));

mkPhaseBoundsDiff = rad2deg(phaseBounds(:,2) - phaseBounds(:,1));
inds = (mkPhaseBoundsDiff < prctile(mkPhaseBoundsDiff,25)) ;

[pac_mk,~, bin_centers]=phase_pac(betaAmp',transpose(lowDeltaAnalytic),num_bins);
[pac_mk_thr]=phase_pac(betaAmp(inds)',transpose(lowDeltaAnalytic(inds)),num_bins);

figure
tmpPAC = detrend(pac_mk(:,1));
tmpPAC_se = pac_mk(:,2);
pacBounds = [tmpPAC - 2*tmpPAC_se, tmpPAC + 2*tmpPAC_se]; 
h1 = plot(bin_centers,tmpPAC,'red','Linewidth',1.5, 'Linestyle','-');
hold on
g3 = shade(bin_centers, pacBounds(:,1),bin_centers, pacBounds(:,2),...
'FillType',[2 1], 'FillAlpha', .4, 'FillColor', 'red', 'Color', 'red');

tmpPAC = detrend(pac_mk_thr(:,1));
tmpPAC_se = pac_mk_thr(:,2);
pacBounds = [tmpPAC - 2*tmpPAC_se, tmpPAC + 2*tmpPAC_se]; 
h2 = plot(bin_centers,tmpPAC,'red','Linewidth',1.5, 'Linestyle','--');
hold on
g3 = shade(bin_centers, pacBounds(:,1),bin_centers, pacBounds(:,2),...
'FillType',[2 1], 'FillAlpha', .4, 'FillColor', 'red', 'Color', 'red');
xlim([-pi,pi]);
grid on
set(gca,'Fontsize', 14)
ylabel('Beta Amp')
%% FIR-Hilbert

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

[confLimits] = hilbConfLimits(useData(:,2),lowDelta(:,2),1000,1.25,.001);
confLimDiff = rad2deg(confLimits(:,2) - confLimits(:,1));
% 
fNQ = Fs/2;
locutoff = 15;                               
hicutoff = 25;
filtorder = 3*fix(Fs/locutoff);
MINFREQ = 0;
transRight  = 0.1; transLeft = 0.1;           % fractional width of transition zones
f=[MINFREQ (1-transRight/10)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+transLeft)*hicutoff/fNQ 1];
m=[0       0                      1            1            0                      0];
filtwts_beta = firls(filtorder,f,m);             % get FIR filter coefficients

beta = filtfilt(filtwts_beta,1,useData);
betaAmp = abs(hilbert(zscore(beta(:,2))));

inds = (confLimDiff<prctile(confLimDiff,25));

[pac_hilb,~, bin_centers]=phase_pac(betaAmp,(hilbert(lowDelta(:,2))),num_bins);
[pac_hilb_thr]=phase_pac(betaAmp(inds),(hilbert(lowDelta(inds,2))),num_bins);

% yyaxis left
tmpPAC = detrend(pac_hilb(:,1));
tmpPAC_se = pac_hilb(:,2);
pacBounds = [tmpPAC - 2*tmpPAC_se, tmpPAC + 2*tmpPAC_se]; 
h3 = plot(bin_centers,tmpPAC,'blue','Linewidth',1.5, 'Linestyle','-');
hold on
g3 = shade(bin_centers, pacBounds(:,1),bin_centers, pacBounds(:,2),...
'FillType',[2 1], 'FillAlpha', .4, 'FillColor', 'blue', 'Color', 'blue');

tmpPAC = detrend(pac_hilb_thr(:,1));
tmpPAC_se = pac_hilb_thr(:,2);
pacBounds = [tmpPAC - 2*tmpPAC_se, tmpPAC + 2*tmpPAC_se]; 
h4 = plot(bin_centers,tmpPAC,'blue','Linewidth',1.5, 'Linestyle','--');
hold on
g3 = shade(bin_centers, pacBounds(:,1),bin_centers, pacBounds(:,2),...
'FillType',[2 1], 'FillAlpha', .4, 'FillColor', 'blue', 'Color', 'blue');
xlim([-pi,pi]);
grid on
set(gca,'Fontsize', 14)
ylabel('Normalized Beta Amplitude')
xlabel('Phase')
legend([h1,h2,h3,h4], {'MK Phase', 'MK Phase (Thr)', 'Hilb Phase' ,'Hilb Phase (Thr)'})

%% Plot data:
figure
betaAmp = abs(zscore(fullStateVec(5,:)) + 1i*zscore(fullStateVec(6,:)));
h  = cline((.001:.001:4e2),fullStateVec(1,:)+fullStateVec(5,:),(betaAmp.^(1/2.5))/max((betaAmp.^1/2.5)))
h.LineWidth = .9;
colormap(flipud(brewermap([],'spectral')))
hh = colorbar;
caxis([-.2,.7])
ylabel(hh,'Beta Amplitude')
set(hh,'YTickLabel',[])
grid on;hold on
i = plot((.001:.001:4e2),fullStateVec(1,:),'Linewidth',2)
j = plot((.001:.001:4e2),fullStateVec(5,:),'Linewidth',1.5, 'Color','blue')
xlim([185,200])
set(gca,'Fontsize', 14)
xlabel('Time (secs)')
legend(j, 'Beta band', 'boxoff')
ylabel('Signal')

