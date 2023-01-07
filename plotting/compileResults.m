ang_var2dev = @(v)rad2deg((sqrt(-2*log(1-v))));
simData = simData_FN;
Fs = 1000;

for i = 1:1000
    hilbPhase = simData(i).hilbPhase;
    mk_phase = simData(i).mk_phase;
    poincarePhase = simData(i).poincarePhase;
    tmpData = simData(i).data;
    
    fNQ = Fs/2;
    locutoff = 4;                               % Low freq passband = [4,7] Hz.
    hicutoff = 8;
    filtorder = 3*fix(Fs/locutoff);
    MINFREQ = 0;
    trans  = 0.15;                      % fractional width of transition zones
    f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
    m=[0       0                      1            1            0                      0];
    filtwts = firls(filtorder,f,m);             % get FIR filter coefficients

    lowAct = filtfilt(filtwts,1,tmpData);
    allAmp(i) = mean(abs(hilbert(lowAct)));
    
    angVarsAll(i,1) = ang_var2dev(1-abs(mean(exp(1i*(hilbPhase - poincarePhase')))));
    angVarsAll(i,2) = ang_var2dev(1-abs(mean(exp(1i*(mk_phase - poincarePhase))))); 
    angVarsAll(i,3) = ang_var2dev(1-abs(mean(exp(1i*(hilbPhase - mk_phase')))));
    
    mk_CI = simData(i).mk_CI;
    hilbConf = simData(i).hilbConf;
    
    tmpHilbConf = rad2deg(hilbConf(:,2) - hilbConf(:,1));
    tmpMKConf = rad2deg(mk_CI(:,2) - mk_CI(:,1));
    allMK_CI(i) = rad2deg(angle(mean(exp(1i*(mk_CI(:,2) - mk_CI(:,1))))));
    
    inds = find((tmpMKConf < prctile(tmpMKConf,25)) | (tmpHilbConf < prctile(tmpHilbConf,25)));
    angVars_CI(i,3) = ang_var2dev(1-abs(mean(exp(1i*(hilbPhase(inds) - mk_phase(inds)')))));
    angVars_CI(i,1) = ang_var2dev(1-abs(mean(exp(1i*(hilbPhase(inds) - poincarePhase(inds)')))));
    angVars_CI(i,2) = ang_var2dev(1-abs(mean(exp(1i*(mk_phase(inds) - poincarePhase(inds))))));
    
%     inds = find((tmpHilbConf < prctile(tmpHilbConf,25)));
%     angVarsAll_hilb(i,3) = ang_var2dev(1-abs(mean(exp(1i*(hilbPhase(inds) - mk_phase(inds)')))));
%     angVarsAll_hilb(i,1) = ang_var2dev(1-abs(mean(exp(1i*(hilbPhase(inds) - poincarePhase(inds)')))));
%     angVarsAll_hilb(i,2) = ang_var2dev(1-abs(mean(exp(1i*(mk_phase(inds) - poincarePhase(inds))))));
end

%% Make violin plot
% figure
% violinplot(angVarsAll,0:2,parula(5))
% hold on
% violinplot(angVars_MK_CI,.5:2.5,winter(5))
% hold on
% violinplot(angVarsAll_hilb,.5:2.5,hot(5))
%%

ang_var2dev = @(v) sqrt(-2*log(1-v));
y = zeros(1000, 2, 3);

y(:,1,:) = angVarsAll;

y(:,2,:) = angVars_CI;

% y(:,3,:) = angVarsAll_hilb;

tmp = brewermap(8, 'RdBu');
tmpCol1 = tmp(1,:);
tmpCol2 = tmp(3,:);
tmpCol3 = tmp(7,:);
% tmpCol4 = tmp(4,:);
x = [1:2];
figure;
h = iosr.statistics.boxPlot(x,y,...
'symbolColor','k',...
'medianColor','k',...
'symbolMarker',{'+','o','d'},...
'boxcolor',{tmpCol1;tmpCol2;tmpCol3},...
'groupLabels',{'Hilbert/Poincaré','State Space/Poincaré','State Space/Hilbert'},...
'showLegend',true);
box on
set(gca,'Fontsize', 16)
set(gca,'XTickLabel', {'All Data', 'High Confidence'})
ylabel('Circular Standard Deviation')