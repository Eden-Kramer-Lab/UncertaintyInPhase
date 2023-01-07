% making a plot of amplitude against credible interval with the difference
% between methods shown in color

ang_var2dev = @(v)rad2deg((sqrt(-2*log(1-v))));
simData = simData_broadband;
Fs = 1000;

for i = 1:50
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
    allAmp(i,:) = (abs(hilbert(lowAct)));
    
    angVarsAll(i,1,:) = (((angle(exp(1i*(hilbPhase - poincarePhase'))))));
    angVarsAll(i,2,:) = (((angle(exp(1i*(mk_phase - poincarePhase)))))); 
    angVarsAll(i,3,:) = (((angle(exp(1i*(hilbPhase - mk_phase'))))));
    
    mk_CI = simData(i).mk_CI;
    hilbConf = simData(i).hilbConf;
    
    tmpHilbConf = rad2deg(hilbConf(:,2) - hilbConf(:,1));
    tmpMKConf = rad2deg(mk_CI(:,2) - mk_CI(:,1));
    allMK_CI(i,:) = (rad2deg(wrapTo2Pi(angle(exp(1i*(mk_CI(:,2) - mk_CI(:,1))))))); 
end

%%
grpLbls = {'Hilbert/Poincaré','State Space/Poincaré','State Space/Hilbert'};
for ii  = 1:3
tmpAngVars = squeeze(angVarsAll(:,ii,:));

cnt = 1; n = 200;
for i = unique(round(allAmp(:)*n,0))'
   newAllAmp(cnt) = i/n;
   newMK_CI(cnt) = rad2deg(angle(mean(exp(1i*deg2rad(allMK_CI(round(allAmp(:)*n,0) == i))))));
   newErrorBw(cnt) = (ang_var2dev(1 - abs(mean(exp(1i*(tmpAngVars(round(allAmp(:)*n,0) == i)))))));
   cnt = cnt +1;
end
figure
tmp1 = 1 - (abs(newErrorBw')/100);
tmp2 = tmp1 - .65;
tmp2 = tmp2/max(tmp2);
tmp2 = 10.^(tmp2)/(max(10.^(tmp2)));
scatter(newAllAmp(:),newMK_CI(:),60,[tmp2,...
            zeros(length(newErrorBw),1),zeros(length(newErrorBw),1)],'filled')
set(gca,'Fontsize', 14)
ylabel('SS Credible Interval')
xlabel('Amplitude')
title(grpLbls{ii})
end
%%
inds = (round(allAmp(:),0) ==5);
hold on
scatter(wrapTo180(allMK_CI(inds)), wrapTo180(tmpAngVars(inds)),40,'filled');