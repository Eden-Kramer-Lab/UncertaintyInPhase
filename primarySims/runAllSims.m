% Code to run the broadband and pink noise and save out data
Fs = 1000;
timeL = 10;
simData_broadband = struct();
parfor i = 1:1000
   [data,snr, mk_phase,phaseBounds, hilb_phase, confLimits, poincarePhase] = ...
            broadBand_sim(Fs,timeL); 
    simData_broadband(i).data = data;
    simData_broadband(i).hilbPhase = hilb_phase;
    simData_broadband(i).hilbConf = confLimits;
    simData_broadband(i).mk_phase = mk_phase;
    simData_broadband(i).mk_CI = phaseBounds;
    simData_broadband(i).poincarePhase = poincarePhase;
    simData_broadband(i).snr = snr;
end

%%
% Code to run AR2 and save out data
Fs = 1000;
timeL = 10;
simData_ar2 = struct();
parfor i = 1:1000
   [data,snr, mk_phase,phaseBounds, hilb_phase, confLimits, poincarePhase] = ...
            AR_2_sim(Fs,timeL); 
    simData_ar2(i).data = data;
    simData_ar2(i).hilbPhase = hilb_phase;
    simData_ar2(i).hilbConf = confLimits;
    simData_ar2(i).mk_phase = mk_phase;
    simData_ar2(i).mk_CI = phaseBounds;
    simData_ar2(i).poincarePhase = poincarePhase;
    simData_ar2(i).snr = snr;
end

%%

% Code to run FN and save out data
Fs = 1000;
timeL = 10;
simData_FN = struct();
parfor i = 1:1000
   [data,snr, mk_phase,phaseBounds, hilb_phase, confLimits, poincarePhase] = ...
            FN_sim(); 
    simData_FN(i).data = data;
    simData_FN(i).hilbPhase = hilb_phase;
    simData_FN(i).hilbConf = confLimits;
    simData_FN(i).mk_phase = mk_phase;
    simData_FN(i).mk_CI = phaseBounds;
    simData_FN(i).poincarePhase = poincarePhase;
    simData_FN(i).snr = snr;
end
