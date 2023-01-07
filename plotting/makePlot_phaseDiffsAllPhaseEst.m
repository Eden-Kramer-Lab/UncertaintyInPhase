%%
% can estimate the relative error with respect to one another

figure
cnt = 1;
for i = [25,0]
tmpHilbConf = rad2deg(confLimits(:,2) - confLimits(:,1));
tmpMKConf = rad2deg(phaseBounds(:,2) - phaseBounds(:,1));
if cnt == 1
    inds = find((tmpHilbConf < prctile(tmpHilbConf,i)) & (tmpMKConf < prctile(tmpMKConf,i)));
else
    inds = find((tmpHilbConf > prctile(tmpHilbConf,i)) & (tmpMKConf > prctile(tmpMKConf,i)));
end

subplot(2,3,1 +(cnt-1)*3)
polarhistogram(poincarePhase(inds) - hilb_phase(inds)',25)
if cnt ==1; title('Poincaré/Hilbert');%rlim([0,800]);
% else;rlim([0,2000]);
end

set(gca,'Fontsize', 14)
subplot(2,3,2+(cnt-1)*3)
polarhistogram(poincarePhase(inds) - mk_phase(inds),25)
if cnt ==1; title('State Space/Poincaré');%rlim([0,800]);
    %else;rlim([0,2000]); 
end
    
set(gca,'Fontsize', 14)
subplot(2,3,3+(cnt-1)*3)
polarhistogram(hilb_phase(inds) - mk_phase(inds)',25)
if cnt ==1; title('Hilbert/State Space');%rlim([0,800]);
    %else; rlim([0,2000]);
end
set(gca,'Fontsize', 14)

cnt = cnt + 1;
end