%%
figure
g(1) = subplot(211)
plot(tspan, data(1:1e4),'Linewidth',2)
grid on
% xlim([0,3])
ylabel('Voltage')
set(gca,'Fontsize', 16)
xlabel('Time (s)')

tmp = brewermap(8, 'RdBu');
tmpCol1 = tmp(1,:);
tmpCol2 = tmp(3,:);
tmpCol3 = tmp(7,:);

g(2) = subplot(212)
h3 = plot(tspan, poincarePhase,'black','Linewidth',2)
hold on

h1 = plot(tspan, mk_phase, 'Color', tmpCol1,'Linewidth',2);
g3 = shade(tspan, phaseBounds(:,1),tspan, phaseBounds(:,2),...
'FillType',[2 1], 'FillAlpha', .4, 'FillColor', 'red', 'Color', tmpCol1);

h2 = plot(tspan, hilb_phase,'Color', tmpCol3,'Linewidth',2);
g2 = shade(tspan, confLimits(:,1),tspan, confLimits(:,2),...
'FillType',[2 1], 'FillAlpha', .4, 'FillColor', 'blue', 'Color', tmpCol3);

grid on
% xlim([0,3])
set(gca,'Fontsize', 16)
xlabel('Time (s)')
ylabel('Phase')
hleg = legend([h1,h2,h3],'location','best');
hleg.String = {'State Space', 'Hilbert', 'Poincaré'};
linkaxes(g, 'x')
