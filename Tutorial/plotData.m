function plotData(time, data, cellName)
% PLOTDATA plots the cell activity in the left subplot, and a histogram of
% the data in a right subplot

f1 = figure;
nbins = 40;

subplot(1,2,1)
plot(time, data, 'k', 'linewidth',1.5)
xlabel('Time (s)')
ylabel('Change in Fluorescence')
title([cellName,' Activity'])

hold on;
subplot(1,2,2)
histogram(data, nbins, 'FaceColor','k')
xlabel('Change in Fluorescence')
ylabel('Number of Bins')
title([cellName, ' Histogram'])

set(f1, 'pos',[10 10 1400 600]);