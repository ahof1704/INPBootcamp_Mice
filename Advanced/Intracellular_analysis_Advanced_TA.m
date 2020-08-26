
% JAC 8-20-2019, updated 8-25-2019
% Modified by CJB 8-24-2020

% In each section, use the suggestions to fill in the variables and the
% rest of the necessary code.

% Data for this task are an intracellular trace from layer 5 cortical neurons in V1 cortex of an awake mouse 
% and several vectors of data about the visual stimuli that were presented
% during the recording.

%Start by loading 'Cell1_data.mat'

% The membrane potential is measured at 2000 Hz with clock precision of 1000 kHz
% The stimulus onsets and offset timestamps have 1000 kHz precision

% intra: membrane voltage (Vm) with spikes (2000 Hz)
% Vm: same trace, with the spikes removed to leave only Vm (2000 Hz)
% timestamp: timestamp vector for the intracellular recording (corresponds
% with the 'Vm' and 'intra' variables, but has 1000 kHz precision)
% on: timestamps of the stimulus onsets (1000 kHz precision)
% off: timestamps of the stimulus offsets (1000 kHz precision)
% stim: orientations of the drifting gratings stimuli that were presented.
% (degrees)

%% Settings
clear, close all


%%  What's in the data?
% Take a look at the intracellular trace.  Plot the intra against the timestamp vector. 

load('../Data/Cell1_data.mat'); % load data
figure('Name','Intracellular trace'); 
plot(timestamp,intra); % plot intra against timestamp
title('Intracellular trace') 
ylabel('Vm')
xlabel('Time (s)')
set(gca,'FontSize',15)

%%  What/when were the stimuli?
% Using the plot you made of the intra trace, add dots to represent the 
% onset and offset times of the drifting grating stimuli. One way to
% accomplish this is to scatter an arbitrary scalar quantity against the
% values found in the 'on' variable. For the matlab scatter function, the
% 'x' and 'y' data to scatter must be the same length. So make a vector
% the same size as the 'on' variable, and then scale it to facilitate
% visualiztion (hint: the maximum of the intra variable works well!)
% Check out matlabs documentation on the 'scatter' function if you get
% stuck, or ask your TA.  
% 
% Note: The 'hold on' command allows you to add plots to an existing figure. 

figure('Name','Intracellular trace and Stimulus'); 
hold on;
plot(timestamp,intra);
scatter(on,ones(size(on))*max(intra),'g');  %plot the onsets (preferably in green)
scatter(off,ones(size(on))*max(intra),'r');  %plot the offsets (preferably in red)
title('Intracellular trace and Stimulus') 
ylabel('Vm')
xlabel('Time (s)')
legend({'Trace', 'Onset', 'Offset'});
set(gca,'FontSize',15)

%%  Next, try this again using subplots.  
% Check the 'subplot' function help for details.
% replot the traces and points as above in the first of two subplots
figure('Name','Intracellular trace, Stimulus and Orientation'); 
ax1 = subplot(211);  % first of two subplots
hold on; % allow ploting of multiple plots
plot(timestamp,intra); % plot the intracellular trace again
scatter(on,ones(size(on))*max(intra),'g');  %plot the onsets (check the scatter function help for symbol details)
scatter(off,ones(size(on))*max(intra),'r');  %plot the offsets the same way
title('Intracellular trace and Stimulus') 
ylabel('Vm')
xlabel('Time (s)')
legend({'Trace', 'Onset', 'Offset'});
set(gca,'FontSize',15)

% now plot the orientations of the stimuli in the second subplot.  remember
% that you can reuse the 'on' timestamps.
ax2 = subplot(212); % second of the two subplots
plot(on,stim);    % plot the stimulus onsets vs the stimulus orientation values
title('Stimulus orienation') 
ylabel('Orientation (degrees)')
xlabel('Time (s)')
legend({'Trace'});
set(gca,'FontSize',15)

% Without zooming in, the plotted orientations are not very useful. But as
% is, zooming in will change the time axis only in one graph, making 
% comparisons difficult. We can use the 'linkaxes' function to yoke the time
% axes across subplots
linkaxes([ax1 ax2],'x') 


%% Find the action potentials and calculate firing rates
% one way to detect spikes is to find local maxima
% use the 'findpeaks' function with 'minPeakHeight' to index the spikes
% using -10 as a threshold 

[pks,locs]= findpeaks(intra,'minPeakHeight',-10);
spikeTimes = timestamp(locs); % get the spike times by indexing the timestamp vector with the locations of found peaks
%now plot the spike times on top of the intra trace to see how well we
%found the spikes
figure('Name','Intracellular trace and detected spikes'); % make a new figure
hold on; % allow ploting of multiple plots
plot(timestamp,intra);%plot the intracellular trace again
scatter(spikeTimes,intra(locs)); % plot the detected spiketimes on top of the trace using your 'spikeTimes' vector and the 'intra' value at the location of each spike
title('Stimulus orienation') 
ylabel('Orientation (degres)')
xlabel('Time (s)')
legend({'Trace', 'Spikes'});
set(gca,'FontSize',15)

% calculate the mean firing rate for the whole trace by using the length of
% the spikeTimes vector (i.e. number of spikes) and the length of the 
% 'intra' variable, taking into account the sampling rate

numSpikes = length(spikeTimes);
wholeTraceMeanFr = numSpikes/(length(intra)/(2000));

% find the baseline firing rate from the initial portion of the trace
% see the 'find' function
baselineSpikes = find(spikeTimes<on(1)); % find all the spikes before the onset of the first stim
baselineNumSpikes = length(baselineSpikes); % find the length of baselineNum, this gives you the number of spikes
lastSpike = spikeTimes(baselineSpikes(end)); % find the time of the last spike in this series
firstSpike = spikeTimes(baselineSpikes(1)); % find the time of the first spike in this series
totalTime = (lastSpike-firstSpike);  %find the number of seconds represented by that interval (take into account the sampling rate of the on/off timestamps [hint: 1,000,000])
baseRate = baselineNumSpikes/totalTime;         % find the mean firing rate


fr = zeros([length(on),1]); % preallocate fr variable
for i =1:length(on) % 1 line challenge! See if you can calculate the firing rate for each stimulus presentation with only one line inside the for loop
    %fr(i)=sum(spikeTimes>on(i) & spikeTimes<off(i))/((off(i)-on(i))/(1000000)); 
    presentationSpikes = spikeTimes>on(i) & spikeTimes<off(i); % get spikes that fall into the presentation window (greater than on(i) and less than off(i), see '&' operator)
    currNumSpikes = sum(presentationSpikes); % count current number of spikes
    presentationTime = off(i)-on(i); % find the duration of current presentation, don't forget sampling rate!
    fr(i) = currNumSpikes /presentationTime; % divide num spikes by time to get firing rate
end


%% How well tuned is this neuron?
% Scatter the firing rates against the orientations to visualize the
% variability of neural responses across trials.  
figure('Name','Stimulus orientation vs Firing rate'); % make a new figure
scatter(stim,fr); %scatter the stimulus orientations against the firing rates
title('Stimulus orientation vs Firing rate') 
xlabel('Orientation (degres)')
ylabel('Firing Rate (Spikes/sec)')
set(gca,'FontSize',15)
%%
% Calculate the mean and std of firing rates for each orientation


uStim = unique(stim); % get unique orientation values
meanFr = zeros([length(uStim),1]); % preallocate for mean firing rate variable
stdFr = zeros([length(uStim),1]); % preallocate for std firing rate variable

for i = 1:length(uStim) % iterate through each unique grading value
    currInd = stim==uStim(i); % get indicies that correspond to current grading
    currFr = fr(currInd); % get current firing rates of interest by indexing 'fr'
    meanFr(i) = mean(currFr); % take the mean
    stdFr(i) = std(currFr); % gake the std
    %concise
    %meanFr(i) = mean(fr(stim==uStim(i))); % take mean of firing rates corresponding to current grading
    %stdFr(i) = std(fr(stim==uStim(i))); % take std of firing rates
end

% What is the neuron's preferred orientation?
% (see 'errorbar' function)
figure('Name','Stimulus orientation vs Firing rate'); % make a new figure
errorbar(uStim,meanFr,stdFr/sqrt(sum(stim==uStim(i)))) % plot the orientations vs the mean firing rates, adding the SEM error bars
title('Stimulus orientation vs Firing rate') 
xlabel('Orientation (degres)')
ylabel('Firing Rate (Spikes/sec)')
set(gca,'FontSize',15)


[maxFr,maxIndex]=max(meanFr) % use max() to get the maximum fr and the index it occured 
prefStim = uStim(maxIndex) % use the maxIndex to get the optimum orientation from uStim

%% %% Part 2: Same task, Vm version
% What about membrane potential tuning?
% Using the Vm vector, calculate the mean membrane potential for each
% stimulus presentation.  You can use code structure similar to that above, only a few lines need to change.

% Calculate the mean Vm response for each orientation

% Plot the Vm tuning curve on the same plot as the spike tuning curve.

trialVmMean = zeros([length(on),1]);
trialVmStd = zeros([length(on),1]);
for i =1:length(on)
    trialVmMean(i) = mean(intra(timestamp>on(i) & timestamp<off(i)));
end

meanMeanVm = zeros([length(uStim),1]); % preallocate for mean firing rate variable
meanStdVm  = zeros([length(uStim),1]); 
for i =1:length(uStim) % iterate through each unique gradient value
    meanMeanVm(i) = mean(trialVmMean(stim==uStim(i))); % take mean of firing rates corresponding to current gradient
    stdMeanVm(i) = std(trialVmMean(stim==uStim(i)));
end
figure('Name','Stimulus orientation vs membrane potential'); % make a new figure
errorbar(uStim,meanMeanVm,stdMeanVm/sqrt(sum(stim==uStim(i))))
title('Stimulus orientation vs membrane potential') 
xlabel('Orientation (degres)')
ylabel('Membrane potential (mV)')
set(gca,'FontSize',15)


% Plot the tuning curves as a polar plot.
% hint: to make it pretty, append the first value to the end to close the
% circle
figure('Name','Stimulus orientation vs Firing rate'); % make a new figure
ppUStim = [uStim uStim(1)];
ppMeanFr = [meanFr; meanFr(1)];
polarplot(ppUStim/180*pi,ppMeanFr)
title('Stimulus orientation vs Firing rate') 
set(gca,'FontSize',15)

%% Polar plot of spike and Vm tuning on same plot
% scale by mean for visualization
figure('Name','Stimulus orientation vs Firing Rate and Membrane Potential'); % make a new figure
ppUStim = [uStim uStim(1)];
ppMeanMeanVm = [meanMeanVm; meanMeanVm(1)];
polarplot(ppUStim/180*pi,ppMeanFr/mean(ppMeanFr))
hold on
polarplot(ppUStim/180*pi,ppMeanMeanVm/mean(ppMeanMeanVm))
title('Stimulus orientation vs Firing Rate and Membrane Potential') 
legend({'Firing rate', 'Membrane potential'});
set(gca,'FontSize',15)

% what do the different shapes tell you about underlying computations?

%% Average Vm to trial onset
figure('Name','Average Membrane Potential to Trail Onset'); % make a new figure
trialVmMean = zeros([length(on),1]);
trialVmStd = zeros([length(on),1]);

indVm = zeros([length(on),5800]);
for i =1:length(on)
    seg = intra(timestamp>on(i) & timestamp<off(i));
    indVm(i,:) = seg(1:5800);
end

meanIndVm = zeros([length(uStim),5800]); % preallocate for mean firing rate variable
for i =1:length(uStim) % iterate through each unique gradient value
    meanIndVm(i,:) = mean(indVm(stim==uStim(i),:)); % take mean of firing rates corresponding to current gradient
end

figure('Name','Mean Stim Triggered Vm')

for i = 1:8
    plotKey = ['42' num2str(i)];
    subplot(plotKey)
    plot(meanIndVm(i,:));
    title(uStim(i))
end

  
