
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


%%  What's in the data?
% Take a look at the intracellular trace.  Plot the intra against the timestamp vector. 

load('Cell1_data.mat'); % load data
plot(timestamp,intra); % plot intra against timestamp

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

figure;
hold on;
plot(timestamp,intra);
scatter(on,ones(size(on))*max(intra),'g');  %plot the onsets (preferably in green)
scatter(off,ones(size(on))*max(intra),'r');  %plot the offsets (preferably in red)

%%  Next, try this again using subplots.  
% Check the 'subplot' function help for details.
% replot the traces and points as above in the first of two subplots
ax1 = subplot(211);  % first of two subplots
hold on; % allow ploting of multiple plots
plot(timestamp,intra); % plot the intracellular trace again
scatter(on,ones(size(on))*max(intra),'g');  %plot the onsets (check the scatter function help for symbol details)
scatter(off,ones(size(on))*max(intra),'r');  %plot the offsets the same way

% now plot the orientations of the stimuli in the second subplot.  remember
% that you can reuse the 'on' timestamps.
ax2 = subplot(212); % second of the two subplots
plot(on,stim);    % plot the stimulus onsets vs the stimulus orientation values

% Without zooming in, the plotted orientations are not very useful. But as
% is, zooming in will change the time axis only in one graph, making 
% comparisons difficult. We can use the 'linkaxes' function to yoke the time
% axes for the subplots
linkaxes([ax1 ax2],'x') 


%% Find the action potentials and calculate firing rates
% one way to detect spikes is to find local maxima
% use the 'findpeaks' function with 'minPeakHeight' to index the spikes
% using -10 as a threshold 

[pks,locs]= findpeaks(intra,'minPeakHeight',-10);
spikeTimes = timestamp(locs); % get the spike times by indexing the timestamp vector with the locations of found peaks
%now plot the spike times on top of the intra trace to see how well we
%found the spikes
figure; % make a new figure
hold on; % allow ploting of multiple plots
plot(timestamp,intra);%plot the intracellular trace again
scatter(spikeTimes,intra(locs)); % plot the detected spiketimes on top of the trace using your 'spikeTimes' vector and the 'intra' value at the location of each spike


% calculate the mean firing rate for the whole trace by using the length of
% the spikeTimes vector (i.e. number of spikes) and the length of the 
% 'intra' variable, taking into account the sampling rate

numSpikes = length(spikeTimes);
wholeTraceMeanFr = numSpikes/(length(intra)/2000);

% find the baseline firing rate from the initial portion of the trace
% see the 'find' function
baselineSpikes = find(spikeTimes<on(1)); % find all the spikes before the onset of the first stim
baselineNumSpikes = length(baselineSpikes); % find the length of baselineNum, this gives you the number of spikes
lastSpike = spikeTimes(baselineSpikes(end)); % find the time of the last spike in this series
firstSpike = spikeTimes(baselineSpikes(1)); % find the time of the first spike in this series
totalTime = (lastSpike-firstSpike)/(1000000);  %find the number of seconds represented by that interval (take into account the sampling rate of the on/off timestamps [hint: 1,000,000])
baseRate = baselineNumSpikes/totalTime;         % find the mean firing rate


fr = zeros([length(on),1]); % preallocate fr variable
for i =1:length(on) % 1 line challenge! See if you can calculate the firing rate for each stimulus presentation with only one line inside the for loop
    %fr(i)=sum(spikeTimes>on(i) & spikeTimes<off(i))/((off(i)-on(i))/(1000000)); 
    presentationSpikes = spikeTimes>on(i) & spikeTimes<off(i); % get spikes that fall into the presentation window (greater than on(i) and less than off(i), see '&' operator)
    currNumSpikes = sum(presentationSpikes); % count current number of spikes
    presentationTime = ((off(i)-on(i))/(1000000)); % find the duration of current presentation, don't forget sampling rate!
    fr(i) = currNumSpikes /presentationTime; % divide num spikes by time to get firing rate
end


%% How well tuned is this neuron?
% Scatter the firing rates against the orientations to visualize the
% variability of neural responses across trials.  

scatter(stim,fr); %scatter the stimulus orientations against the firing rates

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

errorbar(uStim,meanFr,stdFr/sqrt(5)) % plot the orientations vs the mean firing rates, adding the SEM error bars

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
errorbar(uStim,meanMeanVm,stdMeanVm/sqrt(5))



% Bonus points: plot the tuning curves as a polar plot.
% hint: to make it pretty, append the first value to the end to close the
% circle

ppUStim = [uStim uStim(1)];
ppMeanFr = [meanFr; meanFr(1)];
polarplot(ppUStim/180*pi,ppMeanFr)

%% bonus^2 Polar plot of spike and Vm tuning on same plot
% scale by mean for visualization
ppUStim = [uStim uStim(1)];
ppMeanMeanVm = [meanMeanVm; meanMeanVm(1)];
polarplot(ppUStim/180*pi,ppMeanFr/mean(ppMeanFr))
hold on
polarplot(ppUStim/180*pi,ppMeanMeanVm/mean(ppMeanMeanVm))

% what do the different shapes tell you about underlying computations?

%% bonus^3 average Vm to trial onset

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
subplot(421)
plot(meanIndVm(1,:));
title(uStim(1))
subplot(422)
plot(meanIndVm(2,:));
title(uStim(2))
subplot(423)
plot(meanIndVm(3,:));
title(uStim(3))
subplot(424)
plot(meanIndVm(4,:));
title(uStim(4))
subplot(425)
plot(meanIndVm(5,:));
title(uStim(5))
subplot(426)
plot(meanIndVm(6,:));
title(uStim(6))
subplot(427)
plot(meanIndVm(7,:));
title(uStim(7))
subplot(428)
plot(meanIndVm(8,:));
title(uStim(8))

integ = zeros([8,1]);
for i = 1:8
    integ(i)=sum(meanIndVm(i,:)+55.5);
end
%% Part 3: OK, do it again for cell 2.

%% Part 4: extras for those with more coding experience
% Try measuring the spike threshold for each action potential.  
% Is there any apparent tuning of spike threshold?
% How about the rate of Vm change leading up to the action potentials (dV/dt)?
% Measure the Vm fluctuations (SD of the Vm) for each orientation.  Is the
% VmSD tuned?
% Try looking at the power spectra of the Vm trace, is there any tuning of
% specific frequencies?





