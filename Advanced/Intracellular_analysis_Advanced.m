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
fs_trace = 2e3;
fs_stimulus = 1000e3;

%%  What's in the data?
% Take a look at the intracellular trace.  Plot the intra against the timestamp vector. 

load(); % load data
plot(); % plot intra against timestamp
title('Intracellular trace') % Remember to always format your plots.
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

figure;
hold on;
plot(,);
scatter(,,);  %plot the onsets (preferably in green, check the scatter function help for symbol details)
scatter(,,);  %plot the offsets (preferably in red)

%%  Next, try this again using subplots.  
% Check the 'subplot' function help for details.
% replot the traces and points as above in the first of two subplots
ax1 =   % first of two subplots
hold on; % allow ploting of multiple plots
plot(,); % plot the intracellular trace again
scatter(,,);  %plot the onsets 
scatter(,,);  %plot the offsets the same way

% now plot the orientations of the stimuli in the second subplot.  remember
% that you can reuse the 'on' timestamps.
ax2 = ; % second of the two subplots
plot(,);    % plot the stimulus onsets vs the stimulus orientation values

% Without zooming in, the plotted orientations are not very useful. But as
% is, zooming in will change the time axis only in one graph, making 
% comparisons difficult. We can use the 'linkaxes' function to yoke the time
% axes for the subplots
linkaxes([ ],'') 


%% Find the action potentials and calculate firing rates
% one way to detect spikes is to find local maxima
% use the 'findpeaks' function with 'minPeakHeight' to index the spikes
% using -10 as a threshold 

[pks,locs]= findpeaks();
spikeTimes = ; % get the spike times by indexing the timestamp vector with the locations of found peaks
%now plot the spike times on top of the intra trace to see how well we
%found the spikes
figure; % make a new figure
hold on; % allow ploting of multiple plots
plot(,);%plot the intracellular trace again
scatter(,); % plot the detected spiketimes on top of the trace using your 'spikeTimes' vector and the 'intra' value at the location of each spike


% calculate the mean firing rate for the whole trace by using the length of
% the spikeTimes vector (i.e. number of spikes) and the length of the 
% 'intra' variable, taking into account the sampling rate

numSpikes = ;
wholeTraceMeanFr = ;

% find the baseline firing rate from the initial portion of the trace
% see the 'find' function
baselineSpikes = ; % find all the spikes before the onset of the first stim
baselineNumSpikes = ; % find the length of baselineNum, this gives you the number of spikes
lastSpike = ; % find the time of the last spike in this series
firstSpike = ; % find the time of the first spike in this series
totalTime = ;  %find the number of seconds represented by that interval (take into account the sampling rate of the on/off timestamps [hint: 1,000,000])
baseRate = ;         % find the mean firing rate


fr = ; % preallocate fr variable
for i = 
    presentationSpikes = ; % get spikes that fall into the presentation window (greater than on(i) and less than off(i), see '&' operator)
    currNumSpikes = ; % count current number of spikes
    presentationTime = ; % find the duration of current presentation, don't forget sampling rate!
    fr(i) = ; % divide num spikes by time to get firing rate
end


%% How well tuned is this neuron?
% Scatter the firing rates against the orientations to visualize the
% variability of neural responses across trials.  

scatter(,); %scatter the stimulus orientations against the firing rates

%%
% Calculate the mean and std of firing rates for each orientation


uStim = ; % get unique orientation values
meanFr = ; % preallocate for mean firing rate variable
stdFr = ; % preallocate for std firing rate variable

for i = 1:length(uStim) % iterate through each unique grading value
    currInd = ; % get indicies that correspond to current grading
    currFr = ; % get current firing rates of interest by indexing 'fr'
    meanFr(i) = ; % take the mean
    stdFr(i) = ; % gake the std
end

% What is the neuron's preferred orientation?
% (see 'errorbar' function)

errorbar(,) % plot the orientations vs the mean firing rates, adding the SEM error bars

[maxFr,maxIndex]= % use max() to get the maximum fr and the index it occured 
prefStim =  % use the maxIndex to get the optimum orientation from uStim



%% %% Part 2: Same task, Vm version
% What about membrane potential tuning?
% Using the Vm vector, calculate the mean membrane potential for each
% stimulus presentation.  You can use code structure similar to that above, only a few lines need to change.

% Calculate the mean Vm response for each orientation

% Plot the Vm tuning curve on the same plot as the spike tuning curve.

trialVmMean = ;
trialVmStd = ;
for i =1:length(on)
    trialVmMean(i) = ;
end

meanMeanVm = ; % preallocate 
stdMeanVm  = ; 
% can do like you did for spike section, but see if you can make it more
% concise. Ask TA for help!
for i =1: % iterate through each unique gradient value
    meanMeanVm(i) = ; % take mean of firing rates corresponding to current gradient
    stdMeanVm(i) = ;
end

errorbar(,,)



%% Plot the tuning curves as a polar plot. Don't forget to
% convert to radians!
% Hint: to make it pretty, append the first value to the end to close the
% circle

ppUStim = ;
ppMeanFr = ;
polarplot(,)

%% Polar plot of spike, Vm  and membrane fluctuation (std) tuning on same plot
% scale by mean for visualization
ppUStim = ;
ppVmFluct = ;
ppMeanMeanVm =;
polarplot(,)
hold on % must do after calling first polar plot. Otherwise it assumes you want a regular plot, leading to an error
polarplot(,)
polarplot(,)
% what do the different shapes tell you about underlying computations?

