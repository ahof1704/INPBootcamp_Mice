% JAC 8-20-2019, updated 8-25-2019
% Modified by CJB-KAF-AF 8-27-2020
%
% In each section, use the suggestions to fill in the variables and the
% rest of the necessary code.

% Data for this task are an intracellular trace from layer 5 cortical
% neurons in V1 cortex of an awake mouse and several vectors of data about
% the visual stimuli that were presented during the recording.


%% Let's take a look
% load the data by using the load command 
%(or double clicking on the file in the 'current folder' panel of the GUI)
% and then plot the 'intra' variable

load('../Data/Neuron1_data.mat'); % load data
figure('Name','Intracellular Recording');  
plot(intra); % plot 'intra'

% feel free to zoom in to look at the trace!


%%  What happens when a stimulus is presented?
% Let's see what happens when a stimulus is presented. To do this, we will
% plot the 'intra' variable from the begining of the recording to the offset
% of the first stimulus. The 'timestamp' variable records the moment each
% membrane potential value was sampled. To see what indices occured
% before the first offset, we can use the '<' operator within the 'find' 
% function. Knowing the end times of the stimulus presentations are held in
% the 'off' variable, fill in the variable names to get the indices of the 
% 'intra' variable that occur before the first offset.

indices = find(timestamp < off(1));

figure('Name','Intracellular Recording Up to End of 1st Visual Stimulus');  
plot(intra(indices))

% looking at the figure, can you guess when the stimulus was presented?
% check by accessing the first value of the 'on' variable to check

on(1)

% Hmm.. something doesn't seem right. It looks like the stimulus is shown
% at around 6000 on the x axis, but the first 'on' value is ~3!
% If you only pass one vector to the plot function, it plots the index
% number on the x axis, but we want time so we can compare to the first
% stimulus onset. To do this, we can pass the 'timestamp' variable to the
% 'plot' function as well. This tells the function that
% (x,y) = (timestamp,intra). Note: x and y must be the same length!
% (Hint: index 'timestamp' the same way you index 'intra')

plot(timestamp(indices),intra(indices))

% Now compare the stimulus onset time to what we see in the graph. Look
% better?


%% Find the action potentials and calculate firing rates
% An advantage of electrophysiology is obtaining temporally precise spike
% times which provide valuable insignt into neural computation. One way 
% to detect spikes is to find local maxima using the 'findpeaks' function 
% with 'minPeakHeight' set to -10 as a threshold. The use of
% 'minPeakHeight' is an example of passing an option to a function.
% Typically, options are used to pass extra parameters to the function.
% Check the documentation for a function to see what options can be passed!
% The value for an option follows the name of the option, seperated by a comma.

% First, let's see how well we can find each action potential. Call
% findpeaks with the 'minPeakHeight' option without specifying any output 
% variables. Explore the plot to see how well we did! 
% (arrows indicate local maxima that satisfy the option we passed).

findpeaks(intra,'minPeakHeight',-10); 

% Next, we are going to want to record the location of th peaks for
% subsequent analysis. repeat what you did above with the specified output 
% variables 
[pks,locs]= findpeaks(intra,'minPeakHeight',-10);

% Now calculate the mean firing rate for the whole recording by using the 
% length of the 'locs' variable (i.e. number of spikes) and the length of 
% the 'intra' variable (i.e. recording time), taking into account the
% sampling rate (2000 Hz)

numSpikes = length(locs);
recordingMeanFr = numSpikes/(length(intra)/2000);

%% Firing rates per presentation
% Next, we will find the firing rate for each presentation using a for
% loop. 
% Before the for loop, we need to get the spike times and stimulus onset
% and offset times on common ground. To do this, index 'timestamp' with the
% locs variable you found above.

spikeTimes = timestamp(locs);

% Similar to before when we found the values of 'intra' that occured before
% the first offset time, we need to find the spikes that occured between
% the 'on' and 'off' times. We can accomplish this with the '&' operator.
% The operator will return true only if both values are also true. We can
% take advantage of this to only select those spikes that are > the 'on' time
% and < the 'off' time. Matlab is also nice in that it is easy to transfer 
% between boolean logic (true and false) and numbers (1 and 0). 
% By calling sum on a boolean vector, we can count the number of true 
% values, thus counting the number of spikes that fall between the 'on' and
% 'off' times. Lastly, we calculate the total duration of presentation, and
% then divide the number of spikes by this time to get the firing rate.


fr = []; % declare fr as an empty vector so we can dynamically add to it
for i =1:length(on) % for each presentation
    presentationSpikes = spikeTimes>on(i) & spikeTimes<off(i); % get spikes that fall into the presentation window (greater than on(i) and less than off(i), see '&' operator)
    currNumSpikes      = sum(presentationSpikes); % count current number of spikes
    presentationTime   = ((off(i)-on(i))); % find the duration of current presentation
    fr(end+1)          = currNumSpikes/presentationTime; % divide num spikes by time to get firing rate, and append to the end of your vector
end

% let's take a look at a histogram of the firing rates evoked by the
% stimuli using the histogram function
figure('Name','Histogram')
nbins = 15;    % the number of bins we'll use for the histogram
histogram(fr,nbins)
xlabel('Firing Rate (Hz)')
title('Histogram of Stimulus Evoked Firing Rates')

%% Tuning curves!
% Let's make a canonical tuning curve. To do this, we will need to know
% the orientations of gradients we presented. We can do this by calling the
% 'unique' function on 'stim'. This will allow us to go through each 
% stimulus categorically 

uStim = unique(stim);

% Like above, we will use a for loop. This time, it is to find the mean
% firing rate for each stimulus gradient. We will accomplish this by
% finding the trial indices that correspond to a given orientation of stimulus (say, 0
% degrees). Using the indices, we can grab the corresponding firing rates 
% and calculate their mean. 

meanFr = []; % declare meanFr as list so we can dynamically add to it

for i = 1:length(uStim)                % iterate through each unique gradient value (i.e. orientation)
    gradientInd = stim==uStim(i);      % get indicies that correspond to current gradient
    gradientFr = fr(gradientInd);      % get current firing rates of interest by indexing 'fr'
    meanFr(end+1) = mean(gradientFr);  % take the mean
end

%plot the firing rate against the stimulus presentations 
figure('Name','Tuning Curve')
plot(uStim,meanFr)
xlabel('Orientation (degrees)'); 
ylabel('Mean Firing Rate (Hz)'); 
title('Orientation Tuning Curve');
