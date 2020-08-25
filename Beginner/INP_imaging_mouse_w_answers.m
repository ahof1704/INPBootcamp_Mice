%%
%     INP data analysis bootcamp, Yale University
%     Imaging Section: Mouse Data
%
%     8/27/2019 KAF
%     8/24/2020 CJB
%     In each section, use the suggestions to fill in the variables and the
%     rest of the necessary code.

%     Data for this task are traces of the change in Ca2+ fluorescence from 
%     layer 2/3 somatostatin-positive (SOM+) GABAergic interneurons 
%     in V1 cortex of an awake mouse. The file 'sampleCRFdata.mat' contain the 
%     change in fluorescence of individual cells during the presentation of 
%     visual stimuli of varying contrasts. The file 'sampleRFdata.mat' contains
%     responses during the presentation of visual stimuli of varying sizes. 

%     We will start by looking at the data. For each cell, we will 
%     plot a heatmap to visualize the responses to visual stimuli.  Then we
%     will create error bar plots to see if the cells respond differently
%     to distinct contrasts/sizes. 

%     There are some exercises at the end.  Don't do the 'bonus challenges'
%     until you finish the exercises. 

%     Each .mat file will load the following variables:  
%           cellData: a matrix, the change in Ca2+ fluorescence (time points x neurons)
%           stimValue: a vector, the values of the various stimuli 
%           time: a vector, time for the fluorescence trace
%           visOn: a vector, the indices of the visual stimuli onset 
%           visOff: a vector, the indices of the visual stimuli offset 

%% Let's start with the contrast data, 'sampleCRFdata.mat'. 

%  Load the data here

load('/home/katie/Downloads/sampleCRFdata.mat'); 


stimType = 'Contrast (%)';          %  change to 'Size (degrees)' for sampleRFdata.mat

%  Plot your favorite 3 neurons in subplots, link the time axis (linkaxes).  What do you
%  notice about these cells? 
subplot(3,1,1)
plot(time, cellData(:,1)); 
hold on; 
subplot(3,1,2)
plot(time, cellData(:,2)); 
hold on; 
subplot(3,1,3)
plot(time, cellData(:,3)); 


%  Now let's plot a neuron with the times of the visual stimulus as '*'
%  using scatter.  

nrnNum       = 4;                         % pick your favorite neuron to plot first
visOnT       = time(visOn);               % use the vector 'time' to find the time of your visual stim onset  
visOffT      = time(visOff);              % the time of your visual stim offset    


figure; 
plot(time, cellData(:,nrnNum))                                % plot the cell activity
hold on; 
scatter(visOnT    , 10*ones(size(visOnT)),'*g');  %  plot onset times in green
scatter(visOffT   , 10*ones(size(visOffT)),'*r');  %  plot offset times in red


%%  Let's make our heatmap of neuron activity over each visual stimulation trial, with average activity at bottom

% first set some parameters

numCells     =   size(cellData,2);                          % number of neurons

% Let's choose a time window around the visual stimulus onset so that we can
% align our traces.   We'll start with this:
wBeforeT    = 2;                         % in seconds 
wAfterT     = 5;                         % in seconds

% sampling frequency (Hz)
fs           =  length(time)./(time(end)-time(1));     % we know the total amount of time, and the number of points.  What's our sampling frequency?

% given wBeforeT and wAfterT in seconds, how many points will we take before and after stimulus onset?  
% keep in mind that we can't have a fraction of a data point...     (see function 'round')
wBefore      = round(wBeforeT.*fs); 
wAfter       = round(wAfterT.*fs); 

winSz        = wBefore + wAfter +1 ;          % total number of points in our window

timeTrial    = linspace( -wBeforeT , wAfterT  , winSz );          % create a time vector from wBeforeT to wAfterT

nrnNum       = 4;                         % pick your favorite neuron to plot first

% get the data from our chosen neuronNum out of our cellData variable 
data         = cellData(:,nrnNum); 

%% now we will build a matrix 'visResp'. Each row will hold the fluorescence traces in a window around a visual stimuli
%  There will be as many rows as there are visual stimulus presentations

visResp      = zeros(length(visOn), length(timeTrial));       % initialize the matrix using the 'zeros' function.  What will the size of this matrix be? 

% we will align this matrix based on the visual stim onset (i.e. we don't need offset here).  
% loop through the visual stimuli and add the window of activity to the matrix in a new row
for istim = 1:length(visOn) 
    
    % using the stimulus onset (visOn) for that trial,  what is the index of
    % the start/end of the window?
    wStart = visOn(istim) - wBefore; 
    wEnd   = visOn(istim) + wAfter; 
    
    % add data for that window into a new row
    visResp(istim,:) = data(wStart:wEnd);
    
end

h = figure;

%% create the heatmap in the large top subplot
subplot(5,1,1:4)                    %  we're creating 5 x 1 subplots, but using all the first 4 for this heatmap

imagesc(timeTrial, 1:length(visOn), visResp    )                      % create a heatmap with our new matrix! 

caxis([0 3])                        % color axis limits
hold on;
plot( [0 0]  , [1, length(visOn)]   ,'--w','linew',2)      % mark visual stimulation onset with a white dashed line
title(['Neuron ', num2str(nrnNum) ])              % add the neuron number to your title using 'neuronNum'
ylabel('Trial','fontsize',12)

%% find the average response to stimuli over every trial.  Check that you're averaging correctly by confirming the size is right
meanResp = mean(visResp,1); 

%  plot the average response in the the bottom subplot
subplot(5,1,5)
plot(timeTrial,meanResp,'k')
xlabel('Time (s)', 'fontsize',12)
ylabel({'Average ' ; '\Delta F/F_0'}, 'fontsize',12)
xlim([timeTrial(1), timeTrial(end)])

outputFigName = ['Neuron_' , num2str(nrnNum),'_heatmap'];                       % fill in '?' for the correct neuron number.   

% Save the visResp data for each cell.  Make sure it is saving in your
% desired folder (i.e. specify your path)
save(['~/Downloads/Neuron_',num2str(nrnNum)], 'visResp'); 

% Save the plots as .fig, .jpg, .eps.  For .jpg , .eps, open them to make
% sure they've saved as you intended.  
fulloutputFn = fullfile('~/Downloads',outputFigName); 
savefig(fulloutputFn)
saveas(gcf, fulloutputFn, 'epsc')
saveas(gcf, fulloutputFn, 'jpeg')


