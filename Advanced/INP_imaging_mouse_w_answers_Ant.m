%%
%     INP data analysis bootcamp, Yale University
%     Imaging Section: Mouse Data
%
%     8/27/2019 KAF
%     
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

clear all
close all
%  Load the data here
load('../Data/sampleCRFdata.mat'); 
stimType = 'Contrast (%)';          %  change to 'Size (degrees)' for sampleRFdata.mat

%  Plot your favorite 3 neurons in subplots, link the time axis (linkaxes).  What do you
%  notice about these cells? 
idxNeuron = [1 2 4]; % List of neurons
figure('Name','My favority neurons'); 
for i = 1:length(idxNeuron)
    subplot(3,1,i)
    plot(time, cellData(:,idxNeuron(i))); 
    title(['Neuron ' num2str(idxNeuron(i))]); 
end


%  Now let's plot a neuron with the times of the visual stimulus as '*'
%  using scatter.  

nrnNum       = 4;                         % pick your favorite neuron to plot first
visOnT       = time(visOn);               % use the vector 'time' to find the time of your visual stim onset  
visOffT      = time(visOff);              % the time of your visual stim offset    


figure('Name','Cell activity'); 
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
if ~exist('../plots', 'dir')
   mkdir('../plots')
end
save(['../plots/Neuron_',num2str(nrnNum)], 'visResp'); 

% Save the plots as .fig, .jpg, .eps.  For .jpg , .eps, open them to make
% sure they've saved as you intended.  
fulloutputFn = fullfile('../plots',outputFigName); 
savefig(fulloutputFn)
saveas(gcf, fulloutputFn, 'epsc')
saveas(gcf, fulloutputFn, 'jpeg')

%% How long after visual stim onset does the average activity (meanResp) peak?  
[~,mxIdx] = max(meanResp); 
pkTime = timeTrial(mxIdx); 


%% Exercises for the heatmap %% 

% (1) Change the window sizes.  Look at two cycles simultaneously. 
% (2) Change the proportion of subplots... 
% (3) Change the colormap to your favorite

% (4) Create the same plot for all neurons using a for loop
% (5) Check that you've saved all the plots and data for each iteration of the loop

% (6) Make the same plots for our "sampleRFdata".


%% Now we will look at whether the cells respond differently to various contrasts/sizes
%  We will fill a vResp matrix (as our earlier visResp), but this time only
%  use values between visOn and visOff

vResp      =   zeros(length(visOn), max(visOff-visOn)+1);                    % initialize the matrix using the 'zeros' function.  What will the size of this matrix be? 

for istim =  1 : length(visOn)
    
    % add data for that window into a new row
    vResp(istim,:) = data(visOn(istim):visOff(istim));
    
end

%  We will use the average response across each 2s visual presentation
vRespAvg = mean(vResp,2);


stimAxis =  unique(stimValue);            % Find the unique stim values. 
axLen    =   length(stimAxis);               % What is the length of stimAxis?
    
% now we will create vectors for the mean and standard deviation of visual
% stim values across every distinct presentation type

% initialize
meanData = zeros(size(stimAxis)); 
stdData =  zeros(size(stimAxis)); 

% find data for each distinct visual stim, average and take standard deviation

for istim = 1:length(stimAxis)          % iterate through all possible visual stim values       
    vdata       =  vRespAvg(stimValue == stimAxis(istim));                 % what is the vRespAvg for this iteration of stim value ?
    meanData(istim) =  mean(vdata); 
    stdData(istim)  =  std(vdata); 
end
    
        
% error bar plot of mean and standard deviation over observations within each feature
f = figure; 
errorbar( stimAxis , meanData, stdData , 's')                              
title(['Error bar plot for Neuron ', num2str(nrnNum) ])           % replace ?
xlabel(stimType, 'fontsize',12)



%% Antonio's section
% Design a simple spike detector
figure('Name','Spike detection'); 
cleanData = smoothdata(cellData(:,1),1,'movmean',20); %Smooth signal (play with the different methods or implement your own)
yd = diff(cleanData(:,1))./diff(time);
xd = (time(2:end)+time(1:(end-1)))/2;
plot(time,cellData(:,1),time,cleanData,xd,yd)
hold on; 
[pks, visOnDetec] = findpeaks(yd, xd,'MinPeakDistance',2,'MinPeakHeight',max(meanResp)); % Use MinPeakHeight to restrict the type of peak you want (make it a fucntion of the meanResp. We are assuming 2s of distance)
scatter(visOnDetec    , 10*ones(size(visOnDetec)),'+m');  %  plot onset times in green
scatter(visOnT  , 10*ones(size(visOnT)),'*g');  %  plot onset times in green
scatter(visOffT , 10*ones(size(visOffT)),'*r');  %  plot offset times in red
lgd = legend({'Original Data','Smoothed Data','Derivative','Detected Peak', 'VisOn', 'VisOff'},'FontSize',14);

% Compute accuracy of the spike detector
PerformanceSpikeDetec = zeros(length(visOnT),4);
for i = 1:length(visOnT)
    [min_idx, min_idx] = min(abs(visOnT(i) - visOnDetec));
    PerformanceSpikeDetec(i,1:3) = [visOnT(i), locs(min_idx) visOnT(i)-visOnDetec(min_idx)];
    if abs(visOnT(i)-visOnDetec(min_idx)) < 0.8
        PerformanceSpikeDetec(i,4) = 1;
    end
end
TablePerformance = array2table(PerformanceSpikeDetec,...
    'VariableNames',{'visOnT (s)','visOnDetec (s)','Error (s)', 'Detec'});

title(lgd, ['Recall: ' num2str(sum(PerformanceSpikeDetec(:,4)/length(visOnT)))])
xlabel('Time (s)')
ylabel({'Average ' ; '\Delta F/F_0'}, 'fontsize',12)
set(gca,'FontSize',15)


%% Exercises

% (1) Do this for all cells.  Save your plots programatically as "sampleCRFdata_vResp_neuronX.fig", where X = 1,2,3, ... for the number of the neuron you've plotted
% (2) Save your mean and standard deviation data as "sampleCRFdata_vResp_neuronX.mat"
% (3) Repeat this for sampleRFdata
% (4) In the error bar plot, change the squares to filled red circles. Play with the colors and shapes to make it look as you like! 
% (5) For the receptive field size plots, find the stimulus size with the maximum mean 
% (6) Make the error bars show the standard error of the mean instead of standard deviation
% (7) Create a function for your error bar plot. (Hint: look at the plotData function we used yesterday)
% (8) Find the overall average across cells for each stimulus value and plot using your function



%% Additional Exercises for the heatmap %% 

% (1) What could be a problem with our windows winStart and winEnd in the for loop?  How could we fix it?  Hint: think of the edges
% (2) Compile all the data into one big 3D matrix (trial number x window x cell number)
%           Bonus Challenge: Do it without a for loop!  (Hint: check out repmat)
% (3) Plot the average across all cells.
% (4) Make a new plotting function 'makeHeatmap' to make your
%       heatmap/average trace figure.  



%% If you want more of a challenge, ask about the wheel data ! 









