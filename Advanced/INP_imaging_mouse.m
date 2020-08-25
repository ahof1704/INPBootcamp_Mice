%%
%     INP data analysis bootcamp, Yale University
%     Imaging Section: Mouse Data
%
%     8/27/2020 KAF-CJB-AF
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

clear
close all
%  Load the data here
load('../Data/sampleCRFdata_SOMpopln.mat'); %  Use sampleCRFdata (4 neurons) or sampleCRFdata_SOMpopln (31 neurons)
stimType = 'Contrast (%)';         %  change to 'Size (degrees)' for sampleRFdata.mat


%% Task1
%  Plot three random neurons in subplots. Requirements:
%  1) Each subplot should have the number or the neuron as its title.
%  2) Add labels for X and Y axes
%  3) link the time axis (linkaxes)
%  P.S.: What do you notice about these cells? 
%  P.S2: Design a script that would allow you to plot N neurons. So please,
%  don't just call the function "plot" three times :)

numCells     =   size(cellData,2);                          % number of neurons


%% Task2
%  Plot a neuron calcium trace and mark the time stamps of the visual
%  stimulus. Requirement:
%  1) Add title, labels to the axes and legend
%  2) Whatever type of marker/line you decide to use, make sure it has
%  distintc colors for visOn and visOff.

nrnNum       = 4;                         % pick your favorite neuron to plot first
visOnT       = time(visOn);               % use the vector 'time' to find the time of your visual stim onset  
visOffT      = time(visOff);              % the time of your visual stim offset    


%%  Let's make our heatmap of neuron activity over each visual stimulation trial, with average activity at bottom

% first set some parameters

numCells     =   size(cellData,2);                          % number of neurons

% Let's choose a time window around the visual stimulus onset so that we can
% align our traces.   We'll start with this:
wBeforeT    = 2;                         % in seconds 
wAfterT     = 5;                         % in seconds

% sampling frequency (Hz)
fs           =                            % we know the total amount of time, and the number of points.  What's our sampling frequency?

% given wBeforeT and wAfterT in seconds, how many points will we take before and after stimulus onset?  
% keep in mind that we can't have a fraction of a data point...     (see function 'round')
wBefore      = 
wAfter       = 

winSz        =                            % total number of points in our window

timeTrial    = linspace(  ,   ,  );          % create a time vector from wBeforeT to wAfterT

nrnNum       = 4;                         % pick your favorite neuron to plot first

% get the data from our chosen neuronNum out of our cellData variable 
data         = 

%% now we will build a matrix 'visResp'. Each row will hold the fluorescence traces in a window around a visual stimuli
%  There will be as many rows as there are visual stimulus presentations

visResp      =                            % initialize the matrix using the 'zeros' function.  What will the size of this matrix be? 

% we will align this matrix based on the visual stim onset (i.e. we don't need offset here).  
% loop through the visual stimuli and add the window of activity to the matrix in a new row
for i =  
    
    % using the stimulus onset (visOn) for that trial,  what is the index of
    % the start/end of the window?
    wStart = 
    wEnd   = 
    
    % add data for that window into a new row
    visResp(   ) = data(   );
    
end

h = figure;

%% create the heatmap in the large top subplot
subplot(5,1,1:4)                    %  we're creating 5 x 1 subplots, but using all the first 4 for this heatmap

imagesc(     )                      % create a heatmap with our new matrix! 

caxis([0 3])                        % color axis limits
hold on;
plot(   ,    ,'--w','linew',2)      % mark visual stimulation onset with a white dashed line
title(['Neuron ', ? ])              % add the neuron number to your title using 'neuronNum'
ylabel('Trial','fontsize',15)

%% Task 3 
% find the average response to stimuli over every trial.  Check that you're averaging correctly by confirming the size is right
meanResp = ?

%  plot the average response in the the bottom subplot
subplot(5,1,5)
plot(timeTrial,meanResp,'k')
xlabel('Time (s)')
ylabel({'Average ' ; '\Delta F/F_0'})
xlim([timeTrial(1), timeTrial(end)])
set(gca,'FontSize',15)

outputFigName = ['Neuron_' , num2str(nrnNum),'_heatmap'];                       % fill in '?' for the correct neuron number.   

% Save the visResp data for each cell.  Make sure it is saving in your
% desired folder (i.e. specify your path). Obs: create a new folder t osave
% your .mat and plots. Add a command to check if the folder already exist
% before creating the folder.
if ~exist('../plots', 'dir')
   mkdir('../plots')
end
save(['../plots/Neuron_',num2str(nrnNum)], 'visResp'); 

% Save the plots as .fig, .jpg, .eps.  For .jpg , .eps, open them to make
% sure they've saved as you intended.  
fulloutputFn = ?   % Full path to save the image 
savefig(fulloutputFn)
saveas(gcf, fulloutputFn, 'epsc')
saveas(gcf, fulloutputFn, 'jpeg')

%% How long after visual stim onset does the average activity (meanResp) peak?  
[~,mxIdx] = max(meanResp); 
pkTime = timeTrial(mxIdx); 


%% Extra exercises for the heatmap (optional) %% 

% (1) Change the window sizes.  Look at two cycles simultaneously. 
% (2) Change the proportion of subplots... 
% (3) Change the colormap to your favorite

% (4) Create the same plot for all neurons using a for loop
% (5) Check that you've saved all the plots and data for each iteration of the loop

% (6) Make the same plots for our "sampleRFdata".


%% Task 4
%  Now we will look at whether the cells respond differently to various contrasts/sizes
%  We will fill a vResp matrix (as our earlier visResp), but this time only
%  use values between visOn and visOff

vResp      =   zeros(length(visOn), max(visOff-visOn)+1);                    % initialize the matrix using the 'zeros' function.  What will the size of this matrix be? 

for i =  1 : length(visOn)
    
    % add data for that window into a new row
    vResp(i,:) = ;
    
end

%  We will use the average response across each 2s visual presentation
vRespAvg =


stimAxis =                  % Find the unique stim values. 
axLen    =                  % What is the length of stimAxis?

    
% now we will create vectors for the mean and SEM of visual
% stim values across every distinct presentation type

% initialize
meanData = 
semData =  

% find data for each distinct visual stim, average and take the SEM

for istim =          % iterate through all possible visual stim values       
    vdata       =                  % what is the vRespAvg for this iteration of stim value ?
    meanData(istim) =   
    semData(istim)  =  
end
    
        
% error bar plot of mean and SEM over observations within each feature
f = figure; 
errorbar(  , ,  , 's')                              
title(['Error bar plot for Neuron ', ? ])            % replace ?
xlabel(stimType)
ylabel('Response')
set(gca,'FontSize',15)



%% Task 5. Design a spike detector
%% Now things will get interesting. We want you to design a detector for the neuron's response to the stimulus. 
%% As a proxy for the GT, we will use the time stamps of visOn.
%% There are many ways to achieve this. Use your creativity. At the end of the bootcamp, 
%% we will check who evaluate who obtained the best recall rate.

%% Do you want more?

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