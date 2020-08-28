%%
%     INP data analysis bootcamp, Yale University
%     Imaging Section: Mouse Data
%
%     8/27/2019 KAF
%     8/24/2020 CJB, KAF, AOF
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

load('../Data/sampleRFdata.mat'); 


stimType = 'Size (degrees)';          %  change to 'Size (degrees)' for sampleRFdata.mat

%  Plot your favorite 3 neurons in subplots, link the time axis (linkaxes).  What do you
%  notice about these cells? 
figure('Name','My favorite neurons'); 
ax1 = subplot(3,1,1);
plot(time, cellData(:,1));
ylabel('\Delta F/F_0')
hold on; 
ax2 = subplot(3,1,2);
plot(time, cellData(:,2)); 
ylabel('\Delta F/F_0')
hold on; 
ax3 = subplot(3,1,3);
plot(time, cellData(:,3)); 
ylabel('\Delta F/F_0')
xlabel('Time (s)')
linkaxes([ax1 ax2 ax3])

%  Now let's plot a neuron with the times of the visual stimulus as '*'
%  using scatter.  

nrnNum       = 4;                         % pick your favorite neuron to plot first
visOnT       = time(visOn);               % use the vector 'time' to find the time of your visual stim onset  
visOffT      = time(visOff);              % the time of your visual stim offset    


% plot the cell activity.  Zoom in and examine the data
figure('Name','Cell activity'); 
plot(time, cellData(:,nrnNum))                                
hold on; 
scatter(visOnT    , 10*ones(size(visOnT)),'*g');  %  plot onset times in green
scatter(visOffT   , 10*ones(size(visOffT)),'*r');  %  plot offset times in red
ylabel('\Delta F/F_0')
xlabel('Time (s)')
legend({'Ca trace', 'VisOn', 'VisOff'});
set(gca,'FontSize',15)


%% Now we will look at whether the cells respond differently to various stimulus sizes
%  We will fill a vResp matrix, but this time only use values between visOn and visOff

nrnNum       = 4;                         % pick your favorite neuron to plot first

% get the data from our chosen neuronNum out of our cellData variable 
data         = cellData(:,nrnNum); 


vResp      =   zeros(length(visOn), max(visOff-visOn)+1);               % initialize the matrix using the 'zeros' function.  What will the size of this matrix be? 

for istim =  1 : length(visOn)
    
    % add data for that window into a new row
    vResp(istim,:) = data(visOn(istim):visOff(istim));
    
end

%  We will use the average response across each 2s visual presentation
vRespAvg = mean(vResp,2);


stimAxis =  unique(stimValue);            % Find the unique stim values. 
axLen    =   length(stimAxis);               % What is the length of stimAxis?
    
% now we will create vectors for the mean and SEM of visual
% stim values across every distinct presentation type

% initialize
meanData = zeros(size(stimAxis)); 
semData =  zeros(size(stimAxis)); 

% find data for each distinct visual stim, average and take SEM

for istim = 1:length(stimAxis)          % iterate through all possible visual stim values       
    vdata       =  vRespAvg(stimValue == stimAxis(istim));                 % what is the vRespAvg for this iteration of stim value ?
    meanData(istim) =  mean(vdata); 
    n = sum(stimValue == stimAxis(istim));    % how many trials do we have for this stimulus? 
    semData(istim)  =  std(vdata)/sqrt(n);  
end
    
        
% error bar plot of mean and standard deviation over observations within each feature
f = figure; 
errorbar( stimAxis , meanData, semData , 's')                              
title(['Size tuning curve for Neuron ', num2str(nrnNum) ])           % replace ?
xlabel(stimType, 'fontsize',12)
ylabel('\Delta F/F_0')
set(gca,'FontSize',15)

outputFigName = ['Neuron_' , num2str(nrnNum),'_size_tuning_curve'];                       % fill in '?' for the correct neuron number.   

% Save the data visResp , stimAxis, meanData, and semData.  Make sure it is saving in your
% desired folder (i.e. specify your path)
save(['../Data/Neuron_',num2str(nrnNum)], 'vResp', 'meanData','semData');               % change '../Data/' to your own path to save the data

% Save the plots as .fig, .jpg, .eps.  For .jpg , .eps, open them to make
% sure they've saved as you intended.  
fulloutputFn = fullfile('../plots/',outputFigName);    % or ['../plots/',outputFigName]    % save figures in a "plots" folder.   What does all this mean?
if ~exist('../plots','dir')
    mkdir('../plots');
end

savefig(fulloutputFn)                       % use the file name you just generated .   what does this save?
saveas(gcf, fulloutputFn, 'epsc')           % what does gcf mean?   How do you save an .eps figure?  
saveas(gcf, fulloutputFn, 'jpeg')           % save your .jpg figure

