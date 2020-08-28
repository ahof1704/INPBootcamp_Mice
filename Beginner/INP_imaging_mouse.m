%%
%     INP data analysis bootcamp, Yale University
%     Imaging Section: Mouse Data
%
%     8/27/2019 KAF
%     8/24/2020 CJB, KAF, AOF
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
%     plot a the cell activity.  Then we will create error bar plots (tuning curves) 
%     to see if the cells respond differently to distinct contrasts/sizes. 

%     There are some exercises at the end.  Don't do the 'bonus challenges'
%     until you finish the exercises. 

%     Each .mat file will load the following variables:  
%           cellData: a matrix, the change in Ca2+ fluorescence (time points x neurons)
%           stimValue: a vector, the values of the various stimuli 
%           time: a vector, time for the fluorescence trace
%           visOn: a vector, the indices of the visual stimuli onset 
%           visOff: a vector, the indices of the visual stimuli offset 

%% Let's start with the 'sampleRFdata.mat'.   Complete the code as you go


%  Load the data here

load(  );   % add the path to the 'sampleRFdata.mat' file.


stimType = 'Size (degrees)';          %  change to 'Size (degrees)' for sampleRFdata.mat

%  Plot your favorite 3 neurons in subplots, link the time axis (linkaxes).  What do you
%  notice about these cells? 

figure('Name','My favorite neurons'); 
ax1 = subplot(3,1, );
plot( , cellData( , ));   % plot time vs cell activity for your first chosen neuron 
ylabel('\Delta F/F_0')
hold on; 
ax2 = subplot(3,1,  );
plot(, );     % plot time vs cell activity for your second chosen neuron 
ylabel('\Delta F/F_0')
hold on; 
ax3 = subplot(3,1,  );
plot(, );     % plot time vs cell activity for your second chosen neuron 
ylabel('\Delta F/F_0')
xlabel('Time (s)')
linkaxes([ax1 ax2 ax3])




%  Now let's plot a neuron with the times of the visual stimulus as '*'
%  using scatter.  

nrnNum       =  4;                         % pick your favorite neuron to plot first, or choose the given 4th neuron
visOnT       =                            % use the vector 'time' to find the time of your visual stim onset  
visOffT      =                            % the time of your visual stim offset    


% plot the cell activity.  Zoom in and examine the data
figure('Name','Cell activity'); 
plot(    )                                % plot the cell activity over time
hold on; 
scatter(    , 10*ones(size(visOnT)),'*g');  %  plot onset times in green
scatter(    , 10*ones(size(visOffT)),'*r');  %  plot offset times in red
ylabel('\Delta F/F_0')
xlabel('Time (s)')
legend({'Ca trace', 'VisOn', 'VisOff'});
set(gca,'FontSize',15)


%% Now we will look at whether the cells respond differently to various stimulus sizes
%  We will fill a vResp matrix with the fluorescence values between visOn and visOff

nrnNum       = 4;                         % pick your favorite neuron to plot first

% get the data from our chosen nrnNum out of our cellData variable 
data         = cellData( , ); 

% initialize the matrix using the 'zeros' function.  What will the size of
% this matrix be?   We want each row to be a visual stimulus trial, and
% each column to be the fluorescence between the visOn and visOff indices
% for that trial.  
vResp      =   zeros(, );             

for istim =  1 : length(visOn)
    
    % add data for that window into a new row
    vResp(istim,:) = data(visOn(istim):visOff(istim));
    
end

%  Now take the average vResp so we get the average response of each 2s
%  visual presentation for each trial
vRespAvg = mean( , );


stimAxis =  ;            % Find the unique stim values.  Hint: Check out the function 'unique'
axLen    =  ;               % What is the length of stimAxis?  
    
% now we will create vectors for the mean and SEM of visual
% stimuls values across every distinct presentation type

% initialize
meanData = zeros(size(stimAxis));    % do you know why we initialize with zeros of this size? 
semData =  zeros(size(stimAxis)); 

% find data for each distinct visual stim, average and take SEM

for istim = 1:                                   % iterate through all possible visual stim values       
    vdata       =  vRespAvg(  );                 % what are the values of the vRespAvg for this stim value ?
    meanData(istim) =  ; 
    n = ;                                        % how many trials do we have for this stimulus? 
    semData(istim)  =  ;                         % recall the SEM is the standard deviation divided by the square root of the n.  see functions 'std' and 'sqrt'
end
    
        
% Tuning curve (with an error bar plot) using the mean and SEM over
% distinct stimulus sizes
f = figure; 
errorbar(  , ,  , 's')                              
title(['Size tuning curve for Neuron ', ? ])           % replace ?
xlabel(stimType, 'fontsize',12)
ylabel('\Delta F/F_0')
set(gca,'FontSize',15)


outputFigName = ['Neuron_' , ? ,'_size_tuning_curve'];                       % fill in '?' for the correct neuron number.   

% Save the data visResp , stimAxis, meanData, and semData.  Make sure it is saving in your
% desired folder (i.e. specify your path)
save(['../Data/Neuron_',num2str(nrnNum)], 'vResp', 'meanData','semData');      % change '../Data/' to your own path to save the data

% Save the plots as .fig, .jpg, .eps.  For .jpg , .eps, open them to make
% sure they've saved as you intended.  
fulloutputFn = ['../plots/', ? ];         % save figures in a "plots" folder.   What does all this mean?
if ~exist('../plots','dir')
    mkdir('../plots');
end

savefig(  )             % use the file name you just generated .   what does this save?
saveas(gcf,  ,  )       % what does gcf mean?   How do you save an .eps figure?  
saveas(gcf,   ,  )      % save your .jpg figure
