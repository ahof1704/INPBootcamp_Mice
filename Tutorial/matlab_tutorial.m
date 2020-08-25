% clears the existing variables in the workspace 
clear;

% will load as variable cellResps, because I gave it that name when I
% saved the data
load('Data/sampleData.mat');

%% what variables do we have in the workspace?
whos

%% how big is the data? 
help size
size(cellResps)   % a 1920 x 4 matrix (1920 timepoints, 4 cells)

%% plot the data from all the cells
figure;
plot(cellResps);

%% create a time vector and plot
fs = 30; % Hz
numpts = size(cellResps,1); 
time = linspace(0,numpts./fs,numpts);
figure,
plot(time, cellResps);

%% plot the data from the first cell
firstCellData = cellResps(:, 1);

% firstCellData is a 1920 x 1 vector (1920 time points, 1 cell)
size(firstCellData)
figure; plot(time, firstCellData);

%% on the same plot, plot the first half of the data
% subtract 1 so we can see the two data sets

halfLen = length(firstCellData)/2;

hold on;  % to plot of the same figure
plot(time(1:halfLen),firstCellData(1:halfLen)-1);


%% introducing the for loop: plot the data from cells 1, 2, 4
for i = [1 2 4]
   figure;
   plot(time, cellResps(:, i)); 
end

%% plot the data in the same figure
figure;
for i = [1 2 4]
   hold on;
   plot(time, cellResps(:, i)); 
end

%% store names for a legend in a "cell array"
legendLabelsFull = {'Neuron 1','Neuron 2','Neuron 3','Neuron 4'};  
legendLabelsFull{3}   % to access each cell, use curly braces
% this is annoying to write out for every cell. Instead, we make the legendLabels as we go. 

%% label the data
figure;
legendLabels = cell(3, 1);  % initialize 
j = 1;
for i = [1 2 4]
   hold on;
   plot(time, cellResps(:, i));
   legendLabels{j} = ['Neuron ' num2str(i)];   
   j = j+1;
end
legend(legendLabels);
xlabel('Time');
ylabel('Response');

%% we want to take the average of activity for each cell
help mean

%% call a built in matlab function on the data
dataMean = mean(cellResps, 'omitnan')
sizeDataMean = size(dataMean)
% This mean gives 4 values - must be the average across rows. 

%% But what if we want the mean trace (mean of 4 cells over all 1920 values)?

dataMean = mean(cellResps, 2, 'omitnan');
sizeDataMean = size(dataMean)  % Now we have a 1920 x 1 vector
figure; 
plot(time,dataMean)

%% What if we want to plot both the activity and the histogram showing the distribution of activity?  We can make our own plotting function

neuronNum = 2;   % cell we want to plot

plotData(time, cellResps(:,neuronNum), legendLabelsFull{neuronNum})

%% Exercises

% 1) create a loop to plot all 4 cells using your plotData function
% 2) change the color of the plots
% 3) instead of a loop, plot all 4 cells in subplots
% 4) what is the minimum activity in cell 1? 
% 5) what is the maximum activity of all cells?
% 6) find number of points above 1 for each cell


