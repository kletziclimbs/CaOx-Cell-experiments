%% --- L/D quantification of dead cells INS-6015---
%02/10/2018

%quantifies number and area of dead cells - two output files - mean/image
%and mean values per sample (3 images per sample)

%standard values are given for 6 samples, 3 wells/samples, 5 images/well -
%otherwise values need to be adapted accordingly

%input images: red fluorescence images of cells with dead staining (Invitrogen)
%  20x objective LeicaCTR600

%to do:
%create image datastore with all raw tiff files of red (dead cells) channel
%rename output files accordingly

%% ---INITIALIZATION---
clear variables
close all
clc

%% ---create image datastore---
imds = imageDatastore('D:\lab\Experiments\Cell based assays\Live Dead assays\20181001_RPTEC_LD_INS6015_Exp04-4\rawImages\LD_Exp04-4\Ch01', ...
 'FileExtensions',{'.tif'});

%% ---make binary image; count, area---
numberFiles = length(imds.Files);
countData = ones(numberFiles,1);
areaData = ones(numberFiles,1);


for i = 1:length(imds.Files); 
    I = readimage(imds,i);
    sigma = 0.3; % amplitude of edges
    alpha = 2.5; % smoothing of details, typical values are between 0.01 - 10; <1 - increases details; >1smooths details
    I_filter = locallapfilt(I,sigma,alpha);
    bw = im2bw(I_filter , 0.03);%greyscale to binary image

    totalPx = numel(bw); %total pixel number of image
    numWhitePx = sum(bw(:)); %calculate total number of white pixels
    percentWhitePx = numWhitePx/ (totalPx.*0.01); %calculate percentage of area of white pixel - COM adhered


    D = -bwdist(~bw); % distance transform of binary image
    mask = imextendedmin(D,1); % extended minima transform of h-minima transform (suppresses all minima <h)
    
    D2 = imimposemin(D,mask);% % modifies intensity image D; only regional minima where mask is not zero
    Ld2 = watershed(D2); % shape based watershed segmentation to split touching COM
 
    bw(Ld2 == 0) = 0; %overlay shape based watershed on binary image
    
    [B,L,N] = bwboundaries(bw,'noholes');
    
   
    countData(i,1) = N;
    areaData(i,1) = percentWhitePx;
end
%[name] = fileparts(imds.Files(i));
results = table(imds.Files, countData, areaData, 'VariableNames',{ 'filename', 'count', 'area'});

%% ---statistics per well---
R = 5; % specify number of replicate images per well to average over
wells = 3; % number of wells/condition
conditions = 6; % specify number of conditions
results2 = table2array(results(:,2:3));
ResultsStats = zeros(conditions*wells, 4); % stats array for mean values + SD count/area

ResultsStats(1,1) = mean(results2(1:R,1));%calculate mean/ SD of count, column 1, mean/well
ResultsStats(1,2) = std(results2(1:R,1));
ResultsStats(1,3) = mean(results2(1:R,2));%calculate mean/ SD of area, column 2
ResultsStats(1,4) = std(results2(1:R,2));

for i = 1:((conditions*wells)-1) 
    
    ResultsStats(i+1,1) = mean(results2((R*i)+1:(R*i+R),1));
    ResultsStats(i+1,2) = std(results2((R*i)+1:(R*i+R),1));
    
    ResultsStats(i+1,3) = mean(results2((R*i)+1:(R*i+R),2));
    ResultsStats(i+1,4) = std(results2((R*i)+1:(R*i+R),2));
end

filenames2 = table();
filenames2(1,1) = results(1, 'filename');
for i = 1:((conditions*wells)-1)
    filenames2(i+1,1) = results ((i*R)+1, 'filename');
end

ResultsStats2 = array2table (ResultsStats, 'VariableNames' ,{'meanCount', 'SDcount', 'meanArea', 'SDarea'});
ResultsStats2 = [filenames2 ResultsStats2];

writetable(ResultsStats2,'LD_Exp04-3_meanperImage.csv', 'Delimiter', ' ')

%% ---statistics per condition --- 
ResultsStats3 = zeros(conditions, 4); % stats array for mean values + SD count/area
ResultsStats3(1,1) = mean(ResultsStats(1:wells,1));%calculate mean/ SD of count, column 1, mean/well
ResultsStats3(1,2) = std(ResultsStats(1:wells,1));
ResultsStats3(1,3) = mean(ResultsStats(1:wells,2));%calculate mean/ SD of area, column 2
ResultsStats3(1,4) = std(ResultsStats(1:wells,2));

for i = 1:((conditions)-1) 
    
    ResultsStats3(i+1,1) = mean(ResultsStats((wells*i)+1:(wells*i+wells),1));
    ResultsStats3(i+1,2) = std(ResultsStats((wells*i)+1:(wells*i+wells),1));
    
    ResultsStats3(i+1,3) = mean(ResultsStats((wells*i)+1:(wells*i+wells),2));
    ResultsStats3(i+1,4) = std(ResultsStats((wells*i)+1:(wells*i+wells),2));
end

filenames3 = table();
filenames3(1,1) = filenames2(1, 'Var1');
for i = 1:((conditions)-1)
    filenames3(i+1,1) = filenames2 ((i*wells)+1, 'Var1');
end

ResultsStats3 = array2table (ResultsStats3, 'VariableNames' ,{'meanCount', 'SDcount', 'meanArea', 'SDarea'});
ResultsStats3 = [filenames3 ResultsStats3];

writetable(ResultsStats3,'LD_Exp04-3_meanperCondition.csv', 'Delimiter', ' ')

%% --- plot Data --- 
