
function results = myimfcn(im)
%Image Processing Function
%
% IM      - Input image.
% RESULTS - A scalar structure with the processing results.
% 
% input images: cells treated with Calcium oxalate imaged with a
% LeicaCTR6000 microscope, 63x oil objective in DIC mode - 
% crystals appear very bright and can be easily thresholded

%--------------------------------------------------------------------------
% Auto-generated by imageBatchProcessor App. 
%
% When used by the App, this function will be called for every input image
% file automatically. IM contains the input image as a matrix. RESULTS is a
% scalar structure containing the results of this processing function.
%
%--------------------------------------------------------------------------


%%---Thresholding image - binary and calculating area---
Ibw = im2bw(im, 0.08); %thresholding value, make binary image
totalPx = numel(Ibw); %total pixel number of image
numWhitePx = sum(Ibw(:)); %calculate total number of white pixels
percentWhitePx = numWhitePx/ (totalPx.*0.01); %calculate percentage of area of white pixel - COM adhered


results.Ibw = Ibw
results.percentWhitePx = percentWhitePx

%--------------------------------------------------------------------------
