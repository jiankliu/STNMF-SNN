
function [newMap]= colormap_RF(datarange)
%
% function gauss = receptfield(img,...)
%
%	Define the colormap for Receptive field with 0 colored by white,
%   postive by red, negative by blue.
%
%	Input:		data   				image to be fitted
%
%	Output:		gauss				resulting Gaussian best fit parameters
%

N = 256;

% 1 0 0 red
% 1 1 1 white
% 0 0 1 blue

minV = datarange(1);
maxV = datarange(2);
int = (maxV-minV)/N;
NegativeRange = (0:int:abs(minV))/abs(minV);

PostiveRange = (maxV:-int:0)/maxV;

newMap = [NegativeRange'            NegativeRange' ones(size(NegativeRange'));
          ones(size(PostiveRange')) PostiveRange'  PostiveRange'; ];
len = size(newMap,1);
if len>256
    newMap = newMap(1:256,:);
end
newMap(end,:) = [1 0 0];
