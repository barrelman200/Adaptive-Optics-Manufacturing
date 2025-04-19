function [error] = image_error(Intensity1, Intensity2)

% Normalize the two input intensity profiles
Intensity1=Intensity1/sum(Intensity1);
Intensity2=Intensity2/sum(Intensity2);

% Calculate the error image with an L2 norm
errorim = (Intensity1-Intensity2).^2;

% Mean value over the whole image
error = sqrt(sum(errorim(:)));