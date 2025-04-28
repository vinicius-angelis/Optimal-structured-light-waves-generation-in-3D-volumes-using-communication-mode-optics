% Convert image to csv
clear;clc;
mm = 1e-3; um =1e-6; nm = 1e-9;

slmx = 15.36*mm; slmy = 9.6*mm;
slmN = 1920; slmM = 1200;
slmPix = 8*um; 
pixX = 0:1919; pixY = [0:1199]';


x = linspace(-slmx/2,slmx/2,slmN); y = linspace(-slmy/2,slmy/2,slmM);
[X,Y] = meshgrid(x,y);

% Read image

imgname = 'holog_test.png';


a = double(imread(imgname));
phase = round(a/max(max(a))*(2^10-1));

figure; imagesc(x/mm,y/mm,phase); colorbar;
xlabel('mm'); ylabel('mm');

out = cat(2,pixY,phase); % append y-axis label


%% Write data

filename = 'holog_test';

xaxis = {};
for j = 1:slmN+1
    if j == 1
        xaxis(j) = {'Y/X'};
    else
        xaxis(j) = {pixX(j-1)};
    end
end

writecell(xaxis,[filename,'.csv'])

% Append data
dlmwrite([filename,'.csv'],out,'delimiter',',','-append');
%%




