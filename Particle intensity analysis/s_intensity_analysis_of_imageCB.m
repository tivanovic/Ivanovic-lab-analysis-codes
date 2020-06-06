% Sep 9, 2019
% Sep 14, 2019- CB (curved background) version -- version that does not get
% noise from 20x20 pixel areas in the 4 corners, but instead subtracts a
% background image that is acquired without the virus particles present.

clear all; close all;

% 1. Set movie file to load (e.g. '190612_032.tif')
%filename = '190612_032.tif'; moviefile = filename;
[filename,pathname]=uigetfile('*.tif','Please select movie file:');
if isequal(filename,0)
   disp('User selected Cancel')
   return
end
% else
   moviefile = fullfile(pathname, filename);
   disp(['User selected ', moviefile])
% end


% 2. Load movie and average 10 frames
FramesToAverage = [1:10];
aa = zeros(512,512);                                                            % aa = Imported image averaged over specified frames
for kk = FramesToAverage, aa = aa+double(imread(moviefile,'tif',kk)); end       % Load movie frames and add them up...
aa = aa/length(FramesToAverage);                                                % ... and then make aa the average value (and same magnitude as an individual frame).

figure; imagesc(aa); axis image; colormap(hot); colorbar; set(gca,'ydir','normal');
%caxis([16626 50000])
title([filename ': averaged, linear scale (caxis rescaled)']);

figure; imagesc(log10(aa)); axis image; colormap(hot); colorbar; set(gca,'ydir','normal');
title([filename ': averaged, log scale']);


% 1b. Set movie file to load (e.g. '190612_032.tif') that contains
% background without the particles.
%filename = '190612_032.tif'; moviefile = filename;
[filenameBG,pathnameBG]=uigetfile('*.tif','Please select movie file:');
if isequal(filenameBG,0)
   disp('User selected Cancel')
   return
end
% else
   moviefileBG = fullfile(pathnameBG, filenameBG);
   disp(['User selected ', moviefileBG])
% end


% 2b. Load movie and average 10 frames
aaBG = zeros(512,512);                                                            % aa = Imported image averaged over specified frames
for kk = FramesToAverage, aaBG = aaBG+double(imread(moviefileBG,'tif',kk)); end       % Load movie frames and add them up...
aaBG = aaBG/length(FramesToAverage);                                                % ... and then make aa the average value (and same magnitude as an individual frame).

figure; imagesc(aaBG); axis image; colormap(hot); colorbar; set(gca,'ydir','normal');
%caxis([16626 50000])
title([filenameBG ': averaged, linear scale (caxis rescaled)']);

figure; imagesc(log10(aaBG)); axis image; colormap(hot); colorbar; set(gca,'ydir','normal');
title([filenameBG ': averaged, log scale']);


% 3. Remove background
a2 = aa-aaBG;  % No abs value! 
% a2 = aa-aaBG/2;  %TI190914 divide aaBG by 2 because background acquired at 2x the power as the particle samples
% Test
figure; imagesc(abs(a2)); axis image; colormap(hot); colorbar; set(gca,'ydir','normal');
%caxis([16626 50000])
title([filenameBG ': averaged, linear scale (caxis rescaled), curved background (CB) removed']);

figure; imagesc(log10(abs(a2))); axis image; colormap(hot); colorbar; set(gca,'ydir','normal');
title([filenameBG ': averaged, log scale, curved background (CB) removed']);
cvec=caxis; caxis([cvec(2)-2 cvec(2)]);     % Look at top 2 orders of magnitude


save([moviefile(1:end-4) '_intensitytemplate'],'a2');
