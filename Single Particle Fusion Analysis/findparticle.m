%function [centroid avMovie filtImage mergedRGB] = findparticle(moviefile,varargin)
function [centroid avMovie filtImage mergedRGB] = findparticle(movielengths,moviefile,varargin)
% This function detects particles in a multiimage tiff movie.
%
% moviefile = String containing name of multiimage tiff file to be read.
%             Several frames are averaged to create a low-noise template image
%             defined by
%                 initframe  = first frame to be averaged
%                 finalframe = last frame
%                 Xmargin    = margin to leave out on left and right
%                 Ymargin    = margins on top and bottom of each channel
%             The frames should be chosen in an interval where the virus particles
%             are immobilized just prior to the onset of fusion.
%
% Syntax:
%             [...] = findparticle(moviefile)
%             [...] = findparticle(moviefile,'Property',VALUE)
%
% Argument Properties:
% 
% 'Frames'    {[FramesVector],[FramesVector]} or [FramesVector] for one channel.
%             Vector of frames to use for averaged image(s).
%
% 'Box'       {[StartRow EndRow; StartCol EndCol],[StartRow EndRow; StartCol EndCol]} or
%             [StartRow EndRow; StartCol EndCol] for one channel.
%
%----
%
% 'ROWS'      [initROW finalROW]
% 'COLS'      [initCOL finalCOL]
%             This is the column (x) and row (y) area of the image to be
%             read from disk and analyzed. Default setting for ROWS/COLS is
%             [1 512], i.e. all pixels.
%             The image will be read as a single channel.  No channel correlation will be
%             performed.
% 'Channels'  'red' or 'green' 
%             Use this property if you want to specify a single emission channel to
%             be used.  Red uses rows 11:246 and cols 5:502. Green uses rows
%             256:502 and columns 5:502.
%             Use 'ROWS' and 'COLS' if you need to specify a different region
%             for single channel analysis.
%
% 'ChannelOrientation'  'TopBottom' (Default) or 'LeftRight'
%             Specifies how the CCD is divided for 2-color (red/green) images.
%             'TopBottom' = top half contains red emission, bottom contains green
%             'LeftRight' = left half contains red, right half is green
%

% Code updates:
% -------------
% 2009-10-22 - Added the margin to avMovie to make it frame sized (allows
%              simple plotting of detected centroids).
%            - Plot figures showing detected particle centroids on image
%              slice used for search.
% 2009-10-18 - Reoptimized findparticle code.
% 7/9/2007 - modified to include ImageMerge.m (uses normxcorr2)
% 9/3/2007 - modified to generalize for one or two channel particle detection


%% Adjustable parameters
searchSize = 180;               % Cross correlation template size
xmargin = 5;  ymargin = 30;     % used in xcorrel
FiltNoiseSize = 1;

%% Processing optional arguments
% Variables:
% varargin     = variable length input argument list; matlab cell array; each
%                cell contains either the argument property or its value in the order input into the function call
% optargin     = scalar; number of variable arguments input into the function * 2
% mergeImage   = scalar; boolean value telling whether to merge channels or
%                not (1 = merge, 0 = keep channels separate)
% variable     = scalar; counter for reading in variable arguments
% propertyName = string; temporary variable containing argument property names

optargin = size(varargin,2); % Second argument of size is the dimension --> number of cols in vargin
if rem(optargin,2)~=0        % Error checking: optargin must be multiple of 2 b/c need argument property and its values
    error('Invalid syntax. Input property name followed by value');
end

% Loop assigns varargin's to their variable names and values
% --> e.g. Frames = [100 110]
for variable = 1:2:optargin 
    propertyName = varargin{variable};
    eval([propertyName ' = varargin{variable+1};']);    % Remember, propertyName is a string
    if ~ischar(propertyName)
        error('Invalid syntax. Input property name followed by value.');
    end
end

%if ~exist('Frames','var')  Frames = [1];  end                           % If not passed in, use only 1 frame, for only 1 channel.
%if (~iscell(Frames))  Frames = {Frames};  end
if (length(Frames) < 2)  mergeImage = 0;  else  mergeImage = 1;  end    % Merge emission channel images by default unless optional arg's are used

if ~exist('ChannelOrientation','var')
    ChannelOrientation = 'TopBottom';   % Default setting for two
end


%% Image processing
% Variables:
% numFrames = scalar; number of frames of non-moving particles to average
%numFrames = finalframe-initframe+1;
% for kk = 1:length(Frames)
%     ROWS{kk} = Box{kk}(1,:); COLS{kk} = Box{kk}(2,:);
%     avMovie{kk} = uint16( zeros(diff(ROWS{kk})+1, diff(COLS{kk})+1) );
%     for idx = Frames{kk}
%         avMovie{kk} = avMovie{kk} + imread(moviefile,'tif',idx,'PixelRegion',{ROWS{kk},COLS{kk}});
%     end
%     normfac = 1/length(Frames{kk});
%     avMovie{kk} = uint16(avMovie{kk} * normfac);        % Average the frames and bandpass-filter; convert to uint16 - 2 bytes per element: 0 to 65,535
% 
% %    filtImage{kk} = bpass(avMovie{kk},1,5);
%     filtImage{kk} = bpass(avMovie{kk},1.8,5);   % Characteristic noise size around: 1.8 pixels, and real pixel size a bit smaller than: 5 pixels
% end
for kk = 1:length(Frames)
    ROWS{kk} = Box{kk}(1,:); COLS{kk} = Box{kk}(2,:);
    avMovie{kk} = double( zeros(diff(ROWS{kk})+1, diff(COLS{kk})+1) );
    fprintf('              ');
    for idx = Frames{kk}
%        avMovie{kk} = avMovie{kk} + double( imread(moviefile,'tif',idx,'PixelRegion',{ROWS{kk},COLS{kk}}) );
        if(rem(idx,10)==0),  fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b%05d of %05d', idx, Frames{kk}(end));  end;  % Progress on screen
        avMovie{kk} = avMovie{kk} + double( imreadMF(movielengths,moviefile,'tif',idx,'PixelRegion',{ROWS{kk},COLS{kk}}) );
    end
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b%05d of %05d\n', idx, Frames{kk}(end));
%    normfac = 1/length(Frames{kk});
    normfac = 65535/max(avMovie{kk}(:));
    avMovie{kk} = avMovie{kk} * normfac;        % Average the frames and bandpass-filter; convert to uint16 - 2 bytes per element: 0 to 65,535

%    filtImage{kk} = bpass(avMovie{kk},1,5);
%    filtImage{kk} = bpass(avMovie{kk},1.8,5);   % Characteristic noise size around: 1.8 pixels, and real pixel size a bit smaller than: 5 pixels
%    filtImage{kk} = bpass(avMovie{kk},FiltNoiseSize,5);   % Characteristic noise size around: FiltNoiseSize pixels, and real pixel size a bit smaller than: 5 pixels
    filtImage{kk} = bpass(avMovie{kk},FiltNoiseSize,6); % [100614TI] adjusted value for low power movies
    avMovie{kk} = uint16(avMovie{kk});
end

%% Optional 2-color cross correlation and merging
if (mergeImage == 1)
    fprintf('Generating merged two-channel image ... ');   % Merge both channels into a composite image for particle detection

    % If channels are L/R rotate image -90 degrees before calling ImageMerge
    if( strcmp(ChannelOrientation,'LeftRight') ),  for kk = 1:length(Frames), filtImage{kk} = rot90(filtImage{kk},-1); end;  end

    [processedImage ChanOffset mergedRGB] = ImageMerge(filtImage,searchSize,xmargin,ymargin);

    if( strcmp(ChannelOrientation,'LeftRight') )
        processedImage = rot90(processedImage,1);  mergedRGB = imrotate(mergedRGB,90);  ChanOffset = rot90(ChanOffset,1);
    end
    fprintf('done.\n');
else
    processedImage = filtImage{1};
end


%% Main algorithm
% Find high intensity regions and calculate center
[partROI numParticles] = bwlabel(imextendedmax(processedImage,5));
partStats = regionprops(partROI,'Centroid');

centroid = zeros(numParticles,2);
for n = 1:numParticles
    centroid(n,:)= partStats(n).Centroid;
end

if (mergeImage == 1)
    % Shift centroid coordinates to original coordinates in two channel image
    %   ChanOffset --> [row col]
    %   centroid(x_r,y_r x_g,y_g) -->[col row col row]
%     centroid = [centroid(:,1)+ChanOffset(1,2) centroid(:,2)+ChanOffset(1,1) ...
%                 centroid(:,1)+ChanOffset(2,2) centroid(:,2)+ChanOffset(2,1)];
    tempcentroid = centroid;
    for kk = 1:length(Frames)
        centroid(:,2*kk-1) = tempcentroid(:,1) + ChanOffset{kk}(2);    % Shift x coordinate by origin of channel box wrt origin of total image
        centroid(:,2*kk)   = tempcentroid(:,2) + ChanOffset{kk}(1);    % Shift y coordinate by origin of channel box wrt origin of total image
    end
    clear tempcentroid;
end

% Show detected particles
for kk = 1:length(Frames)
    figure; imagesc(avMovie{kk}); set(gca,'ydir','normal'); axis image; colorbar; colormap jet;
    hold on; plot (centroid(:,2*kk-1),centroid(:,2*kk),'or'); hold off;
    title(['avMovie{' num2str(kk) '}: Detected particle centroids on average movie intensity']);
    figure; imagesc(log10(single(avMovie{kk}-min(avMovie{kk}(:))))); set(gca,'ydir','normal'); axis image; colorbar; colormap hsv;
    hold on; plot (centroid(:,2*kk-1),centroid(:,2*kk),'or'); hold off;
    title(['avMovie{' num2str(kk) '}: Detected particle centroids on log10(average movie intensity)']);
    figure; imagesc(log10(single(filtImage{kk}))); set(gca,'ydir','normal'); axis image; colorbar; colormap gray;
    hold on; plot (centroid(:,2*kk-1),centroid(:,2*kk),'or'); hold off;
    title(['filtImage{' num2str(kk) '}: Detected particle centroids on log10(filtered avMovie intensity)']);
end

for kk = 1:length(Frames)
    centroid(:,2*kk-1) = centroid(:,2*kk-1) + COLS{kk}(1)-1;    % Shift x coordinate by origin of channel box wrt origin of total image
    centroid(:,2*kk)   = centroid(:,2*kk)   + ROWS{kk}(1)-1;    % Shift y coordinate by origin of channel box wrt origin of total image
end
