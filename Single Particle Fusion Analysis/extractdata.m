%function trajectories = extractdata(moviefile, centroid, frames, track_toggle)
function trajectories = extractdata(movielengths, moviefile, centroid, frames, track_toggle)
% This function takes the output from findparticle. First a square ROI is
% defined around each centroid position. The corresponding multi image
% tiff file is then read frame by frame and pixels within each ROI are
% integrated and saved in a structure array called trajectory. The red
% and/or green channels are saved as separate fields within the structure
% (trajectories.redTraj or trajectories.greenTraj). The channels to be
% extracted are determined by the number of xy columns in centroid array
% and where in the CCD image they belong.  Red is assumed to be the top
% half (ROWS 1-256) and green is the bottom half (ROWS 257-512).

% track_toggle determines whether stage drift will be tracked:
%   track_toggle == 0 skip this function
%   track_toggle == 1 carry out drift analysis


% Determine from which channels to extract data
[numParticles, temp] = size(centroid);  numChannels = temp/2;  clear temp;

%% Track microscope stage drift
% %if track_toggle==1
% %   offset = driftTracker(moviefile,numFrames);
% %else
%     offset = zeros(numFrames,2);
% %end

for nchan = 1:numChannels   %length(frames)        % Go through frame set for each channel
    % Extract trajectories from images
    % if (do_green),  trajectories.greenTraj = zeros(numFrames,numParticles);  end
    % if (do_red),    trajectories.redTraj = zeros(numFrames,numParticles);    end
    trajectories{nchan} = double( zeros(length(frames{nchan}), numParticles) );

    fprintf('              ');
    for nframe = 1:length(frames{nchan})
%        curframe = imread(moviefile,'tif',frames{nchan}(nframe));
        curframe = imreadMF(movielengths, moviefile,'tif',frames{nchan}(nframe));
        if(rem(nframe,10)==0),  fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b%05d of %05d', frames{nchan}(nframe), frames{nchan}(end));  end;  % Progress on screen

        for np = 1:numParticles
            columnvect = round(centroid(np,2*nchan-1)) + [-2 1];   % Define ROI box around particle centroid for intensity integration
            rowvect    = round(centroid(np,2*nchan  )) + [-2 1];
            trajectories{nchan}(nframe,np) = ...
                sum(sum( curframe([rowvect(1):rowvect(2)]       , ...   % + offset(nframe,1), ...
                                  [columnvect(1):columnvect(2)] ) ));   % + offset(nframe,2)) ));
        end
    end    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b%05d of %05d\n', frames{nchan}(nframe), frames{nchan}(end));
end


% do_red = true; do_green = true;
% if (numChannels == 1)
%     if centroid(1,2) < 256 % Red channel is top half of CCD (centroid is(xCOL,yROW))
%         do_green = false;
%     else
%         do_red = false;
%     end
% end
% 
% % Generate x and y boundaries (in matrix indices) of ROI for each particle - in green and red channels
% if (do_green)   % Build green particle ROI (4x4 pixels)
%     greenColumnvect = zeros(numParticles,2); % Green ROI column vector  % [20091018MP] Contains [low high] column index for ROI "box" for each particle
%     greenRowvect    = zeros(numParticles,2); %       ROI row vector     % [20091018MP] Contains [low high] row index for ROI "box" for each particle
%     for n = 1:numParticles
%         greenColumnvect(n,1:2) = round(centroid(n,end-1)) + [-2 1];    % [20091018MP] end-1 is used because x is either 1st of 2 columns (green channel only) or 3rd of 4 columns (red+green)
%         greenRowvect(n,1:2)    = round(centroid(n,end))   + [-2 1];    % [20091018MP] end is used because y is either 2nd of 2 columns (green channel only) or 4th of 4 columns (red+green)
%     end
% end
% 
% if (do_red)     % Build red particle ROI (4x4 pixels)
%     redColumnvect = zeros(numParticles,2);  % Poreformation  ROI column vector
%     redRowvect    = zeros(numParticles,2);  %                ROI row vector
%     for n = 1:numParticles
%         redColumnvect(n,1:2) = round(centroid(n,1)) + [-2 1];
%         redRowvect(n,1:2)    = round(centroid(n,2)) + [-2 1];
%     end
% end
% 
% %% Track microscope stage drift
% %if track_toggle==1
%  %   offset = driftTracker(moviefile,numFrames);
% %else
%     offset = zeros(numFrames,2);
% %end
% 
% %% Extract trajectories from images
% if (do_green),  trajectories.greenTraj = zeros(numFrames,numParticles);  end
% if (do_red),    trajectories.redTraj = zeros(numFrames,numParticles);    end
% 
% for nframe = 1:numFrames
%     curframe = imread(moviefile,'tif',nframe);
%     for n = 1:numParticles
%         if (do_green)
%             trajectories.greenTraj(nframe,n) = ...
%                 sum(sum( curframe([greenRowvect(n,1):greenRowvect(n,2)] + offset(nframe,1), ...
%                                   [greenColumnvect(n,1):greenColumnvect(n,2)] + offset(nframe,2)) ));
%         end
%         if (do_red)
%             trajectories.redTraj(nframe,n) = ...
%                 sum(sum( curframe([redRowvect(n,1):redRowvect(n,2)] + offset(nframe,1), ...
%                                   [redColumnvect(n,1):redColumnvect(n,2)] + offset(nframe,2)) ));
%         end
%     end
% end
