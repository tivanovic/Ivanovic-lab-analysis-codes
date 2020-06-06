function [idxgreenchannel, idxredchannel] = parseinterleavedmovieframes(filename, StartFrame, EndFrame, SaveRedGreenMovies, movielengths)

% This function sorts movie with alternating channel images
% NUMSKIPS defines a maximum number of consecutive same channel images to be expected
%
% Syntax:  [idxgreenchannel, idxredchannel] = parseinterleavedmovieframes(filename, StartFrame, EndFrame, SaveRedGreenMovies)
%
% idxgreenchannel/idxredchannel are green/red channel frame numbers
% filename is the name of the multiimage tif file, or a cell array of multiple files.
% StartFrame is the number of the first frame to be analyzed
% EndFrame is the number of the last frame to be analyzed
%
% Authors: T. Ivanovic (ivanovic AT crystal DOT harvard DOT edu) and M. Popovic - 090924

NUMSKIPS = 2;    % MUST BE EVEN NUMBER; up to how many consecutive frames in the same channel to expect
if(nargin < 4)  SaveRedGreenMovies = 0;  end    % By default don't save individual red/green movies

for k = StartFrame:EndFrame
    frameimage = imreadMF(movielengths,filename,[],k);
    totalintensity(k) = sum(frameimage(:));
end
clear frameimage;

for k = 1:(NUMSKIPS+1)
    shiftedmatrix(k,:) = totalintensity(k:end-(NUMSKIPS+1)+k);
end

amin = min(shiftedmatrix);
amax = max(shiftedmatrix);
athreshold = 0.5*(amin + amax);
clear shiftedmatrix;
% Extrapolate threshold to edges
athreshold = [ones(1,((NUMSKIPS+1)-1)/2) * athreshold(1), athreshold, ones(1,((NUMSKIPS+1)-1)/2) * athreshold(end)];

isgreenchannel  = (totalintensity > athreshold);
idxgreenchannel = find( isgreenchannel );
idxredchannel   = find( ~isgreenchannel );

if(SaveRedGreenMovies)
    for k=1:numFrames
        frameimage = imread(filename,[],k);
        [pathstr,name,ext,versn]=fileparts(filename);
        if isgreenchannel(k)  %%totalintensity(k)>athreshold(k)
            filename2 = fullfile(pathstr,[[name '_greenchannel'] ext versn]);
            imwrite(frameimage, filename2, 'tif', 'Compression', 'none', 'WriteMode', 'append');
    %imwrite(myImage(:,:,i), ‘10-pk_img.tif’,'tif’, ‘Compression’, ‘none’, ‘WriteMode’, ‘append’);
        else
            filename2 = fullfile(pathstr,[[name '_redchannel'] ext versn]);
            imwrite(frameimage, filename2, 'tif', 'Compression', 'none', 'WriteMode', 'append');
        end
    end
end
