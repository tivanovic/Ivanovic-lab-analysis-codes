function phprobe = getphprobefromtiff(filename,frameset) %,boundbox)
% GETPHPROBEFROMTIFF - Extract ph probe data from TIFF image movie
%
% Syntax: phprobe = getphprobefromtiff(tiff_filename, frameRange)
%
% frameRange    - [startframe endframe] Frames from which to get phprobe,
%                 inclusive.

% Code updates:
% -------------
% 2009-10-25 - Changed code to use frameset passed in.  Considered using
%              box for each channel as bounding box but decided against for
%              now.
%            - If topleft/right clicks are in top half plane, then using
%              whole image for background, if they are in bottom half
%              plane, then using only bottom half.
% 2009-07-03 - Added help text header.
%            - Deleting phprobe box selection figure right away after clicks
%            - Changed numFrames (end frame) to frameRange to control both
%              start and end frame of extracted data.
%
% TODO:      - Integrate this with trajectory extraction in extractdata.m,
% -----        so that the movie is read only once for both operations, not once for each.

%if (length(frameRange) < 2) frameRange = [1 frameRange]; end  % If length 1, then assume frameRange was end frame

% Get boxes
%A = double(imread(filename,[],frameRange(1))).';  % [20090704] Use first frame in range for box selection.
A = double(imread(filename,[],frameset(1))).';  % [20090704] Use first frame in range for box selection.
for kk = 1:4    % Average 1st five frames
    A = A + double(imread(filename,[],frameset(1+kk))).';  % [20090704] Use first frame in range for box selection.
end

figure; imagesc(A.'); axis image;
%rectangle('Position', [boundbox(1,1), boundbox(2,1), boundbox(1,2)-boundbox(1,1), boundbox(2,2)-boundbox(2,1)], 'EdgeColor', [1 1 0]);

title('Click on inner corner of top-left box');
xyTL = ginput(1);
xyTLint = round(xyTL);
title('Click on inner corner of top-right box');
xyTR = ginput(1);
xyTRint = round(xyTR);
title('Click on inner corner of bottom-left box');
xyBL = ginput(1);
xyBLint = round(xyBL);
title('Click on inner corner of bottom-right box');
xyBR = ginput(1);
xyBRint = round(xyBR);
title('Click on top-left corner of central box');
xyTLC = ginput(1);
xyTLCint = round(xyTLC);
title('Click on bottom-right of central box');
xyBRC = ginput(1);
xyBRCint = round(xyBRC);

close; drawnow;
fprintf('Processing frames for phprobe data...\n');

% Find background
if (xyTLint(2) <= 256 || xyTRint(2) <= 256)  ytop = 1;  else  ytop = 257;  end

%for k = frameRange(1):frameRange(2)   % 1:numFrames    % [20090704MP] Make possible providing arbitrary frame range
for k = 1:length(frameset)
   A = double(imread(filename,[],frameset(k))).';

   Atopleft = A(1:xyTLint(1), ytop:xyTLint(2));
   Atopright = A(xyTRint(1):end, ytop:xyTRint(2));
   Abotleft = A(1:xyBLint(1), xyBLint(2):end);
   Abotright= A(xyBRint(1):end, xyBRint(2):end);
   Acenter = A(xyTLCint(1):xyBRCint(1), xyTLCint(2):xyBRCint(2)); 
   
   Ibackground(k,1) = sum([Atopleft(:); Atopright(:); Abotleft(:); Abotright(:)]) / length([Atopleft(:); Atopright(:); ...
                      Abotleft(:); Abotright(:)]);
   Isignal(k,1) = sum(Acenter(:)) / length(Acenter(:));
end

Inormalized = Isignal - Ibackground;

phprobe = Inormalized;
