function [time_ms,info,IsHamamatsuMovie] = getandortiffinfo(filename)
% GETANDORTIFFINFO - Extracts proprietary data structures from Andor Camera
%                    software saved multi-image TIFF files.
%
% Syntax:  [time_ms, info, IsHamamatsuMovie] = getandortiffinfo(filename)
%
% Inputs:       filename - name of single or multi-image TIFF file
%
% Outputs:      time_ms  - vector with number of elements equal to number
%                          of frames in multi-image TIFF, containing time
%                          (milliseconds) of each frame.
%
%               info     - [OPTIONAL, size NFRAMES] complete data structure
%                          with all proprietary Andor tag data in each frame.
%
%               IsHamamatsuMovie - 0 = no, 1 = yes; whether the movie
%                          consists of interleaved red channel and green
%                          channel frames (rather than both channels on every frame).
%
%
% Author: M. Popovic (milos AT mit DOT edu) - Feb 22, 2009
%

% Andor Camera TIFF Tags:
%               34361    - Andor tag (first slide only) with some header
%                          information and file name.
%               34362    - Andor tag (every slide) with time (ms) info.
%               34363    - Andor tag (first slide only) with unknown info.
%
% Hamamatsu C9100-13 Camera TIFF Tags:
%               ImageDescription - standard TIFF-format tag, stores time in
%               hh:mm:ss.0000 format following 'Time_From_Start = ' string.
%
% Uses TIFF processing code based on imtifinfo.m (Revision 1.1.6.5)
% function provided in Matlab R2008b.
%
% -------------------------------------------------------------------------

% Code updates:
% -------------
% Oct 24, 2009 - Updated code to output a flag indicating interleaved movies.
% Jun 12, 2009 - Updated code to add recognition of Hamamatsu C9100-13
%                camera type TIF tags storing time data records.
% Feb 22, 2009 - Initial code (supporting Andor camera TIF tags that store
%                frame time data records).

TIFFTAG_IMAGEDESCRIPTION = 270;            % Define some standard records
FLAG_CAMERADETECTED = 0;                   % Used to track whether camera
                                           % type has been detected, for
                                           % purposes of diagnostic output.

% 1) Open TIFF file assuming little-endian ordering
%
% TIFF files might be little-endian or big-endian.  Start with
% little-endian.  If we're wrong, we'll catch it down below and
% reopen the file.
[fid, msg] = fopen(filename, 'r', 'ieee-le');
if (fid == -1)
    error('MATLAB:imtifinfo:fileOpen', ...
          'Unable to open file "%s" for reading: %s.', filename, msg);
end


% 2) Read TIFF signature in first 4 bytes, reopen if big-endian
sig = fread(fid, 4, 'uint8')';
if (~isequal(sig, [73 73 42 0]) && ...     % Not a TIFF file?
    ~isequal(sig, [77 77 0 42]))
    fclose(fid);
    error('MATLAB:imtifinfo:notTIFF', ...
          'Not a TIFF file');
end

if (sig(1) ~= 73)                          % Opened with wrong endianness?
%    byteOrder = 'big-endian';             % Whoops!  Must reopen the file.
    pos = ftell(fid);
    fclose(fid);
    fid = fopen(filename, 'r', 'ieee-be');
    fseek(fid, pos, 'bof');
%else
%    byteOrder = 'little-endian';
end


% 3) Start going through IFD (Image File Directory), which is a "file
%    allocation table" of the contained slides in a multiimage TIFF.
nextIFDOffset = fread(fid, 1, 'uint32');   % [MP] Get offset of first slide
k = 0;  % number of images

while (nextIFDOffset ~= 0)
    status = fseek(fid, nextIFDOffset, 'bof');
    if (status ~= 0)
        % The seek to find the next IFD just failed.  This is probably
        % because of a nonconforming TIFF file that didn't store the
        % offset to the next IDF correctly at the end of the last IFD.
        % The best we can do here is assume there aren't any more IFD's.
        break;
    end
    
    k = k + 1;                             % [MP] Found the next slide

    tagCount = fread(fid, 1, 'uint16');
    tagPos = ftell(fid);
    %
    % Process the tags
    %
    for p = 1:tagCount
        fseek(fid, tagPos, 'bof');
        tagID = fread(fid, 1, 'uint16');

        switch tagID

            % Ignore most tags, process only Andor tags 34361, 34362, 34363

            % Andor TIFF tag with generic camera data (first frame only)
            case 34361          % Andor Camera proprietary tag in each TIFF image frame
%                type = fread(fid, 1, 'uint16');   % Should be type == 1, for BYTE type
                fseek(fid, 2, 'cof');              % Skip 2 bytes, they simply give type which is 1
                count = fread(fid, 1, 'uint32');   % How many bytes long is the tag data?
                if (count <= 4)
                    data1 = fread(fid, count, 'uchar=>char');
                else
                    offset = fread(fid, 1, 'uint32');
                    fseek(fid, offset, 'bof');
                    data1 = fread(fid, count, 'uchar=>char');
                end
                info(k).Andor1 = data1.';

            % Andor TIFF tag containing frame time (ms) in each frame, as
            % the first element (other 15 numbers are zeros)
            case 34362          % Andor Camera proprietary tag in each TIFF image frame
%                type = fread(fid, 1, 'uint16');   % Should be type == 1, for BYTE type
                fseek(fid, 2, 'cof');              % Skip 2 bytes, they simply give type which is 1
                count = fread(fid, 1, 'uint32');   % How many bytes long is the tag data?
                if (count <= 4)
                    data2 = fread(fid, count/4, 'double=>double');
                else
                    offset = fread(fid, 1, 'uint32');
                    fseek(fid, offset, 'bof');
                    data2 = fread(fid, count/4, 'double=>double');
                end
                info(k).Andor2 = data2.';

            % Andor TIFF tag (first frame only)
            case 34363          % Andor Camera proprietary tag in each TIFF image frame
%                type = fread(fid, 1, 'uint16');   % Should be type == 1, for BYTE type
                fseek(fid, 2, 'cof');              % Skip 2 bytes, they simply give type which is 1
                count = fread(fid, 1, 'uint32');   % How many bytes long is the tag data?
                if (count <= 4)
                    data3 = fread(fid, count, 'uchar=>char');
                else
                    offset = fread(fid, 1, 'uint32');
                    fseek(fid, offset, 'bof');
                    data3 = fread(fid, count, 'uchar=>char');
                end
                info(k).Andor3 = data3.';


            % For second type of camera: Hamamatsu C9100-13
            case TIFFTAG_IMAGEDESCRIPTION
                % ImageDescription tag
                fseek(fid, 2, 'cof');
                count = fread(fid, 1, 'uint32');
                if (count <= 4)
                    info(k).ImageDescription = char(fread(fid, count, 'uint8')');
                else
                    offset = fread(fid, 1, 'uint32');
                    fseek(fid, offset, 'bof');
                    info(k).ImageDescription = char(fread(fid, count, 'uint8')');
                end

                % don't use null-terminated strings in MATLAB
                if (info(k).ImageDescription(end) == 0)
                    info(k).ImageDescription(end) = [];
                end

        end  %%%% switch
        
        tagPos = tagPos + 12;
        
    end  %%%% for
    
    fseek(fid, tagPos, 'bof');
    nextIFDOffset = fread(fid, 1, 'uint32');
    
end  %%%% while

numImages = k;                              % [MP] Error if no images found
if (numImages == 0)
    fclose(fid);
    error('MATLAB:imgifinfo:noImages', ...
        'No images found in TIFF file');
end

% Extract time from collected records
for k = 1:numImages
    if isfield(info(k),'Andor2')                % Andor camera
        if(~FLAG_CAMERADETECTED)  FLAG_CAMERADETECTED = 1; fprintf('Detected Andor camera image file.\n'); IsHamamatsuMovie = 0;  end
        time_ms(k,1) = info(k).Andor2(1);       % [MP] Fill time column vector
    elseif isfield(info(k),'ImageDescription')
        idx = strfind(info(k).ImageDescription, 'Time_From_Start');
        if ~isempty(idx)                        % Hamamatsu C9100-13 camera
            if(~FLAG_CAMERADETECTED)  FLAG_CAMERADETECTED = 1; fprintf('Detected Hamamatsu C9100-13 camera image file.\n'); IsHamamatsuMovie = 1;  end
            hh = str2double(info(k).ImageDescription(idx+18:idx+19));
            mm = str2double(info(k).ImageDescription(idx+21:idx+22));
            ss = str2double(info(k).ImageDescription(idx+24:idx+25));
            ms = str2double(info(k).ImageDescription(idx+27:idx+30));
            time_ms(k,1) = ((hh*60 + mm)*60 + ss)*1000 + ms*0.1;
        end
    end
end
if(~FLAG_CAMERADETECTED)  fprintf('No known camera image file detected for time data extraction.  Time variable set to empty vector.\n'); time_ms = [];  end

fclose(fid);
