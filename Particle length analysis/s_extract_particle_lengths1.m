% Script to extract particle lengths from TIF images
clear all; close all;

% Get all TIF files in folder
ss = dir('*.tif');

for k = 1:length(ss)
    [~,fnameroot] = fileparts(ss(k).name);
    if( ~exist([fnameroot '.mat'], 'file') )
        % Read in and plot image
        A = imread(ss(k).name);
        fullscreen = get(0,'ScreenSize');
        figure('Position',[0 -50 fullscreen(3) fullscreen(4)]);
        colormap(gray); imagesc(A); axis image;
        disp(['Working on file ' ss(k).name ' ...']);

        % First, get scale
        title([ss(k).name ': Please click on two sides of scale bar so the connecting line runs through it.']);
        clickdata = ginput;
        jj = round( mean(clickdata(:,2)) );
        % Extract 1D line image through the scale bar
        scalebarimage = A( jj, round(clickdata(1,1)):round(clickdata(2,1)), 1 );
    %    plot(scaleimage);
        % Find min/max extents of scalebar
        scalebarindices = find(scalebarimage == 0);
        imin = min(scalebarindices); imax = max(scalebarindices);  % Pixels at which scalebar starts and ends

%        title('Please enter scalebar length in microns:');
        Lbar = input('Please enter scalebar length in um : ');
        Lpixel = Lbar / (imax - imin + 1);

        % Second, get particle data
        title([ss(k).name ': For each particle, click 2+ times to create length contour and hit enter; when done all particle hit enter again.']);

        clickdata = 0;
        particles = [];
        while ~isempty(clickdata)
            clickdata = ginput;
            if ~isempty(clickdata)
                particles(end+1).vertices = clickdata;

                % Compute lengths in pixels
                for mm = 1:length( particles(end).vertices )-1
                    particles(end).edgelen(mm) = ...
                        sqrt( (particles(end).vertices(mm+1,1) - particles(end).vertices(mm,1))^2 + ...
                              (particles(end).vertices(mm+1,2) - particles(end).vertices(mm,2))^2 );
                          line([particles(end).vertices(mm,1),particles(end).vertices(mm+1,1)],[particles(end).vertices(mm,2),particles(end).vertices(mm+1,2)])
                end
                particles(end).totallen = sum(particles(end).edgelen);

                % Put in length data in microns
                particles(end).Lbar_um = Lbar;
                particles(end).Lpixel_um = Lpixel;
                particles(end).totallen_um = particles(end).totallen * Lpixel;
            end
        end
        
        save(fnameroot, 'particles');
        saveas(gcf, fnameroot, 'fig' );
        print( gcf, '-dpng', fnameroot );
    else
        disp(['File ' fnameroot '.mat already exists; delete if you want to redo.']);
    end
end

disp('Done.  All data in ''particles'' structure.');
