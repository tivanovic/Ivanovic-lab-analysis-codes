% Experimental particle-selector UI
% Mar 12, 2012

% Code updates:
% -------------
% 2012-06-11    - Added Gaussian step-fit saving/loading.
%
% To do:
% ------
% - which kind of fit / how long a movie to consider
% - collect particles
% - export fit data to global workspace
%
% - change trace plotting to not replot but rather change xdata and ydata
%   only (will be faster)
% - when you do smaller timescale fits, then load only that much of the
%   movie
% - when dragging movie time indicator line (red), don't let it go left of
%   start or right of end
%
% - set color scale strategy in drawMovieZoom (constant color axis,
%   variable, log scale constant,???)
% - fix bug that clicking fast on slider calls callback many times
%   simultaneously so they all work at once and overwrite each other
% - make main figure window "modal", i.e. if you run it a second time, it
%   just brings up the same one, not makes a second copy.
% - don't redraw page if slider dragged slightly so stays on same page, but
%   callback still called b/c of floating point slider value slightly changed

function partselgui2
%clear all; close all;
global centroid timems hemiindex hemifusion hemifirstderiv hemiEventStartIntensity ...
       pHdropHemi NumptsAvgStartIntensity ndrop tdrop maxHemiDerivIndex ...
       poreindex porefirstderiv poreEventStartIntensity ...
       pHdropPore maxPoreDerivIndex idxgreenframes pHdrop dataset % poreformation

% fGaussFit = fittype('yofs + h*erfc((x-xo)/w)/2 * exp(-m*(x-xo))','coefficients',{'h','m','w','xo','yofs'},'independent','x');   % Drop step
fGaussFit = fittype('yofs + h*(erf((x-xo)/w)+1)/2 * exp(-m*(x-xo))','coefficients',{'h','m','w','xo','yofs'},'independent','x'); % Rise step


% Parameters
nx = 3; ny = 3;     % # of subplots in x and y
xsf = 0.75;         % Subplots are this (left) fraction of the x width of figure
ysf = 0.95;         % Subplots are this (top) fraction of the y height of figure
%npages = 20;
wsize = 50;         % Window size for movie zoom-in

% Create main figure
%f = figure;
scrsz = get(0,'ScreenSize');
hfig = figure('Position',[1 1 scrsz(3) scrsz(4)],'NumberTitle','off', ...
              'Name','Particle selection GUI'); %,'Resize','off');

% Create subplots
hsubplot = zeros(nx*ny,1);
for mm = 1:ny
    for nn = 1:nx
%         h(1+(nn-1)+(mm-1)*nx) = axes('OuterPosition',[(nn-1)/nx (mm-1)/ny 1/nx 1/ny]);
%         h(end+1) = axes('OuterPosition',[(nn-1)/nx (mm-1)/ny 1/nx 1/ny]);
%         h(end+1) = axes('OuterPosition',[(nn-1)*xsf/nx (1-ysf)+(mm-1)*ysf/ny xsf/nx ysf/ny], ...
%                         'Position',[(nn-1)*xsf/nx+0.02 (1-ysf)+(mm-1)*ysf/ny+0.03 xsf/nx-0.03 ysf/ny-0.04]);
        hsubplot((mm-1)*nx+nn) = axes('OuterPosition',[(nn-1)*xsf/nx (1-ysf)+(ny-mm)*ysf/ny xsf/nx ysf/ny], ...
                                    'Position',[(nn-1)*xsf/nx+0.02 (1-ysf)+(ny-mm)*ysf/ny+0.03 xsf/nx-0.03 ysf/ny-0.04]);
%        title(sprintf('Figure %d', nn+(mm-1)*nx));       % Test of subplot arrangement order
%        set(h(mm),'OuterPosition',[0 0 1/nx 1/ny]);
    end
end

% Initialize data set parameters
numParticles = size(centroid,1);
numPages     = ceil(numParticles/(nx*ny));
curPage      = 1;
flagSelect   = false(numParticles,2);    % Has the particle been left/right-click selected?  All start unselected.
nwstart = 30;   % We set here the start and end points of window subsection of particle trace to do fitting on (currently only for Gaussian integral fit)
%nwstart = pHdropHemi;
nwend = 150;
particledata = [];          % Starting a data structure to contain *all* input and output data.

% Create controls
h1=uicontrol('Parent', hfig, 'Style','slider','Units','normalized', ...
             'Position',[0.01 0.02 0.3 0.02], 'Min', 1, 'Value', curPage, 'Max', numPages, ...
             'SliderStep', [1/(numPages-1) 1/(numPages-1)], ...
             'TooltipString', 'Page through particle trace plots', ...
             'Callback',@sliderPages_Callback);

hpagelabel = uicontrol('Parent', hfig, 'Style','text','Units','normalized', ...
             'Position',[0.31 0.02 0.05 0.02],'String','Page','HorizontalAlignment','Left','BackgroundColor',get(hfig,'Color'));

hmovie = axes('OuterPosition',[xsf 0.5 (1-xsf) 0.5],'XTickLabel',[],'YTickLabel',[], ...
              'Position',[xsf 0.6 (1-xsf)-0.01 0.4-0.01]);
hmovieparticlemarker = []; hframeline = []; isWholeZoomMovieLoaded = 0; hzoomimage = [];

hmoviezoom = axes('OuterPosition',[xsf 0 (1-xsf) 0.5],'XTickLabel',[],'YTickLabel',[], ...
                  'Position',[xsf 0.1 (1-xsf)-0.01 0.4-0.01]);

haxisscaleselectbuttongroup = uibuttongroup('Parent', hfig, 'Units', 'normalized', ...
              'Position', [xsf 0.0625 0.3*(1-xsf) 0.035], 'BackgroundColor',get(hfig,'Color'), 'Title', 'Scale');
haxisscaleselectmax = uicontrol('Parent', haxisscaleselectbuttongroup, 'Style','radiobutton', 'Units', 'normalized', ...
              'Position', [0 0 0.5 1], 'String', 'Global', 'Value', 1);
haxisscaleselect1 = uicontrol('Parent', haxisscaleselectbuttongroup, 'Style','radiobutton', 'Units', 'normalized', ...
              'Position', [0.5 0 0.5 1], 'String', 'Region', 'Value', 0);

% hstepfittype = uicontrol('Parent', hfig, 'Style', 'popup', 'String', 'Max deriv|Gaussian integral|none', ...
%                          'Units','normalized','Position', [xsf 0 0.25*(1-xsf) 0.05], 'Callback', @popupStepfittype_Callback);
hstepfitmaxderiv = uicontrol('Parent', hfig, 'Style', 'checkbox', 'String', 'Max deriv', ...
                         'Units','normalized','Position', [xsf 0 0.25*(1-xsf) 0.05], 'Value', 1, 'Callback', @checkboxStepfitMaxderiv_Callback);
hstepfitgausstep = uicontrol('Parent', hfig, 'Style', 'checkbox', 'String', 'Gaussian integral', ...
                         'Units','normalized','Position', [xsf+0.2*(1-xsf) 0 0.25*(1-xsf) 0.05], 'Value', 1, 'Callback', @checkboxStepfitGausstep_Callback);

hcollectbutton = uicontrol('Parent', hfig, 'Style', 'pushbutton', 'String', 'Collect particles', ...
                           'Units','normalized','Position', [1-0.3*(1-xsf) 0 0.3*(1-xsf) 0.05], 'Callback', @pushbuttonCollect_Callback);

% Pre-color previously selected particles using saved data
selectHemiWinners = []; selectPoreWinners = []; selectHemiWinners_old = []; selectPoreWinners_old = [];
loadPreselectInfo
flagSelect(selectHemiWinners,1) = true;
flagSelect(selectPoreWinners,2) = true;

% Initialize variables for fit data to be exported at end
hemiEventFittedMean = []; hemiEventFittedWidth = [];

% Fit all the data
%fitGaussianStep
loadPrefitInfo

% Get movie dataset info
[~,moviename,~] = fileparts(pwd);
moviepointer  = ['../' moviename '_datasetconfig']; % e.g. ../120309_003_datasetconfig.mat
Qmov = load(moviepointer);
%    imread([dirname '.tif'],'Index','Info','PixelRegion');
tiffinfo = imfinfo(['../' moviename '.tif']);
Azoom = [];         % Start with blank matrix for zoomin movie
curplot = 1; curparticle = 1;    % Default start value.
curzoom = 1;                     % Particle currently zoomed
curframe = round(mean(Qmov.frames)); % Starting frame for time

% Draw particle trace data
drawMoviePlot
plotsRedraw
drawMovieZoom;


function loadPreselectInfo
    %% 1) If available, load previously selected events from a saved datafile
    disp([mfilename ': Plotting all event trajectories to allow user to select actual events'])
    [~,dirname] = fileparts(pwd);     % Find name of current folder
    aa = dir([dirname 'events*']);
    if(~isempty(aa))
        for kk = 1:length(aa),  aadates(kk) = datenum(aa(kk).date);  end;  idx = find(aadates == max(aadates));     % Find newest file
        [filename,pathname] = uigetfile({[dirname 'events*.mat'],'MITHEMI event data file ([directoryname]events*.mat)'}, ...
                                        'Select saved event file to load as starting point, or click CANCEL to start blank', aa(idx).name);
        if ~(isequal(filename,0) || isequal(pathname,0))  % If there's saved event data...
            QQ=load(filename);
            if(isfield(QQ,'selectHemiWinners')), selectHemiWinners = QQ.selectHemiWinners; selectHemiWinners_old = selectHemiWinners; end;
            if(isfield(QQ,'selectPoreWinners')), selectPoreWinners = QQ.selectPoreWinners; selectPoreWinners_old = selectPoreWinners; end;
        end
    else
        fprintf('No saved event data found in current folder.  Starting with blank selection.\n');
    end
%     if ~exist('selectHemiWinners','var'), selectHemiWinners = []; end;
%     if ~exist('selectPoreWinners','var'), selectPoreWinners = []; end;
end

function loadPrefitInfo
    %% If available, load previously done fits of Gaussian step
    [~,dirname] = fileparts(pwd);     % Find name of current folder
    aa = dir([dirname 'fits*']);
    if(~isempty(aa))
        for kk = 1:length(aa),  aadates(kk) = datenum(aa(kk).date);  end;  idx = find(aadates == max(aadates));     % Find newest file
        [filename,pathname] = uigetfile({[dirname 'fits*.mat'],'MITHEMI fits data file ([directoryname]fits*.mat)'}, ...
                                        'Select saved fits file to load, or click CANCEL to (re)run the fitting procedure', aa(idx).name);
        if ~(isequal(filename,0) || isequal(pathname,0))  % If there's saved event data...
            QQ = load(filename);
            hemiEventFittedMean = QQ.hemiEventFittedMean;
            hemiEventFittedWidth = QQ.hemiEventFittedWidth;
            particledata = QQ.particledata;
        else
            fprintf('No saved fits datafile was selected.  Starting pre-fitting procedure.\n');
            fitGaussianStep
        end
    else
        fprintf('No saved fits data found in current folder.  Starting pre-fitting procedure.\n');
        fitGaussianStep
    end
end

function sliderPages_Callback(hObject,eventdata)
    curPage = round(get(hObject,'Value'));
    set(hObject,'Value',curPage);
%     fprintf('Clicked! Page %d...\n', curPage);

    % TODO: if we changed the page, only then...
    plotsRedraw
    drawMovieZoom;
%    updateMoviePlot     % Only update positions of particles.

    guidata(hObject,eventdata);
end

function plotsRedraw
%     disp('Redrawing..');
    set(hpagelabel,'String',sprintf('Page %d of %d',curPage,numPages));
    isWholeZoomMovieLoaded = 0;     % When you redraw traces, clear the "movie cache"


% %    fGaussFit = fittype('yofs + h*erfc((x-xo)/w)/2 * exp(-m*(x-xo))','coefficients',{'h','m','w','xo','yofs'},'independent','x');   % Drop step
%     fGaussFit = fittype('yofs + h*(erf((x-xo)/w)+1)/2 * exp(-m*(x-xo))','coefficients',{'h','m','w','xo','yofs'},'independent','x'); % Rise step

    if ~isempty(hmovieparticlemarker),  delete(hmovieparticlemarker); hmovieparticlemarker = [];  end
    timesec = timems/1000;
    
    curzoom = nx*ny*(curPage-1)+1;                          % Zoom on first particle on page (top-left)
%    for curplot = 1:28
    for m = 1:ny
        for k = 1:nx
            curplot = (m-1)*nx + k;                         % Subplot #
            curparticle = nx*ny*(curPage-1) + curplot;
            
            if (curparticle > numParticles)
                if( strcmp(get(hsubplot(curplot),'Visible'),'on')  ),  set(hsubplot(curplot),'Visible','off');  end
            else
                if( strcmp(get(hsubplot(curplot),'Visible'),'off') ),  set(hsubplot(curplot),'Visible','on');   end
                axes(hsubplot(curplot));

                % Plot particle traces in the subplot
                hold off;
%                set(hsubplot(curplot), 'NextPlot', 'replace');  % Same as "hold on" but can specify for which axes (not current axes)

                if exist('hemifusion','var')
                     plot(timesec(hemiindex), hemifusion(:,curparticle), 'Color', [0 .5 0]);
%                    plot(timesec(hemiindex(1:nwend)), hemifusion((1:nwend),curparticle), 'Color', [0 .5 0]);
% %                    plot(hsubplot(curplot), timesec(hemiindex(1:nwend)), hemifusion((1:nwend),curparticle), 'Color', [0 .5 0]);
                    hold on;
%                    set(hsubplot(curplot), 'NextPlot', 'add');  % Same as "hold on" but can specify for which axes (not current axes)
                     hemifusionsmoothed = cumsum(hemifirstderiv(:,curparticle));
                    hemifusionsmoothed = hemiEventStartIntensity(curparticle)+hemifusionsmoothed-hemifusionsmoothed(pHdropHemi-NumptsAvgStartIntensity/2);
                    plot(timesec(hemiindex), hemifusionsmoothed, 'Color', [0 1 1])     % [20090705TIMP] Show also smoothed curve
                    ndrop = maxHemiDerivIndex(curparticle); tdrop = timesec(hemiindex(ndrop));
                    plot(timesec(hemiindex(ndrop)), hemifusionsmoothed(ndrop), ...
                          'r+','LineWidth',2,'MarkerSize',8)     % [20090705TIMP] Show also smoothed curve

%                     % Gaussian-integral step fit to event trace
%         %            yfit = fit(t(pHdrop:end), ydata(pHdrop:end), f, 'StartPoint', [(mean(ydata(pHdrop:ndrop))-mean(ydata(ndrop:end))) 0 (max(t)-min(t))/100 tdrop mean(ydata(ndrop:end))], 'Lower',[-inf 0 -inf -inf -inf]);
% %                     yfit = fit(timesec(hemiindex(pHdropHemi:end)), hemifusion(pHdropHemi:end,curparticle), fGaussFit, ...
% %                           'StartPoint', [(mean(hemifusion(pHdropHemi:ndrop,curparticle))-mean(hemifusion(ndrop:end,curparticle))) ...
% %                           0 (max(timesec(hemiindex))-min(timesec(hemiindex)))/100 tdrop mean(hemifusion(ndrop:end,curparticle))], ...
% %                           'Lower',[-inf 0 -inf -inf -inf]);
% %                     plot(timesec(hemiindex(pHdropHemi:end)), yfit(timesec(hemiindex(pHdropHemi:end))), 'k');
% %                     hemiEventFittedMean(curparticle) = yfit.xo-timesec(idxgreenframes(pHdrop));
% %                     hemiEventFittedWidth(curparticle) = 2*yfit.w;
%                     
% %                     yfit = fit(timesec(hemiindex(pHdropHemi:nwend)), hemifusion(pHdropHemi:nwend,curparticle), fGaussFit, ...
% %                           'StartPoint', [(max(hemifusion(pHdropHemi:nwend,curparticle))-min(hemifusion(pHdropHemi:nwend,curparticle))) ...
% %                           0 (max(timesec(hemiindex))-min(timesec(hemiindex)))/100 timesec(hemiindex(nwend))/2 min(hemifusion(pHdropHemi:nwend,curparticle))], ...
% %                           'Lower',[-inf 0 -inf -inf -inf]);
% 
% %                     yydata = hemifusion(pHdropHemi:nwend,curparticle);
% %                     yymin = min(yydata); yymax = max(yydata);
% %                     yfit = fit(timesec(hemiindex(pHdropHemi:nwend)), yydata, fGaussFit, ...
% %                           'StartPoint', [(max(yydata)-min(yydata)) ...
% %                           0 (max(timesec(hemiindex))-min(timesec(hemiindex)))/100 timesec(hemiindex(nwend))/2 min(hemifusion(pHdropHemi:nwend,curparticle))], ...
% %                           'Lower',[-inf 0 -inf -inf -inf]);
%                     yydata = hemifusion(pHdropHemi:nwend,curparticle);
%                     yymin = min(yydata); yymax = max(yydata);
%                     yynorm = (yydata-yymin)/(yymax-yymin);
%                     ttdata = timesec(hemiindex(pHdropHemi:nwend));
%                     ttmin = min(ttdata); ttmax = max(ttdata);
%                     ttnorm = (ttdata - ttmin)/(ttmax-ttmin);
%                     yfit = fit(ttnorm, yynorm, fGaussFit, ...
%                           'StartPoint', [1 0 1/100 0.5 0], ...
%                           'Lower',[-inf 0 -inf -inf -inf]);
% 
%                     yfit.h    = yfit.h  * (yymax-yymin);
%                     yfit.xo   = yfit.xo * (ttmax-ttmin) + ttmin;
%                     yfit.w    = yfit.w  * (ttmax-ttmin);
%                     yfit.yofs = yfit.yofs * (yymax-yymin) + yymin;
%                     yfit.m    = yfit.m  / (ttmax-ttmin);
%                     clear yydata yymin yymax yynorm ttdata ttmin ttmax ttnorm;


                    plot(timesec(hemiindex(nwstart:nwend)), particledata(curparticle).yfit(timesec(hemiindex(nwstart:nwend))), 'k');
% %                    plot(hsubplot(curplot), timesec(hemiindex(pHdropHemi:nwend)), yfit(timesec(hemiindex(pHdropHemi:nwend))), 'k');
%                     hemiEventFittedMean(curparticle) = yfit.xo-timesec(idxgreenframes(pHdrop));
%                     hemiEventFittedWidth(curparticle) = 2*yfit.w;
                end
                if exist('poreformation','var')  % This needs to be fixed to be in line with the hemifusion code above, right now this is old/deprecated.
                    plot(timesec(poreindex), poreformation(:,curparticle), 'Color', [0 0 1]);
                    hold on;
                    poreformationsmoothed = cumsum(porefirstderiv(:,curparticle));
                    poreformationsmoothed = poreEventStartIntensity(curparticle)+poreformationsmoothed-poreformationsmoothed(pHdropPore-NumptsAvgStartIntensity/2);
                    plot(timesec(poreindex), poreformationsmoothed, 'Color', [0 1 1]);     % [20100614TIMP] added smoothed line to poreformation too (note: this smoothed line is NOT the one used for detecting slope)
                    ndrop = maxPoreDerivIndex(curparticle); tdrop = timesec(poreindex(ndrop));
                    plot(timesec(poreindex(ndrop)), poreformationsmoothed(ndrop), ...
                          'rx','LineWidth',2,'MarkerSize',8);

                    % Gaussian-integral step fit to event trace
                    yfit = fit(timesec(poreindex(pHdropPore:end)), poreformation(pHdropPore:end,curparticle), fGaussFit, ...
                           'StartPoint', [(mean(poreformation(pHdropPore:ndrop,curparticle))-mean(poreformation(ndrop:end,curparticle))) ...
                           0 (max(timesec(poreindex))-min(timesec(poreindex)))/100 tdrop mean(poreformation(ndrop:end,curparticle))], ...
                           'Lower',[-inf 0 -inf -inf -inf]);
                    plot(timesec(poreindex(pHdropPore:end)), yfit(timesec(poreindex(pHdropPore:end))), 'k');
                    poreEventFittedMean(curparticle) = yfit.xo-timesec(idxgreenframes(pHdrop));
                    poreEventFittedWidth(curparticle) = 2*yfit.w; % 2*sigma*sqrt(2)
                end

        %         % Make uniform y axes
        %         ylim([hemiymin hemiymax]);    % [20090705TIMP] Show also smoothed curve

                % Plot vertical line at pHdrop time
                yminmax = ylim;
                plot(timesec(idxgreenframes(pHdrop*[1 1])), yminmax, 'Color', [1 1 1]*0.5);
%                plot(hsubplot(curplot), timesec(idxgreenframes(pHdrop*[1 1])), yminmax, 'Color', [1 1 1]*0.5);

                %ymin = min([min(hemifusion(pHdrop:end,curparticle)) min(poreformation(pHdrop:end,curparticle))]);
                %ymax = max([max(hemifusion(pHdrop:end,curparticle)) max(poreformation(pHdrop:end,curparticle))]);
                %axis(ylim([ymin-0.05*ymin ymax+0.05*ymax]))


                % Plot vertical line at particle detection frame center
%                hframeline(curplot) = plot(hsubplot(curplot), timesec(idxgreenframes(curframe*[1 1])), yminmax, 'Color', [1 0 0]*0.75);
                hframeline(curplot) = plot(timesec(idxgreenframes(curframe*[1 1])), yminmax, 'Color', [1 0 0]*0.75);
                set(hframeline(curplot),'ButtonDownFcn',{@frameline_ButtonDownFcn,curplot});
                %get(hframeline(curplot))


                % Color background according to left/right click selection
                % state of current particle.
                if (flagSelect(curparticle,1)), hhred  = 0.8; else hhred  = 1; end   % Left click
                if (flagSelect(curparticle,2)), hhblue = 0.8; else hhblue = 1; end   % Right click
                set(hsubplot(curplot),'Color',[hhred 1 hhblue]);

                grid on;
%                set(hsubplot(curplot), 'XGrid', 'on', 'YGrid', 'on');   % Turn on grid for specified (rather than default) axes


                % Set up ButtonDownFcn for selecting hemi/pore by
                % left/right clicking background of each subplot:
        %         ClickFcn = ['hh = get(gcbo,''Color''); set(gcbo,''Color'',[1 1 1-hh(3)])'];      % [MP-2009-02-21]
%                 ClickFcn = ['hh = get(gcbo,''Color''); kk = get(get(gcbo,''Parent''),''SelectionType''); ' ...
%                             'if(strcmp(kk,''normal'')) set(gcbo,''Color'',[1.8-hh(1) hh(2) hh(3)]); ' ...
%                             'elseif(strcmp(kk,''alt'')) set(gcbo,''Color'',[hh(1) hh(2) 1.8-hh(3)]); end'];      % [MP-2009-02-21]
%                ClickFcn = @subplotBackgroundClick_Callback;      % [MP-2009-02-21]
                set(hsubplot(curplot),'ButtonDownFcn',@subplotBackgroundClick_Callback);               % [MP-2009-02-21]
                % Selection colors are: white [1 1 1], cyan [0.8 1 1], yellow [1 1 0.8] and green [0.8 1 0.8]

%                 %set (gca, 'FontSize',3)
%                 set (gca, 'FontSize',10,'Tag',num2str(curparticle));   % [MP-2009-02-21]
                set (hsubplot(curplot), 'FontSize',10,'Tag',num2str(curparticle));   % [MP-2009-02-21]
                title(hsubplot(curplot), ['particle # ' num2str(curparticle)],'FontSize',10);

%                set(hsubplot(curplot), 'XLimMode', 'auto', 'YLimMode', 'auto'); drawnow;  % Scale axes and update right away.

                 axes(hmovie); hold on;
%                 %curplot = 1;
%                 %curparticle = nx*ny*(curPage-1) + curplot;
%                set(hmovie,'NextPlot','add');
                hmovieparticlemarker(curplot) = plot(centroid(curparticle,1), centroid(curparticle,2), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 0.6*[1 1 1]);
%                hmovieparticlemarker(curplot) = plot(hmovie, centroid(curparticle,1), centroid(curparticle,2), 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 0.6*[1 1 1]);
%                get(hmovieparticlemarker)
            end
        end
    end
end

% ButtonDownFcn for clicking on background of a subplot
function subplotBackgroundClick_Callback(hObject,eventdata)
%     hh = get(gcbo,''Color'');
%     kk = get(get(gcbo,''Parent''),''SelectionType'');
%     if(strcmp(kk,''normal''))
%         set(gcbo,''Color'',[1.8-hh(1) hh(2) hh(3)]);
%     elseif(strcmp(kk,''alt''))
%         set(gcbo,''Color'',[hh(1) hh(2) 1.8-hh(3)]);
%     end
    disp('Clicked on figure');
%    hObject = gcbo;
    hh = get(hObject,'Color');
    kk = get(get(hObject,'Parent'),'SelectionType');
    curparticle = str2num( get(hObject,'Tag') );
    if(strcmp(kk,'normal'))
        set(hObject,'Color',[1.8-hh(1) hh(2) hh(3)]);
        flagSelect(curparticle,1) = ~flagSelect(curparticle,1);  % Toggle particle on/off
    elseif(strcmp(kk,'alt'))
        set(hObject,'Color',[hh(1) hh(2) 1.8-hh(3)]);
        flagSelect(curparticle,2) = ~flagSelect(curparticle,2);  % Toggle particle on/off
    end
    guidata(hObject,eventdata);
    
    % Set the hemi and pore selection value for that particle
end

function drawMoviePlot    
    A = zeros(512,512);
    for kk = Qmov.frames(1):Qmov.frames(2)
%        A = A + double( imread(['../' moviename '.tif'],'Index',kk) );
        A = A + double( imread(['../' moviename '.tif'],'Index',kk,'Info',tiffinfo) );
    end
    A = A / (Qmov.frames(2)-Qmov.frames(1)+1);
    
    axes(hmovie);
    imagesc(A); axis image; colormap(hot);
end

% function updateMoviePlot
%     curplot = 1;
%     curparticle = nx*ny*(curPage-1) + curplot;
%     set(hmovieparticlemarker,'XData',centroid(curparticle,1));
%     set(hmovieparticlemarker,'YData',centroid(curparticle,2));
% end

function [rvec,cvec] = drawMovieZoom
    % curparticle?? kk (which frame)??
    set(hmovieparticlemarker,'MarkerEdgeColor',0.6*[1 1 1]);                % Color all particles gray except selected one white
    set(hmovieparticlemarker(curzoom-nx*ny*(curPage-1)),'MarkerEdgeColor',[1 1 1]);
    
    rvec = int16( round( centroid(curzoom,2) ) + [-1/2 1/2]*wsize );    % Set local window extents
    cvec = int16( round( centroid(curzoom,1) ) + [-1/2 1/2]*wsize );
    if (rvec(1) < 0),   rvec = rvec - rvec(1);     end                  % Manage particles near edges
    if (rvec(2) > 511), rvec = rvec - rvec(2)+511; end
    if (cvec(1) < 0),   cvec = cvec - cvec(1);     end
    if (cvec(2) > 511), cvec = cvec - cvec(2)+511; end

    for kk = curframe                                                       % Load only current frame
%        Azoom(:,:,kk) = double( imread(['../' moviename '.tif'],'Index',kk,'PixelRegion',{rvec,cvec}) );
        Azoom(:,:,kk) = double( imread(['../' moviename '.tif'],'Index',kk,'Info',tiffinfo,'PixelRegion',{rvec+1,cvec+1}) );
    end

    axes(hmoviezoom);
    hzoomimage = imagesc(cvec(1):cvec(2), rvec(1):rvec(2), Azoom(:,:,curframe)); axis image; colormap(hot);
    colorbar;
    title(sprintf('Particle #%d',curzoom));
%    imagesc(Azoom(:,:,55)); axis image; colormap(hot);
    % TODO: Set color scale strategy here.
    
    % Draw crosshairs at particle position
    hold on;
    hxhairs = plot(centroid(curzoom,1)*[1 1], rvec, '--', cvec, centroid(curzoom,2)*[1 1], '--', 'Color', 0.6*[1 1 1]);
    hold off;
    
    % Get caxis configuration and
    %if(auto is selected)
    %    set(hmoviezoom,'CLimMode','auto');
    %elseif(fixed manual)
    set(hmoviezoom,'CLimMode','manual');
    switch get(haxisscaleselectbuttongroup,'SelectedObject')
        case haxisscaleselectmax
            cc(1) = min(Azoom(:)); cc(2) = max(Azoom(:));
        case haxisscaleselect1
            cc(1) = min(Azoom(nwstart:nwend)); cc(2) = max(Azoom(nwstart:nwend));
        otherwise
            error('Invalid selection in ''Scale'' radio button group.  Terminating GUI.');
    end
    set(hmoviezoom,'CLim',cc);
    %elseif(variable hand-entered manual:   
    %end
end

function frameline_ButtonDownFcn(hObject,eventdata,plotnum)
    set(hfig,'WindowButtonUpFcn',@fig_WindowButtonUpFcn);
    set(hfig,'WindowButtonMotionFcn',{@fig_WindowButtonMotionFcn,plotnum});
%    curzoom = plotnum;
    curzoom = nx*ny*(curPage-1) + plotnum;
    [rvec,cvec] = drawMovieZoom;
    drawnow

    % Load remainder of zoomed-in movie
    if(isWholeZoomMovieLoaded ~= plotnum)
        for kk = [1:curframe-1, curframe+1:Qmov.numFrames]                      % Load all numFrames of movie for local window region
    %        Azoom(:,:,kk) = double( imread(['../' moviename '.tif'],'Index',kk,'PixelRegion',{rvec,cvec}) );
            Azoom(:,:,kk) = double( imread(['../' moviename '.tif'],'Index',kk,'Info',tiffinfo,'PixelRegion',{rvec+1,cvec+1}) );
    %        x = (kk-1)/(Qmov.numFrames-1);
    %        set(hframeline(plotnum),'Color',[cos(x*pi)^2 0 1-cos(x*pi)^2]); drawnow;
            if(mod(kk,10)==0),
                set(get(hmoviezoom,'Title'),'String',sprintf('Particle #%d: Loading frame %d of %d',curzoom,kk,Qmov.numFrames)); drawnow;
            end
        end
        set(get(hmoviezoom,'Title'),'String',sprintf('Particle #%d',curzoom));

        isWholeZoomMovieLoaded = plotnum;
    end
end

function fig_WindowButtonUpFcn(hObject,eventdata)
    set(hfig,'WindowButtonMotionFcn',[]);
    set(hfig,'WindowButtonUpFcn',[]);
end

function fig_WindowButtonMotionFcn(hObject,eventdata,plotnum)
% Get x position in axes and move red line there.
    %xpos = get(hfig,'CurrentPoint'); xpos = xpos(1);
    xpos = get(hsubplot(plotnum),'CurrentPoint'); xpos = xpos(1);
%    fprintf('Mouse is at x = %d\n', xpos);
    set(hframeline(plotnum),'XData',[1 1]*xpos);

% Get corresponding frame # and put it into moviezoom axes.
    frametime = round(xpos);
%    nframe = frametime;
    timesec = timems/1000;
% TODO: later change hemiindex to whatever is actually being shown (hemi/pore/other channel, combo, etc.)
    idx = find(abs(timesec(hemiindex)-frametime) == min(abs(timesec(hemiindex)-frametime))); idx = idx(1);
    nframe = hemiindex(idx);
    set(hzoomimage,'CData',Azoom(:,:,nframe));
end

function popupStepfittype_Callback(hObject,eventdata)
    % empty for now
end

function pushbuttonCollect_Callback(hObject,eventdata)
    % empty for now
    selectHemiWinners = find(flagSelect(:,1));
    selectPoreWinners = find(flagSelect(:,2));
    assignin('base','selectHemiWinners',selectHemiWinners);
    assignin('base','selectPoreWinners',selectPoreWinners);
    disp('Exported to base workspace the selected particle information (selectHemiWinners, selectPoreWinners).');

    % From plotall3b.m:
    %% 3) Save selected data in event file (choose filename to go after last file saved)
    if( isempty(setxor(selectHemiWinners_old,selectHemiWinners)) && ...
        isempty(setxor(selectPoreWinners_old,selectPoreWinners)) )     % [20090704TIMP] If no event selections changed...
        FLAG_SAVE = 0;
    else                                                    % ... otherwise:
        FLAG_SAVE = 1;
    end

    clear selectHemiWinners_old selectPoreWinners_old


    if(FLAG_SAVE)
        [~,dirname] = fileparts(pwd);     % Find name of current folder
        aa = dir([dirname 'events*']);
        if (isempty(aa))
            filename = [dirname 'events000'];
        else
            for kk = 1:length(aa),  aadates(kk) = datenum(aa(kk).date);  end
            idx = find(aadates == max(aadates));
            fprintf('Most recent event file is: %s\n', aa(idx).name);
            filename = aa(idx).name; [~, filename] = fileparts(filename);
            filename = [filename(1:end-3) num2str(str2double(filename(end-2:end))+1,'%03d')];
        end
        fprintf('Saving new events data to file: %s\n', [filename '.mat']);
        if(~exist([filename '.mat'],'file'))
    %         if(exist('selectHemiWinners','var')),  save(filename,'selectHemiWinners','selectHemiWinners','-append');  end
    %         if(exist('selectPoreWinners','var')),  save(filename,'selectPoreWinners','selectPoreWinners','-append');  end
            save(filename,'selectHemiWinners','selectHemiWinners','selectPoreWinners','selectPoreWinners');
        else
            fprintf([mfilename ': Can''t save event data to file because it already exists.']);
        end
        clear pathstr dirname aa kk aadates idx dummy filename ext
    else
        fprintf('No event selections changed.  Not saving results to new event file.\n');
    end

    % Export fit data to workspace:
    % Get gaussian integral fit step detection data
    hemiSelectedEventFittedMean     = hemiEventFittedMean(selectHemiWinners);
    hemiSelectedEventFittedWidth    = hemiEventFittedWidth(selectHemiWinners);
    hemiSelectedEventFittedLeftEnd  = hemiSelectedEventFittedMean - hemiSelectedEventFittedWidth/2;
    hemiSelectedEventFittedRightEnd = hemiSelectedEventFittedMean + hemiSelectedEventFittedWidth/2;
    assignin('base','hemiSelectedEventFittedMean',hemiSelectedEventFittedMean);
    assignin('base','hemiSelectedEventFittedWidth',hemiSelectedEventFittedMean);
    assignin('base','hemiSelectedEventFittedLeftEnd',hemiSelectedEventFittedMean);
    assignin('base','hemiSelectedEventFittedRightEnd',hemiSelectedEventFittedMean);

    % Get plain maxDeriv step detection data
    pHtime=timems(idxgreenframes(pHdrop))/1000;
    hemiEvent = timems(hemiindex(maxHemiDerivIndex))/1000; % get the event times for each index
    hemiEvent = hemiEvent(selectHemiWinners) - pHtime; % ti correct
    assignin('base','hemiEvent',hemiEvent);
    
    %save analysis to disk
    disp('save analysis results to');
    disp(pwd);
    saveinput = input('(y/n)?','s');
    if saveinput=='y'
       save([dataset 'auto_analysis']);
    end
end

function fitGaussianStep
    timesec = timems/1000;
    warning off;
    for kk = 1:numParticles
        if (mod(kk,10) == 0),  fprintf('\rPrefitting Gaussian step to data of particle # %d... ', kk);  end
       % yydata = hemifusion(pHdropHemi:nwend,kk);
        yydata = hemifusion(nwstart:nwend,kk);
        yymin = min(yydata); yymax = max(yydata);
        yynorm = (yydata-yymin)/(yymax-yymin);
       % ttdata = timesec(hemiindex(pHdropHemi:nwend));
        ttdata = timesec(hemiindex(nwstart:nwend));
        ttmin = min(ttdata); ttmax = max(ttdata);
        ttnorm = (ttdata - ttmin)/(ttmax-ttmin);
        yfit = fit(ttnorm, yynorm, fGaussFit, ...
              'StartPoint', [1 0 1/100 0.5 0], ...
              'Lower',[-inf 0 -inf -inf -inf]);

%         particledata(kk).yfit.h    = yfit.h  * (yymax-yymin);
%         particledata(kk).yfit.xo   = yfit.xo * (ttmax-ttmin) + ttmin;
%         particledata(kk).yfit.w    = yfit.w  * (ttmax-ttmin);
%         particledata(kk).yfit.yofs = yfit.yofs * (yymax-yymin) + yymin;
%         particledata(kk).yfit.m    = yfit.m  / (ttmax-ttmin);
% %        clear yydata yymin yymax yynorm ttdata ttmin ttmax ttnorm;
%         hemiEventFittedMean(kk) = particledata(kk).yfit.xo-timesec(idxgreenframes(pHdrop));
%         hemiEventFittedWidth(kk) = 2*particledata(kk).yfit.w;
        yfit.h    = yfit.h  * (yymax-yymin);
        yfit.xo   = yfit.xo * (ttmax-ttmin) + ttmin;
        yfit.w    = yfit.w  * (ttmax-ttmin);
        yfit.yofs = yfit.yofs * (yymax-yymin) + yymin;
        yfit.m    = yfit.m  / (ttmax-ttmin);
%        clear yydata yymin yymax yynorm ttdata ttmin ttmax ttnorm;
        hemiEventFittedMean(kk) = yfit.xo-timesec(idxgreenframes(pHdrop));
        hemiEventFittedWidth(kk) = 2*yfit.w;
        particledata(kk).yfit = yfit;
    end
    warning on;
    fprintf('\nDone pre-fitting.\n');

    % Save fit data to file [20120611]
    [~,dirname] = fileparts(pwd);     % Find name of current folder
    aa = dir([dirname 'fits*']);
    if (isempty(aa))
        filename = [dirname 'fits00'];
    else
        for kk = 1:length(aa),  aadates(kk) = datenum(aa(kk).date);  end
        idx = find(aadates == max(aadates));
        fprintf('Most recent fits file is: %s\n', aa(idx).name);
        filename = aa(idx).name; [~, filename] = fileparts(filename);
        filename = [filename(1:end-2) num2str(str2double(filename(end-1:end))+1,'%02d')];
    end
    fprintf('Saving new fits data to file: %s\n', [filename '.mat']);
    if(~exist([filename '.mat'],'file'))
        save(filename,'hemiEventFittedMean','hemiEventFittedWidth','particledata');
    else
        fprintf([mfilename ': Can''t save fits data to file because it already exists.']);
    end
end

end
