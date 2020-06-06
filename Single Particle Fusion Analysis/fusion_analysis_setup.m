% viral fusion analysis configuration script
%savefolder = 'C:\Documents and Settings\Irene Kim\My Documents\MATLAB\fusion_data_analysis';
%savefolder = 'D:\Tijana\posao\steve harrison\matlabcodes\fusion_data_analysis';
savefolder = fileparts(which('fusion_analysis_setup'));     % [090218 Milos]

update_config = 'n';
while update_config =='n';
    DataID = {[],'poreformation','hemifusion'};
    signalID = {'dequenching','dissipation'};
    % setup red channel
    disp('what signal is contained in the RED channel?')
    disp(' 0 = no signal   1 = content dye    2 = lipid dye')
    channel.red.ID = DataID{input('--> ')+1};

    if channel.red.ID ~=0
        disp('event is a dequenching spike or signal dissipation?')
        disp('dequenching = 1    dissipation = 2')
        channel.red.signal = signalID{input('-->')};
    else channel.red.signal=[];
    end
    %setup green channel
    disp('what signal is contained in the GREEN channel?')
    disp(' 0 = no signal   1 = content dye    2 = lipid dye')
    channel.green.ID =DataID{input('--> ')+1};

    if channel.green.ID ~=0
        disp('event is a dequenching spike or signal dissipation?')
        disp('dequenching = 1    dissipation = 2')
        channel.green.signal = signalID{input('-->')};
        else channel.green.signal=[];
    end
    %%
    disp('analysis settings:')
    disp('red channel')
    disp(channel.red)
    disp(' ')
    disp('green channel')
    disp(channel.green)
    disp(' ')
    if strcmp(channel.red.ID,channel.green.ID)==1;
        warning('content/lipid signals cannot be in the same channel!') %#ok<WNTAG>
    end
    disp(' ')
    disp('is this correct?')
    update_config= input('y/n -->','s');
    if update_config=='n'
        disp('restarting setup...')
        disp('.')
        disp('.')
        disp('.')
    end
end

save([savefolder '\' 'fusion_config_settings'],'channel')
clear DataID signal ID update_congfig


