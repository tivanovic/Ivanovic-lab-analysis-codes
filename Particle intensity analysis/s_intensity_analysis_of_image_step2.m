% 4. Select intensity template file to load (e.g. '190612_032_intensitytemplate.mat')
%filename = '190612_032_intensitytemplate.tif'; moviefile = filename;
[filename,pathname]=uigetfile('*.mat','Please select processed movie file still image:');
if isequal(filename,0)
   disp('User selected Cancel')
   return
end
% else
   moviefile = fullfile(pathname, filename);
   disp(['User selected ', moviefile])
   Q=load(moviefile);
   a2=Q.a2; clear Q;
% end


% 5. Use ROI editor to generate particle masks

roiwindow = croieditor2(log10(abs(a2)));
caxis([2.4 4.4]);
