% 6. Post-processing of the selected ROIs

roidata.image = roiwindow.image;
roidata.roi   = roiwindow.roi;
roidata.labels= roiwindow.labels;
roidata.number= roiwindow.number;

for kk = 1:roidata.number
    bwmask = (roidata.labels == kk);
    a2masked = a2 .* double(bwmask);

    particle(kk).intensity = sum(a2masked(:));
    particle(kk).stats = regionprops(bwmask,'Centroid');

%     figure; imagesc(a2masked); axis image; colormap hot; colorbar;
%     title(sprintf('Particle %d',kk));
end

for kk = 1:roidata.number
    partintvector(kk) = particle(kk).intensity;                                 % Particle intensity-data vector
end
figure; hist(partintvector)
partintvector

save([moviefile(1:end-4) '_processedIntensityData'],'particle','partintvector');
