% Save and plot the histogram of inner axonal diameter in corpus callosum
% based on Fig. 4 in (Aboitiz, et al., Brain Research, 1992, 598:143-153).

filenames = {'genu.mat','body_anterior.mat','body_middle.mat','body_posterior.mat','splenium.mat'};
roinames = {'GENU','ANT.BODY','MID.BODY','POST.BODY','SPLENIUM'};
cc = struct([]);

figure; set(gcf,'unit','inch','position',[0,0,4,12]);
for i = 1:numel(filenames)
    % load file
    filename = filenames{i};
    data = load(filename);
    
    cc(i).name = roinames{i};
    
    % inner axonal diameter
    di = data.dist(:,1);
    cc(i).diameter = 0.2*round(di/0.2);
    
    % frequency/normalized count
    ci = data.dist(:,2);
    cc(i).frequency = ci/sum(ci);
    
    subplot(5,1,i)
    h = bar(cc(i).diameter,cc(i).frequency*100,1);
    title(cc(i).name)
    xlim([0 9]); ylim([0 30])
    box on; pbaspect([2 1 1])
    set(gca,'xtick',0:9,'ytick',0:5:30);
    if i == 5, xlabel('Inner Diameter (µm)','fontsize',16); end
    if i == 3, ylabel('Frequency (%)','fontsize',16); end
end

save(fullfile('CC_diameter_histogram.mat'),'cc')